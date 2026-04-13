/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2020 NKI/AVL, Netherlands Cancer Institute
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if __has_include("config.hpp")
# include "config.hpp"
#endif

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <cif++/cif++.hpp>
#include <mcfp/mcfp.hpp>

#include "dssp.hpp"
#include "revision.hpp"

namespace fs = std::filesystem;

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what(const std::exception &e)
{
	std::cerr << e.what() << '\n';
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception &nested)
	{
		std::cerr << " >> ";
		print_what(nested);
	}
}

// --------------------------------------------------------------------

int d_main(int argc, const char *argv[])
{
	using namespace std::literals;

	auto &config = mcfp::config::instance();

	config.init("Usage: mkdssp [options] input-file [output-file]",
		mcfp::make_option<std::string>("output-format", "Output format, can be either 'dssp' for classic DSSP or 'mmcif' for annotated mmCIF. The default is chosen based on the extension of the output file, if any."),
		mcfp::make_option<short>("min-pp-stretch", 3, "Minimal number of residues having PSI/PHI in range for a PP helix, default is 3"),
		mcfp::make_option("write-other", "If set, write the type OTHER for loops, default is to leave this out"),
		mcfp::make_option("no-dssp-categories", "If set, will suppress output of new DSSP output in mmCIF format"),

		mcfp::make_option("calculate-accessibility", "Default is to not calculate the surface accessibility when the output format is mmCIF"),

		mcfp::make_option<std::string>("mmcif-dictionary", "Path to the mmcif_pdbx.dic file to use instead of default"),

		mcfp::make_option("help,h", "Display help message"),
		mcfp::make_option("version", "Print version"),
		mcfp::make_option("verbose,v", "verbose output"),
		mcfp::make_option("quiet", "Reduce verbose output to a minimum"),

		mcfp::make_hidden_option<int>("debug,d", "Debug level (for even more verbose output)"));

	config.parse(argc, argv);

	// --------------------------------------------------------------------

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help") or config.operands().empty())
	{
		std::cerr << config << '\n';
		exit(config.has("help") ? 0 : 1);
	}

	if (config.has("output-format") and config.get<std::string>("output-format") != "dssp" and config.get<std::string>("output-format") != "mmcif")
	{
		std::cerr << "Output format should be one of 'dssp' or 'mmcif'\n";
		exit(1);
	}

	if (config.count("quiet"))
		cif::VERBOSE = -1;
	else
		cif::VERBOSE = config.count("verbose");

	// --------------------------------------------------------------------

	// private mmcif_pdbx dictionary?
	if (config.has("mmcif-dictionary"))
	{
		fs::path mmcif_dict = config.get<std::string>("mmcif-dictionary");

		cif::add_file_resource("mmcif_pdbx.dic", mmcif_dict);

		// Try to be smart, maybe dssp-extension.dic is at that location as well?
		auto dir = fs::canonical(mmcif_dict.parent_path());
		if (auto dssp_dict = cif::load_resource("dssp-extension.dic"); dssp_dict == nullptr and fs::exists(dir / "dssp-extension.dic"))
			cif::add_data_directory(dir);
	}

	cif::file f;

	try
	{
		cif::gzio::ifstream in(config.operands().front());
		if (not in.is_open())
		{
			std::cerr << "Could not open file\n";
			exit(1);
		}

		if (cif::VERBOSE > 0)
			std::cerr << "Loading file...";

		f.load(in);

		if (cif::VERBOSE > 0)
			std::cerr << " fixup file...";

		cif::pdb::fixup_pdbx(f);

		if (cif::VERBOSE > 0)
			std::cerr << " done\n";
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';

		f = cif::pdb::read(config.operands().front());
	}

	// --------------------------------------------------------------------

	short pp_stretch = 3;
	if (config.has("min-pp-stretch"))
		pp_stretch = config.get<short>("min-pp-stretch");

	bool writeOther = config.has("write-other");

	std::string fmt;
	if (config.has("output-format"))
		fmt = config.get<std::string>("output-format");

	fs::path output;
	if (config.operands().size() > 1)
		output = config.operands()[1];

	if (fmt.empty() and not output.empty())
	{
		if (output.extension() == ".gz" or output.extension() == ".xz")
		{
			if (output.stem().extension() == ".dssp")
				fmt = "dssp";
			else
				fmt = "cif";
		}
		else if (output.extension() == ".dssp")
			fmt = "dssp";
		else
			fmt = "cif";
	}

	if (fmt == "dssp")
	{
		// See if the data will fit at all
		auto &db = f.front();
		for (const auto &[chain_id, seq_nr] : db["pdbx_poly_seq_scheme"].rows<std::string, int>("pdb_strand_id", "pdb_seq_num"))
		{
			if (chain_id.length() > 1 or seq_nr > 99999)
			{
				std::cerr << "The data in this file won't fit in the old DSSP format, please use the mmCIF format instead.\n";
				exit(2);
			}
		}
	}

	dssp dssp(f.front(), 1, pp_stretch, fmt == "dssp" or config.has("calculate-accessibility"));

	if (not output.empty())
	{
		cif::gzio::ofstream out(output);

		if (not out.is_open())
		{
			std::cerr << "Could not open output file\n";
			exit(1);
		}

		if (fmt == "dssp")
			dssp.write_legacy_output(out);
		else
		{
			dssp.annotate(f.front(), writeOther, not config.has("no-dssp-categories"));
			out << f.front();
		}
	}
	else
	{
		if (fmt == "dssp")
			dssp.write_legacy_output(std::cout);
		else
		{
			dssp.annotate(f.front(), writeOther, not config.has("no-dssp-categories"));
			std::cout << f.front();
		}
	}

	return 0;
}

// --------------------------------------------------------------------

int main(int argc, const char *argv[])
{
	int result = 0;

	try
	{
#if defined(DATA_DIR)
		cif::add_data_directory(DATA_DIR);
#endif
		result = d_main(argc, argv);
	}
	catch (const std::exception &ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
