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
#include "config.hpp"
#endif

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <cfp/cfp.hpp>
#include <gxrio.hpp>

#include <cif++/pdb/io.hpp>

#include "DSSP.hpp"

#include "dssp_wrapper.hpp"
#include "revision.hpp"

namespace fs = std::filesystem;

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what(const std::exception &e)
{
	std::cerr << e.what() << std::endl;
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

	auto &config = cfp::config::instance();

	config.init("Usage: mkdssp [options] input-file [output-file]",
		cfp::make_option<std::string>("output-format", "Output format, can be either 'dssp' for classic DSSP or 'mmcif' for annotated mmCIF. The default is chosen based on the extension of the output file, if any."),
		cfp::make_option<short>("min-pp-stretch", 3, "Minimal number of residues having PSI/PHI in range for a PP helix, default is 3"),
		cfp::make_option("write-other", "If set, write the type OTHER for loops, default is to leave this out"),

		// cfp::make_option("components",			po::value<std::string,	"Location of the components.cif file from CCD")
	    // cfp::make_option("extra-compounds",		po::value<std::string,	"File containing residue information for extra compounds in this specific target, should be either in CCD format or a CCP4 restraints file")
		cfp::make_option<std::string>("mmcif-dictionary", "Path to the mmcif_pdbx.dic file to use instead of default"),

		cfp::make_option("help,h", "Display help message"),
		cfp::make_option("version", "Print version"),
		cfp::make_option("verbose,v", "verbose output"),

		cfp::make_hidden_option<int>("debug,d", "Debug level (for even more verbose output)"));

	config.parse(argc, argv);

	// --------------------------------------------------------------------

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help"))
	{
		std::cerr << config << std::endl;
		exit(0);
	}

	if (config.operands().empty())
	{
		std::cerr << "Input file not specified" << std::endl;
		exit(1);
	}

	if (config.has("output-format") and config.get<std::string>("output-format") != "dssp" and config.get<std::string>("output-format") != "mmcif")
	{
		std::cerr << "Output format should be one of 'dssp' or 'mmcif'" << std::endl;
		exit(1);
	}

	cif::VERBOSE = config.count("verbose");
	if (config.has("debug"))
		cif::VERBOSE = config.get<int>("debug");

	// --------------------------------------------------------------------

	// Load extra CCD definitions, if any

	// if (config.has("compounds"))
	// 	cif::add_file_resource("components.cif", config.get<std::string>("compounds"));
	// else if (config.has("components"))
	// 	cif::add_file_resource("components.cif", config.get<std::string>("components"));

	// if (config.has("extra-compounds"))
	// 	mmcif::CompoundFactory::instance().pushDictionary(config.get<std::string>("extra-compounds"));

	// And perhaps a private mmcif_pdbx dictionary
	if (config.has("mmcif-dictionary"))
		cif::add_file_resource("mmcif_pdbx.dic", config.get<std::string>("mmcif-dictionary"));

	gxrio::ifstream in(config.operands().front());
	if (not in.is_open())
	{
		std::cerr << "Could not open file" << std::endl;
		exit(1);
	}

	cif::file f = cif::pdb::read(in);
	if (not f.is_valid())
	{
		std::cerr << "Could not validate file" << std::endl;
		exit(1);
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
			output = output.stem();

		if (output.extension() == ".dssp")
			fmt = "dssp";
		else
			fmt = "cif";
	}

	dssp::DSSP dssp(f.front(), 1, pp_stretch, fmt == "dssp");

	if (not output.empty())
	{
		gxrio::ofstream out(output);

		if (not out.is_open())
		{
			std::cerr << "Could not open output file" << std::endl;
			exit(1);
		}

		if (fmt == "dssp")
			writeDSSP(dssp, out);
		else
			annotateDSSP(f.front(), dssp, writeOther, out);
	}
	else
	{
		if (fmt == "dssp")
			writeDSSP(dssp, std::cout);
		else
			annotateDSSP(f.front(), dssp, writeOther, std::cout);
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
