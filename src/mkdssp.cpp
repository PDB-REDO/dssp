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
#include <iostream>
#include <filesystem>
#include <fstream>

#include <boost/format.hpp>
#include <boost/date_time/gregorian/formatters.hpp>

#include <cif++/Structure.hpp>
#include <cif++/Secondary.hpp>
#include <cif++/CifUtils.hpp>
#include <cif++/Cif2PDB.hpp>

#include <boost/program_options.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "dssp.hpp"

namespace fs = std::filesystem;
namespace io = boost::iostreams;
namespace po = boost::program_options;

// --------------------------------------------------------------------

// recursively print exception whats:
void print_what (const std::exception& e)
{
	std::cerr << e.what() << std::endl;
	try
	{
		std::rethrow_if_nested(e);
	}
	catch (const std::exception& nested)
	{
		std::cerr << " >> ";
		print_what(nested);
	}
}

// --------------------------------------------------------------------

int d_main(int argc, const char* argv[])
{
	using namespace std::literals;

	po::options_description visible_options(argv[0] + " [options] input-file [output-file]"s);
	visible_options.add_options()
		("dict",				po::value<std::vector<std::string>>(),
															"Dictionary file containing restraints for residues in this specific target, can be specified multiple times.")
		("output-format",		po::value<std::string>(),	"Output format, can be either 'dssp' for classic DSSP or 'mmcif' for annotated mmCIF. The default is chosen based on the extension of the output file, if any.")
		("min-pp-stretch",		po::value<short>(),			"Minimal number of residues having PSI/PHI in range for a PP helix, default is 3")
		("write-other",										"If set, write the type OTHER for loops, default is to leave this out")

		("help,h",											"Display help message")
		("version",											"Print version")
		("verbose,v",										"verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("xyzin,i",				po::value<std::string>(),	"coordinates file")
		("output,o",            po::value<std::string>(),	"Output to this file")
		("debug,d",				po::value<int>(),			"Debug level (for even more verbose output)")

		("compounds",			po::value<std::string>(),	"Location of the components.cif file from CCD")
		("components",			po::value<std::string>(),	"Location of the components.cif file from CCD, alias")
		("extra-compounds",		po::value<std::string>(),	"File containing residue information for extra compounds in this specific target, should be either in CCD format or a CCP4 restraints file")
		("mmcif-dictionary",	po::value<std::string>(),	"Path to the mmcif_pdbx.dic file to use instead of default")
		;

	po::options_description cmdline_options;
	cmdline_options.add(visible_options).add(hidden_options);

	po::positional_options_description p;
	p.add("xyzin", 1);
	p.add("output", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(p).run(), vm);

	po::notify(vm);

	// --------------------------------------------------------------------

	if (vm.count("version"))
	{
		std::cout << argv[0] << " version " << get_version_string() << std::endl;
		exit(0);
	}

	if (vm.count("help"))
	{
		std::cerr << visible_options << std::endl;
		exit(0);
	}
	
	if (vm.count("xyzin") == 0)
	{
		std::cerr << "Input file not specified" << std::endl;
		exit(1);
	}

	if (vm.count("output-format") and vm["output-format"].as<std::string>() != "dssp" and vm["output-format"].as<std::string>() != "mmcif")
	{
		std::cerr << "Output format should be one of 'dssp' or 'mmcif'" << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	// --------------------------------------------------------------------

	// Load extra CCD definitions, if any

	if (vm.count("compounds"))
		cif::addFileResource("components.cif", vm["compounds"].as<std::string>());
	else if (vm.count("components"))
		cif::addFileResource("components.cif", vm["components"].as<std::string>());
	
	if (vm.count("extra-compounds"))
		mmcif::CompoundFactory::instance().pushDictionary(vm["extra-compounds"].as<std::string>());
	
	// And perhaps a private mmcif_pdbx dictionary

	if (vm.count("mmcif-dictionary"))
		cif::addFileResource("mmcif_pdbx_v50.dic", vm["mmcif-dictionary"].as<std::string>());

	if (vm.count("dict"))
	{
		for (auto dict: vm["dict"].as<std::vector<std::string>>())
			mmcif::CompoundFactory::instance().pushDictionary(dict);
	}

	mmcif::File f(vm["xyzin"].as<std::string>());
	mmcif::Structure structure(f, 1, mmcif::StructureOpenOptions::SkipHydrogen);

	// --------------------------------------------------------------------

	short pp_stretch = 3;
	if (vm.count("min-pp-stretch"))
		pp_stretch = vm["min-pp-stretch"].as<short>();

	bool writeOther = vm.count("write-other");

	std::string fmt;
	if (vm.count("output-format"))
		fmt = vm["output-format"].as<std::string>();

	if (fmt.empty() and vm.count("output"))
	{
		fs::path output = vm["output"].as<std::string>();
	
		if (output.extension() == ".gz")
			output = output.stem();
		else if (output.extension() == ".bz2")
			output = output.stem();

		if (output.extension() == ".dssp")
			fmt = "dssp";
		else
			fmt = "cif";
	}


	mmcif::DSSP dssp(structure, pp_stretch, fmt == "dssp");

	if (vm.count("output"))
	{
		fs::path output = vm["output"].as<std::string>();
		std::ofstream of(output, std::ios_base::out | std::ios_base::binary);

		if (not of.is_open())
		{
			std::cerr << "Could not open output file" << std::endl;
			exit(1);
		}

		io::filtering_stream<io::output> out;
		
		if (output.extension() == ".gz")
			out.push(io::gzip_compressor());

		out.push(of);

		if (fmt == "dssp")
			writeDSSP(structure, dssp, out);
		else
			annotateDSSP(structure, dssp, writeOther, out);
	}
	else
	{
		if (fmt == "dssp")
			writeDSSP(structure, dssp, std::cout);
		else
			annotateDSSP(structure, dssp, writeOther, std::cout);
	}
	
	return 0;
}

// --------------------------------------------------------------------

int main(int argc, const char* argv[])
{
	int result = 0;

	try
	{
#if defined(DATA_DIR)
		cif::addDataDirectory(DATA_DIR);
#endif
		load_version_info();

		result = d_main(argc, argv);
	}
	catch (const std::exception& ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
