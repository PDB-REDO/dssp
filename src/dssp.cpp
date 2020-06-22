//           Copyright Maarten L. Hekkelman 2020
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <exception>
#include <iostream>
#include <fstream>

#include <cif++/Config.h>
#include <cif++/Structure.h>
#include <cif++/Secondary.h>
#include <cif++/CifUtils.h>

#include <boost/program_options.hpp>

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

void writeDSSP(const mmcif::Structure& structure, std::ostream& os)
{
	mmcif::DSSP dssp(structure);

}

void annotateDSSP(const mmcif::Structure& structure, std::ostream& os)
{

}

// --------------------------------------------------------------------

int main(int argc, char* argv[])
{
	using namespace std::literals;

	po::options_description visible_options(argv[0] + " input-file [output-file] [options]"s);
	visible_options.add_options()
		("xyzin",				po::value<std::string>(),	"coordinates file")
		("output",              po::value<std::string>(),	"Output to this file")

		("dict",				po::value<std::vector<std::string>>(),
															"Dictionary file containing restraints for residues in this specific target, can be specified multiple times.")

		("help,h",											"Display help message")
		("version",											"Print version")

		("output-format",		po::value<std::string>(),	"Output format, can be either 'dssp' for classic DSSP or 'mmcif' for annotated mmCIF. The default is chosen based on the extension of the output file, if any.")

#if not USE_RSRC
		("rsrc-dir",			po::value<std::string>(),	"Directory containing the 'resources' used by this application")
#endif

		("verbose,v",										"verbose output")
		;
	
	po::options_description hidden_options("hidden options");
	hidden_options.add_options()
		("debug,d",				po::value<int>(),			"Debug level (for even more verbose output)")
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
		std::cout << argv[0] << " version " << VERSION << std::endl;
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

	if (vm.count("output-format") and vm["output-format"].as<std::string>() != "dsps" and vm["output-format"].as<std::string>() != "mmcif")
	{
		std::cerr << "Output format should be one of 'dssp' or 'mmcif'" << std::endl;
		exit(1);
	}

	cif::VERBOSE = vm.count("verbose") != 0;
	if (vm.count("debug"))
		cif::VERBOSE = vm["debug"].as<int>();

	// --------------------------------------------------------------------
	
#if USE_RSRC
	cif::rsrc_loader::init();
#else
	std::string rsrc_dir = ".";
	cif::rsrc_loader::init({ { cif::rsrc_loader_type::file, rsrc_dir }});
#endif

	// --------------------------------------------------------------------
	
	if (vm.count("dict"))
	{
		for (auto dict: vm["dict"].as<std::vector<std::string>>())
			mmcif::CompoundFactory::instance().pushDictionary(dict);
	}

	mmcif::File f(vm["xyzin"].as<std::string>());
	mmcif::Structure structure(f);
	
	// --------------------------------------------------------------------

	auto fmt = vm["output-format"].as<std::string>();
	
	if (vm.count("output"))
	{
		std::ofstream of(vm["output"].as<std::string>());
		if (not of.is_open())
		{
			std::cerr << "Could not open output file" << std::endl;
			exit(1);
		}
		
		if (fmt == "dssp")
			writeDSSP(structure, of);
		else
			annotateDSSP(structure, of);
	}
	else
	{
		if (fmt == "dssp")
			writeDSSP(structure, std::cout);
		else
			annotateDSSP(structure, std::cout);
	}
	
	return 0;
}
