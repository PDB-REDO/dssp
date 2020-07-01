//           Copyright Maarten L. Hekkelman 2020
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <exception>
#include <iostream>
#include <fstream>

#include <boost/format.hpp>
#include <boost/date_time/gregorian/formatters.hpp>

#include <cif++/Config.hpp>
#include <cif++/Structure.hpp>
#include <cif++/Secondary.hpp>
#include <cif++/CifUtils.hpp>
#include <cif++/Cif2PDB.hpp>
#include <cif++/FixDMC.hpp>

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

std::string ResidueToDSSPLine(const mmcif::DSSP::ResidueInfo& info)
{
/*
	This is the header line for the residue lines in a DSSP file:

	#  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
 */
	boost::format kDSSPResidueLine(
		"%5.5d%5.5d%1.1s%1.1s %c  %c %c%c%c%c%c%c%c%4.4d%4.4d%c%4.4d %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f");

	auto& residue = info.residue();

	if (residue.asymID().length() > 1)
		throw std::runtime_error("This file contains data that won't fit in the original DSSP format");

	char code = 'X';
	if (mmcif::kAAMap.find(residue.compoundID()) != mmcif::kAAMap.end())
		code = mmcif::kAAMap.at(residue.compoundID());

	if (code == 'C')	// a cysteine
	{
		auto ssbridgenr = info.ssBridgeNr();
		if (ssbridgenr)
			code = 'a' + ((ssbridgenr - 1) % 26);
	}

	char ss;
	switch (info.ss())
	{
		case mmcif::ssAlphahelix:	ss = 'H'; break;
		case mmcif::ssBetabridge:	ss = 'B'; break;
		case mmcif::ssStrand:		ss = 'E'; break;
		case mmcif::ssHelix_3:		ss = 'G'; break;
		case mmcif::ssHelix_5:		ss = 'I'; break;
		case mmcif::ssTurn:			ss = 'T'; break;
		case mmcif::ssBend:			ss = 'S'; break;
		case mmcif::ssLoop:			ss = ' '; break;
	}

	char helix[3] = { ' ', ' ', ' ' };
	for (auto stride: { 3, 4, 5})
	{
		switch (info.helix(stride))
		{
			case mmcif::Helix::None:		helix[stride - 3] = ' '; break;
			case mmcif::Helix::Start:		helix[stride - 3] = '>'; break;
			case mmcif::Helix::End:			helix[stride - 3] = '<'; break;
			case mmcif::Helix::StartAndEnd:	helix[stride - 3] = 'X'; break;
			case mmcif::Helix::Middle:		helix[stride - 3] = '0' + stride; break;
		}
	}

	char bend = ' ';
	if (info.bend())
		bend = 'S';

	double alpha = residue.alpha();
	char chirality = alpha == 360 ? ' ' : (alpha < 0 ? '-' : '+');
	
	uint32_t bp[2] = {};
	char bridgelabel[2] = { ' ', ' ' };
	for (uint32_t i: { 0, 1 })
	{
		const auto& [p, ladder, parallel] = info.bridgePartner(i);
		if (not p)
			continue;

		bp[i] = p.nr() % 10000;  // won't fit otherwise...
		bridgelabel[i] = (parallel ? 'a' : 'A') + ladder % 26;
	}

	char sheet = ' ';
	if (info.sheet() != 0)
		sheet = 'A' + (info.sheet() - 1) % 26;

	std::string NHO[2], ONH[2];
	for (int i: { 0, 1 })
	{
		const auto& [donor, donorE] = info.donor(i);
		const auto& [acceptor, acceptorE] = info.acceptor(i);
		
		NHO[i] = ONH[i] = "0, 0.0";

		if (acceptor)
		{
			auto d = acceptor.nr() - info.nr();
			NHO[i] = (boost::format("%d,%3.1f") % d % acceptorE).str();
		}

		if (donor)
		{
			auto d = donor.nr() - info.nr();
			ONH[i] = (boost::format("%d,%3.1f") % d % donorE).str();
		}
	}

	auto ca = residue.atomByID("CA");
	auto const& [cax, cay, caz] = ca.location();

	return (kDSSPResidueLine % info.nr() % ca.authSeqID() % ca.pdbxAuthInsCode() % ca.authAsymID() % code %
		ss % helix[0] % helix[1] % helix[2] % bend % chirality % bridgelabel[0] % bridgelabel[1] %
		bp[0] % bp[1] % sheet % floor(info.accessibility() + 0.5) %
		NHO[0] % ONH[0] % NHO[1] % ONH[1] %
		residue.tco() % residue.kappa() % alpha % residue.phi() % residue.psi() %
		cax % cay % caz).str();
}

void writeDSSP(const mmcif::Structure& structure, const mmcif::DSSP& dssp, std::ostream& os)
{
	const std::string kFirstLine("==== Secondary Structure Definition by the program DSSP, NKI version 3.0                           ==== ");
	boost::format kHeaderLine("%1% %|127t|%2%");

	using namespace boost::gregorian;

	auto stats = dssp.GetStatistics();

	date today = day_clock::local_day();

	auto& cf = structure.getFile().file();

	os << kHeaderLine % (kFirstLine + "DATE=" + to_iso_extended_string(today)) % '.' << std::endl
	   << kHeaderLine % "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637" % '.' << std::endl
	   << GetPDBHEADERLine(cf, 127) << '.' << std::endl
	   << GetPDBCOMPNDLine(cf, 127) << '.' << std::endl
	   << GetPDBSOURCELine(cf, 127) << '.' << std::endl
	   << GetPDBAUTHORLine(cf, 127) << '.' << std::endl;

	os << boost::format("%5.5d%3.3d%3.3d%3.3d%3.3d TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN) %|127t|%c") %
			 stats.nrOfResidues % stats.nrOfChains % stats.nrOfSSBridges % stats.nrOfIntraChainSSBridges % (stats.nrOfSSBridges - stats.nrOfIntraChainSSBridges) % '.' << std::endl;
		 os << kHeaderLine % (boost::format("%8.1f   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)") % stats.accessibleSurface) % '.' << std::endl;

	// hydrogenbond summary

	os << kHeaderLine % (
		boost::format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES")
			% stats.nrOfHBonds % (stats.nrOfHBonds * 100.0 / stats.nrOfResidues)) % '.' << std::endl;

	os << kHeaderLine % (
		boost::format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES")
			% stats.nrOfHBondsInParallelBridges % (stats.nrOfHBondsInParallelBridges * 100.0 / stats.nrOfResidues)) % '.' << std::endl;

	os << kHeaderLine % (
		boost::format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES")
			% stats.nrOfHBondsInAntiparallelBridges % (stats.nrOfHBondsInAntiparallelBridges * 100.0 / stats.nrOfResidues)) % '.' << std::endl;

	boost::format kHBondsLine("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I%c%1.1d), SAME NUMBER PER 100 RESIDUES");
	for (int k = 0; k < 11; ++k)
		os << kHeaderLine % (kHBondsLine % stats.nrOfHBondsPerDistance[k] % (stats.nrOfHBondsPerDistance[k] * 100.0 / stats.nrOfResidues) % (k - 5 < 0 ? '-' : '+') % abs(k - 5)) % '.' << std::endl;

	// histograms...
	os << "  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           ." << std::endl;

	for (auto hi: stats.residuesPerAlphaHelixHistogram)
		os << boost::format("%3.3d") % hi;
	os << "    RESIDUES PER ALPHA HELIX         ." << std::endl;

	for (auto hi: stats.parallelBridgesPerLadderHistogram)
		os << boost::format("%3.3d") % hi;
	os << "    PARALLEL BRIDGES PER LADDER      ." << std::endl;

	for (auto hi: stats.antiparallelBridgesPerLadderHistogram)
		os << boost::format("%3.3d") % hi;
	os << "    ANTIPARALLEL BRIDGES PER LADDER  ." << std::endl;

	for (auto hi: stats.laddersPerSheetHistogram)
		os << boost::format("%3.3d") % hi;
	os << "    LADDERS PER SHEET                ." << std::endl;

	// per residue information

	os << "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA         CHAIN" << std::endl;
	boost::format kDSSPResidueLine(
		"%5.5d        !%c             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0");

	int last = 0;
	for (auto ri: dssp)
	{
		// insert a break line whenever we detect missing residues
		// can be the transition to a different chain, or missing residues in the current chain

		if (ri.nr() != last + 1)
			os << (kDSSPResidueLine % (last + 1) % (ri.chainBreak() == mmcif::ChainBreak::NewChain ? '*' : ' ')) << std::endl;
		
		os << ResidueToDSSPLine(ri) << std::endl;
		last = ri.nr();
	}

	// std::vector<const MResidue*> residues;

	// foreach (const MChain* chain, protein.GetChains())
	// {
	// 	foreach (const MResidue* residue, chain->GetResidues())
	// 		residues.push_back(residue);
	// }

	// // keep residues sorted by residue number as assigned during reading the PDB file
	// sort(residues.begin(), residues.end(), boost::bind(&MResidue::GetNumber, _1) < boost::bind(&MResidue::GetNumber, _2));

	// const MResidue* last = nullptr;
	// foreach (const MResidue* residue, residues)
	// {
	// 	// insert a break line whenever we detect missing residues
	// 	// can be the transition to a different chain, or missing residues in the current chain
	// 	if (last != nullptr and last->GetNumber() + 1 != residue->GetNumber())
	// 	{
	// 		char breaktype = ' ';
	// 		if (last->GetChainID() != residue->GetChainID())
	// 			breaktype = '*';
	// 		os << (kDSSPResidueLine % (last->GetNumber() + 1) % breaktype) << std::endl;
	// 	}
	// 	os << ResidueToDSSPLine(*residue) << std::endl;
	// 	last = residue;
	// }
}

void annotateDSSP(const mmcif::Structure& structure, const mmcif::DSSP& dssp, std::ostream& os)
{

}

// --------------------------------------------------------------------

int d_main(int argc, const char* argv[])
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

		("create-missing",									"Create missing backbone atoms")

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

	if (vm.count("output-format") and vm["output-format"].as<std::string>() != "dssp" and vm["output-format"].as<std::string>() != "mmcif")
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
	mmcif::Structure structure(f, mmcif::StructureOpenOptions::SkipHydrogen);

	// --------------------------------------------------------------------
	
	if (vm.count("create-missing"))
		mmcif::CreateMissingBackboneAtoms(structure, true);

	// --------------------------------------------------------------------

	mmcif::DSSP dssp(structure);

	auto fmt = vm.count("output-format") ? vm["output-format"].as<std::string>() : "";
	
	if (vm.count("output"))
	{
		std::ofstream of(vm["output"].as<std::string>());
		if (not of.is_open())
		{
			std::cerr << "Could not open output file" << std::endl;
			exit(1);
		}
		
		if (fmt == "dssp")
			writeDSSP(structure, dssp, of);
		else
			annotateDSSP(structure, dssp, of);
	}
	else
	{
		if (fmt == "dssp")
			writeDSSP(structure, dssp, std::cout);
		else
			annotateDSSP(structure, dssp, std::cout);
	}
	
	return 0;
}

// --------------------------------------------------------------------

int main(int argc, const char* argv[])
{
	int result = 0;

	try
	{
		result = d_main(argc, argv);
	}
	catch (const std::exception& ex)
	{
		print_what(ex);
		exit(1);
	}

	return result;
}
