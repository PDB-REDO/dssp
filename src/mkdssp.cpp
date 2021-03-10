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
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>

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

namespace {
	std::string gVersionNr, gVersionDate, VERSION_STRING;
}

void load_version_info()
{
	const std::regex
		rxVersionNr(R"(build-(\d+)-g[0-9a-f]{7}(-dirty)?)"),
		rxVersionDate(R"(Date: +(\d{4}-\d{2}-\d{2}).*)"),
		rxVersionNr2(R"(mkdssp-version: (\d+(?:\.\d+)+))");

#if __has_include("revision.hpp")
#include "revision.hpp"
#else
	const char* kRevision = "";
#endif

	struct membuf : public std::streambuf
	{
		membuf(char* data, size_t length)       { this->setg(data, data, data + length); }
	} buffer(const_cast<char*>(kRevision), sizeof(kRevision));

	std::istream is(&buffer);

	std::string line;

	while (getline(is, line))
	{
		std::smatch m;

		if (std::regex_match(line, m, rxVersionNr))
		{
			gVersionNr = m[1];
			if (m[2].matched)
				gVersionNr += '*';
			continue;
		}

		if (std::regex_match(line, m, rxVersionDate))
		{
			gVersionDate = m[1];
			continue;
		}

		// always the first, replace with more specific if followed by the other info
		if (std::regex_match(line, m, rxVersionNr2))
		{
			gVersionNr = m[1];
			continue;
		}
	}

	if (not VERSION_STRING.empty())
		VERSION_STRING += "\n";
	VERSION_STRING += gVersionNr + " " + gVersionDate;
}

std::string get_version_nr()
{
	return gVersionNr/* + '/' + cif::get_version_nr()*/;
}

std::string get_version_date()
{
	return gVersionDate;
}

// --------------------------------------------------------------------

std::string ResidueToDSSPLine(const mmcif::DSSP::ResidueInfo& info)
{
/*
	This is the header line for the residue lines in a DSSP file:

	#  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
 */
	boost::format kDSSPResidueLine(
		"%5.5d%5.5d%1.1s%1.1s %c  %c%c%c%c%c%c%c%c%c%4.4d%4.4d%c%4.4d %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f");

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
		case mmcif::ssHelix_PPII:	ss = 'P'; break;
		case mmcif::ssTurn:			ss = 'T'; break;
		case mmcif::ssBend:			ss = 'S'; break;
		case mmcif::ssLoop:			ss = ' '; break;
	}

	char helix[4] = { ' ', ' ', ' ', ' ' };
	for (mmcif::HelixType helixType: { mmcif::HelixType::rh_3_10, mmcif::HelixType::rh_alpha, mmcif::HelixType::rh_pi, mmcif::HelixType::rh_pp })
	{
		switch (info.helix(helixType))
		{
			case mmcif::Helix::None:		helix[static_cast<int>(helixType)] = ' '; break;
			case mmcif::Helix::Start:		helix[static_cast<int>(helixType)] = '>'; break;
			case mmcif::Helix::End:			helix[static_cast<int>(helixType)] = '<'; break;
			case mmcif::Helix::StartAndEnd:	helix[static_cast<int>(helixType)] = 'X'; break;
			case mmcif::Helix::Middle:		helix[static_cast<int>(helixType)] = (helixType == mmcif::HelixType::rh_pp ? 'P' : ('3' + static_cast<int>(helixType))); break;
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
		ss % helix[3] % helix[0] % helix[1] % helix[2] % bend % chirality % bridgelabel[0] % bridgelabel[1] %
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

	os << "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA" << std::endl;
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
}

void annotateDSSP(mmcif::Structure& structure, const mmcif::DSSP& dssp, bool writeOther, std::ostream& os)
{
	auto& db = structure.getFile().data();

	if (dssp.empty())
	{
		if (cif::VERBOSE)
			std::cout << "No secondary structure information found" << std::endl;
	}
	else
	{
		// replace all struct_conf and struct_conf_type records
		auto& structConfType = db["struct_conf_type"];
		structConfType.clear();

		auto& structConf = db["struct_conf"];
		structConf.clear();

		std::map<std::string,int> foundTypes;

		auto st = dssp.begin(), lt = st;
		auto lastSS = st->ss();

		for (auto t = dssp.begin(); ; lt = t, ++t)
		{
			bool stop = t == dssp.end();

			bool flush = (stop or t->ss() != lastSS);

			if (flush and (writeOther or lastSS != mmcif::SecondaryStructureType::ssLoop))
			{
				auto& rb = st->residue();
				auto& re = lt->residue();

				std::string id;
				switch (lastSS)
				{
					case mmcif::SecondaryStructureType::ssHelix_3:
						id = "HELX_RH_3T_P";
						break;

					case mmcif::SecondaryStructureType::ssAlphahelix:
						id = "HELX_RH_AL_P";
						break;

					case mmcif::SecondaryStructureType::ssHelix_5:
						id = "HELX_RH_PI_P";
						break;

					case mmcif::SecondaryStructureType::ssHelix_PPII:
						id = "HELX_LH_PP_P";
						break;

					case mmcif::SecondaryStructureType::ssTurn:
						id = "TURN_TY1_P";
						break;

					case mmcif::SecondaryStructureType::ssBend:
						id = "BEND";
						break;

					case mmcif::SecondaryStructureType::ssBetabridge:
					case mmcif::SecondaryStructureType::ssStrand:
						id = "STRN";
						break;

					case mmcif::SecondaryStructureType::ssLoop:
						id = "OTHER";
						break;
				}

				if (foundTypes.count(id) == 0)
				{
					structConfType.emplace({
						{ "id", id },
						{ "criteria", "DSSP" }
					});
					foundTypes[id] = 1;
				}

				structConf.emplace({
					{ "conf_type_id", id },
					{ "id", id + std::to_string(foundTypes[id]++) },
					// { "pdbx_PDB_helix_id", vS(12, 14) },
					{ "beg_label_comp_id", rb.compoundID() },
					{ "beg_label_asym_id", rb.asymID() },
					{ "beg_label_seq_id", rb.seqID() },
					{ "pdbx_beg_PDB_ins_code", rb.authInsCode() },
					{ "end_label_comp_id", re.compoundID() },
					{ "end_label_asym_id", re.asymID() },
					{ "end_label_seq_id", re.seqID() },
					{ "pdbx_end_PDB_ins_code", re.authInsCode() },
		
					{ "beg_auth_comp_id", rb.compoundID() },
					{ "beg_auth_asym_id", rb.authAsymID() },
					{ "beg_auth_seq_id", rb.authSeqID() },
					{ "end_auth_comp_id", re.compoundID() },
					{ "end_auth_asym_id", re.authAsymID() },
					{ "end_auth_seq_id", re.authSeqID() }

					// { "pdbx_PDB_helix_class", vS(39, 40) },
					// { "details", vS(41, 70) },
					// { "pdbx_PDB_helix_length", vI(72, 76) }
				});

				st = t;
			}

			if (lastSS != t->ss())
			{
				st = t;
				lastSS = t->ss();
			}

			if (stop)
				break;
		}
	}

	db.add_software("dssp", "other", get_version_nr(), get_version_date());

	db.write(os);
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
		std::cout << argv[0] << " version " << VERSION_STRING << std::endl;
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
		fs::path p = vm["output"].as<std::string>();
	
		if (p.extension() == ".gz")
			p = p.stem();
		else if (p.extension() == ".bz2")
			p = p.stem();

		if (p.extension() == ".dssp")
			fmt = "dssp";
		else
			fmt = "cif";
	}


	mmcif::DSSP dssp(structure, pp_stretch, fmt == "dssp");

	if (vm.count("output"))
	{
		fs::path p = vm["output"].as<std::string>();
		std::ofstream of(p, std::ios_base::out | std::ios_base::binary);

		if (not of.is_open())
		{
			std::cerr << "Could not open output file" << std::endl;
			exit(1);
		}

		io::filtering_stream<io::output> out;
		
		if (p.extension() == ".gz")
			out.push(io::gzip_compressor());
		else if (p.extension() == ".bz2")
			out.push(io::bzip2_compressor());

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
