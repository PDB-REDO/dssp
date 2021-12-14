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

#include <boost/date_time/gregorian/formatters.hpp>
#include <boost/format.hpp>

#include <cif++/Cif2PDB.hpp>
#include <cif++/CifUtils.hpp>

#include "dssp.hpp"

using ResInfo = mmcif::DSSP::ResidueInfo;

// --------------------------------------------------------------------

bool operator<(const ResInfo &a, const ResInfo &b)
{
	return a.nr() < b.nr();
}

struct ResInfoLess
{
	bool operator()(ResInfo const &a, ResInfo const &b) const
	{
		return a.nr() < b.nr();
	}
};

// --------------------------------------------------------------------

namespace
{
std::string gVersionNr, gVersionDate;
}

void load_version_info()
{
	const std::regex
		rxVersionNr(R"(build-(\d+)-g[0-9a-f]{7}(-dirty)?)"),
		rxVersionDate(R"(Date: +(\d{4}-\d{2}-\d{2}).*)"),
		rxVersionNr2(R"(mkdssp-version: (\d+(?:\.\d+)+))");

#include "revision.hpp"

	struct membuf : public std::streambuf
	{
		membuf(char *data, size_t length) { this->setg(data, data, data + length); }
	} buffer(const_cast<char *>(kRevision), sizeof(kRevision));

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
}

std::string get_version_nr()
{
	return gVersionNr /* + '/' + cif::get_version_nr()*/;
}

std::string get_version_date()
{
	return gVersionDate;
}

std::string get_version_string()
{
	return gVersionNr + " " + gVersionDate;
}

// --------------------------------------------------------------------

std::string ResidueToDSSPLine(const ResInfo &info)
{
	/*
	    This is the header line for the residue lines in a DSSP file:

	    #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
	 */
	boost::format kDSSPResidueLine(
		"%5.5d%5.5d%1.1s%1.1s %c  %c%c%c%c%c%c%c%c%c%4.4d%4.4d%c%4.4d %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f");

	auto &residue = info.residue();

	if (residue.asymID().length() > 1)
		throw std::runtime_error("This file contains data that won't fit in the original DSSP format");

	char code = 'X';
	if (mmcif::kAAMap.find(residue.compoundID()) != mmcif::kAAMap.end())
		code = mmcif::kAAMap.at(residue.compoundID());

	if (code == 'C') // a cysteine
	{
		auto ssbridgenr = info.ssBridgeNr();
		if (ssbridgenr)
			code = 'a' + ((ssbridgenr - 1) % 26);
	}

	char ss;
	switch (info.ss())
	{
		case mmcif::ssAlphahelix: ss = 'H'; break;
		case mmcif::ssBetabridge: ss = 'B'; break;
		case mmcif::ssStrand: ss = 'E'; break;
		case mmcif::ssHelix_3: ss = 'G'; break;
		case mmcif::ssHelix_5: ss = 'I'; break;
		case mmcif::ssHelix_PPII: ss = 'P'; break;
		case mmcif::ssTurn: ss = 'T'; break;
		case mmcif::ssBend: ss = 'S'; break;
		case mmcif::ssLoop: ss = ' '; break;
	}

	char helix[4] = {' ', ' ', ' ', ' '};
	for (mmcif::HelixType helixType : {mmcif::HelixType::rh_3_10, mmcif::HelixType::rh_alpha, mmcif::HelixType::rh_pi, mmcif::HelixType::rh_pp})
	{
		switch (info.helix(helixType))
		{
			case mmcif::Helix::None: helix[static_cast<int>(helixType)] = ' '; break;
			case mmcif::Helix::Start: helix[static_cast<int>(helixType)] = '>'; break;
			case mmcif::Helix::End: helix[static_cast<int>(helixType)] = '<'; break;
			case mmcif::Helix::StartAndEnd: helix[static_cast<int>(helixType)] = 'X'; break;
			case mmcif::Helix::Middle: helix[static_cast<int>(helixType)] = (helixType == mmcif::HelixType::rh_pp ? 'P' : static_cast<char>('3' + static_cast<int>(helixType))); break;
		}
	}

	char bend = ' ';
	if (info.bend())
		bend = 'S';

	double alpha = residue.alpha();
	char chirality = alpha == 360 ? ' ' : (alpha < 0 ? '-' : '+');

	uint32_t bp[2] = {};
	char bridgelabel[2] = {' ', ' '};
	for (uint32_t i : {0, 1})
	{
		const auto &[p, ladder, parallel] = info.bridgePartner(i);
		if (not p)
			continue;

		bp[i] = p.nr() % 10000; // won't fit otherwise...
		bridgelabel[i] = (parallel ? 'a' : 'A') + ladder % 26;
	}

	char sheet = ' ';
	if (info.sheet() != 0)
		sheet = 'A' + (info.sheet() - 1) % 26;

	std::string NHO[2], ONH[2];
	for (int i : {0, 1})
	{
		const auto &[donor, donorE] = info.donor(i);
		const auto &[acceptor, acceptorE] = info.acceptor(i);

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
	auto const &[cax, cay, caz] = ca.location();

	return (kDSSPResidueLine % info.nr() % ca.authSeqID() % ca.pdbxAuthInsCode() % ca.authAsymID() % code %
			ss % helix[3] % helix[0] % helix[1] % helix[2] % bend % chirality % bridgelabel[0] % bridgelabel[1] %
			bp[0] % bp[1] % sheet % floor(info.accessibility() + 0.5) %
			NHO[0] % ONH[0] % NHO[1] % ONH[1] %
			residue.tco() % residue.kappa() % alpha % residue.phi() % residue.psi() %
			cax % cay % caz)
	    .str();
}

void writeDSSP(const mmcif::Structure &structure, const mmcif::DSSP &dssp, std::ostream &os)
{
	const std::string kFirstLine("==== Secondary Structure Definition by the program DSSP, NKI version 4.0                           ==== ");
	boost::format kHeaderLine("%1% %|127t|%2%");

	using namespace boost::gregorian;

	auto stats = dssp.GetStatistics();

	date today = day_clock::local_day();

	auto &cf = structure.getFile().file();

	os << kHeaderLine % (kFirstLine + "DATE=" + to_iso_extended_string(today)) % '.' << std::endl
	   << kHeaderLine % "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637" % '.' << std::endl
	   << GetPDBHEADERLine(cf, 127) << '.' << std::endl
	   << GetPDBCOMPNDLine(cf, 127) << '.' << std::endl
	   << GetPDBSOURCELine(cf, 127) << '.' << std::endl
	   << GetPDBAUTHORLine(cf, 127) << '.' << std::endl;

	os << boost::format("%5.5d%3.3d%3.3d%3.3d%3.3d TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN) %|127t|%c") %
			  stats.nrOfResidues % stats.nrOfChains % stats.nrOfSSBridges % stats.nrOfIntraChainSSBridges % (stats.nrOfSSBridges - stats.nrOfIntraChainSSBridges) % '.'
	   << std::endl;
	os << kHeaderLine % (boost::format("%8.1f   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)") % stats.accessibleSurface) % '.' << std::endl;

	// hydrogenbond summary

	os << kHeaderLine % (boost::format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES") % stats.nrOfHBonds % (stats.nrOfHBonds * 100.0 / stats.nrOfResidues)) % '.' << std::endl;

	os << kHeaderLine % (boost::format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES") % stats.nrOfHBondsInParallelBridges % (stats.nrOfHBondsInParallelBridges * 100.0 / stats.nrOfResidues)) % '.' << std::endl;

	os << kHeaderLine % (boost::format("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES") % stats.nrOfHBondsInAntiparallelBridges % (stats.nrOfHBondsInAntiparallelBridges * 100.0 / stats.nrOfResidues)) % '.' << std::endl;

	boost::format kHBondsLine("%5.5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I%c%1.1d), SAME NUMBER PER 100 RESIDUES");
	for (int k = 0; k < 11; ++k)
		os << kHeaderLine % (kHBondsLine % stats.nrOfHBondsPerDistance[k] % (stats.nrOfHBondsPerDistance[k] * 100.0 / stats.nrOfResidues) % (k - 5 < 0 ? '-' : '+') % abs(k - 5)) % '.' << std::endl;

	// histograms...
	os << "  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           ." << std::endl;

	for (auto hi : stats.residuesPerAlphaHelixHistogram)
		os << boost::format("%3.3d") % hi;
	os << "    RESIDUES PER ALPHA HELIX         ." << std::endl;

	for (auto hi : stats.parallelBridgesPerLadderHistogram)
		os << boost::format("%3.3d") % hi;
	os << "    PARALLEL BRIDGES PER LADDER      ." << std::endl;

	for (auto hi : stats.antiparallelBridgesPerLadderHistogram)
		os << boost::format("%3.3d") % hi;
	os << "    ANTIPARALLEL BRIDGES PER LADDER  ." << std::endl;

	for (auto hi : stats.laddersPerSheetHistogram)
		os << boost::format("%3.3d") % hi;
	os << "    LADDERS PER SHEET                ." << std::endl;

	// per residue information

	os << "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA" << std::endl;
	boost::format kDSSPResidueLine(
		"%5.5d        !%c             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0");

	int last = 0;
	for (auto ri : dssp)
	{
		// insert a break line whenever we detect missing residues
		// can be the transition to a different chain, or missing residues in the current chain

		if (ri.nr() != last + 1)
			os << (kDSSPResidueLine % (last + 1) % (ri.chainBreak() == mmcif::ChainBreak::NewChain ? '*' : ' ')) << std::endl;

		os << ResidueToDSSPLine(ri) << std::endl;
		last = ri.nr();
	}
}

// --------------------------------------------------------------------

using resinfo_set = std::set<ResInfo,ResInfoLess>;

// class Sheet
// {
//   public:
// 	Sheet(int sheetNr)
// 		: mSheetNr(sheetNr)
// 	{
// 	}

// 	static void write(cif::Datablock &db, const std::vector<Sheet> &sheets);
// 	static std::vector<Sheet> getSheets(mmcif::DSSP const &dssp);

//   private:

// 	void write(cif::Datablock &db) const;
// 	void createStrands();

// 	int mSheetNr;

// 	struct Ladder
// 	{
// 		int nr;
// 		bool parallel;
// 		resinfo_set bridges[2];

// 		resinfo_set strand(int i)
// 		{
// 			assert(i == 0 or i == 1);
// 			return bridges[i];
// 		}
// 	};

// 	std::vector<Ladder> mLadders;

// 	struct Strand
// 	{
// 		int nr;
// 		resinfo_set residues;

// 		bool overlaps(const Strand &s) const
// 		{
// assert(
// 	(find_first_of(residues.begin(), residues.end(), s.residues.begin(), s.residues.end()) != residues.end()) ==
// 	(find_first_of(s.residues.begin(), s.residues.end(), residues.begin(), residues.end()) != s.residues.end())
// );

// 			return find_first_of(residues.begin(), residues.end(), s.residues.begin(), s.residues.end()) != residues.end() or
// 				   find_first_of(s.residues.begin(), s.residues.end(), residues.begin(), residues.end()) != s.residues.end();
// 		}

// 		bool overlaps(const resinfo_set &s) const
// 		{
// assert(
// 	(find_first_of(residues.begin(), residues.end(), s.begin(), s.end()) != residues.end()) ==
// 	(find_first_of(s.begin(), s.end(), residues.begin(), residues.end()) != s.end())
// );

// 			return find_first_of(residues.begin(), residues.end(), s.begin(), s.end()) != residues.end() or
// 				   find_first_of(s.begin(), s.end(), residues.begin(), residues.end()) != s.end();
// 		}
// 	};

// 	std::vector<Strand> mStrands;
// };

// void Sheet::createStrands()
// {
// 	// for (auto &ladder : mLadders)
// 	// {
// 	// 	for (int i : { 0, 1 })
// 	// 	{
// 	// 		auto s1 = ladder.strand(i);

// 	// 		if (std::find_if(mStrands.begin(), mStrands.end(), [&s1](Strand &s) { return s.residues == s1; }) == mStrands.end())
// 	// 			mStrands.emplace_back(Strand{{ -1, -1 }, std::move(s1) });
// 	// 	}
// 	// }

// 	// auto findStrand = [this](ResInfo const &res)
// 	// {
// 	// 	auto si = std::find_if(mStrands.begin(), mStrands.end(), [&res](Strand const &s) { return s.residues.count(res); });
		
// 	// 	assert(si != mStrands.end());
		
// 	// 	return si - mStrands.begin();
// 	// };

// 	// for (auto &ladder : mLadders)
// 	// {
// 	// 	for (int i : { 0, 1 })
// 	// 	{
// 	// 		auto sIx1 = findStrand(bp1);
// 	// 		auto sIx2 = findStrand(bp2);

// 	// 		auto &s1 = mStrands[sIx1];
// 	// 		auto &s2 = mStrands[sIx2];

// 	// 		if (s1.ladder[0] != ladder.nr)
// 	// 		{
// 	// 			assert(s1.ladder[0] == -1);
// 	// 			s1.ladder[0] = ladder.nr;
// 	// 		}

// 	// 		if (s2.ladder[1] != ladder.nr)
// 	// 		{
// 	// 			assert(s2.ladder[1] == -1);
// 	// 			s2.ladder[1] = ladder.nr;
// 	// 		}
// 	// 	}
// 	// }

// }

// void Sheet::write(cif::Datablock &db) const
// {
// 	auto &struct_sheet = db["struct_sheet"];


// }

// void Sheet::write(cif::Datablock &db, const std::vector<Sheet> &sheets)
// {
// 	for (auto sheet_cat : { "struct_sheet", "struct_sheet_order", "struct_sheet_range", "pdbx_struct_sheet_hbond" })
// 	{
// 		auto &cat = db[sheet_cat];
// 		cat.clear();
// 	}

// 	for (auto sheet : sheets)
// 		sheet.write(db);
// }

// std::vector<Sheet> Sheet::getSheets(mmcif::DSSP const &dssp)
// {
// 	std::vector<Sheet> result;

// 	// std::map<int,std::pair<resinfo_set,resinfo_set>> ladders;

// 	for (auto &info : dssp)
// 	{
// 		if (info.sheet() == 0)
// 			continue;
		
// 		auto si = std::find_if(result.begin(), result.end(), [nr = info.sheet()](Sheet &s) { return s.mSheetNr == nr; });

// 		Sheet &sheet = si != result.end() ? *si : result.emplace_back(info.sheet());
// 		auto &ladders = sheet.mLadders;

// 		for (uint32_t i : {0, 1})
// 		{
// 			const auto &[p, ladder, parallel] = info.bridgePartner(i);
// 			if (not p)
// 				continue;

// 			auto li = std::find_if(ladders.begin(), ladders.end(), [nr=ladder](auto &l) { return l.nr == nr; });
// 			if (li == ladders.end())
// 				li = ladders.insert(ladders.end(), { ladder, parallel });

// 			if (li->bridges[0].count(info) or li->bridges[1].count(info))
// 				continue;

// 			li->bridges[0].insert(info);
// 			li->bridges[1].insert(p);
// 		}
// 	}

// 	int strand_nr = 0;
// 	for (auto &sheet : result)
// 	{
// 		for (auto &ladder : sheet.mLadders)
// 		{
// 			for (int i : { 0, 1 })
// 			{
// 				auto si = std::find_if(sheet.mStrands.begin(), sheet.mStrands.end(), [&ladder,i](Strand &str){ return str.overlaps(ladder.bridges[i]); });
// 				if (si != sheet.mStrands.end())
// 					continue;
				
// 				sheet.mStrands.emplace_back(Strand{strand_nr++, ladder.bridges[i]});
// 			}
// 		}

// 		assert(sheet.mLadders.size() + 1 == sheet.mStrands.size());
// 	}


// 	for (auto &sheet : result)
// 		sheet.createStrands();

// 	return result;
// }

// --------------------------------------------------------------------

void writeHBonds(cif::Datablock &db, const mmcif::DSSP &dssp)
{
	using ResidueInfo = mmcif::DSSP::ResidueInfo;

	auto &hb = db["dssp_hbond"];

	for (auto &info : dssp)
	{
		auto write_res = [&](const std::string &prefix, ResidueInfo const &info, double energy)
		{
			auto &res = info.residue();

			hb.back().assign({
				{ prefix + "label_comp_id", res.compoundID() },
				{ prefix + "label_seq_id", res.seqID() },
				{ prefix + "label_asym_id", res.asymID() },
				// { prefix + "auth_comp_id", res.compoundID() },
				{ prefix + "auth_seq_id", res.authSeqID() },
				{ prefix + "auth_asym_id", res.authAsymID() },
				{ prefix + "pdbx_PDB_ins_code", res.authInsCode() }
			});

			if (not prefix.empty())
				hb.back().assign({
					{ prefix + "energy", energy, "%.1f" }
				});
		};

		hb.emplace({
			{ "id", hb.getUniqueID("") }
		});

		write_res("", info, 0);

		for (int i : {0, 1})
		{
			const auto &&[ donor, donorEnergy ] = info.donor(i);

			if (donor)
				write_res(i ? "donor_2_" : "donor_1_", donor, donorEnergy);

			const auto &&[ acceptor, acceptorEnergy ] = info.acceptor(i);

			if (acceptor)
				write_res(i ? "acceptor_2_" : "acceptor_1_", acceptor, acceptorEnergy);
		}
	}
}

using res_list = std::vector<ResInfo>;
using strand_list = std::vector<std::tuple<int,res_list>>;

int strandNrForResidue(const strand_list &strands, ResInfo const &res)
{
	int result = 0;
	for (auto &&[sheet, strand] : strands)
	{
		if (std::find(strand.begin(), strand.end(), res) != strand.end())
			break;
		++result;
	}
	assert(result < strands.size());
	return result + 1;
};

strand_list writeSheets(cif::Datablock &db, const mmcif::DSSP &dssp)
{
	using ss_type = mmcif::SecondaryStructureType;

	for (auto sheet_cat : { "struct_sheet", "struct_sheet_order", "struct_sheet_range", "pdbx_struct_sheet_hbond" })
	{
		auto &cat = db[sheet_cat];
		cat.clear();
	}

	std::vector<std::tuple<int,res_list>> strands;

	ss_type ss = ss_type::ssLoop;

	for (auto &info : dssp)
	{
		ss_type iss = info.ss();

		if (iss == ss)
		{
			if (ss == ss_type::ssStrand)
				std::get<1>(strands.back()).emplace_back(info);
			continue;
		}

		ss = iss;

		if (ss != ss_type::ssStrand)
			continue;

		int sheetNr = info.sheet();
		assert(sheetNr);

		strands.emplace_back(std::make_tuple(sheetNr, res_list{ info }));
	}

	auto &struct_sheet = db["struct_sheet"];

	int lastSheet = -1;
	for (const auto &[ sheetNr, strand ] : strands)
	{
		if (sheetNr != lastSheet)
		{
			struct_sheet.emplace({
				{ "id", cif::cifIdForNumber(sheetNr - 1) },
				{ "number_strands", std::count_if(strands.begin(), strands.end(), [nr=sheetNr](std::tuple<int,res_list> const &s) { return std::get<0>(s) == nr; }) }
			});

			lastSheet = sheetNr;
		}
	}

	for (auto &&[sheet, strand] : strands)
		std::sort(strand.begin(), strand.end(), ResInfoLess());

	std::set<std::tuple<int,int>> ladderSeen;
	auto &struct_sheet_order = db["struct_sheet_order"];

	for (auto &info : dssp)
	{
		if (info.ss() != ss_type::ssStrand)
			continue;

		int s1 = strandNrForResidue(strands, info);

		for (int i : { 0, 1 })
		{
			const auto &[p, ladder, parallel] = info.bridgePartner(i);
			if (not p or p.sheet() != info.sheet())
				continue;

			int s2 = strandNrForResidue(strands, p);
			assert(s1 != s2);

			auto k = s1 > s2 ? std::make_tuple(s2, s1) : std::make_tuple(s1, s2);
			if (ladderSeen.count(k))
				continue;
			
			struct_sheet_order.emplace({
				{ "sheet_id", cif::cifIdForNumber(info.sheet() - 1) },
				{ "range_id_1", s1 },
				{ "range_id_2", s2 },
				{ "sense", parallel ? "parallel" : "anti-parallel" }
			});

			ladderSeen.insert(k);
		}
	}

	auto &struct_sheet_range = db["struct_sheet_range"];

	int strandID = 0;
	for (auto &&[sheet, strand] : strands)
	{
		auto &beg = strand.front().residue();
		auto &end = strand.back().residue();

		struct_sheet_range.emplace({
			{ "sheet_id", cif::cifIdForNumber(sheet - 1) },
			{ "id", ++strandID },
			{ "beg_label_comp_id", beg.compoundID() },
			{ "beg_label_asym_id", beg.asymID() },
			{ "beg_label_seq_id", beg.seqID() },
			{ "pdbx_beg_PDB_ins_code", beg.authInsCode() },
			{ "end_label_comp_id", end.compoundID() },
			{ "end_label_asym_id", end.asymID() },
			{ "end_label_seq_id", end.seqID() },
			{ "pdbx_end_PDB_ins_code", end.authInsCode() },
			{ "beg_auth_comp_id", beg.compoundID() },
			{ "beg_auth_asym_id", beg.authAsymID() },
			{ "beg_auth_seq_id", beg.authSeqID() },
			{ "end_auth_comp_id", end.compoundID() },
			{ "end_auth_asym_id", end.authAsymID() },
			{ "end_auth_seq_id", end.authSeqID() }
		});
	}

	return strands;
}

void annotateDSSP(mmcif::Structure &structure, const mmcif::DSSP &dssp, bool writeOther, std::ostream &os)
{
	auto &db = structure.getFile().data();

	if (dssp.empty())
	{
		if (cif::VERBOSE)
			std::cout << "No secondary structure information found" << std::endl;
	}
	else
	{
		// start by writing out all hbonds
		writeHBonds(db, dssp);

		// replace all struct_conf and struct_conf_type records
		auto &structConfType = db["struct_conf_type"];
		structConfType.clear();

		auto &structConf = db["struct_conf"];
		structConf.clear();

		std::map<std::string, int> foundTypes;

		auto st = dssp.begin(), lt = st;
		auto lastSS = st->ss();

		for (auto t = dssp.begin();; lt = t, ++t)
		{
			bool stop = t == dssp.end();

			bool flush = (stop or t->ss() != lastSS);

			if (flush and (writeOther or lastSS != mmcif::SecondaryStructureType::ssLoop))
			{
				auto &rb = st->residue();
				auto &re = lt->residue();

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
					structConfType.emplace({{"id", id},
						{"criteria", "DSSP"}});
					foundTypes[id] = 1;
				}

				structConf.emplace({
					{"conf_type_id", id},
					{"id", id + std::to_string(foundTypes[id]++)},
					// { "pdbx_PDB_helix_id", vS(12, 14) },
					{"beg_label_comp_id", rb.compoundID()},
					{"beg_label_asym_id", rb.asymID()},
					{"beg_label_seq_id", rb.seqID()},
					{"end_label_comp_id", re.compoundID()},
					{"end_label_asym_id", re.asymID()},
					{"end_label_seq_id", re.seqID()},

					{"beg_auth_comp_id", rb.compoundID()},
					{"beg_auth_asym_id", rb.authAsymID()},
					{"beg_auth_seq_id", rb.authSeqID()},
					{"pdbx_beg_PDB_ins_code", rb.authInsCode()},
					{"end_auth_comp_id", re.compoundID()},
					{"end_auth_asym_id", re.authAsymID()},
					{"end_auth_seq_id", re.authSeqID()},
					{"pdbx_end_PDB_ins_code", re.authInsCode()},

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

		// write out our concept of sheets
		auto strands = writeSheets(db, dssp);

		// extra info
		auto &dssp_info = db["dssp_res_info"];

		for (auto info : dssp)
		{
			auto &res = info.residue();

			char ss;
			switch (info.ss())
			{
				case mmcif::ssAlphahelix: ss = 'H'; break;
				case mmcif::ssBetabridge: ss = 'B'; break;
				case mmcif::ssStrand: ss = 'E'; break;
				case mmcif::ssHelix_3: ss = 'G'; break;
				case mmcif::ssHelix_5: ss = 'I'; break;
				case mmcif::ssHelix_PPII: ss = 'P'; break;
				case mmcif::ssTurn: ss = 'T'; break;
				case mmcif::ssBend: ss = 'S'; break;
				case mmcif::ssLoop: ss = '.'; break;
			}

			char helix[4] = {'.', '.', '.', '.'};
			for (mmcif::HelixType helixType : {mmcif::HelixType::rh_3_10, mmcif::HelixType::rh_alpha, mmcif::HelixType::rh_pi, mmcif::HelixType::rh_pp})
			{
				switch (info.helix(helixType))
				{
					case mmcif::Helix::None: helix[static_cast<int>(helixType)] = '.'; break;
					case mmcif::Helix::Start: helix[static_cast<int>(helixType)] = '>'; break;
					case mmcif::Helix::End: helix[static_cast<int>(helixType)] = '<'; break;
					case mmcif::Helix::StartAndEnd: helix[static_cast<int>(helixType)] = 'X'; break;
					case mmcif::Helix::Middle: helix[static_cast<int>(helixType)] = (helixType == mmcif::HelixType::rh_pp ? 'P' : static_cast<char>('3' + static_cast<int>(helixType))); break;
				}
			}

			double alpha = res.alpha();

			const auto &[ca_x, ca_y, ca_z] = res.atomByID("CA").location();

			std::string ladderID[2];
			if (info.sheet())
			{
				for (uint32_t i : {0, 1})
				{
					const auto &[p, ladder, parallel] = info.bridgePartner(i);
					if (not p)
						continue;

					ladderID[i] = cif::cifIdForNumber(ladder);
				}
			}

			dssp_info.emplace({
				{"id", info.nr()},
				{"label_comp_id", res.compoundID()},
				{"label_seq_id", res.seqID()},
				{"label_asym_id", res.asymID()},
				{"auth_comp_id", res.compoundID()},
				{"auth_seq_id", res.authSeqID()},
				{"auth_asym_id", res.authAsymID()},
				{"pdbx_PDB_ins_code", res.authInsCode()},
				{"cysteine_bridge_nr", info.ssBridgeNr() ? std::to_string(info.ssBridgeNr()) : "."},
				{"secondary_structure_code", ss},
				{"helix_3_10", helix[0]},
				{"helix_alpha", helix[1]},
				{"helix_pi", helix[2]},
				{"helix_pp", helix[3]},
				{"bend", info.bend() ? 'S' : '.'},
				{"chirality", alpha == 360 ? '.' : (alpha < 0 ? '-' : '+')},
				{"sheet_id", info.sheet() ? cif::cifIdForNumber(info.sheet() - 1) : ""},
				{"sheet_range_id", info.ss() == mmcif::SecondaryStructureType::ssStrand ? std::to_string(strandNrForResidue(strands, info)) : "" },
				{"tco", res.tco(), "%.3f"},
				{"kappa", res.kappa(), "%.1f"},
				{"alpha", res.alpha(), "%.1f"},
				{"phi", res.phi(), "%.1f"},
				{"psi", res.psi(), "%.1f"},
				{"ca_cartn_x", ca_x, "%.3f"},
				{"ca_cartn_y", ca_y, "%.3f"},
				{"ca_cartn_z", ca_z, "%.3f"}
			});
		}
	}

	db.add_software("dssp", "model annotation", get_version_nr(), get_version_date());

	db.write(os);
}
