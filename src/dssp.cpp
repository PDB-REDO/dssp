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
#include <boost/algorithm/string.hpp>

#include <cif++/Cif2PDB.hpp>
#include <cif++/CifUtils.hpp>

#include "dssp.hpp"

namespace ba = boost::algorithm;

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

std::map<std::tuple<std::string,int>,int> writeSheets(cif::Datablock &db, const mmcif::DSSP &dssp)
{
	using res_list = std::vector<ResInfo>;
	using ss_type = mmcif::SecondaryStructureType;

	// clean up old info first

	for (auto sheet_cat : { "struct_sheet", "struct_sheet_order", "struct_sheet_range", "struct_sheet_hbond", "pdbx_struct_sheet_hbond" })
	{
		auto &cat = db[sheet_cat];
		cat.clear();
	}

	// create a list of strands, based on the SS info in DSSP. Store sheet number along with the strand.

	std::vector<std::tuple<int,res_list>> strands;
	std::map<std::tuple<std::string,int>,int> sheetMap;		// mapping from bridge number (=info.sheet()) in DSSP to sheet ID

	ss_type ss = ss_type::ssLoop;

	for (auto &info : dssp)
	{
		std::tuple<std::string,int> sheetID{ info.residue().asymID(), info.sheet() };
		ss_type iss = info.ss();

		if (iss == ss and ss == ss_type::ssStrand and sheetMap[sheetID] == std::get<0>(strands.back()))
		{
			std::get<1>(strands.back()).emplace_back(info);
			continue;
		}

		ss = iss;

		if (ss != ss_type::ssStrand)
			continue;

		if (not sheetMap.count(sheetID))
			sheetMap[sheetID] = sheetMap.size();

		strands.emplace_back(std::make_tuple(sheetMap[sheetID], res_list{ info }));
	}

	// sort the strands vector
	std::sort(strands.begin(), strands.end(), [](auto &a, auto &b)
	{
		const auto &[sheetA, strandsA] = a;
		const auto &[sheetB, strandsB] = b;

		int d = sheetA - sheetB;
		if (d == 0)
			d = strandsA.front().nr() - strandsB.front().nr();
		return d < 0;
	});

	// write out the struct_sheet, since all info is available now
	auto &struct_sheet = db["struct_sheet"];

	int lastSheet = -1;
	for (const auto &[ sheetNr, strand ] : strands)
	{
		if (sheetNr != lastSheet)
		{
			struct_sheet.emplace({
				{ "id", cif::cifIdForNumber(sheetNr) },
				{ "number_strands", std::count_if(strands.begin(), strands.end(), [nr=sheetNr](std::tuple<int,res_list> const &s) { return std::get<0>(s) == nr; }) }
			});

			lastSheet = sheetNr;
		}
	}

	// Each residue resides in a single strand which is part of a single sheet
	// this function returns the sequence number inside the sheet for the strand
	// containing res
	auto strandNrForResidue = [&strands,&sheetMap](ResInfo const &res)
	{
		for (const auto &[k, iSheet] : sheetMap)
		{
			int result = 0;
			for (auto &&[sheet, strand] : strands)
			{
				if (sheet != iSheet)
					continue;

				if (std::find(strand.begin(), strand.end(), res) != strand.end())
					return result;

				++result;
			}
		}

		assert(false);
		return -1;
	};

	// This map is used to record the sense of a ladder, and can be used
	// to detect ladders already seen.
	std::map<std::tuple<int,int,int>, std::tuple<int,bool>> ladderSense;

	for (auto &info : dssp)
	{
		if (info.ss() != ss_type::ssStrand)
			continue;

		int s1 = strandNrForResidue(info);

		for (int i : { 0, 1 })
		{
			const auto &[p, ladder, parallel] = info.bridgePartner(i);
			if (not p or p.residue().asymID() != info.residue().asymID() or p.sheet() != info.sheet() or p.ss() != ss_type::ssStrand)
				continue;

			int s2 = strandNrForResidue(p);
			// assert(s1 != s2);
			if (s2 == s1)
				continue;

			std::tuple<std::string,int> sheetID{info.residue().asymID(), info.sheet() };
			int sheet = sheetMap[sheetID];

			auto k = s1 > s2 ? std::make_tuple(sheet, s2, s1) : std::make_tuple(sheet, s1, s2);
			if (ladderSense.count(k))
			{
				// assert(ladderSense[k] == std::make_tuple(ladder, parallel));
				assert(std::get<1>(ladderSense[k]) == parallel);
				continue;
			}

			ladderSense.emplace(k, std::make_tuple(ladder, parallel));
		}
	}

	auto &struct_sheet_order = db["struct_sheet_order"];
	auto &struct_sheet_hbond = db["struct_sheet_hbond"];

	for (const auto&[key, value] : ladderSense)
	{
		const auto &[sheet, s1, s2] = key;
		const auto &[ladder, parallel] = value;

		struct_sheet_order.emplace({
			{ "sheet_id", cif::cifIdForNumber(sheet) },
			// { "dssp_ladder_id", cif::cifIdForNumber(ladder) },
			{ "range_id_1", s1 + 1 },
			{ "range_id_2", s2 + 1 },
			{ "sense", parallel ? "parallel" : "anti-parallel" }
		});

		res_list strand1, strand2;

		int strandIx = 0;
		for (auto const &s : strands)
		{
			const auto &[sSheet, strand] = s;
			if (sSheet != sheet)
				continue;
			
			if (strandIx == s1)
				strand1 = strand;
			else if (strandIx == s2)
			{
				strand2 = strand;
				break;
			}

			++strandIx;
		}

		assert(not (strand1.empty() or strand2.empty()));

		int beg1SeqID = 0, beg2SeqID = 0, end1SeqID = 0, end2SeqID = 0;
		std::string beg1AtomID, beg2AtomID, end1AtomID, end2AtomID;

		if (parallel)
		{
			// I.	a	d	II.	a	d		parallel
			//		  \			  /
			//		b	e		b	e								<= the residues forming the bridge
			// 		  /			  \                      ..
			//		c	f		c	f

			for (auto b : strand1)
			{
				for (int i : { 0, 1 })
				{
					const auto &[e, ladder1, parallel1] = b.bridgePartner(i);
					auto esi = std::find(strand2.begin(), strand2.end(), e);

					if (esi == strand2.end())
						continue;
					
					auto bi = std::find(dssp.begin(), dssp.end(), b);
					assert(bi != dssp.end() and bi != dssp.begin());

					auto a = *std::prev(bi);
					auto c = *std::next(bi);

					auto ei = std::find(dssp.begin(), dssp.end(), e);
					assert(ei != dssp.end() and ei != dssp.begin());
					auto d = *std::prev(ei);
					auto f = *std::next(ei);

					if (TestBond(e, a) and TestBond(c, e))					// case I.
					{
						beg1SeqID = a.residue().seqID();
						beg2SeqID = e.residue().seqID();

						beg1AtomID = "O";
						beg2AtomID = "N";
					}
					else if (TestBond(b, d) and TestBond(f, b))				// case II.
					{
						beg1SeqID = b.residue().seqID();
						beg2SeqID = d.residue().seqID();

						beg1AtomID = "N";
						beg2AtomID = "O";
					}

					break;
				}

				if (beg1SeqID)
					break;
			}

			std::reverse(strand1.begin(), strand1.end());
			std::reverse(strand2.begin(), strand2.end());

			for (auto b : strand1)
			{
				for (int i : { 0, 1 })
				{
					const auto &[e, ladder1, parallel1] = b.bridgePartner(i);
					auto esi = std::find(strand2.begin(), strand2.end(), e);

					if (esi == strand2.end())
						continue;
					
					auto bi = std::find(dssp.begin(), dssp.end(), b);
					assert(bi != dssp.end() and bi != dssp.begin());

					auto a = *std::next(bi);
					auto c = *std::prev(bi);

					auto ei = std::find(dssp.begin(), dssp.end(), e);
					assert(ei != dssp.end() and ei != dssp.begin());
					auto d = *std::next(ei);
					auto f = *std::prev(ei);

					if (TestBond(a, e) and TestBond(e, c))					// case I.
					{
						end1SeqID = a.residue().seqID();
						end2SeqID = e.residue().seqID();

						end1AtomID = "N";
						end2AtomID = "O";
					}
					else if (TestBond(d, b) and TestBond(b, f))				// case II.
					{
						end1SeqID = b.residue().seqID();
						end2SeqID = d.residue().seqID();

						end1AtomID = "O";
						end2AtomID = "N";
					}

					break;
				}

				if (end1SeqID)
					break;
			}
		}
		else
		{
			// III.	a <- f	IV. a	  f		antiparallel
			//
			//		b	 e      b <-> e							<= the residues forming the bridge
			//
			//		c -> d		c     d

			std::reverse(strand2.begin(), strand2.end());

			for (auto b : strand1)
			{
				for (int i : { 0, 1 })
				{
					const auto &[e, ladder1, parallel1] = b.bridgePartner(i);
					auto esi = std::find(strand2.begin(), strand2.end(), e);

					if (esi == strand2.end())
						continue;
					
					auto bi = std::find(dssp.begin(), dssp.end(), b);
					assert(bi != dssp.end() and bi != dssp.begin());

					auto a = *std::prev(bi);
					auto c = *std::next(bi);

					auto ei = std::find(dssp.begin(), dssp.end(), e);
					assert(ei != dssp.end() and ei != dssp.begin());
					auto d = *std::prev(ei);
					auto f = *std::next(ei);

					if (TestBond(f, a) and TestBond(c, d))				// case III.
					{
						beg1SeqID = a.residue().seqID();
						beg2SeqID = f.residue().seqID();

						beg1AtomID = "O";
						beg2AtomID = "N";
					}
					else if (TestBond(b, e) and TestBond(e, b))			// case IV.
					{
						beg1SeqID = b.residue().seqID();
						beg2SeqID = e.residue().seqID();

						beg1AtomID = "N";
						beg2AtomID = "O";
					}

					break;
				}

				if (beg1SeqID)
					break;
			}

			std::reverse(strand1.begin(), strand1.end());
			std::reverse(strand2.begin(), strand2.end());

			for (auto b : strand1)
			{
				for (int i : { 0, 1 })
				{
					const auto &[e, ladder1, parallel1] = b.bridgePartner(i);
					auto esi = std::find(strand2.begin(), strand2.end(), e);

					if (esi == strand2.end())
						continue;
					
					auto bi = std::find(dssp.begin(), dssp.end(), b);
					assert(bi != dssp.end() and bi != dssp.begin());

					auto a = *std::next(bi);
					auto c = *std::prev(bi);

					auto ei = std::find(dssp.begin(), dssp.end(), e);
					assert(ei != dssp.end() and ei != dssp.begin());
					auto d = *std::next(ei);
					auto f = *std::prev(ei);

					if (TestBond(a, f) and TestBond(d, c))				// case III.
					{
						end1SeqID = a.residue().seqID();
						end2SeqID = f.residue().seqID();

						end1AtomID = "N";
						end2AtomID = "O";
					}
					else if (TestBond(b, e) and TestBond(e, b))			// case IV.
					{
						end1SeqID = b.residue().seqID();
						end2SeqID = e.residue().seqID();

						end1AtomID = "N";
						end2AtomID = "O";
					}

					break;
				}

				if (end1SeqID)
					break;
			}
		}

		struct_sheet_hbond.emplace({
			{ "sheet_id", cif::cifIdForNumber(sheet) },
			{ "range_id_1", s1 + 1 },
			{ "range_id_2", s2 + 1 },
			{ "range_1_beg_label_seq_id", beg1SeqID },	
			{ "range_1_beg_label_atom_id", beg1AtomID },	
			{ "range_2_beg_label_seq_id", beg2SeqID },	
			{ "range_2_beg_label_atom_id", beg2AtomID },	
			{ "range_1_end_label_seq_id", end1SeqID },	
			{ "range_1_end_label_atom_id", end1AtomID },	
			{ "range_2_end_label_seq_id", end2SeqID },	
			{ "range_2_end_label_atom_id", end2AtomID }	
		});


		// if (parallel)
		// {
		// 	if (TestBond(c, e) and TestBond(e, a))
		// 		row.emplace({
		// 			{ "sheet_id", cif::cifIdForNumber(sheetMap[info.sheet()]) },
		// 			{ "range_id_1", s1 + 1 },
		// 			{ "range_id_2", s2 + 1 }
		// 			{ "range_1_beg_label_seq_id", "" },
		// 			{ "range_1_beg_label_atom_id", "" },
		// 			{ "range_2_beg_label_seq_id", "" },
		// 			{ "range_2_beg_label_atom_id", "" },
		// 			{ "range_1_end_label_seq_id", "" },
		// 			{ "range_1_end_label_atom_id", "" },
		// 			{ "range_2_end_label_seq_id", "" },
		// 			{ "range_2_end_label_atom_id", "" }
		// 		});
		// }



	}

	auto &struct_sheet_range = db["struct_sheet_range"];

	for (const auto &[key, iSheet] : sheetMap)
	{
		for (auto &&[sheet, strand] : strands)
		{
			if (sheet != iSheet)
				continue;

			std::sort(strand.begin(), strand.end(), ResInfoLess());

			auto &beg = strand.front().residue();
			auto &end = strand.back().residue();

			struct_sheet_range.emplace({
				{ "sheet_id", cif::cifIdForNumber(sheet) },
				{ "id", strandNrForResidue(strand.front()) + 1 },
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
	}


	return sheetMap;
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
		auto sheetMap = writeSheets(db, dssp);

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
			std::string bridge[2];

			if (info.sheet())
			{
				for (uint32_t i : {0, 1})
				{
					const auto &[p, ladder, parallel] = info.bridgePartner(i);
					if (not p)
						continue;

					ladderID[i] = cif::cifIdForNumber(ladder);

					if (parallel)
						ba::to_lower(ladderID[i]);

					bridge[i] = std::to_string(p.nr());
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
				{"sheet_id", info.sheet() ? cif::cifIdForNumber(sheetMap[{info.residue().asymID(), info.sheet()}]) : ""},
				{"dssp_ladder_id_1", ladderID[0]},
				{"dssp_ladder_id_2", ladderID[1]},
				{"bridge_partner_id_1", bridge[0]},
				{"bridge_partner_id_2", bridge[1]},
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
