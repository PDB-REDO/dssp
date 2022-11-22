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

#include <cif++/pdb/io.hpp>

#include "dssp-io.hpp"
#include "revision.hpp"

// --------------------------------------------------------------------

std::string ResidueToDSSPLine(const dssp::residue_info &info)
{
	/*
	    This is the header line for the residue lines in a DSSP file:

	    #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
	 */

	// auto& residue = info.residue();
	auto &residue = info;

	if (residue.asym_id().length() > 1)
		throw std::runtime_error("This file contains data that won't fit in the original DSSP format");

	char code = residue.compound_letter();

	if (code == 'C') // a cysteine
	{
		auto ssbridgenr = info.ssBridgeNr();
		if (ssbridgenr)
			code = 'a' + ((ssbridgenr - 1) % 26);
	}

	char ss = ' ';
	switch (info.type())
	{
		case dssp::structure_type::Alphahelix: ss = 'H'; break;
		case dssp::structure_type::Betabridge: ss = 'B'; break;
		case dssp::structure_type::Strand: ss = 'E'; break;
		case dssp::structure_type::Helix_3: ss = 'G'; break;
		case dssp::structure_type::Helix_5: ss = 'I'; break;
		case dssp::structure_type::Helix_PPII: ss = 'P'; break;
		case dssp::structure_type::Turn: ss = 'T'; break;
		case dssp::structure_type::Bend: ss = 'S'; break;
		case dssp::structure_type::Loop: ss = ' '; break;
	}

	char helix[4] = {' ', ' ', ' ', ' '};
	for (dssp::helix_type helixType : {dssp::helix_type::_3_10, dssp::helix_type::alpha, dssp::helix_type::pi, dssp::helix_type::pp})
	{
		switch (info.helix(helixType))
		{
			case dssp::helix_position_type::None: helix[static_cast<int>(helixType)] = ' '; break;
			case dssp::helix_position_type::Start: helix[static_cast<int>(helixType)] = '>'; break;
			case dssp::helix_position_type::End: helix[static_cast<int>(helixType)] = '<'; break;
			case dssp::helix_position_type::StartAndEnd: helix[static_cast<int>(helixType)] = 'X'; break;
			case dssp::helix_position_type::Middle: helix[static_cast<int>(helixType)] = (helixType == dssp::helix_type::pp ? 'P' : static_cast<char>('3' + static_cast<int>(helixType))); break;
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
		const auto &[p, ladder, parallel] = info.bridge_partner(i);
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
			NHO[i] = cif::format("%d,%3.1f", d, acceptorE).str();
		}

		if (donor)
		{
			auto d = donor.nr() - info.nr();
			ONH[i] = cif::format("%d,%3.1f", d, donorE).str();
		}
	}

	// auto ca = residue.atomByID("CA");
	auto const &[cax, cay, caz] = residue.ca_location();

	return cif::format("%5d%5d%1.1s%1.1s %c  %c%c%c%c%c%c%c%c%c%4d%4d%c%4.0f %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f",
		info.nr(), residue.pdb_seq_num(), residue.pdb_ins_code(), residue.pdb_strand_id(), code,
			ss, helix[3], helix[0], helix[1], helix[2], bend, chirality, bridgelabel[0], bridgelabel[1],
			bp[0], bp[1], sheet, floor(info.accessibility() + 0.5),
			NHO[0], ONH[0], NHO[1], ONH[1],
			residue.tco(), residue.kappa(), alpha, residue.phi(), residue.psi(),
			cax, cay, caz)
	    .str();
}

void writeDSSP(const dssp &dssp, std::ostream &os)
{
	using namespace std::chrono;

	auto stats = dssp.get_statistics();

	std::time_t today = system_clock::to_time_t(system_clock::now());
	std::tm *tm = std::gmtime(&today);

	os << "==== Secondary Structure Definition by the program DSSP, NKI version 4.0                           ==== DATE=" << std::put_time(tm, "%F") << "        ." << std::endl
	   << "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637                                                              ." << std::endl
	   << dssp.get_pdb_header_line(dssp::pdb_record_type::HEADER) << '.' << std::endl
	   << dssp.get_pdb_header_line(dssp::pdb_record_type::COMPND) << '.' << std::endl
	   << dssp.get_pdb_header_line(dssp::pdb_record_type::SOURCE) << '.' << std::endl
	   << dssp.get_pdb_header_line(dssp::pdb_record_type::AUTHOR) << '.' << std::endl;

	os << cif::format("%5d%3d%3d%3d%3d TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .",
			  stats.count.residues, stats.count.chains, stats.count.SS_bridges, stats.count.intra_chain_SS_bridges, (stats.count.SS_bridges - stats.count.intra_chain_SS_bridges))
	   << std::endl;

	os << cif::format("%8.1f   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .", stats.accessible_surface) << std::endl;

	// hydrogenbond summary

	os << cif::format("%5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .", stats.count.H_bonds, (stats.count.H_bonds * 100.0 / stats.count.residues)) << std::endl;

	os << cif::format("%5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .", stats.count.H_bonds_in_parallel_bridges, (stats.count.H_bonds_in_parallel_bridges * 100.0 / stats.count.residues)) << std::endl;

	os << cif::format("%5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .", stats.count.H_bonds_in_antiparallel_bridges, (stats.count.H_bonds_in_antiparallel_bridges * 100.0 / stats.count.residues)) << std::endl;

	for (int k = 0; k < 11; ++k)
		os << cif::format("%5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I%c%1d), SAME NUMBER PER 100 RESIDUES                              .", stats.count.H_Bonds_per_distance[k], (stats.count.H_Bonds_per_distance[k] * 100.0 / stats.count.residues), (k - 5 < 0 ? '-' : '+'), abs(k - 5)) << std::endl;

	// histograms...
	os << "  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           ." << std::endl;

	for (auto hi : stats.histogram.residues_per_alpha_helix)
		os << cif::format("%3d", hi);
	os << "    RESIDUES PER ALPHA HELIX         ." << std::endl;

	for (auto hi : stats.histogram.parallel_bridges_per_ladder)
		os << cif::format("%3d", hi);
	os << "    PARALLEL BRIDGES PER LADDER      ." << std::endl;

	for (auto hi : stats.histogram.antiparallel_bridges_per_ladder)
		os << cif::format("%3d", hi);
	os << "    ANTIPARALLEL BRIDGES PER LADDER  ." << std::endl;

	for (auto hi : stats.histogram.ladders_per_sheet)
		os << cif::format("%3d", hi);
	os << "    LADDERS PER SHEET                ." << std::endl;

	// per residue information

	os << "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA" << std::endl;

	int last = 0;
	for (auto ri : dssp)
	{
		// insert a break line whenever we detect missing residues
		// can be the transition to a different chain, or missing residues in the current chain

		if (ri.nr() != last + 1)
			os << cif::format("%5d        !%c             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0",
				(last + 1), (ri.chain_break() == dssp::chain_break_type::NewChain ? '*' : ' ')) << std::endl;

		os << ResidueToDSSPLine(ri) << std::endl;
		last = ri.nr();
	}
}

// --------------------------------------------------------------------

void writeHBonds(cif::datablock &db, const dssp &dssp)
{
	using ResidueInfo = dssp::residue_info;

	auto &hb = db["dssp_hbond"];

	for (auto &res : dssp)
	{
		auto write_res = [&](const std::string &prefix, ResidueInfo const &info, double energy)
		{
			hb.back().assign({
				{ prefix + "label_comp_id", res.compound_id() },
				{ prefix + "label_seq_id", res.seq_id() },
				{ prefix + "label_asym_id", res.asym_id() },
				// { prefix + "auth_comp_id", res.compound_id() },
				{ prefix + "auth_seq_id", res.auth_seq_id() },
				{ prefix + "auth_asym_id", res.auth_asym_id() },
				{ prefix + "pdbx_PDB_ins_code", res.pdb_ins_code() }
			});

			if (not prefix.empty())
				hb.back().assign({
					{ prefix + "energy", energy, 1 }
				});
		};

		hb.emplace({
			{ "id", hb.get_unique_id("") }
		});

		write_res("", res, 0);

		for (int i : {0, 1})
		{
			const auto &&[ donor, donorEnergy ] = res.donor(i);

			if (donor)
				write_res(i ? "donor_2_" : "donor_1_", donor, donorEnergy);

			const auto &&[ acceptor, acceptorEnergy ] = res.acceptor(i);

			if (acceptor)
				write_res(i ? "acceptor_2_" : "acceptor_1_", acceptor, acceptorEnergy);
		}
	}
}

std::map<std::tuple<std::string,int>,int> writeSheets(cif::datablock &db, const dssp &dssp)
{
	using res_list = std::vector<dssp::residue_info>;
	using ss_type = dssp::structure_type;

	// clean up old info first

	for (auto sheet_cat : { "struct_sheet", "struct_sheet_order", "struct_sheet_range", "struct_sheet_hbond", "pdbx_struct_sheet_hbond" })
	{
		auto &cat = db[sheet_cat];
		cat.clear();
	}

	// create a list of strands, based on the SS info in DSSP. Store sheet number along with the strand.

	std::vector<std::tuple<int,res_list>> strands;
	std::map<std::tuple<std::string,int>,int> sheetMap;		// mapping from bridge number (=info.sheet()) in DSSP to sheet ID

	ss_type ss = ss_type::Loop;

	for (auto &res : dssp)
	{
		std::tuple<std::string,int> sheetID{ res.asym_id(), res.sheet() };
		ss_type iss = res.type();

		if (iss == ss and ss == ss_type::Strand and sheetMap[sheetID] == std::get<0>(strands.back()))
		{
			std::get<1>(strands.back()).emplace_back(res);
			continue;
		}

		ss = iss;

		if (ss != ss_type::Strand)
			continue;

		if (not sheetMap.count(sheetID))
			sheetMap[sheetID] = sheetMap.size();

		strands.emplace_back(std::make_tuple(sheetMap[sheetID], res_list{ res }));
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
				{ "id", cif::cif_id_for_number(sheetNr) },
				{ "number_strands", std::count_if(strands.begin(), strands.end(), [nr=sheetNr](std::tuple<int,res_list> const &s) { return std::get<0>(s) == nr; }) }
			});

			lastSheet = sheetNr;
		}
	}

	// Each residue resides in a single strand which is part of a single sheet
	// this function returns the sequence number inside the sheet for the strand
	// containing res
	auto strandNrForResidue = [&strands,&sheetMap](dssp::residue_info const &res)
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

	for (auto &res : dssp)
	{
		if (res.type() != ss_type::Strand)
			continue;

		int s1 = strandNrForResidue(res);

		for (int i : { 0, 1 })
		{
			const auto &[p, ladder, parallel] = res.bridge_partner(i);
			if (not p or p.asym_id() != res.asym_id() or p.sheet() != res.sheet() or p.type() != ss_type::Strand)
				continue;

			int s2 = strandNrForResidue(p);
			// assert(s1 != s2);
			if (s2 == s1)
				continue;

			std::tuple<std::string,int> sheetID{res.asym_id(), res.sheet() };
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
			{ "sheet_id", cif::cif_id_for_number(sheet) },
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
					const auto &[e, ladder1, parallel1] = b.bridge_partner(i);
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

					if (test_bond(e, a) and test_bond(c, e))					// case I.
					{
						beg1SeqID = a.seq_id();
						beg2SeqID = e.seq_id();

						beg1AtomID = "O";
						beg2AtomID = "N";
					}
					else if (test_bond(b, d) and test_bond(f, b))				// case II.
					{
						beg1SeqID = b.seq_id();
						beg2SeqID = d.seq_id();

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
					const auto &[e, ladder1, parallel1] = b.bridge_partner(i);
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

					if (test_bond(a, e) and test_bond(e, c))					// case I.
					{
						end1SeqID = a.seq_id();
						end2SeqID = e.seq_id();

						end1AtomID = "N";
						end2AtomID = "O";
					}
					else if (test_bond(d, b) and test_bond(b, f))				// case II.
					{
						end1SeqID = b.seq_id();
						end2SeqID = d.seq_id();

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
					const auto &[e, ladder1, parallel1] = b.bridge_partner(i);
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

					if (test_bond(f, a) and test_bond(c, d))				// case III.
					{
						beg1SeqID = a.seq_id();
						beg2SeqID = f.seq_id();

						beg1AtomID = "O";
						beg2AtomID = "N";
					}
					else if (test_bond(b, e) and test_bond(e, b))			// case IV.
					{
						beg1SeqID = b.seq_id();
						beg2SeqID = e.seq_id();

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
					const auto &[e, ladder1, parallel1] = b.bridge_partner(i);
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

					if (test_bond(a, f) and test_bond(d, c))				// case III.
					{
						end1SeqID = a.seq_id();
						end2SeqID = f.seq_id();

						end1AtomID = "N";
						end2AtomID = "O";
					}
					else if (test_bond(b, e) and test_bond(e, b))			// case IV.
					{
						end1SeqID = b.seq_id();
						end2SeqID = e.seq_id();

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
			{ "sheet_id", cif::cif_id_for_number(sheet) },
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
		// 	if (test_bond(c, e) and test_bond(e, a))
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

			std::sort(strand.begin(), strand.end());

			auto &beg = strand.front();
			auto &end = strand.back();

			struct_sheet_range.emplace({
				{ "sheet_id", cif::cif_id_for_number(sheet) },
				{ "id", strandNrForResidue(strand.front()) + 1 },
				{ "beg_label_comp_id", beg.compound_id() },
				{ "beg_label_asym_id", beg.asym_id() },
				{ "beg_label_seq_id", beg.seq_id() },
				{ "pdbx_beg_PDB_ins_code", beg.pdb_ins_code() },
				{ "end_label_comp_id", end.compound_id() },
				{ "end_label_asym_id", end.asym_id() },
				{ "end_label_seq_id", end.seq_id() },
				{ "pdbx_end_PDB_ins_code", end.pdb_ins_code() },
				{ "beg_auth_comp_id", beg.compound_id() },
				{ "beg_auth_asym_id", beg.auth_asym_id() },
				{ "beg_auth_seq_id", beg.auth_seq_id() },
				{ "end_auth_comp_id", end.compound_id() },
				{ "end_auth_asym_id", end.auth_asym_id() },
				{ "end_auth_seq_id", end.auth_seq_id() }
			});
		}
	}


	return sheetMap;
}


void annotateDSSP(cif::datablock &db, const dssp &dssp, bool writeOther, std::ostream &os)
{
	using namespace std::literals;

	if (dssp.empty())
	{
		if (cif::VERBOSE > 0)
			std::cout << "No secondary structure information found" << std::endl;
	}
	else
	{
		auto stats = dssp.get_statistics();

		auto &dssp_statistics = db["dssp_statistics"];

		dssp_statistics.emplace({
			{ "entry_id", db.name() },
			{ "nr_of_residues", stats.count.residues },
			{ "nr_of_chains", stats.count.chains },
			{ "nr_of_ss_bridges_total", stats.count.SS_bridges },
			{ "nr_of_ss_bridges_intra_chain", stats.count.intra_chain_SS_bridges },
			{ "nr_of_ss_bridges_inter_chain", stats.count.SS_bridges - stats.count.intra_chain_SS_bridges },

			{ "accessible_surface_of_protein", stats.accessible_surface }
		});

		auto &dssp_hbonds = db["dssp_hbond_statistics"];

		dssp_hbonds.emplace({
			{ "entry_id", db.name() },
			{ "type", "O(I)-->H-N(J)" },
			{ "count", stats.count.H_bonds },
			{ "count_per_100", stats.count.H_bonds * 100.0 / stats.count.residues, 1 }
		});

		dssp_hbonds.emplace({
			{ "entry_id", db.name() },
			{ "type", "PARALLEL BRIDGES" },
			{ "count", stats.count.H_bonds_in_parallel_bridges },
			{ "count_per_100", stats.count.H_bonds_in_parallel_bridges * 100.0 / stats.count.residues, 1 }
		});

		dssp_hbonds.emplace({
			{ "entry_id", db.name() },
			{ "type", "ANTIPARALLEL BRIDGES" },
			{ "count", stats.count.H_bonds_in_antiparallel_bridges },
			{ "count_per_100", stats.count.H_bonds_in_antiparallel_bridges * 100.0 / stats.count.residues, 1 }
		});

		for (int k = 0; k < 11; ++k)
			dssp_hbonds.emplace({
				{ "entry_id", db.name() },
				{ "type", "O(I)-->H-N(I"s + char(k - 5 < 0 ? '-' : '+') + std::to_string(abs(k - 5)) + ")" },
				{ "count", stats.count.H_Bonds_per_distance[k] },
				{ "count_per_100", stats.count.H_Bonds_per_distance[k] * 100.0 / stats.count.residues, 1 }
			});

		auto &dssp_histograms = db["dssp_histograms"];
		
		using histogram_data_type = std::tuple<const char *, const uint32_t *>;
		for (const auto &[type, values] : std::vector<histogram_data_type>{
				{ "parallel_bridges_per_ladder", stats.histogram.residues_per_alpha_helix },
				{ "parallel_bridges_per_ladder", stats.histogram.parallel_bridges_per_ladder },
				{ "antiparallel_bridges_per_ladder", stats.histogram.antiparallel_bridges_per_ladder },
				{ "ladders_per_sheet", stats.histogram.ladders_per_sheet }
			})
		{
			auto hi = dssp_histograms.emplace({
				{ "entry_id", db.name() },
				{ "type", type },
				{ "1", values[0] },
				{ "2", values[1] },
				{ "3", values[2] },
				{ "4", values[3] },
				{ "5", values[4] },
				{ "6", values[5] },
				{ "7", values[6] },
				{ "8", values[7] },
				{ "9", values[8] },
				{ "10", values[9] },
				{ "11", values[10] },
				{ "12", values[11] },
				{ "13", values[12] },
				{ "14", values[13] },
				{ "15", values[14] },
				{ "16", values[15] },
				{ "17", values[16] },
				{ "18", values[17] },
				{ "19", values[18] },
				{ "20", values[19] },
				{ "21", values[20] },
				{ "22", values[21] },
				{ "23", values[22] },
				{ "24", values[23] },
				{ "25", values[24] },
				{ "26", values[25] },
				{ "27", values[26] },
				{ "28", values[27] },
				{ "29", values[28] },
				{ "30", values[29] },
			});
		}

		// A approximation of the old format

		auto &dssp_residue = db["dssp_residue"];

		for (auto res : dssp)
		{

			/*
				This is the header line for the residue lines in a DSSP file:

				#  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
			*/

			std::string ss_bridge = ".";
			if (res.ssBridgeNr() != 0)
				ss_bridge = std::to_string(res.ssBridgeNr());

			std::string ss = { res.type() == dssp::structure_type::Loop ? '.' : static_cast<char>(res.type()) };

			std::string helix[4] = { ".", ".", ".", "."};
			for (dssp::helix_type helixType : {dssp::helix_type::_3_10, dssp::helix_type::alpha, dssp::helix_type::pi, dssp::helix_type::pp})
			{
				switch (res.helix(helixType))
				{
					// case dssp::helix_position_type::None: helix[static_cast<int>(helixType)] = ' '; break;
					case dssp::helix_position_type::Start:
						helix[static_cast<int>(helixType)] = { '>' };
						break;

					case dssp::helix_position_type::End:
						helix[static_cast<int>(helixType)] = { '<' };
						break;

					case dssp::helix_position_type::StartAndEnd:
						helix[static_cast<int>(helixType)] = { 'X' };
						break;

					case dssp::helix_position_type::Middle:
						if (helixType == dssp::helix_type::pp)
							helix[static_cast<int>(helixType)] = { 'P' };
						else
							helix[static_cast<int>(helixType)] = { static_cast<char>('3' + static_cast<int>(helixType)) };
						break;
					
					default:
						break;
				}
			}

			std::string bend = ".";
			if (res.bend())
				bend = "S";

			double alpha = res.alpha();
			std::string chirality = ".";
			if (alpha != 360)
				chirality = alpha < 0 ? "-" : "+";

			std::string bridgepair[2] = { ".", "." };
			std::string bridgelabel[2] = { ".", "." };

			for (uint32_t i : {0, 1})
			{
				const auto &[p, ladder, parallel] = res.bridge_partner(i);
				if (not p)
					continue;

				bridgepair[i] = std::to_string(p.nr());
				bridgelabel[i] = (parallel ? "p-" : "a-") + std::to_string(ladder);
			}

			std::string sheet = ".";
			if (res.sheet() != 0)
				sheet = std::to_string(res.sheet());

			// std::string NHO[2], ONH[2];
			// for (int i : {0, 1})
			// {
			// 	const auto &[donor, donorE] = res.donor(i);
			// 	const auto &[acceptor, acceptorE] = res.acceptor(i);

			// 	NHO[i] = ONH[i] = "0, 0.0";

			// 	if (acceptor)
			// 	{
			// 		auto d = acceptor.nr() - res.nr();
			// 		NHO[i] = cif::format("%d,%3.1f", d, acceptorE).str();
			// 	}

			// 	if (donor)
			// 	{
			// 		auto d = donor.nr() - res.nr();
			// 		ONH[i] = cif::format("%d,%3.1f", d, donorE).str();
			// 	}
			// }

			auto const &[cax, cay, caz] = res.ca_location();

			dssp_residue.emplace({
				{ "entry_id", db.name() },
				{ "label_asym_id", res.asym_id() },
				{ "label_seq_id", res.seq_id() },
				{ "label_comp_id", res.compound_id() },

				{ "secondary_structure", ss },

				{ "ss_bridge", ss_bridge },

				{ "helix_3_10", helix[0] },
				{ "helix_alpha", helix[1] },
				{ "helix_pi", helix[2] },
				{ "helix_pp", helix[3] },

				{ "bend", bend },
				{ "chirality", chirality },

				{ "bridge_pair_1", bridgepair[0] },
				{ "bridge_label_1", bridgelabel[0] },

				{ "bridge_pair_2", bridgepair[1] },
				{ "bridge_label_2", bridgelabel[1] },

				{ "sheet", sheet },

				{ "accesssibility", res.accessibility(), 1 },

				{ "TCO", res.tco(), 3 },
				{ "kappa", res.kappa(), 1 },
				{ "alpha", res.alpha(), 1 },
				{ "phi", res.phi(), 1 },
				{ "psi", res.psi(), 1 },

				{ "x_ca", cax, 1 },
				{ "y_ca", cay, 1 },
				{ "z_ca", caz, 1 },


				
				
			});

		}










		// replace all struct_conf and struct_conf_type records
		auto &structConfType = db["struct_conf_type"];
		structConfType.clear();

		auto &structConf = db["struct_conf"];
		structConf.clear();

		std::map<std::string, int> foundTypes;

		auto st = dssp.begin(), lt = st;
		auto lastSS = st->type();

		for (auto t = dssp.begin();; lt = t, ++t)
		{
			bool stop = t == dssp.end();

			bool flush = (stop or t->type() != lastSS);

			if (flush and (writeOther or lastSS != dssp::structure_type::Loop))
			{
				auto &rb = *st;
				auto &re = *lt;

				std::string id;
				switch (lastSS)
				{
					case dssp::structure_type::Helix_3:
						id = "HELX_RH_3T_P";
						break;

					case dssp::structure_type::Alphahelix:
						id = "HELX_RH_AL_P";
						break;

					case dssp::structure_type::Helix_5:
						id = "HELX_RH_PI_P";
						break;

					case dssp::structure_type::Helix_PPII:
						id = "HELX_LH_PP_P";
						break;

					case dssp::structure_type::Turn:
						id = "TURN_TY1_P";
						break;

					case dssp::structure_type::Bend:
						id = "BEND";
						break;

					case dssp::structure_type::Betabridge:
					case dssp::structure_type::Strand:
						id = "STRN";
						break;

					case dssp::structure_type::Loop:
						id = "OTHER";
						break;
				}

				if (foundTypes.count(id) == 0)
				{
					structConfType.emplace({ { "id", id },
						{ "criteria", "DSSP" } });
					foundTypes[id] = 1;
				}

				structConf.emplace({
					{ "conf_type_id", id },
					{ "id", id + std::to_string(foundTypes[id]++) },
					// { "pdbx_PDB_helix_id", vS(12, 14) },
					{ "beg_label_comp_id", rb.compound_id() },
					{ "beg_label_asym_id", rb.asym_id() },
					{ "beg_label_seq_id", rb.seq_id() },
					{ "pdbx_beg_PDB_ins_code", rb.pdb_ins_code() },
					{ "end_label_comp_id", re.compound_id() },
					{ "end_label_asym_id", re.asym_id() },
					{ "end_label_seq_id", re.seq_id() },
					{ "pdbx_end_PDB_ins_code", re.pdb_ins_code() },

					{ "beg_auth_comp_id", rb.compound_id() },
					{ "beg_auth_asym_id", rb.auth_asym_id() },
					{ "beg_auth_seq_id", rb.auth_seq_id() },
					{ "end_auth_comp_id", re.compound_id() },
					{ "end_auth_asym_id", re.auth_asym_id() },
					{ "end_auth_seq_id", re.auth_seq_id() }

					// { "pdbx_PDB_helix_class", vS(39, 40) },
				    // { "details", vS(41, 70) },
				    // { "pdbx_PDB_helix_length", vI(72, 76) }
				});

				st = t;
			}

			if (stop)
				break;

			if (lastSS != t->type())
			{
				st = t;
				lastSS = t->type();
			}
		}
	}

	auto &software = db["software"];
	software.emplace({
		{ "pdbx_ordinal", software.get_unique_id("") },
		{ "name", "dssp" },
		{ "version", kVersionNumber },
		{ "date", kBuildDate },
		{ "classification", "model annotation" } });

	db.write(os);
}
