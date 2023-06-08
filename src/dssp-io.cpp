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

#include "dssp-io.hpp"
#include "revision.hpp"

#include <cif++/pdb/io.hpp>
#include <cif++/dictionary_parser.hpp>

#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>


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

	char helix[4] = { ' ', ' ', ' ', ' ' };
	for (dssp::helix_type helixType : { dssp::helix_type::_3_10, dssp::helix_type::alpha, dssp::helix_type::pi, dssp::helix_type::pp })
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

	double alpha = residue.alpha().value_or(360);
	char chirality = alpha == 360 ? ' ' : (alpha < 0 ? '-' : '+');

	uint32_t bp[2] = {};
	char bridgelabel[2] = { ' ', ' ' };
	for (uint32_t i : { 0, 1 })
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
	for (int i : { 0, 1 })
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
		residue.tco().value_or(0),
		residue.kappa().value_or(360),
		residue.alpha().value_or(360),
		residue.phi().value_or(360),
		residue.psi().value_or(360),
		cax, cay, caz)
	    .str();
}

void writeDSSP(const dssp &dssp, std::ostream &os)
{
	using namespace std::chrono;

	auto stats = dssp.get_statistics();

	std::time_t today = system_clock::to_time_t(system_clock::now());
	std::tm *tm = std::gmtime(&today);

	os << "==== Secondary Structure Definition by the program DSSP, NKI version 4.3                           ==== DATE=" << std::put_time(tm, "%F") << "        ." << std::endl
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
					  (last + 1), (ri.chain_break() == dssp::chain_break_type::NewChain ? '*' : ' '))
			   << std::endl;

		os << ResidueToDSSPLine(ri) << std::endl;
		last = ri.nr();
	}
}

// --------------------------------------------------------------------

void writeBridgePairs(cif::datablock &db, const dssp &dssp)
{
	auto &hb = db["dssp_struct_bridge_pairs"];

	hb.add_column("id");
	hb.add_column("label_comp_id");
	hb.add_column("label_seq_id");
	hb.add_column("label_asym_id");
	hb.add_column("auth_seq_id");
	hb.add_column("auth_asym_id");
	hb.add_column("pdbx_PDB_ins_code");

	// force right order
	for (std::string da : { "acceptor_", "donor_" })
	{
		for (std::string i : { "1_", "2_" })
		{
			for (std::string n : { "label_comp_id", "label_seq_id", "label_asym_id", "auth_seq_id", "auth_asym_id", "pdbx_PDB_ins_code", "energy" })
				hb.add_column(da + i + n);
		}
	}

	for (auto &res : dssp)
	{
		cif::row_initializer data({
			{ "id", hb.get_unique_id("") },
			{ "label_comp_id", res.compound_id() },
			{ "label_seq_id", res.seq_id() },
			{ "label_asym_id", res.asym_id() },
			// { "auth_comp_id", res.compound_id() },
			{ "auth_seq_id", res.auth_seq_id() },
			{ "auth_asym_id", res.auth_asym_id() },
			{ "pdbx_PDB_ins_code", res.pdb_ins_code() }
		});

		for (int i : { 0, 1 })
		{
			const auto &&[acceptor, acceptorEnergy] = res.acceptor(i);
			const auto &&[donor, donorEnergy] = res.donor(i);

			if (acceptor)
			{
				if (i == 0)
				{
					data.emplace_back("acceptor_1_label_comp_id", acceptor.compound_id());
					data.emplace_back("acceptor_1_label_seq_id", acceptor.seq_id());
					data.emplace_back("acceptor_1_label_asym_id", acceptor.asym_id());
					// data.emplace_back("acceptor_1_auth_comp_id", acceptor.compound_id());
					data.emplace_back("acceptor_1_auth_seq_id", acceptor.auth_seq_id());
					data.emplace_back("acceptor_1_auth_asym_id", acceptor.auth_asym_id());
					data.emplace_back("acceptor_1_pdbx_PDB_ins_code", acceptor.pdb_ins_code());
					data.emplace_back("acceptor_1_energy", acceptorEnergy, 1);
				}
				else
				{
					data.emplace_back("acceptor_2_label_comp_id", acceptor.compound_id());
					data.emplace_back("acceptor_2_label_seq_id", acceptor.seq_id());
					data.emplace_back("acceptor_2_label_asym_id", acceptor.asym_id());
					// data.emplace_back("acceptor_2_auth_comp_id", acceptor.compound_id());
					data.emplace_back("acceptor_2_auth_seq_id", acceptor.auth_seq_id());
					data.emplace_back("acceptor_2_auth_asym_id", acceptor.auth_asym_id());
					data.emplace_back("acceptor_2_pdbx_PDB_ins_code", acceptor.pdb_ins_code());
					data.emplace_back("acceptor_2_energy", acceptorEnergy, 1);
				}
			}

			if (donor)
			{
				if (i == 0)
				{
					data.emplace_back("donor_1_label_comp_id", donor.compound_id());
					data.emplace_back("donor_1_label_seq_id", donor.seq_id());
					data.emplace_back("donor_1_label_asym_id", donor.asym_id());
					// data.emplace_back("donor_1_auth_comp_id", donor.compound_id());
					data.emplace_back("donor_1_auth_seq_id", donor.auth_seq_id());
					data.emplace_back("donor_1_auth_asym_id", donor.auth_asym_id());
					data.emplace_back("donor_1_pdbx_PDB_ins_code", donor.pdb_ins_code());
					data.emplace_back("donor_1_energy", donorEnergy, 1);
				}
				else
				{
					data.emplace_back("donor_2_label_comp_id", donor.compound_id());
					data.emplace_back("donor_2_label_seq_id", donor.seq_id());
					data.emplace_back("donor_2_label_asym_id", donor.asym_id());
					// data.emplace_back("donor_2_auth_comp_id", donor.compound_id());
					data.emplace_back("donor_2_auth_seq_id", donor.auth_seq_id());
					data.emplace_back("donor_2_auth_asym_id", donor.auth_asym_id());
					data.emplace_back("donor_2_pdbx_PDB_ins_code", donor.pdb_ins_code());
					data.emplace_back("donor_2_energy", donorEnergy, 1);
				}
			}
		}

		hb.emplace(std::move(data));
	}
}

void writeSheets(cif::datablock &db, const dssp &dssp)
{
	using res_list = std::vector<dssp::residue_info>;
	using ss_type = dssp::structure_type;

	// clean up old info first

	for (auto sheet_cat : { "struct_sheet", "struct_sheet_order", "struct_sheet_range", "struct_sheet_hbond", "pdbx_struct_sheet_hbond" })
		db.erase(remove_if(db.begin(), db.end(), [sheet_cat](const cif::category &cat)
					 { return cat.name() == sheet_cat; }),
			db.end());

	// create a list of strands, based on the SS info in DSSP. Store sheet number along with the strand.

	std::vector<std::tuple<int, res_list>> strands;
	std::map<std::tuple<std::string, int>, int> sheetMap; // mapping from bridge number (=info.sheet()) in DSSP to sheet ID

	ss_type ss = ss_type::Loop;

	for (auto &res : dssp)
	{
		std::tuple<std::string, int> sheetID{ res.asym_id(), res.sheet() };
		ss_type iss = res.type();

		if (iss == ss and ss == ss_type::Strand and
			sheetMap.count(sheetID) > 0 and sheetMap[sheetID] == std::get<0>(strands.back()))
		{
			std::get<1>(strands.back()).emplace_back(res);
			continue;
		}

		ss = iss;

		if (ss != ss_type::Strand)
			continue;

		if (not sheetMap.count(sheetID))
			sheetMap[sheetID] = static_cast<int>(sheetMap.size());

		strands.emplace_back(sheetMap[sheetID], res_list{ res });
	}

	// sort the strands vector
	std::sort(strands.begin(), strands.end(), [](auto &a, auto &b)
		{
		const auto &[sheetA, strandsA] = a;
		const auto &[sheetB, strandsB] = b;

		int d = sheetA - sheetB;
		if (d == 0)
			d = strandsA.front().nr() - strandsB.front().nr();
		return d < 0; });

	// write out the struct_sheet, since all info is available now
	auto &struct_sheet = db["struct_sheet"];

	int lastSheet = -1;
	for (const auto &[sheetNr, strand] : strands)
	{
		if (sheetNr != lastSheet)
		{
			struct_sheet.emplace({
				{ "id", cif::cif_id_for_number(sheetNr) },
				{ "number_strands",
					std::count_if(strands.begin(), strands.end(), [nr = sheetNr](std::tuple<int, res_list> const &s)
						{ return std::get<0>(s) == nr; })
				}
			});

			lastSheet = sheetNr;
		}
	}

	// Each residue resides in a single strand which is part of a single sheet
	// this function returns the sequence number inside the sheet for the strand
	// containing res
	std::map<int,int> strandMap;
	int strandNr = 0;

	for (auto &&[sheet, strand] : strands)
	{
		for (auto &r : strand)
		{
			if (r.type() != ss_type::Strand)
				continue;
			strandMap[r.nr()] = strandNr;
		}

		++strandNr;
	}

	// This map is used to record the sense of a ladder, and can be used
	// to detect ladders already seen.
	std::map<std::tuple<int, int, int>, std::tuple<int, bool>> ladderSense;

	for (auto &res : dssp)
	{
		if (res.type() != ss_type::Strand)
			continue;

		int s1 = strandMap[res.nr()];

		for (int i : { 0, 1 })
		{
			const auto &[p, ladder, parallel] = res.bridge_partner(i);
			if (not p or p.asym_id() != res.asym_id() or p.sheet() != res.sheet() or p.type() != ss_type::Strand)
				continue;

			int s2 = strandMap[p.nr()];
			// assert(s1 != s2);
			if (s2 == s1)
				continue;

			std::tuple<std::string, int> sheetID{ res.asym_id(), res.sheet() };
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
	auto &pdbx_struct_sheet_hbond = db["pdbx_struct_sheet_hbond"];
	auto &pdbx_poly_seq_scheme = db["pdbx_poly_seq_scheme"];

	for (const auto &[key, value] : ladderSense)
	{
		const auto &[sheet, s1, s2] = key;
		const auto &[ladder, parallel] = value;

		struct_sheet_order.emplace({ { "sheet_id", cif::cif_id_for_number(sheet) },
			// { "dssp_struct_ladder_id", cif::cifIdForNumber(ladder) },
			{ "range_id_1", s1 + 1 },
			{ "range_id_2", s2 + 1 },
			{ "sense", parallel ? "parallel" : "anti-parallel" } });

		res_list strand1, strand2;

		for (int strandIx = 0; static_cast<size_t>(strandIx) < strands.size(); ++strandIx)
		{
			const auto &[sSheet, strand] = strands[strandIx];
			if (sSheet != sheet)
				continue;

			if (strandIx == s1)
				strand1 = strand;
			else if (strandIx == s2)
			{
				strand2 = strand;
				break;
			}
		}

		assert(not(strand1.empty() or strand2.empty()));

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

					if (test_bond(e, a) and test_bond(c, e)) // case I.
					{
						beg1SeqID = a.seq_id();
						beg2SeqID = e.seq_id();

						beg1AtomID = "O";
						beg2AtomID = "N";
					}
					else if (test_bond(b, d) and test_bond(f, b)) // case II.
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

					if (test_bond(a, e) and test_bond(e, c)) // case I.
					{
						end1SeqID = a.seq_id();
						end2SeqID = e.seq_id();

						end1AtomID = "N";
						end2AtomID = "O";
					}
					else if (test_bond(d, b) and test_bond(b, f)) // case II.
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

					if (test_bond(f, a) and test_bond(c, d)) // case III.
					{
						beg1SeqID = a.seq_id();
						beg2SeqID = f.seq_id();

						beg1AtomID = "O";
						beg2AtomID = "N";
					}
					else if (test_bond(b, e) and test_bond(e, b)) // case IV.
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

					if (test_bond(a, f) and test_bond(d, c)) // case III.
					{
						end1SeqID = a.seq_id();
						end2SeqID = f.seq_id();

						end1AtomID = "N";
						end2AtomID = "O";
					}
					else if (test_bond(b, e) and test_bond(e, b)) // case IV.
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

		struct_sheet_hbond.emplace({ { "sheet_id", cif::cif_id_for_number(sheet) },
			{ "range_id_1", s1 + 1 },
			{ "range_id_2", s2 + 1 },
			{ "range_1_beg_label_seq_id", beg1SeqID },
			{ "range_1_beg_label_atom_id", beg1AtomID },
			{ "range_2_beg_label_seq_id", beg2SeqID },
			{ "range_2_beg_label_atom_id", beg2AtomID },
			{ "range_1_end_label_seq_id", end1SeqID },
			{ "range_1_end_label_atom_id", end1AtomID },
			{ "range_2_end_label_seq_id", end2SeqID },
			{ "range_2_end_label_atom_id", end2AtomID } });

		using namespace cif::literals;

		auto b1 = pdbx_poly_seq_scheme.find_first("asym_id"_key == strand1.front().asym_id() and "seq_id"_key == beg1SeqID);
		auto e1 = pdbx_poly_seq_scheme.find_first("asym_id"_key == strand1.front().asym_id() and "seq_id"_key == end1SeqID);
		auto b2 = pdbx_poly_seq_scheme.find_first("asym_id"_key == strand2.front().asym_id() and "seq_id"_key == beg2SeqID);
		auto e2 = pdbx_poly_seq_scheme.find_first("asym_id"_key == strand2.front().asym_id() and "seq_id"_key == end2SeqID);

		if (not(b1 and e1 and b2 and e2))
		{
			if (cif::VERBOSE > 0)
				std::cerr << "error looking up strand ends" << std::endl;
			continue;
		}

		pdbx_struct_sheet_hbond.emplace({
			{ "sheet_id", cif::cif_id_for_number(sheet) },
			{ "range_id_1", s1 + 1 },
			{ "range_id_2", s2 + 1 },

			{ "range_1_label_atom_id", beg1AtomID },
			{ "range_1_label_comp_id", b1["mon_id"].as<std::string>() },
			{ "range_1_label_asym_id", b1["asym_id"].as<std::string>() },
			{ "range_1_label_seq_id", b1["seq_id"].as<std::string>() },
			{ "range_1_PDB_ins_code", b1["pdb_ins_code"].as<std::string>() },
			{ "range_1_auth_atom_id", beg1AtomID },
			{ "range_1_auth_comp_id", b1["auth_mon_id"].as<std::string>() },
			{ "range_1_auth_asym_id", b1["pdb_strand_id"].as<std::string>() },
			{ "range_1_auth_seq_id", b1["auth_seq_num"].as<std::string>() },
			{ "range_2_label_atom_id", beg2AtomID },
			{ "range_2_label_comp_id", b2["mon_id"].as<std::string>() },
			{ "range_2_label_asym_id", b2["asym_id"].as<std::string>() },
			{ "range_2_label_seq_id", b2["seq_id"].as<std::string>() },
			{ "range_2_PDB_ins_code", b2["pdb_ins_code"].as<std::string>() },
			{ "range_2_auth_atom_id", beg2AtomID },
			{ "range_2_auth_comp_id", b2["auth_mon_id"].as<std::string>() },
			{ "range_2_auth_asym_id", b2["pdb_strand_id"].as<std::string>() },
			{ "range_2_auth_seq_id", b2["auth_seq_num"].as<std::string>() },
		});
	}

	auto &struct_sheet_range = db["struct_sheet_range"];

	for (const auto &[key, iSheet] : sheetMap)
	{
		for (auto &&[sheet, strand] : strands)
		{
			if (sheet != iSheet)
				continue;

			std::sort(strand.begin(), strand.end(), [](dssp::residue_info const &a, dssp::residue_info const &b)
				{ return a.nr() < b.nr(); });

			auto &beg = strand.front();
			auto &end = strand.back();

			struct_sheet_range.emplace({
				{ "sheet_id", cif::cif_id_for_number(sheet) },
				{ "id", strandMap[strand.front().nr()] + 1 },
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
				{ "end_auth_seq_id", end.auth_seq_id() } });
		}
	}
}

void writeLadders(cif::datablock &db, const dssp &dssp)
{
	// Write out the DSSP ladders
	struct ladder_info
	{
		ladder_info(const std::string &label, bool parallel, const dssp::residue_info &a, const dssp::residue_info &b)
			: label(label)
			, parallel(parallel)
			, pairs({ { a, b } })
		{
		}

		std::string label;
		bool parallel;
		std::vector<std::pair<dssp::residue_info, dssp::residue_info>> pairs;
	};

	std::vector<ladder_info> ladders;

	for (auto res : dssp)
	{
		for (int i : { 0, 1 })
		{
			const auto &[p, ladder, parallel] = res.bridge_partner(i);

			if (not p)
				continue;

			auto label = cif::cif_id_for_number(ladder);

			bool is_new = true;
			for (auto &l : ladders)
			{
				if (l.label != label)
					continue;

				assert(l.parallel == parallel);

				if (find_if(l.pairs.begin(), l.pairs.end(), [na = p.nr(), nb = res.nr()](const auto &p)
						{ return p.first.nr() == na and p.second.nr() == nb; }) != l.pairs.end())
				{
					is_new = false;
					break;
				}

				l.pairs.emplace_back(res, p);
				is_new = false;
				break;
			}

			if (not is_new)
				continue;

			ladders.emplace_back(label, parallel, res, p);
		}
	}

	auto &dssp_struct_ladder = db["dssp_struct_ladder"];

	for (const auto &l : ladders)
	{
		const auto &[beg1, beg2] = l.pairs.front();
		const auto &[end1, end2] = l.pairs.back();

		dssp_struct_ladder.emplace({ { "id", l.label },
			{ "type", l.parallel ? "parallel" : "anti-parallel" },

			{ "beg_1_label_comp_id", beg1.compound_id() },
			{ "beg_1_label_asym_id", beg1.asym_id() },
			{ "beg_1_label_seq_id", beg1.seq_id() },
			{ "pdbx_beg_1_PDB_ins_code", beg1.pdb_ins_code() },
			{ "end_1_label_comp_id", end1.compound_id() },
			{ "end_1_label_asym_id", end1.asym_id() },
			{ "end_1_label_seq_id", end1.seq_id() },
			{ "pdbx_end_1_PDB_ins_code", end1.pdb_ins_code() },
			{ "beg_1_auth_comp_id", beg1.compound_id() },
			{ "beg_1_auth_asym_id", beg1.auth_asym_id() },
			{ "beg_1_auth_seq_id", beg1.auth_seq_id() },
			{ "end_1_auth_comp_id", end1.compound_id() },
			{ "end_1_auth_asym_id", end1.auth_asym_id() },
			{ "end_1_auth_seq_id", end1.auth_seq_id() },

			{ "beg_2_label_comp_id", beg2.compound_id() },
			{ "beg_2_label_asym_id", beg2.asym_id() },
			{ "beg_2_label_seq_id", beg2.seq_id() },
			{ "pdbx_beg_2_PDB_ins_code", beg2.pdb_ins_code() },
			{ "end_2_label_comp_id", end2.compound_id() },
			{ "end_2_label_asym_id", end2.asym_id() },
			{ "end_2_label_seq_id", end2.seq_id() },
			{ "pdbx_end_2_PDB_ins_code", end2.pdb_ins_code() },
			{ "beg_2_auth_comp_id", beg2.compound_id() },
			{ "beg_2_auth_asym_id", beg2.auth_asym_id() },
			{ "beg_2_auth_seq_id", beg2.auth_seq_id() },
			{ "end_2_auth_comp_id", end2.compound_id() },
			{ "end_2_auth_asym_id", end2.auth_asym_id() },
			{ "end_2_auth_seq_id", end2.auth_seq_id() } });
	}
}

void writeStatistics(cif::datablock &db, const dssp &dssp)
{
	using namespace std::literals;

	auto stats = dssp.get_statistics();

	auto &dssp_statistics = db["dssp_statistics"];

	auto stats_i = dssp_statistics.emplace({ { "entry_id", db.name() },
		{ "nr_of_residues", stats.count.residues },
		{ "nr_of_chains", stats.count.chains },
		{ "nr_of_ss_bridges_total", stats.count.SS_bridges },
		{ "nr_of_ss_bridges_intra_chain", stats.count.intra_chain_SS_bridges },
		{ "nr_of_ss_bridges_inter_chain", stats.count.SS_bridges - stats.count.intra_chain_SS_bridges } });
	
	if (stats.accessible_surface > 0)
		(*stats_i)["accessible_surface_of_protein"] = stats.accessible_surface;

	auto &dssp_struct_hbonds = db["dssp_statistics_hbond"];

	dssp_struct_hbonds.emplace({ { "entry_id", db.name() },
		{ "type", "O(I)-->H-N(J)" },
		{ "count", stats.count.H_bonds },
		{ "count_per_100", stats.count.H_bonds * 100.0 / stats.count.residues, 1 } });

	dssp_struct_hbonds.emplace({ { "entry_id", db.name() },
		{ "type", "PARALLEL BRIDGES" },
		{ "count", stats.count.H_bonds_in_parallel_bridges },
		{ "count_per_100", stats.count.H_bonds_in_parallel_bridges * 100.0 / stats.count.residues, 1 } });

	dssp_struct_hbonds.emplace({ { "entry_id", db.name() },
		{ "type", "ANTIPARALLEL BRIDGES" },
		{ "count", stats.count.H_bonds_in_antiparallel_bridges },
		{ "count_per_100", stats.count.H_bonds_in_antiparallel_bridges * 100.0 / stats.count.residues, 1 } });

	for (int k = 0; k < 11; ++k)
		dssp_struct_hbonds.emplace({ { "entry_id", db.name() },
			{ "type", "O(I)-->H-N(I"s + char(k - 5 < 0 ? '-' : '+') + std::to_string(abs(k - 5)) + ")" },
			{ "count", stats.count.H_Bonds_per_distance[k] },
			{ "count_per_100", stats.count.H_Bonds_per_distance[k] * 100.0 / stats.count.residues, 1 } });

	auto &dssp_statistics_histogram = db["dssp_statistics_histogram"];

	using histogram_data_type = std::tuple<const char *, const uint32_t *>;
	for (const auto &[type, values] : std::vector<histogram_data_type>{
			 { "residues_per_alpha_helix", stats.histogram.residues_per_alpha_helix },
			 { "parallel_bridges_per_ladder", stats.histogram.parallel_bridges_per_ladder },
			 { "antiparallel_bridges_per_ladder", stats.histogram.antiparallel_bridges_per_ladder },
			 { "ladders_per_sheet", stats.histogram.ladders_per_sheet } })
	{
		auto hi = dssp_statistics_histogram.emplace({
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
}

void writeSummary(cif::datablock &db, const dssp &dssp)
{
	bool writeAccessibility = dssp.get_statistics().accessible_surface > 0;

	// A approximation of the old format

	auto &dssp_struct_summary = db["dssp_struct_summary"];

	// prime the category with the field labels we need, this is to ensure proper order in writing out the data.

	for (auto label : { "entry_id", "label_comp_id", "label_asym_id", "label_seq_id", "secondary_structure", "ss_bridge", "helix_3_10", "helix_alpha", "helix_pi", "helix_pp", "bend", "chirality", "ladder_1", "ladder_2", "sheet", "accessibility", "TCO", "kappa", "alpha", "phi", "psi", "x_ca", "y_ca", "z_ca"})
		dssp_struct_summary.add_column(label);

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

		std::string helix[4] = { ".", ".", ".", "." };
		for (dssp::helix_type helixType : { dssp::helix_type::_3_10, dssp::helix_type::alpha, dssp::helix_type::pi, dssp::helix_type::pp })
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

		std::string chirality = ".";
		if (res.alpha().has_value())
			chirality = *res.alpha() < 0 ? "-" : "+";

		std::string ladders[2] = { ".", "." };

		for (uint32_t i : { 0, 1 })
		{
			const auto &[p, ladder, parallel] = res.bridge_partner(i);
			if (not p)
				continue;

			ladders[i] = cif::cif_id_for_number(ladder);
		}

		auto const &[cax, cay, caz] = res.ca_location();

		cif::row_initializer data{
			{ "entry_id", db.name() },
			{ "label_comp_id", res.compound_id() },
			{ "label_asym_id", res.asym_id() },
			{ "label_seq_id", res.seq_id() },

			{ "secondary_structure", ss },

			{ "ss_bridge", ss_bridge },

			{ "helix_3_10", helix[0] },
			{ "helix_alpha", helix[1] },
			{ "helix_pi", helix[2] },
			{ "helix_pp", helix[3] },

			{ "bend", bend },
			{ "chirality", chirality },

			{ "ladder_1", ladders[0] },
			{ "ladder_2", ladders[1] },

			{ "sheet", res.sheet() ? cif::cif_id_for_number(res.sheet() - 1) : "." },

			{ "x_ca", cax, 1 },
			{ "y_ca", cay, 1 },
			{ "z_ca", caz, 1 },
		};

		if (writeAccessibility)
			data.emplace_back("accessibility", res.accessibility(), 1);

		if (res.tco().has_value())
			data.emplace_back("TCO", *res.tco(), 3);
		else
			data.emplace_back("TCO", ".");

		if (res.kappa().has_value())
			data.emplace_back("kappa", *res.kappa(), 1);
		else
			data.emplace_back("kappa", ".");

		if (res.alpha().has_value())
			data.emplace_back("alpha", *res.alpha(), 1);
		else
			data.emplace_back("alpha", ".");

		if (res.phi().has_value())
			data.emplace_back("phi", *res.phi(), 1);
		else
			data.emplace_back("phi", ".");

		if (res.psi().has_value())
			data.emplace_back("psi", *res.psi(), 1);
		else
			data.emplace_back("psi", ".");
		
		dssp_struct_summary.emplace(std::move(data));
	}
}

void annotateDSSP(cif::datablock &db, const dssp &dssp, bool writeOther, bool writeExperimental, std::ostream &os)
{
	using namespace std::literals;

	auto &validator = const_cast<cif::validator &>(*db.get_validator());
	if (validator.get_validator_for_category("dssp_struct_summary") == nullptr)
	{
		auto dssp_extension = cif::load_resource("dssp-extension.dic");
		if (dssp_extension)
			cif::extend_dictionary(validator, *dssp_extension);
	}

	if (dssp.empty())
	{
		if (cif::VERBOSE > 0)
			std::cout << "No secondary structure information found" << std::endl;
	}
	else
	{
		if (writeExperimental)
		{
			writeBridgePairs(db, dssp);
			writeSheets(db, dssp);
			writeLadders(db, dssp);
			writeStatistics(db, dssp);
			writeSummary(db, dssp);
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
	software.emplace({ { "pdbx_ordinal", software.get_unique_id("") },
		{ "name", "dssp" },
		{ "version", kVersionNumber },
		{ "date", kBuildDate },
		{ "classification", "model annotation" } });

	db.write(os);
}
