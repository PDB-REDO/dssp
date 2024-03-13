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

#include <cif++.hpp>
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

	if (residue.pdb_strand_id().length() > 1)
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

	std::string version = kVersionNumber;
	if (version.length() < 10)
		version.insert(version.end(), 10 - version.length(), ' ');

	os << "==== Secondary Structure Definition by the program DSSP, NKI version " << version << "                    ==== DATE=" << std::put_time(tm, "%F") << "        ." << std::endl
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

	hb.add_item("id");
	hb.add_item("label_comp_id");
	hb.add_item("label_seq_id");
	hb.add_item("label_asym_id");
	hb.add_item("auth_seq_id");
	hb.add_item("auth_asym_id");
	hb.add_item("pdbx_PDB_ins_code");

	// force right order
	for (std::string da : { "acceptor_", "donor_" })
	{
		for (std::string i : { "1_", "2_" })
		{
			for (std::string n : { "label_comp_id", "label_seq_id", "label_asym_id", "auth_seq_id", "auth_asym_id", "pdbx_PDB_ins_code", "energy" })
				hb.add_item(da + i + n);
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
	using namespace cif::literals;

	using res_list = std::vector<dssp::residue_info>;

	// clean up old info first

	for (auto sheet_cat : { "struct_sheet", "struct_sheet_order", "struct_sheet_range", "struct_sheet_hbond", "pdbx_struct_sheet_hbond" })
		db.erase(remove_if(db.begin(), db.end(), [sheet_cat](const cif::category &cat)
					 { return cat.name() == sheet_cat; }),
			db.end());

	// create a list of strands, based on the SS info in DSSP. Store sheet number along with the strand.

	std::map<std::tuple<int,int>, res_list> strands;
	std::set<int> sheetNrs;

	for (auto &res : dssp)
	{
		if (res.type() != dssp::structure_type::Strand and res.type() != dssp::structure_type::Betabridge)
			continue;

		strands[{res.sheet(), res.strand()}].emplace_back(res);
		sheetNrs.insert(res.sheet());
	}

	// --------------------------------------------------------------------

	auto &struct_sheet = db["struct_sheet"];
	auto &struct_sheet_range = db["struct_sheet_range"];

	for (int sheetNr : sheetNrs)
	{
		auto sheetID = cif::cif_id_for_number(sheetNr - 1);

		struct_sheet.emplace({
			{ "id", sheetID },
			{ "number_strands",
				std::count_if(strands.begin(), strands.end(), [nr = sheetNr](std::tuple<std::tuple<int,int>, res_list> const &s)
				{
					const auto &[strandID, strand] = s;
					return strand.front().sheet() == nr;
				})
			}
		});

		for (auto &&[strandTuple, strand] : strands)
		{
			if (strand.front().sheet() != sheetNr)
				continue;
			
			std::string strandID = cif::cif_id_for_number(strand.front().strand() - 1);

			std::sort(strand.begin(), strand.end(), [](dssp::residue_info const &a, dssp::residue_info const &b)
				{ return a.nr() < b.nr(); });

			auto &beg = strand.front();
			auto &end = strand.back();

			struct_sheet_range.emplace({
				{ "sheet_id", sheetID },
				{ "id", strandID },
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
		ladder_info(int label, int sheet, bool parallel, const dssp::residue_info &a, const dssp::residue_info &b)
			: ladder(label)
			, sheet(sheet)
			, parallel(parallel)
			, pairs({ { a, b } })
		{
		}

		bool operator<(const ladder_info &rhs) const
		{
			return ladder < rhs.ladder;
		}

		int ladder;
		int sheet;
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

			bool is_new = true;
			for (auto &l : ladders)
			{
				if (l.ladder != ladder)
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

			ladders.emplace_back(ladder, res.sheet() - 1, parallel, res, p);
		}
	}

	std::sort(ladders.begin(), ladders.end());

	auto &dssp_struct_ladder = db["dssp_struct_ladder"];

	for (const auto &l : ladders)
	{
		const auto &[beg1, beg2] = l.pairs.front();
		const auto &[end1, end2] = l.pairs.back();

		dssp_struct_ladder.emplace({
			{ "id", cif::cif_id_for_number(l.ladder) },
			{ "sheet_id", cif::cif_id_for_number(l.sheet) },
			{ "range_id_1", cif::cif_id_for_number(beg1.strand() - 1) },
			{ "range_id_2", cif::cif_id_for_number(beg2.strand() - 1) },
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

	for (auto label : { "entry_id", "label_comp_id", "label_asym_id", "label_seq_id", "secondary_structure",
			"ss_bridge", "helix_3_10", "helix_alpha", "helix_pi", "helix_pp", "bend", "chirality", "sheet",
			"strand", "ladder_1", "ladder_2", "accessibility", "TCO", "kappa", "alpha", "phi", "psi",
			"x_ca", "y_ca", "z_ca"})
		dssp_struct_summary.add_item(label);

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

			{ "sheet", res.sheet() ? cif::cif_id_for_number(res.sheet() - 1) : "." },
			{ "strand", res.strand() ? cif::cif_id_for_number(res.strand() - 1) : "." },
			{ "ladder_1", ladders[0] },
			{ "ladder_2", ladders[1] },

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

void annotateDSSP(cif::datablock &db, const dssp &dssp, bool writeOther, bool writeExperimental)
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
		{ "date", kRevisionDate },
		{ "classification", "model annotation" }
	});
}
