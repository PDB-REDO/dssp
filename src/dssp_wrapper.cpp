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

#include <cif++.hpp>
#include <date/date.h>
#include <pdbx++.hpp>

#include "dssp_wrapper.hpp"
#include "revision.hpp"

// --------------------------------------------------------------------

std::string ResidueToDSSPLine(const dssp::DSSP::residue_info &info)
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

	char ss;
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
			NHO[i] = pdbx::format("%d,%3.1f", d, acceptorE).str();
		}

		if (donor)
		{
			auto d = donor.nr() - info.nr();
			ONH[i] = pdbx::format("%d,%3.1f", d, donorE).str();
		}
	}

	// auto ca = residue.atomByID("CA");
	auto const &[cax, cay, caz] = residue.ca_location();

	return pdbx::format("%5d%5d%1.1s%1.1s %c  %c%c%c%c%c%c%c%c%c%4d%4d%c%4.0f %11s%11s%11s%11s  %6.3f%6.1f%6.1f%6.1f%6.1f %6.1f %6.1f %6.1f",
		info.nr(), residue.pdb_seq_num(), residue.pdb_ins_code(), residue.pdb_strand_id(), code,
			ss, helix[3], helix[0], helix[1], helix[2], bend, chirality, bridgelabel[0], bridgelabel[1],
			bp[0], bp[1], sheet, floor(info.accessibility() + 0.5),
			NHO[0], ONH[0], NHO[1], ONH[1],
			residue.tco(), residue.kappa(), alpha, residue.phi(), residue.psi(),
			cax, cay, caz)
	    .str();
}

void writeDSSP(const dssp::DSSP &dssp, std::ostream &os)
{
	using namespace date;
	using namespace std::chrono;

	auto stats = dssp.get_statistics();

	auto today = system_clock::now();

	os << "==== Secondary Structure Definition by the program DSSP, NKI version 4.0                           ==== DATE=" << format("%F", today) << "        ." << std::endl
	   << "REFERENCE W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637                                                              ." << std::endl
	   << dssp.get_pdb_header_line(dssp::DSSP::pdb_record_type::HEADER) << '.' << std::endl
	   << dssp.get_pdb_header_line(dssp::DSSP::pdb_record_type::COMPND) << '.' << std::endl
	   << dssp.get_pdb_header_line(dssp::DSSP::pdb_record_type::SOURCE) << '.' << std::endl
	   << dssp.get_pdb_header_line(dssp::DSSP::pdb_record_type::AUTHOR) << '.' << std::endl;

	os << pdbx::format("%5d%3d%3d%3d%3d TOTAL NUMBER OF RESIDUES, NUMBER OF CHAINS, NUMBER OF SS-BRIDGES(TOTAL,INTRACHAIN,INTERCHAIN)                .",
			  stats.count.residues, stats.count.chains, stats.count.SS_bridges, stats.count.intra_chain_SS_bridges, (stats.count.SS_bridges - stats.count.intra_chain_SS_bridges))
	   << std::endl;

	os << pdbx::format("%8.1f   ACCESSIBLE SURFACE OF PROTEIN (ANGSTROM**2)                                                                         .", stats.accessible_surface) << std::endl;

	// hydrogenbond summary

	os << pdbx::format("%5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(J)  , SAME NUMBER PER 100 RESIDUES                              .", stats.count.H_bonds, (stats.count.H_bonds * 100.0 / stats.count.residues)) << std::endl;

	os << pdbx::format("%5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN     PARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .", stats.count.H_bonds_in_parallel_bridges, (stats.count.H_bonds_in_parallel_bridges * 100.0 / stats.count.residues)) << std::endl;

	os << pdbx::format("%5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS IN ANTIPARALLEL BRIDGES, SAME NUMBER PER 100 RESIDUES                              .", stats.count.H_bonds_in_antiparallel_bridges, (stats.count.H_bonds_in_antiparallel_bridges * 100.0 / stats.count.residues)) << std::endl;

	for (int k = 0; k < 11; ++k)
		os << pdbx::format("%5d%5.1f   TOTAL NUMBER OF HYDROGEN BONDS OF TYPE O(I)-->H-N(I%c%1d), SAME NUMBER PER 100 RESIDUES                              .", stats.count.H_Bonds_per_distance[k], (stats.count.H_Bonds_per_distance[k] * 100.0 / stats.count.residues), (k - 5 < 0 ? '-' : '+'), abs(k - 5)) << std::endl;

	// histograms...
	os << "  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30     *** HISTOGRAMS OF ***           ." << std::endl;

	for (auto hi : stats.histogram.residues_per_alpha_helix)
		os << pdbx::format("%3d", hi);
	os << "    RESIDUES PER ALPHA HELIX         ." << std::endl;

	for (auto hi : stats.histogram.parallel_bridges_per_ladder)
		os << pdbx::format("%3d", hi);
	os << "    PARALLEL BRIDGES PER LADDER      ." << std::endl;

	for (auto hi : stats.histogram.antiparallel_bridges_per_ladder)
		os << pdbx::format("%3d", hi);
	os << "    ANTIPARALLEL BRIDGES PER LADDER  ." << std::endl;

	for (auto hi : stats.histogram.ladders_per_sheet)
		os << pdbx::format("%3d", hi);
	os << "    LADDERS PER SHEET                ." << std::endl;

	// per residue information

	os << "  #  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA" << std::endl;

	int last = 0;
	for (auto ri : dssp)
	{
		// insert a break line whenever we detect missing residues
		// can be the transition to a different chain, or missing residues in the current chain

		if (ri.nr() != last + 1)
			os << pdbx::format("%5d        !%c             0   0    0      0, 0.0     0, 0.0     0, 0.0     0, 0.0   0.000 360.0 360.0 360.0 360.0    0.0    0.0    0.0",
				(last + 1), (ri.chain_break() == dssp::chain_break_type::NewChain ? '*' : ' ')) << std::endl;

		os << ResidueToDSSPLine(ri) << std::endl;
		last = ri.nr();
	}
}

void annotateDSSP(cif::datablock &db, const dssp::DSSP &dssp, bool writeOther, std::ostream &os)
{
	if (dssp.empty())
	{
		if (cif::VERBOSE)
			std::cout << "No secondary structure information found" << std::endl;
	}
	else
	{
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

	// db.add_software("dssp", "other", kVersionNumber, kBuildDate);
	// db["software"].emplace({

	// });

	db.write(os);
}
