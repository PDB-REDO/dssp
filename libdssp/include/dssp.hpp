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

#pragma once

/// \file DSSP.hpp
/// Calculate DSSP-like secondary structure information.

#include <cif++.hpp>

#include <filesystem>
#include <iterator>
#include <ranges>

class dssp
{
  public:
	struct residue;

	enum class structure_type : char
	{
		Loop = ' ',
		Alphahelix = 'H',
		Betabridge = 'B',
		Strand = 'E',
		Helix_3 = 'G',
		Helix_5 = 'I',
		Helix_PPII = 'P',
		Turn = 'T',
		Bend = 'S'
	};

	enum class helix_type
	{
		_3_10,
		alpha,
		pi,
		pp
	};

	enum class helix_position_type
	{
		None,
		Start,
		End,
		StartAndEnd,
		Middle
	};

	static constexpr size_t kHistogramSize = 30;

	struct statistics
	{
		struct
		{
			uint32_t residues, chains, SS_bridges, intra_chain_SS_bridges, H_bonds;
			uint32_t H_bonds_in_antiparallel_bridges, H_bonds_in_parallel_bridges;
			uint32_t H_Bonds_per_distance[11];
		} count;

		double accessible_surface;

		struct
		{
			uint32_t residues_per_alpha_helix[kHistogramSize];
			uint32_t parallel_bridges_per_ladder[kHistogramSize];
			uint32_t antiparallel_bridges_per_ladder[kHistogramSize];
			uint32_t ladders_per_sheet[kHistogramSize];
		} histogram;
	};

	enum class chain_break_type
	{
		None,
		NewChain,
		Gap
	};

	dssp(const cif::datablock &db, int model_nr, int min_poly_proline_stretch_length, bool calculateSurfaceAccessibility);
	dssp(const cif::mm::structure &s, int min_poly_proline_stretch_length, bool calculateSurfaceAccessibility);

	~dssp();

	dssp(const dssp &) = delete;
	dssp &operator=(const dssp &) = delete;

	[[nodiscard]] statistics get_statistics() const;

	class iterator;
	using res_iterator = typename std::vector<residue>::iterator;

	class residue_info
	{
	  public:
		friend class iterator;

		residue_info() = default;
		residue_info(const residue_info &rhs) = default;
		residue_info &operator=(const residue_info &rhs) = default;

		explicit operator bool() const { return not empty(); }
		[[nodiscard]] bool empty() const { return m_impl == nullptr; }

		[[nodiscard]] std::string asym_id() const;
		[[nodiscard]] int seq_id() const;
		[[nodiscard]] std::string alt_id() const;
		[[nodiscard]] std::string compound_id() const;
		[[nodiscard]] char compound_letter() const; // Single letter for residue compound type, or 'X' in case it is not known

		[[nodiscard]] std::string auth_asym_id() const;
		[[nodiscard]] int auth_seq_id() const;

		[[nodiscard]] std::string pdb_strand_id() const;
		[[nodiscard]] int pdb_seq_num() const;
		[[nodiscard]] std::string pdb_ins_code() const;

		[[nodiscard]] std::optional<float> alpha() const;
		[[nodiscard]] std::optional<float> kappa() const;
		[[nodiscard]] std::optional<float> phi() const;
		[[nodiscard]] std::optional<float> psi() const;
		[[nodiscard]] std::optional<float> tco() const;
		[[nodiscard]] std::optional<float> omega() const;

		[[nodiscard]] bool is_pre_pro() const;
		[[nodiscard]] bool is_cis() const { return std::abs(omega().value_or(360)) < 30.0f; }

		[[nodiscard]] float chiral_volume() const;

		[[nodiscard]] std::size_t nr_of_chis() const;
		[[nodiscard]] float chi(std::size_t index) const;

		[[nodiscard]] std::vector<float> chis() const
		{
			std::vector<float> result;
			for (size_t i = 0; i < nr_of_chis(); ++i)
				result.push_back(chi(i));
			return result;
		}

		[[nodiscard]] std::tuple<float, float, float> ca_location() const;

		[[nodiscard]] chain_break_type chain_break() const;

		/// \brief the internal number in DSSP
		[[nodiscard]] int nr() const;

		[[nodiscard]] structure_type type() const;

		[[nodiscard]] int ssBridgeNr() const;

		[[nodiscard]] helix_position_type helix(helix_type helixType) const;

		[[nodiscard]] bool is_alpha_helix_end_before_start() const;

		[[nodiscard]] bool bend() const;

		[[nodiscard]] double accessibility() const;

		/// \brief returns resinfo, ladder and parallel
		[[nodiscard]] std::tuple<residue_info, int, bool> bridge_partner(int i) const;

		[[nodiscard]] int sheet() const;
		[[nodiscard]] int strand() const;

		/// \brief return resinfo and the energy of the bond
		[[nodiscard]] std::tuple<residue_info, double> acceptor(int i) const;
		[[nodiscard]] std::tuple<residue_info, double> donor(int i) const;

		/// \brief Simple compare equals
		bool operator==(const residue_info &rhs) const
		{
			return m_impl == rhs.m_impl;
		}

		/// \brief Returns \result true if there is a bond between two residues
		friend bool test_bond(residue_info const &a, residue_info const &b);

		[[nodiscard]] residue_info next() const;

	  private:
		residue_info(residue *res)
			: m_impl(res)
		{
		}

		residue *m_impl = nullptr;
	};

	class iterator
	{
	  public:
		using iterator_category = std::bidirectional_iterator_tag;
		using value_type = const residue_info;
		using difference_type = std::ptrdiff_t;
		using pointer = value_type *;
		using reference = value_type &;

		iterator() = default;
		iterator(const iterator &i) = default;
		iterator(residue *res);
		iterator &operator=(const iterator &i) = default;

		reference operator*() const { return m_current; }
		pointer operator->() const { return &m_current; }

		iterator &operator++();
		iterator operator++(int)
		{
			auto tmp(*this);
			this->operator++();
			return tmp;
		}

		iterator &operator--();
		iterator operator--(int)
		{
			auto tmp(*this);
			this->operator--();
			return tmp;
		}

		bool operator==(const iterator &rhs) const { return m_current.m_impl == rhs.m_current.m_impl; }
		bool operator!=(const iterator &rhs) const { return m_current.m_impl != rhs.m_current.m_impl; }

	  private:
		residue_info m_current;
	};

	static_assert(std::input_iterator<iterator>);

	using value_type = residue_info;

	// To access residue info by key, i.e. LabelAsymID and LabelSeqID
	using key_type = std::tuple<std::string, int>;

	[[nodiscard]] iterator begin() const;
	[[nodiscard]] iterator end() const;

	[[nodiscard]] iterator cbegin() const;
	[[nodiscard]] iterator cend() const;

	residue_info operator[](const key_type &key) const;

	[[nodiscard]] bool empty() const { return begin() == end(); }

	// --------------------------------------------------------------------
	// Writing out the data, either in legacy format...

	void write_legacy_output(std::ostream &os) const;

	// ... or as annotation in the cif::datablock
	void annotate(cif::datablock &db, bool writeOther, bool writeDSSPCategories) const;

	// convenience method, when creating old style DSSP files

	enum class pdb_record_type
	{
		HEADER,
		COMPND,
		SOURCE,
		AUTHOR
	};

	[[nodiscard]] std::string get_pdb_header_line(pdb_record_type pdb_record) const;

  private:
	struct DSSP_impl *m_impl;
};

static_assert(std::ranges::input_range<dssp>);