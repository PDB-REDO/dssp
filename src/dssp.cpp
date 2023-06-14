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

// Calculate DSSP-like secondary structure information

#include "dssp.hpp"
#include "queue.hpp"

#include "dssp-io.hpp"

#include <deque>
#include <iomanip>
#include <numeric>
#include <thread>

using residue = dssp::residue;
using statistics = dssp::statistics;
using structure_type = dssp::structure_type;
using helix_type = dssp::helix_type;
using helix_position_type = dssp::helix_position_type;
using chain_break_type = dssp::chain_break_type;

using queue_type = blocking_queue<std::tuple<uint32_t,uint32_t>>;

// --------------------------------------------------------------------

const double
	kPI = 3.141592653589793238462643383279502884;

struct point
{
	float mX, mY, mZ;
};

inline point operator-(const point &lhs, const point &rhs)
{
	return { lhs.mX - rhs.mX, lhs.mY - rhs.mY, lhs.mZ - rhs.mZ };
}

inline point operator*(const point &lhs, float rhs)
{
	return { lhs.mX * rhs, lhs.mY * rhs, lhs.mZ * rhs };
}

inline float distance_sq(const point &a, const point &b)
{
	return (a.mX - b.mX) * (a.mX - b.mX) +
	       (a.mY - b.mY) * (a.mY - b.mY) +
	       (a.mZ - b.mZ) * (a.mZ - b.mZ);
}

inline float distance(const point &a, const point &b)
{
	return std::sqrt(distance_sq(a, b));
}

inline float dot_product(const point &a, const point &b)
{
	return a.mX * b.mX + a.mY * b.mY + a.mZ * b.mZ;
}

inline point cross_product(const point &a, const point &b)
{
	return {
		a.mY * b.mZ - b.mY * a.mZ,
		a.mZ * b.mX - b.mZ * a.mX,
		a.mX * b.mY - b.mX * a.mY
	};
}

float dihedral_angle(const point &p1, const point &p2, const point &p3, const point &p4)
{
	point v12 = p1 - p2; // vector from p2 to p1
	point v43 = p4 - p3; // vector from p3 to p4

	point z = p2 - p3; // vector from p3 to p2

	point p = cross_product(z, v12);
	point x = cross_product(z, v43);
	point y = cross_product(z, x);

	float u = dot_product(x, x);
	float v = dot_product(y, y);

	float result = 360;
	if (u > 0 and v > 0)
	{
		u = dot_product(p, x) / std::sqrt(u);
		v = dot_product(p, y) / std::sqrt(v);
		if (u != 0 or v != 0)
			result = std::atan2(v, u) * 180.f / static_cast<float>(kPI);
	}

	return result;
}

float cosinus_angle(const point &p1, const point &p2, const point &p3, const point &p4)
{
	point v12 = p1 - p2;
	point v34 = p3 - p4;

	float result = 0;

	float x = dot_product(v12, v12) * dot_product(v34, v34);
	if (x > 0)
		result = dot_product(v12, v34) / std::sqrt(x);

	return result;
}

enum residue_type : char
{
	kUnknownResidue = 'X',

	//
	kAlanine = 'A',       //	ala
	kArginine = 'R',      //	arg
	kAsparagine = 'N',    //	asn
	kAsparticAcid = 'D',  //	asp
	kCysteine = 'C',      //	cys
	kGlutamicAcid = 'E',  //	glu
	kGlutamine = 'Q',     //	gln
	kGlycine = 'G',       //	gly
	kHistidine = 'H',     //	his
	kIsoleucine = 'I',    //	ile
	kLeucine = 'L',       //	leu
	kLysine = 'K',        //	lys
	kMethionine = 'M',    //	met
	kPhenylalanine = 'F', //	phe
	kProline = 'P',       //	pro
	kSerine = 'S',        //	ser
	kThreonine = 'T',     //	thr
	kTryptophan = 'W',    //	trp
	kTyrosine = 'Y',      //	tyr
	kValine = 'V',        //	val
};

struct
{
	residue_type type;
	char name[4];
} const kResidueInfo[] = {
	{ kUnknownResidue, "UNK" },
	{ kAlanine, "ALA" },
	{ kArginine, "ARG" },
	{ kAsparagine, "ASN" },
	{ kAsparticAcid, "ASP" },
	{ kCysteine, "CYS" },
	{ kGlutamicAcid, "GLU" },
	{ kGlutamine, "GLN" },
	{ kGlycine, "GLY" },
	{ kHistidine, "HIS" },
	{ kIsoleucine, "ILE" },
	{ kLeucine, "LEU" },
	{ kLysine, "LYS" },
	{ kMethionine, "MET" },
	{ kPhenylalanine, "PHE" },
	{ kProline, "PRO" },
	{ kSerine, "SER" },
	{ kThreonine, "THR" },
	{ kTryptophan, "TRP" },
	{ kTyrosine, "TYR" },
	{ kValine, "VAL" }
};

constexpr residue_type MapResidue(std::string_view inName)
{
	residue_type result = kUnknownResidue;

	for (auto &ri : kResidueInfo)
	{
		if (inName == ri.name)
		{
			result = ri.type;
			break;
		}
	}

	return result;
}

struct HBond
{
	residue *res;
	double energy;
};

enum class bridge_type
{
	None,
	Parallel,
	AntiParallel
};

struct bridge
{
	bridge_type type;
	uint32_t sheet, ladder;
	std::set<bridge *> link;
	std::deque<uint32_t> i, j;
	std::string chainI, chainJ;

	bool operator<(const bridge &b) const { return chainI < b.chainI or (chainI == b.chainI and i.front() < b.i.front()); }
};

struct bridge_partner
{
	residue *m_residue;
	uint32_t ladder;
	bool parallel;
};

// --------------------------------------------------------------------

const float
	// kSSBridgeDistance = 3.0f,
	kMinimalDistance = 0.5f,
	kMinimalCADistance = 9.0f,
	kMinHBondEnergy = -9.9f,
	kMaxHBondEnergy = -0.5f,
	kCouplingConstant = -27.888f, //	= -332 * 0.42 * 0.2
	kMaxPeptideBondLength = 2.5f;

const float
	kRadiusN = 1.65f,
	kRadiusCA = 1.87f,
	kRadiusC = 1.76f,
	kRadiusO = 1.4f,
	kRadiusSideAtom = 1.8f,
	kRadiusWater = 1.4f;

struct dssp::residue
{
	residue(residue &&) = default;
	residue &operator=(residue &&) = default;

	residue(int model_nr,
		std::string_view pdb_strand_id, int pdb_seq_num, std::string_view pdb_ins_code)
		: mPDBStrandID(pdb_strand_id)
		, mPDBSeqNum(pdb_seq_num)
		, mPDBInsCode(pdb_ins_code)
		, mChainBreak(chain_break_type::None)
		, m_model_nr(model_nr)
	{
		// update the box containing all atoms
		mBox[0].mX = mBox[0].mY = mBox[0].mZ = std::numeric_limits<float>::max();
		mBox[1].mX = mBox[1].mY = mBox[1].mZ = -std::numeric_limits<float>::max();

		mH = point{ std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max() };

		for (auto &p : m_chiralAtoms)
			p = {};
	}

	void addAtom(cif::row_handle atom)
	{
		std::string asymID, compID, atomID, type, authAsymID;
		std::optional<std::string> altID;
		int seqID, authSeqID;
		std::optional<int> model;
		float x, y, z;

		cif::tie(asymID, compID, atomID, altID, type, seqID, model, x, y, z, authAsymID, authSeqID) =
			atom.get("label_asym_id", "label_comp_id", "label_atom_id", "label_alt_id", "type_symbol", "label_seq_id",
				"pdbx_PDB_model_num", "Cartn_x", "Cartn_y", "Cartn_z",
				"auth_asym_id", "auth_seq_id");

		if (model and model != m_model_nr)
			return;

		if (m_seen == 0)
		{
			mAsymID = asymID;
			mCompoundID = compID;
			mSeqID = seqID;

			mAuthSeqID = authSeqID;
			mAuthAsymID = authAsymID;

			mType = MapResidue(mCompoundID);

			if (altID)
				mAltID = *altID;
		}

		if (atomID == "CA")
		{
			m_seen |= 1;
			mCAlpha = { x, y, z };
			ExtendBox(mCAlpha, kRadiusCA + 2 * kRadiusWater);

			if (mType == kValine)
				m_chiralAtoms[1] = mCAlpha;
		}
		else if (atomID == "C")
		{
			m_seen |= 2;
			mC = { x, y, z };
			ExtendBox(mC, kRadiusC + 2 * kRadiusWater);
		}
		else if (atomID == "N")
		{
			m_seen |= 4;
			mH = mN = { x, y, z };
			ExtendBox(mN, kRadiusN + 2 * kRadiusWater);
		}
		else if (atomID == "O")
		{
			m_seen |= 8;
			mO = { x, y, z };
			ExtendBox(mO, kRadiusO + 2 * kRadiusWater);
		}
		else if (type != "H")
		{
			m_seen |= 16;
			mSideChain.emplace_back(atomID, point{ x, y, z });
			ExtendBox(point{ x, y, z }, kRadiusSideAtom + 2 * kRadiusWater);

			if (mType == kLeucine)
			{
				if (atomID == "CG")
					m_chiralAtoms[0] = { x, y, z };
				else if (atomID == "CB")
					m_chiralAtoms[1] = { x, y, z };
				else if (atomID == "CD1")
					m_chiralAtoms[2] = { x, y, z };
				else if (atomID == "CD2")
					m_chiralAtoms[3] = { x, y, z };
			}
			else if (mType == kValine)
			{
				if (atomID == "CB")
					m_chiralAtoms[0] = { x, y, z };
				else if (atomID == "CG1")
					m_chiralAtoms[2] = { x, y, z };
				else if (atomID == "CG2")
					m_chiralAtoms[3] = { x, y, z };
			}
		}
	}

	void finish()
	{
		const int kSeenAll = (1 bitor 2 bitor 4 bitor 8);
		mComplete = (m_seen bitand kSeenAll) == kSeenAll;

		if (mType == kValine or mType == kLeucine)
			mChiralVolume = dot_product(m_chiralAtoms[1] - m_chiralAtoms[0],
				cross_product(m_chiralAtoms[2] - m_chiralAtoms[0], m_chiralAtoms[3] - m_chiralAtoms[0]));

		mRadius = mBox[1].mX - mBox[0].mX;
		if (mRadius < mBox[1].mY - mBox[0].mY)
			mRadius = mBox[1].mY - mBox[0].mY;
		if (mRadius < mBox[1].mZ - mBox[0].mZ)
			mRadius = mBox[1].mZ - mBox[0].mZ;

		mCenter.mX = (mBox[0].mX + mBox[1].mX) / 2;
		mCenter.mY = (mBox[0].mY + mBox[1].mY) / 2;
		mCenter.mZ = (mBox[0].mZ + mBox[1].mZ) / 2;
	}

	void assignHydrogen()
	{
		// assign the Hydrogen
		mH = mN;

		if (mType != kProline and mPrev != nullptr)
		{
			auto pc = mPrev->mC;
			auto po = mPrev->mO;

			float CODistance = static_cast<float>(distance(pc, po));

			mH.mX += (pc.mX - po.mX) / CODistance;
			mH.mY += (pc.mY - po.mY) / CODistance;
			mH.mZ += (pc.mZ - po.mZ) / CODistance;
		}
	}

	void SetSecondaryStructure(structure_type inSS) { mSecondaryStructure = inSS; }
	structure_type GetSecondaryStructure() const { return mSecondaryStructure; }

	void SetBetaPartner(uint32_t n, residue &inResidue, uint32_t inLadder, bool inParallel)
	{
		assert(n == 0 or n == 1);

		mBetaPartner[n].m_residue = &inResidue;
		mBetaPartner[n].ladder = inLadder;
		mBetaPartner[n].parallel = inParallel;
	}

	bridge_partner GetBetaPartner(uint32_t n) const
	{
		assert(n == 0 or n == 1);
		return mBetaPartner[n];
	}

	void SetSheet(uint32_t inSheet) { mSheet = inSheet; }
	uint32_t GetSheet() const { return mSheet; }

	void SetStrand(uint32_t inStrand) { mStrand = inStrand; }
	uint32_t GetStrand() const { return mStrand; }

	bool IsBend() const { return mBend; }
	void SetBend(bool inBend) { mBend = inBend; }

	helix_position_type GetHelixFlag(helix_type helixType) const
	{
		size_t stride = static_cast<size_t>(helixType);
		assert(stride < 4);
		return mHelixFlags[stride];
	}

	bool IsHelixStart(helix_type helixType) const
	{
		size_t stride = static_cast<size_t>(helixType);
		assert(stride < 4);
		return mHelixFlags[stride] == helix_position_type::Start or mHelixFlags[stride] == helix_position_type::StartAndEnd;
	}

	void SetHelixFlag(helix_type helixType, helix_position_type inHelixFlag)
	{
		size_t stride = static_cast<size_t>(helixType);
		assert(stride < 4);
		mHelixFlags[stride] = inHelixFlag;
	}

	void SetSSBridgeNr(uint8_t inBridgeNr)
	{
		if (mType != kCysteine)
			throw std::runtime_error("Only cysteine residues can form sulphur bridges");
		mSSBridgeNr = inBridgeNr;
	}

	uint8_t GetSSBridgeNr() const
	{
		if (mType != kCysteine)
			throw std::runtime_error("Only cysteine residues can form sulphur bridges");
		return mSSBridgeNr;
	}

	float CalculateSurface(const std::vector<residue> &inResidues);
	float CalculateSurface(const point &inAtom, float inRadius, const std::vector<residue *> &inNeighbours);

	bool AtomIntersectsBox(const point &atom, float inRadius) const
	{
		return atom.mX + inRadius >= mBox[0].mX and
		       atom.mX - inRadius <= mBox[1].mX and
		       atom.mY + inRadius >= mBox[0].mY and
		       atom.mY - inRadius <= mBox[1].mY and
		       atom.mZ + inRadius >= mBox[0].mZ and
		       atom.mZ - inRadius <= mBox[1].mZ;
	}

	void ExtendBox(const point &atom, float inRadius)
	{
		if (mBox[0].mX > atom.mX - inRadius)
			mBox[0].mX = atom.mX - inRadius;
		if (mBox[0].mY > atom.mY - inRadius)
			mBox[0].mY = atom.mY - inRadius;
		if (mBox[0].mZ > atom.mZ - inRadius)
			mBox[0].mZ = atom.mZ - inRadius;
		if (mBox[1].mX < atom.mX + inRadius)
			mBox[1].mX = atom.mX + inRadius;
		if (mBox[1].mY < atom.mY + inRadius)
			mBox[1].mY = atom.mY + inRadius;
		if (mBox[1].mZ < atom.mZ + inRadius)
			mBox[1].mZ = atom.mZ + inRadius;
	}

	point get_atom(std::string_view name)
	{
		if (name == "CA")
			return mCAlpha;
		else if (name == "C")
			return mC;
		else if (name == "N")
			return mN;
		else if (name == "O")
			return mO;
		else if (name == "H")
			return mH;
		else
		{
			for (const auto &[n, p] : mSideChain)
			{
				if (n == name)
					return p;
			}
		}

		return {};
	}

	residue *mNext = nullptr;
	residue *mPrev = nullptr;

	// const Monomer &mM;

	std::string mAsymID;
	int mSeqID;
	std::string mCompoundID;
	std::string mAltID;
	int mNumber;
	bool mComplete = false;

	std::string mAuthAsymID;
	int mAuthSeqID;

	std::string mPDBStrandID;
	int mPDBSeqNum;
	std::string mPDBInsCode;

	point mCAlpha, mC, mN, mO, mH;
	point mBox[2] = {};
	float mRadius;
	point mCenter;
	std::vector<std::tuple<std::string, point>> mSideChain;
	float mAccessibility = 0;
	float mChiralVolume = 0;

	// float mAlpha = 360, mKappa = 360, mPhi = 360, mPsi = 360, mTCO = 0, mOmega = 360;
	std::optional<float> mAlpha, mKappa, mPhi, mPsi, mTCO, mOmega;

	residue_type mType;
	uint8_t mSSBridgeNr = 0;
	structure_type mSecondaryStructure = structure_type::Loop;
	HBond mHBondDonor[2] = {}, mHBondAcceptor[2] = {};
	bridge_partner mBetaPartner[2] = {};
	uint32_t mSheet = 0;
	uint32_t mStrand = 0;	// Added to ease the writing of mmCIF's struct_sheet and friends
	helix_position_type mHelixFlags[4] = { helix_position_type::None, helix_position_type::None, helix_position_type::None, helix_position_type::None }; //
	bool mBend = false;
	chain_break_type mChainBreak = chain_break_type::None;

	int m_seen = 0, m_model_nr = 0;
	std::array<point, 4> m_chiralAtoms;
};

// --------------------------------------------------------------------

class accumulator
{
  public:
	struct candidate
	{
		point location;
		double radius;
		double distance;

		bool operator<(const candidate &rhs) const
		{
			return distance < rhs.distance;
		}
	};

	void operator()(const point &a, const point &b, double d, double r)
	{
		double distance = distance_sq(a, b);

		d += kRadiusWater;
		r += kRadiusWater;

		double test = d + r;
		test *= test;

		if (distance < test and distance > 0.0001)
		{
			candidate c = { b - a, r * r, distance };

			m_x.push_back(c);
			push_heap(m_x.begin(), m_x.end());
		}
	}

	void sort()
	{
		sort_heap(m_x.begin(), m_x.end());
	}

	std::vector<candidate> m_x;
};

// we use a fibonacci sphere to calculate the even distribution of the dots
class MSurfaceDots
{
  public:
	static MSurfaceDots &Instance();

	size_t size() const { return mPoints.size(); }
	const point &operator[](size_t inIx) const { return mPoints[inIx]; }
	double weight() const { return mWeight; }

  private:
	MSurfaceDots(int32_t inN);

	std::vector<point> mPoints;
	double mWeight;
};

MSurfaceDots &MSurfaceDots::Instance()
{
	const int32_t kN = 200;

	static MSurfaceDots sInstance(kN);
	return sInstance;
}

MSurfaceDots::MSurfaceDots(int32_t N)
{
	auto P = 2 * N + 1;

	const float kGoldenRatio = (1 + std::sqrt(5.0f)) / 2;

	mWeight = (4 * kPI) / P;

	for (auto i = -N; i <= N; ++i)
	{
		float lat = std::asin((2.0f * i) / P);
		float lon = static_cast<float>(std::fmod(i, kGoldenRatio) * 2 * kPI / kGoldenRatio);

		mPoints.emplace_back(point{ std::sin(lon) * std::cos(lat), std::cos(lon) * std::cos(lat), std::sin(lat) });
	}
}

float residue::CalculateSurface(const point &inAtom, float inRadius, const std::vector<residue *> &inNeighbours)
{
	accumulator accumulate;

	for (auto r : inNeighbours)
	{
		if (r->AtomIntersectsBox(inAtom, inRadius))
		{
			accumulate(inAtom, r->mN, inRadius, kRadiusN);
			accumulate(inAtom, r->mCAlpha, inRadius, kRadiusCA);
			accumulate(inAtom, r->mC, inRadius, kRadiusC);
			accumulate(inAtom, r->mO, inRadius, kRadiusO);

			for (const auto &[name, atom] : r->mSideChain)
				accumulate(inAtom, atom, inRadius, kRadiusSideAtom);
		}
	}

	accumulate.sort();

	float radius = inRadius + kRadiusWater;
	float surface = 0;

	MSurfaceDots &surfaceDots = MSurfaceDots::Instance();

	for (size_t i = 0; i < surfaceDots.size(); ++i)
	{
		point xx = surfaceDots[i] * radius;

		bool free = true;
		for (size_t k = 0; free and k < accumulate.m_x.size(); ++k)
			free = accumulate.m_x[k].radius < distance_sq(xx, accumulate.m_x[k].location);

		if (free)
			surface += static_cast<float>(surfaceDots.weight());
	}

	return surface * radius * radius;
}

float residue::CalculateSurface(const std::vector<residue> &inResidues)
{
	std::vector<residue *> neighbours;

	for (auto &r : inResidues)
	{
		point center = r.mCenter;
		float radius = r.mRadius;

		if (distance_sq(mCenter, center) < (mRadius + radius) * (mRadius + radius))
			neighbours.push_back(const_cast<residue *>(&r));
	}

	mAccessibility = CalculateSurface(mN, kRadiusN, neighbours) +
	                 CalculateSurface(mCAlpha, kRadiusCA, neighbours) +
	                 CalculateSurface(mC, kRadiusC, neighbours) +
	                 CalculateSurface(mO, kRadiusO, neighbours);

	for (const auto &[name, atom] : mSideChain)
		mAccessibility += CalculateSurface(atom, kRadiusSideAtom, neighbours);

	return mAccessibility;
}

void CalculateAccessibilities(std::vector<residue> &inResidues, statistics &stats)
{
	stats.accessible_surface = 0;
	for (auto &residue : inResidues)
		stats.accessible_surface += residue.CalculateSurface(inResidues);
}

// --------------------------------------------------------------------
// TODO: use the angle to improve bond energy calculation.

double CalculateHBondEnergy(residue &inDonor, residue &inAcceptor)
{
	double result = 0;

	if (inDonor.mType != kProline)
	{
		double distanceHO = distance(inDonor.mH, inAcceptor.mO);
		double distanceHC = distance(inDonor.mH, inAcceptor.mC);
		double distanceNC = distance(inDonor.mN, inAcceptor.mC);
		double distanceNO = distance(inDonor.mN, inAcceptor.mO);

		if (distanceHO < kMinimalDistance or distanceHC < kMinimalDistance or distanceNC < kMinimalDistance or distanceNO < kMinimalDistance)
			result = kMinHBondEnergy;
		else
			result = kCouplingConstant / distanceHO - kCouplingConstant / distanceHC + kCouplingConstant / distanceNC - kCouplingConstant / distanceNO;

		// DSSP compatibility mode:
		result = std::round(result * 1000) / 1000;

		if (result < kMinHBondEnergy)
			result = kMinHBondEnergy;
	}

	// update donor
	if (result < inDonor.mHBondAcceptor[0].energy)
	{
		inDonor.mHBondAcceptor[1] = inDonor.mHBondAcceptor[0];
		inDonor.mHBondAcceptor[0].res = &inAcceptor;
		inDonor.mHBondAcceptor[0].energy = result;
	}
	else if (result < inDonor.mHBondAcceptor[1].energy)
	{
		inDonor.mHBondAcceptor[1].res = &inAcceptor;
		inDonor.mHBondAcceptor[1].energy = result;
	}

	// and acceptor
	if (result < inAcceptor.mHBondDonor[0].energy)
	{
		inAcceptor.mHBondDonor[1] = inAcceptor.mHBondDonor[0];
		inAcceptor.mHBondDonor[0].res = &inDonor;
		inAcceptor.mHBondDonor[0].energy = result;
	}
	else if (result < inAcceptor.mHBondDonor[1].energy)
	{
		inAcceptor.mHBondDonor[1].res = &inDonor;
		inAcceptor.mHBondDonor[1].energy = result;
	}

	return result;
}

// --------------------------------------------------------------------

void CalculateHBondEnergies(std::vector<residue> &inResidues, std::vector<std::tuple<uint32_t, uint32_t>> &q)
{
	std::unique_ptr<cif::progress_bar> progress;
	if (cif::VERBOSE == 0 or cif::VERBOSE == 1)
		progress.reset(new cif::progress_bar(q.size(), "calculate hbond energies"));

	for (const auto &[i, j] : q)
	{
		auto &ri = inResidues[i];
		auto &rj = inResidues[j];

		CalculateHBondEnergy(ri, rj);
		if (j != i + 1)
			CalculateHBondEnergy(rj, ri);
		
		if (progress)
			progress->consumed(1);
	}
}

// --------------------------------------------------------------------

bool NoChainBreak(const residue *a, const residue *b)
{
	bool result = a->mAsymID == b->mAsymID;
	for (auto r = a; result and r != b; r = r->mNext)
	{
		auto next = r->mNext;
		if (next == nullptr)
			result = false;
		else
			result = next->mNumber == r->mNumber + 1;
	}
	return result;
}

bool NoChainBreak(const residue &a, const residue &b)
{
	return NoChainBreak(&a, &b);
}

// --------------------------------------------------------------------

bool TestBond(const residue *a, const residue *b)
{
	return (a->mHBondAcceptor[0].res == b and a->mHBondAcceptor[0].energy < kMaxHBondEnergy) or
	       (a->mHBondAcceptor[1].res == b and a->mHBondAcceptor[1].energy < kMaxHBondEnergy);
}

bool test_bond(dssp::residue_info const &a, dssp::residue_info const &b)
{
	return a and b and TestBond(a.m_impl, b.m_impl);
}

// --------------------------------------------------------------------

bridge_type TestBridge(const residue &r1, const residue &r2)
{                                           // I.	a	d	II.	a	d		parallel
	auto a = r1.mPrev;                      //		  \			  /
	auto b = &r1;                           //		b	e		b	e
	auto c = r1.mNext;                      // 		  /			  \                      ..
	auto d = r2.mPrev;                      //		c	f		c	f
	auto e = &r2;                           //
	auto f = r2.mNext;                      // III.	a <- f	IV. a	  f		antiparallel
	                                        //
	bridge_type result = bridge_type::None; //		b	 e      b <-> e
	                                        //
	                                        //		c -> d		c     d

	if (a and c and NoChainBreak(a, c) and d and f and NoChainBreak(d, f))
	{
		if ((TestBond(c, e) and TestBond(e, a)) or (TestBond(f, b) and TestBond(b, d)))
			result = bridge_type::Parallel;
		else if ((TestBond(c, d) and TestBond(f, a)) or (TestBond(e, b) and TestBond(b, e)))
			result = bridge_type::AntiParallel;
	}

	return result;
}

// --------------------------------------------------------------------
// return true if any of the residues in bridge a is identical to any of the residues in bridge b
bool Linked(const bridge &a, const bridge &b)
{
	return find_first_of(a.i.begin(), a.i.end(), b.i.begin(), b.i.end()) != a.i.end() or
	       find_first_of(a.i.begin(), a.i.end(), b.j.begin(), b.j.end()) != a.i.end() or
	       find_first_of(a.j.begin(), a.j.end(), b.i.begin(), b.i.end()) != a.j.end() or
	       find_first_of(a.j.begin(), a.j.end(), b.j.begin(), b.j.end()) != a.j.end();
}

// --------------------------------------------------------------------

void CalculateBetaSheets(std::vector<residue> &inResidues, statistics &stats, std::vector<std::tuple<uint32_t, uint32_t>> &q)
{
	// if (cif::VERBOSE)
	// 	std::cerr << "calculating beta sheets" << std::endl;

	std::unique_ptr<cif::progress_bar> progress;
	if (cif::VERBOSE == 0 or cif::VERBOSE == 1)
		progress.reset(new cif::progress_bar(q.size(), "calculate beta sheets"));

	// Calculate Bridges
	std::vector<bridge> bridges;

	for (const auto &[i, j] : q)
	{
		if (progress)
			progress->consumed(1);

		auto &ri = inResidues[i];
		auto &rj = inResidues[j];

		bridge_type type = TestBridge(ri, rj);
		if (type == bridge_type::None)
			continue;

		bool found = false;
		for (bridge &bridge : bridges)
		{
			if (type != bridge.type or i != bridge.i.back() + 1)
				continue;

			if (type == bridge_type::Parallel and bridge.j.back() + 1 == j)
			{
				bridge.i.push_back(i);
				bridge.j.push_back(j);
				found = true;
				break;
			}

			if (type == bridge_type::AntiParallel and bridge.j.front() - 1 == j)
			{
				bridge.i.push_back(i);
				bridge.j.push_front(j);
				found = true;
				break;
			}
		}

		if (not found)
		{
			bridge bridge = {};

			bridge.type = type;
			bridge.i.push_back(i);
			bridge.chainI = ri.mAsymID;
			bridge.j.push_back(j);
			bridge.chainJ = rj.mAsymID;

			bridges.push_back(bridge);
		}
	}

	// extend ladders
	std::sort(bridges.begin(), bridges.end());

	for (uint32_t i = 0; i < bridges.size(); ++i)
	{
		for (uint32_t j = i + 1; j < bridges.size(); ++j)
		{
			uint32_t ibi = bridges[i].i.front();
			uint32_t iei = bridges[i].i.back();
			uint32_t jbi = bridges[i].j.front();
			uint32_t jei = bridges[i].j.back();
			uint32_t ibj = bridges[j].i.front();
			uint32_t iej = bridges[j].i.back();
			uint32_t jbj = bridges[j].j.front();
			uint32_t jej = bridges[j].j.back();

			if (bridges[i].type != bridges[j].type or
				NoChainBreak(inResidues[std::min(ibi, ibj)], inResidues[std::max(iei, iej)]) == false or
				NoChainBreak(inResidues[std::min(jbi, jbj)], inResidues[std::max(jei, jej)]) == false or
				ibj - iei >= 6 or
				(iei >= ibj and ibi <= iej))
			{
				continue;
			}

			bool bulge;
			if (bridges[i].type == bridge_type::Parallel)
				bulge = ((jbj - jei < 6 and ibj - iei < 3) or (jbj - jei < 3));
			else
				bulge = ((jbi - jej < 6 and ibj - iei < 3) or (jbi - jej < 3));

			if (bulge)
			{
				bridges[i].i.insert(bridges[i].i.end(), bridges[j].i.begin(), bridges[j].i.end());
				if (bridges[i].type == bridge_type::Parallel)
					bridges[i].j.insert(bridges[i].j.end(), bridges[j].j.begin(), bridges[j].j.end());
				else
					bridges[i].j.insert(bridges[i].j.begin(), bridges[j].j.begin(), bridges[j].j.end());
				bridges.erase(bridges.begin() + j);
				--j;
			}
		}
	}

	// Sheet
	std::set<bridge *> ladderset;
	for (bridge &bridge : bridges)
	{
		ladderset.insert(&bridge);

		size_t n = bridge.i.size();
		if (n > dssp::kHistogramSize)
			n = dssp::kHistogramSize;

		if (bridge.type == bridge_type::Parallel)
			stats.histogram.parallel_bridges_per_ladder[n - 1] += 1;
		else
			stats.histogram.antiparallel_bridges_per_ladder[n - 1] += 1;
	}

	uint32_t sheet = 1, ladder = 0;
	while (not ladderset.empty())
	{
		std::set<bridge *> sheetset;
		sheetset.insert(*ladderset.begin());
		ladderset.erase(ladderset.begin());

		bool done = false;
		while (not done)
		{
			done = true;
			for (bridge *a : sheetset)
			{
				for (bridge *b : ladderset)
				{
					if (Linked(*a, *b))
					{
						sheetset.insert(b);
						ladderset.erase(b);
						done = false;
						break;
					}
				}
				if (not done)
					break;
			}
		}

		for (bridge *bridge : sheetset)
		{
			bridge->ladder = ladder;
			bridge->sheet = sheet;
			bridge->link = sheetset;

			++ladder;
		}

		size_t nrOfLaddersPerSheet = sheetset.size();
		if (nrOfLaddersPerSheet > dssp::kHistogramSize)
			nrOfLaddersPerSheet = dssp::kHistogramSize;
		if (nrOfLaddersPerSheet == 1 and (*sheetset.begin())->i.size() > 1)
			stats.histogram.ladders_per_sheet[0] += 1;
		else if (nrOfLaddersPerSheet > 1)
			stats.histogram.ladders_per_sheet[nrOfLaddersPerSheet - 1] += 1;

		++sheet;
	}

	for (bridge &bridge : bridges)
	{
		// find out if any of the i and j set members already have
		// a bridge assigned, if so, we're assigning bridge 2

		uint32_t betai = 0, betaj = 0;

		for (uint32_t l : bridge.i)
		{
			if (inResidues[l].GetBetaPartner(0).m_residue != nullptr)
			{
				betai = 1;
				break;
			}
		}

		for (uint32_t l : bridge.j)
		{
			if (inResidues[l].GetBetaPartner(0).m_residue != nullptr)
			{
				betaj = 1;
				break;
			}
		}

		structure_type ss = structure_type::Betabridge;
		if (bridge.i.size() > 1)
			ss = structure_type::Strand;

		if (bridge.type == bridge_type::Parallel)
		{
			stats.count.H_bonds_in_parallel_bridges += bridge.i.back() - bridge.i.front() + 2;

			std::deque<uint32_t>::iterator j = bridge.j.begin();
			for (uint32_t i : bridge.i)
				inResidues[i].SetBetaPartner(betai, inResidues[*j++], bridge.ladder, true);

			j = bridge.i.begin();
			for (uint32_t i : bridge.j)
				inResidues[i].SetBetaPartner(betaj, inResidues[*j++], bridge.ladder, true);
		}
		else
		{
			stats.count.H_bonds_in_antiparallel_bridges += bridge.i.back() - bridge.i.front() + 2;

			std::deque<uint32_t>::reverse_iterator j = bridge.j.rbegin();
			for (uint32_t i : bridge.i)
				inResidues[i].SetBetaPartner(betai, inResidues[*j++], bridge.ladder, false);

			j = bridge.i.rbegin();
			for (uint32_t i : bridge.j)
				inResidues[i].SetBetaPartner(betaj, inResidues[*j++], bridge.ladder, false);
		}

		for (uint32_t i = bridge.i.front(); i <= bridge.i.back(); ++i)
		{
			if (inResidues[i].GetSecondaryStructure() != structure_type::Strand)
				inResidues[i].SetSecondaryStructure(ss);
			inResidues[i].SetSheet(bridge.sheet);
		}

		for (uint32_t i = bridge.j.front(); i <= bridge.j.back(); ++i)
		{
			if (inResidues[i].GetSecondaryStructure() != structure_type::Strand)
				inResidues[i].SetSecondaryStructure(ss);
			inResidues[i].SetSheet(bridge.sheet);
		}
	}

	// Construct the 'strands'
	// For mmCIF output, this is needed and since we now have the information available
	// it is best to do the calculation here.

	// strands are ranges of residues of length > 1 that form beta bridges in a sheet.

	for (uint32_t iSheet = 1; iSheet < sheet; ++iSheet)
	{
		std::vector<std::tuple<uint32_t,uint32_t>> strands;
		for (auto &bridge : bridges)
		{
			if (bridge.sheet != iSheet)
				continue;

			for (auto &range : { bridge.i, bridge.j})
			{
				auto imin = range.front();
				auto imax = range.back();

				if (imin == imax)
					continue;

				if (imin > imax)
					std::swap(imin, imax);

				auto ii = find_if(strands.begin(), strands.end(), [a = imin, b = imax] (std::tuple<uint32_t,uint32_t> &t)
				{
					auto &&[start, end] = t;

					bool result = false;
					if (start <= b and end >= a)
					{
						result = true;
						if (start > a)
							start = a;
						if (end < b)
							end = b;
					}

					return result;
				});

				if (ii == strands.end())
					strands.emplace_back(imin, imax);
			}
		}

		std::sort(strands.begin(), strands.end());

		// collapse ranges that overlap
		if (strands.size() > 1)
		{
			auto si = strands.begin();
			while (std::next(si) != strands.end())
			{
				auto &&[afirst, alast] = *si;
				auto &&[bfirst, blast] = *(std::next(si));

				if (alast >= bfirst)
				{
					bfirst = afirst;
					si = strands.erase(si);
					continue;
				}

				++si;
			}
		}

		for (size_t i = 0; i < strands.size(); ++i)
		{
			const auto &[first, last] = strands[i];
			for (auto nr = first; nr <= last; ++nr)
			{
				assert(inResidues[nr].mStrand == 0);
				inResidues[nr].SetStrand(i + 1);
			}
		}
	}
}

// --------------------------------------------------------------------

void CalculateAlphaHelices(std::vector<residue> &inResidues, statistics &stats, bool inPreferPiHelices = true)
{
	if (cif::VERBOSE)
		std::cerr << "calculating alpha helices" << std::endl;

	// Helix and Turn
	for (helix_type helixType : { helix_type::_3_10, helix_type::alpha, helix_type::pi })
	{
		uint32_t stride = static_cast<uint32_t>(helixType) + 3;

		for (uint32_t i = 0; i + stride < inResidues.size(); ++i)
		{
			if (NoChainBreak(inResidues[i], inResidues[i + stride]) and TestBond(&inResidues[i + stride], &inResidues[i]))
			{
				inResidues[i + stride].SetHelixFlag(helixType, helix_position_type::End);
				for (uint32_t j = i + 1; j < i + stride; ++j)
				{
					if (inResidues[j].GetHelixFlag(helixType) == helix_position_type::None)
						inResidues[j].SetHelixFlag(helixType, helix_position_type::Middle);
				}

				if (inResidues[i].GetHelixFlag(helixType) == helix_position_type::End)
					inResidues[i].SetHelixFlag(helixType, helix_position_type::StartAndEnd);
				else
					inResidues[i].SetHelixFlag(helixType, helix_position_type::Start);
			}
		}
	}

	for (auto &r : inResidues)
	{
		if (r.mKappa.has_value())
			r.SetBend(*r.mKappa > 70);
	}

	for (uint32_t i = 1; i + 4 < inResidues.size(); ++i)
	{
		if (inResidues[i].IsHelixStart(helix_type::alpha) and inResidues[i - 1].IsHelixStart(helix_type::alpha))
		{
			for (uint32_t j = i; j <= i + 3; ++j)
				inResidues[j].SetSecondaryStructure(structure_type::Alphahelix);
		}
	}

	for (uint32_t i = 1; i + 3 < inResidues.size(); ++i)
	{
		if (inResidues[i].IsHelixStart(helix_type::_3_10) and inResidues[i - 1].IsHelixStart(helix_type::_3_10))
		{
			bool empty = true;
			for (uint32_t j = i; empty and j <= i + 2; ++j)
				empty = inResidues[j].GetSecondaryStructure() == structure_type::Loop or inResidues[j].GetSecondaryStructure() == structure_type::Helix_3;
			if (empty)
			{
				for (uint32_t j = i; j <= i + 2; ++j)
					inResidues[j].SetSecondaryStructure(structure_type::Helix_3);
			}
		}
	}

	for (uint32_t i = 1; i + 5 < inResidues.size(); ++i)
	{
		if (inResidues[i].IsHelixStart(helix_type::pi) and inResidues[i - 1].IsHelixStart(helix_type::pi))
		{
			bool empty = true;
			for (uint32_t j = i; empty and j <= i + 4; ++j)
				empty = inResidues[j].GetSecondaryStructure() == structure_type::Loop or inResidues[j].GetSecondaryStructure() == structure_type::Helix_5 or
				        (inPreferPiHelices and inResidues[j].GetSecondaryStructure() == structure_type::Alphahelix);
			if (empty)
			{
				for (uint32_t j = i; j <= i + 4; ++j)
					inResidues[j].SetSecondaryStructure(structure_type::Helix_5);
			}
		}
	}

	for (uint32_t i = 1; i + 1 < inResidues.size(); ++i)
	{
		if (inResidues[i].GetSecondaryStructure() == structure_type::Loop)
		{
			bool isTurn = false;
			for (helix_type helixType : { helix_type::_3_10, helix_type::alpha, helix_type::pi })
			{
				uint32_t stride = 3 + static_cast<uint32_t>(helixType);
				for (uint32_t k = 1; k < stride and not isTurn; ++k)
					isTurn = (i >= k) and inResidues[i - k].IsHelixStart(helixType);
			}

			if (isTurn)
				inResidues[i].SetSecondaryStructure(structure_type::Turn);
			else if (inResidues[i].IsBend())
				inResidues[i].SetSecondaryStructure(structure_type::Bend);
		}
	}

	std::string asym;
	size_t helixLength = 0;
	for (auto &r : inResidues)
	{
		if (r.mAsymID != asym)
		{
			helixLength = 0;
			asym = r.mAsymID;
		}

		if (r.GetSecondaryStructure() == structure_type::Alphahelix)
			++helixLength;
		else if (helixLength > 0)
		{
			if (helixLength > dssp::kHistogramSize)
				helixLength = dssp::kHistogramSize;

			stats.histogram.residues_per_alpha_helix[helixLength - 1] += 1;
			helixLength = 0;
		}
	}
}

// --------------------------------------------------------------------

void CalculatePPHelices(std::vector<residue> &inResidues, statistics &stats, int stretch_length)
{
	if (cif::VERBOSE)
		std::cerr << "calculating pp helices" << std::endl;

	size_t N = inResidues.size();

	const float epsilon = 29;
	const float phi_min = -75 - epsilon;
	const float phi_max = -75 + epsilon;
	const float psi_min = 145 - epsilon;
	const float psi_max = 145 + epsilon;

	std::vector<float> phi(N), psi(N);

	for (uint32_t i = 1; i + 1 < inResidues.size(); ++i)
	{
		phi[i] = static_cast<float>(inResidues[i].mPhi.value_or(360));
		psi[i] = static_cast<float>(inResidues[i].mPsi.value_or(360));
	}

	for (uint32_t i = 1; i + 3 < inResidues.size(); ++i)
	{
		switch (stretch_length)
		{
			case 2:
			{
				if (phi_min > phi[i + 0] or phi[i + 0] > phi_max or
					phi_min > phi[i + 1] or phi[i + 1] > phi_max)
					continue;

				if (psi_min > psi[i + 0] or psi[i + 0] > psi_max or
					psi_min > psi[i + 1] or psi[i + 1] > psi_max)
					continue;

				// auto phi_avg = (phi[i + 0] + phi[i + 1]) / 2;
				// auto phi_sq = (phi[i + 0] - phi_avg) * (phi[i + 0] - phi_avg) +
				// 			  (phi[i + 1] - phi_avg) * (phi[i + 1] - phi_avg);

				// if (phi_sq >= 200)
				// 	continue;

				// auto psi_avg = (psi[i + 0] + psi[i + 1]) / 2;
				// auto psi_sq = (psi[i + 0] - psi_avg) * (psi[i + 0] - psi_avg) +
				// 			  (psi[i + 1] - psi_avg) * (psi[i + 1] - psi_avg);

				// if (psi_sq >= 200)
				// 	continue;

				switch (inResidues[i].GetHelixFlag(helix_type::pp))
				{
					case helix_position_type::None:
						inResidues[i].SetHelixFlag(helix_type::pp, helix_position_type::Start);
						break;

					case helix_position_type::End:
						inResidues[i].SetHelixFlag(helix_type::pp, helix_position_type::Middle);
						break;

					default:
						break;
				}

				inResidues[i + 1].SetHelixFlag(helix_type::pp, helix_position_type::End);

				if (inResidues[i].GetSecondaryStructure() == structure_type::Loop)
					inResidues[i].SetSecondaryStructure(structure_type::Helix_PPII);
				if (inResidues[i + 1].GetSecondaryStructure() == structure_type::Loop)
					inResidues[i + 1].SetSecondaryStructure(structure_type::Helix_PPII);
			}
			break;

			case 3:
			{
				if (phi_min > phi[i + 0] or phi[i + 0] > phi_max or
					phi_min > phi[i + 1] or phi[i + 1] > phi_max or
					phi_min > phi[i + 2] or phi[i + 2] > phi_max)
					continue;

				if (psi_min > psi[i + 0] or psi[i + 0] > psi_max or
					psi_min > psi[i + 1] or psi[i + 1] > psi_max or
					psi_min > psi[i + 2] or psi[i + 2] > psi_max)
					continue;

				// auto phi_avg = (phi[i + 0] + phi[i + 1] + phi[i + 2]) / 3;
				// auto phi_sq = (phi[i + 0] - phi_avg) * (phi[i + 0] - phi_avg) +
				// 			  (phi[i + 1] - phi_avg) * (phi[i + 1] - phi_avg) +
				// 			  (phi[i + 2] - phi_avg) * (phi[i + 2] - phi_avg);

				// if (phi_sq >= 300)
				// 	continue;

				// auto psi_avg = (psi[i + 0] + psi[i + 1] + psi[i + 2]) / 3;
				// auto psi_sq = (psi[i + 0] - psi_avg) * (psi[i + 0] - psi_avg) +
				// 			  (psi[i + 1] - psi_avg) * (psi[i + 1] - psi_avg) +
				// 			  (psi[i + 2] - psi_avg) * (psi[i + 2] - psi_avg);

				// if (psi_sq >= 300)
				// 	continue;

				switch (inResidues[i].GetHelixFlag(helix_type::pp))
				{
					case helix_position_type::None:
						inResidues[i].SetHelixFlag(helix_type::pp, helix_position_type::Start);
						break;

					case helix_position_type::End:
						inResidues[i].SetHelixFlag(helix_type::pp, helix_position_type::StartAndEnd);
						break;

					default:
						break;
				}

				inResidues[i + 1].SetHelixFlag(helix_type::pp, helix_position_type::Middle);
				inResidues[i + 2].SetHelixFlag(helix_type::pp, helix_position_type::End);

				if (inResidues[i + 0].GetSecondaryStructure() == structure_type::Loop)
					inResidues[i + 0].SetSecondaryStructure(structure_type::Helix_PPII);

				if (inResidues[i + 1].GetSecondaryStructure() == structure_type::Loop)
					inResidues[i + 1].SetSecondaryStructure(structure_type::Helix_PPII);

				if (inResidues[i + 2].GetSecondaryStructure() == structure_type::Loop)
					inResidues[i + 2].SetSecondaryStructure(structure_type::Helix_PPII);

				break;
			}

			default:
				throw std::runtime_error("Unsupported stretch length");
		}
	}
}

// --------------------------------------------------------------------

struct DSSP_impl
{
	DSSP_impl(const cif::datablock &db, int model_nr, int min_poly_proline_stretch_length);

	auto findRes(const std::string &asymID, int seqID)
	{
		return std::find_if(mResidues.begin(), mResidues.end(), [&](auto &r)
			{ return r.mAsymID == asymID and r.mSeqID == seqID; });
	}

	void calculateSurface();
	void calculateSecondaryStructure();

	std::string GetPDBHEADERLine();
	std::string GetPDBCOMPNDLine();
	std::string GetPDBSOURCELine();
	std::string GetPDBAUTHORLine();

	const cif::datablock &mDB;
	std::vector<residue> mResidues;
	std::vector<std::pair<residue *, residue *>> mSSBonds;
	int m_min_poly_proline_stretch_length;
	statistics mStats = {};
};

// --------------------------------------------------------------------

DSSP_impl::DSSP_impl(const cif::datablock &db, int model_nr, int min_poly_proline_stretch_length)
	: mDB(db)
	, m_min_poly_proline_stretch_length(min_poly_proline_stretch_length)
{
	using namespace cif::literals;

	if (cif::VERBOSE)
		std::cerr << "loading residues" << std::endl;

	int resNumber = 0;

	auto &pdbx_poly_seq_scheme = mDB["pdbx_poly_seq_scheme"];
	auto &atom_site = mDB["atom_site"];

	using key_type = std::tuple<std::string,int>;
	using index_type = std::map<key_type, size_t>;

	index_type index;

	mResidues.reserve(pdbx_poly_seq_scheme.size());

	for (const auto &[asym_id, seq_id, pdb_strand_id, pdb_seq_num, pdb_ins_code]
		: pdbx_poly_seq_scheme.rows<std::string,int, std::string, int, std::string>("asym_id", "seq_id", "pdb_strand_id", "pdb_seq_num", "pdb_ins_code"))
	{
		index[{asym_id, seq_id}] = mResidues.size();
		mResidues.emplace_back(model_nr, pdb_strand_id, pdb_seq_num, pdb_ins_code);
	}

	for (auto atom : atom_site)
	{
		std::string asym_id;
		int seq_id;

		cif::tie(asym_id, seq_id) = atom.get("label_asym_id", "label_seq_id");
		auto i = index.find({asym_id, seq_id});
		if (i == index.end())
			continue;
		
		mResidues[i->second].addAtom(atom);
	}

	for (auto &residue : mResidues)
		residue.finish();
	
	mResidues.erase(std::remove_if(mResidues.begin(), mResidues.end(), [](const dssp::residue &r) { return not r.mComplete; }), mResidues.end());
	mStats.count.chains = 1;

	chain_break_type brk = chain_break_type::NewChain;
	for (size_t i = 0; i < mResidues.size(); ++i)
	{
		auto &residue = mResidues[i];
		++resNumber;

		if (i > 0)
		{
			if (distance(mResidues[i - 1].mC, mResidues[i].mN) > kMaxPeptideBondLength)
			{
				++mStats.count.chains;
				if (mResidues[i - 1].mAsymID == mResidues[i].mAsymID)
					brk = chain_break_type::Gap;
				else
					brk = chain_break_type::NewChain;

				++resNumber;
			}
		}

		residue.mChainBreak = brk;
		residue.mNumber = resNumber;

		brk = chain_break_type::None;
	}

	mStats.count.residues = static_cast<uint32_t>(mResidues.size());

	for (size_t i = 0; i + 1 < mResidues.size(); ++i)
	{
		auto &cur = mResidues[i];

		auto &next = mResidues[i + 1];
		next.mPrev = &cur;
		cur.mNext = &next;
	}

	for (size_t i = 0; i < mResidues.size(); ++i)
	{
		auto &cur = mResidues[i];

		if (i >= 2 and i + 2 < mResidues.size())
		{
			auto &prevPrev = mResidues[i - 2];
			auto &nextNext = mResidues[i + 2];

			if (NoChainBreak(prevPrev, nextNext) and prevPrev.mSeqID + 4 == nextNext.mSeqID)
			{
				float ckap = cosinus_angle(
					cur.mCAlpha,
					prevPrev.mCAlpha,
					nextNext.mCAlpha,
					cur.mCAlpha);
				float skap = std::sqrt(1 - ckap * ckap);
				cur.mKappa = std::atan2(skap, ckap) * static_cast<float>(180 / kPI);
			}
		}

		if (i + 1 < mResidues.size())
		{
			auto &next = mResidues[i + 1];
			next.assignHydrogen();
			if (NoChainBreak(cur, next))
			{
				cur.mPsi = dihedral_angle(cur.mN, cur.mCAlpha, cur.mC, next.mN);
				cur.mOmega = dihedral_angle(cur.mCAlpha, cur.mC, next.mN, next.mCAlpha);
			}
		}

		if (i > 0)
		{
			auto &prev = mResidues[i - 1];
			if (NoChainBreak(prev, cur))
			{
				cur.mTCO = cosinus_angle(cur.mC, cur.mO, prev.mC, prev.mO);
				cur.mPhi = dihedral_angle(prev.mC, cur.mN, cur.mCAlpha, cur.mC);
			}
		}

		if (i >= 1 and i + 2 < mResidues.size())
		{
			auto &prev = mResidues[i - 1];
			auto &next = mResidues[i + 1];
			auto &nextNext = mResidues[i + 2];

			if (NoChainBreak(prev, nextNext))
				cur.mAlpha = dihedral_angle(prev.mCAlpha, cur.mCAlpha, next.mCAlpha, nextNext.mCAlpha);
		}
	}
}

void DSSP_impl::calculateSecondaryStructure()
{
	if (cif::VERBOSE)
		std::cerr << "calculating secondary structure" << std::endl;

	using namespace cif::literals;

	for (auto [asym1, seq1, asym2, seq2] : mDB["struct_conn"].find<std::string, int, std::string, int>("conn_type_id"_key == "disulf",
			 "ptnr1_label_asym_id", "ptnr1_label_seq_id", "ptnr2_label_asym_id", "ptnr2_label_seq_id"))
	{
		auto r1 = findRes(asym1, seq1);
		if (r1 == mResidues.end())
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Missing (incomplete?) residue for SS bond when trying to find " << asym1 << '/' << seq1 << std::endl;
			continue;
			// throw std::runtime_error("Invalid file, missing residue for SS bond");
		}

		auto r2 = findRes(asym2, seq2);
		if (r2 == mResidues.end())
		{
			if (cif::VERBOSE > 0)
				std::cerr << "Missing (incomplete?) residue for SS bond when trying to find " << asym2 << '/' << seq2 << std::endl;
			continue;
			// throw std::runtime_error("Invalid file, missing residue for SS bond");
		}

		mSSBonds.emplace_back(&*r1, &*r2);
	}

	// Prefetch the c-alpha positions. No, really, that might be the trick

	std::vector<point> cAlphas;
	cAlphas.reserve(mResidues.size());
	for (auto &r : mResidues)
		cAlphas.emplace_back(r.mCAlpha);

	std::unique_ptr<cif::progress_bar> progress;
	if (cif::VERBOSE == 0 or cif::VERBOSE == 1)
		progress.reset(new cif::progress_bar((mResidues.size() * (mResidues.size() - 1)) / 2, "calculate distances"));

	// Calculate the HBond energies
	std::vector<std::tuple<uint32_t,uint32_t>> near;

	for (uint32_t i = 0; i + 1 < mResidues.size(); ++i)
	{
		auto cai = cAlphas[i];

		for (uint32_t j = i + 1; j < mResidues.size(); ++j)
		{
			auto caj = cAlphas[j];

			if (distance_sq(cai, caj) > (kMinimalCADistance * kMinimalCADistance))
				continue;

			near.emplace_back(i, j);
		}

		if (progress)
			progress->consumed(mResidues.size() - i - 1);
	}

	if (cif::VERBOSE > 0)
		std::cerr << "Considering " << near.size() << " pairs of residues" << std::endl;

	progress.reset(nullptr);

	CalculateHBondEnergies(mResidues, near);
	CalculateBetaSheets(mResidues, mStats, near);
	CalculateAlphaHelices(mResidues, mStats);
	CalculatePPHelices(mResidues, mStats, m_min_poly_proline_stretch_length);

	if (cif::VERBOSE > 1)
	{
		for (auto &r : mResidues)
		{
			char helix[5] = {};
			for (helix_type helixType : { helix_type::_3_10, helix_type::alpha, helix_type::pi, helix_type::pp })
			{
				switch (r.GetHelixFlag(helixType))
				{
					case helix_position_type::Start: helix[static_cast<int>(helixType)] = '>'; break;
					case helix_position_type::Middle: helix[static_cast<int>(helixType)] = helixType == helix_type::pp ? 'P' : '3' + static_cast<char>(helixType); break;
					case helix_position_type::StartAndEnd: helix[static_cast<int>(helixType)] = 'X'; break;
					case helix_position_type::End: helix[static_cast<int>(helixType)] = '<'; break;
					case helix_position_type::None: helix[static_cast<int>(helixType)] = ' '; break;
				}
			}

			auto id = r.mAsymID + ':' + std::to_string(r.mSeqID) + '/' + r.mCompoundID;

			std::cerr << id << std::string(12 - id.length(), ' ')
					  << char(r.mSecondaryStructure) << ' '
					  << helix
					  << std::endl;
		}
	}

	// finish statistics
	mStats.count.SS_bridges = static_cast<uint32_t>(mSSBonds.size());

	mStats.count.intra_chain_SS_bridges = 0;
	uint8_t ssBondNr = 0;
	for (const auto &[a, b] : mSSBonds)
	{
		if (a == b)
		{
			if (cif::VERBOSE > 0)
				std::cerr << "In the SS bonds list, the residue " << a->mAsymID << ':' << a->mSeqID << " is bonded to itself" << std::endl;
			continue;
		}

		if (a->mAsymID == b->mAsymID and NoChainBreak(a, b))
			++mStats.count.intra_chain_SS_bridges;

		a->mSSBridgeNr = b->mSSBridgeNr = ++ssBondNr;
	}

	mStats.count.H_bonds = 0;
	for (auto &r : mResidues)
	{
		auto donor = r.mHBondDonor;

		for (int i = 0; i < 2; ++i)
		{
			if (donor[i].res != nullptr and donor[i].energy < kMaxHBondEnergy)
			{
				++mStats.count.H_bonds;
				auto k = donor[i].res->mNumber - r.mNumber;
				if (k >= -5 and k <= 5)
					mStats.count.H_Bonds_per_distance[k + 5] += 1;
			}
		}
	}
}

void DSSP_impl::calculateSurface()
{
	CalculateAccessibilities(mResidues, mStats);
}

// --------------------------------------------------------------------

// Truncate lines in pseudo PDB format to this length
const int kTruncateAt = 127;

std::string FixStringLength(std::string s, std::string::size_type l = kTruncateAt)
{
	if (s.length() > l)
		s = s.substr(0, l - 4) + "... ";
	else if (s.length() < l)
		s.append(l - s.length(), ' ');

	return s;
}

std::string cif2pdbDate(const std::string &d)
{
	const std::regex rx(R"((\d{4})-(\d{2})(?:-(\d{2}))?)");
	const char *kMonths[12] = {
		"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
	};

	std::smatch m;
	std::ostringstream os;

	if (std::regex_match(d, m, rx))
	{
		int year = std::stoi(m[1].str());
		int month = std::stoi(m[2].str());

		if (m[3].matched)
			os << std::setw(2) << std::setfill('0') << stoi(m[3].str()) << '-';
		os << kMonths[month - 1] << '-' << std::setw(2) << std::setfill('0') << (year % 100);
	}

	return os.str();
}

std::string cif2pdbAuth(std::string name)
{
	const std::regex rx(R"(([^,]+), (\S+))");

	std::smatch m;
	if (std::regex_match(name, m, rx))
		name = m[2].str() + m[1].str();

	return name;
}

std::string DSSP_impl::GetPDBHEADERLine()
{
	std::string keywords;
	auto &cat1 = mDB["struct_keywords"];

	for (auto r : cat1)
	{
		keywords = FixStringLength(r["pdbx_keywords"].as<std::string>(), 40);
		break;
	}

	std::string date;
	for (auto r : mDB["pdbx_database_status"])
	{
		date = r["recvd_initial_deposition_date"].as<std::string>();
		if (date.empty())
			continue;
		date = cif2pdbDate(date);
		break;
	}

	if (date.empty())
	{
		for (auto r : mDB["database_PDB_rev"])
		{
			date = r["date_original"].as<std::string>();
			if (date.empty())
				continue;
			date = cif2pdbDate(date);
			break;
		}
	}

	date = FixStringLength(date, 9);

	//   0         1         2         3         4         5         6         7         8
	//   HEADER    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxDDDDDDDDD   IIII
	char header[] =
		"HEADER    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxDDDDDDDDD   IIII";

	std::copy(keywords.begin(), keywords.end(), header + 10);
	std::copy(date.begin(), date.end(), header + 50);

	std::string id = mDB.name();
	if (id.length() < 4)
		id.insert(id.end(), 4 - id.length(), ' ');
	else if (id.length() > 4)
		id.erase(id.begin() + 4, id.end());

	std::copy(id.begin(), id.end(), header + 62);

	return FixStringLength(header);
}

std::string DSSP_impl::GetPDBCOMPNDLine()
{
	// COMPND
	using namespace std::placeholders;
	using namespace cif::literals;

	int molID = 0;
	std::vector<std::string> cmpnd;

	for (auto r : mDB["entity"].find("type"_key == "polymer"))
	{
		std::string entityID = r["id"].as<std::string>();

		++molID;
		cmpnd.push_back("MOL_ID: " + std::to_string(molID));

		std::string molecule = r["pdbx_description"].as<std::string>();
		cmpnd.push_back("MOLECULE: " + molecule);

		auto poly = mDB["entity_poly"].find("entity_id"_key == entityID);
		if (not poly.empty())
		{
			std::string chains = poly.front()["pdbx_strand_id"].as<std::string>();
			cif::replace_all(chains, ",", ", ");
			cmpnd.push_back("CHAIN: " + chains);
		}

		std::string fragment = r["pdbx_fragment"].as<std::string>();
		if (not fragment.empty())
			cmpnd.push_back("FRAGMENT: " + fragment);

		for (auto sr : mDB["entity_name_com"].find("entity_id"_key == entityID))
		{
			std::string syn = sr["name"].as<std::string>();
			if (not syn.empty())
				cmpnd.push_back("SYNONYM: " + syn);
		}

		std::string mutation = r["pdbx_mutation"].as<std::string>();
		if (not mutation.empty())
			cmpnd.push_back("MUTATION: " + mutation);

		std::string ec = r["pdbx_ec"].as<std::string>();
		if (not ec.empty())
			cmpnd.push_back("EC: " + ec);

		if (r["src_method"] == "man" or r["src_method"] == "syn")
			cmpnd.push_back("ENGINEERED: YES");

		std::string details = r["details"].as<std::string>();
		if (not details.empty())
			cmpnd.push_back("OTHER_DETAILS: " + details);
	}

	return FixStringLength("COMPND    " + cif::join(cmpnd, "; "), kTruncateAt);
}

std::string DSSP_impl::GetPDBSOURCELine()
{
	// SOURCE

	using namespace cif::literals;

	int molID = 0;
	std::vector<std::string> source;

	for (auto r : mDB["entity"])
	{
		if (r["type"] != "polymer")
			continue;

		std::string entityID = r["id"].as<std::string>();

		++molID;
		source.push_back("MOL_ID: " + std::to_string(molID));

		if (r["src_method"] == "syn")
			source.push_back("SYNTHETIC: YES");

		auto &gen = mDB["entity_src_gen"];
		const std::pair<const char *, const char *> kGenSourceMapping[] = {
			{ "gene_src_common_name", "ORGANISM_COMMON" },
			{ "pdbx_gene_src_gene", "GENE" },
			{ "gene_src_strain", "STRAIN" },
			{ "pdbx_gene_src_cell_line", "CELL_LINE" },
			{ "pdbx_gene_src_organelle", "ORGANELLE" },
			{ "pdbx_gene_src_cellular_location", "CELLULAR_LOCATION" },
			{ "pdbx_gene_src_scientific_name", "ORGANISM_SCIENTIFIC" },
			{ "pdbx_gene_src_ncbi_taxonomy_id", "ORGANISM_TAXID" },
			{ "pdbx_host_org_scientific_name", "EXPRESSION_SYSTEM" },
			{ "pdbx_host_org_ncbi_taxonomy_id", "EXPRESSION_SYSTEM_TAXID" },
			{ "pdbx_host_org_strain", "EXPRESSION_SYSTEM_STRAIN" },
			{ "pdbx_host_org_variant", "EXPRESSION_SYSTEM_VARIANT" },
			{ "pdbx_host_org_cellular_location", "EXPRESSION_SYSTEM_CELLULAR_LOCATION" },
			{ "pdbx_host_org_vector_type", "EXPRESSION_SYSTEM_VECTOR_TYPE" },
			{ "pdbx_host_org_vector", "EXPRESSION_SYSTEM_VECTOR" },
			{ "pdbx_host_org_gene", "EXPRESSION_SYSTEM_GENE" },
			{ "plasmid_name", "EXPRESSION_SYSTEM_PLASMID" }
		};

		for (auto gr : gen.find("entity_id"_key == entityID))
		{
			for (auto m : kGenSourceMapping)
			{
				std::string cname, sname;
				tie(cname, sname) = m;

				std::string s = gr[cname].as<std::string>();
				if (not s.empty())
					source.push_back(sname + ": " + s);
			}
		}

		auto &nat = mDB["entity_src_nat"];
		const std::pair<const char *, const char *> kNatSourceMapping[] = {
			{ "common_name", "ORGANISM_COMMON" },
			{ "strain", "STRAIN" },
			{ "pdbx_organism_scientific", "ORGANISM_SCIENTIFIC" },
			{ "pdbx_ncbi_taxonomy_id", "ORGANISM_TAXID" },
			{ "pdbx_cellular_location", "CELLULAR_LOCATION" },
			{ "pdbx_plasmid_name", "PLASMID" },
			{ "pdbx_organ", "ORGAN" },
			{ "details", "OTHER_DETAILS" }
		};

		for (auto nr : nat.find("entity_id"_key == entityID))
		{
			for (auto m : kNatSourceMapping)
			{
				std::string cname, sname;
				tie(cname, sname) = m;

				std::string s = nr[cname].as<std::string>();
				if (not s.empty())
					source.push_back(sname + ": " + s);
			}
		}
	}

	return FixStringLength("SOURCE    " + cif::join(source, "; "), kTruncateAt);
}

std::string DSSP_impl::GetPDBAUTHORLine()
{
	// AUTHOR
	std::vector<std::string> author;
	for (auto r : mDB["audit_author"])
		author.push_back(cif2pdbAuth(r["name"].as<std::string>()));

	return FixStringLength("AUTHOR    " + cif::join(author, "; "), kTruncateAt);
}

// --------------------------------------------------------------------

std::string dssp::residue_info::asym_id() const
{
	return m_impl->mAsymID;
}

std::string dssp::residue_info::compound_id() const
{
	return m_impl->mCompoundID;
}

char dssp::residue_info::compound_letter() const
{
	return MapResidue(compound_id());
}

int dssp::residue_info::seq_id() const
{
	return m_impl->mSeqID;
}

std::string dssp::residue_info::alt_id() const
{
	return m_impl->mAltID;
}

std::string dssp::residue_info::auth_asym_id() const
{
	return m_impl->mAuthAsymID;
}

int dssp::residue_info::auth_seq_id() const
{
	return m_impl->mAuthSeqID;
}

std::string dssp::residue_info::pdb_strand_id() const
{
	return m_impl->mPDBStrandID;
}

int dssp::residue_info::pdb_seq_num() const
{
	return m_impl->mPDBSeqNum;
}

std::string dssp::residue_info::pdb_ins_code() const
{
	return m_impl->mPDBInsCode;
}

std::optional<float> dssp::residue_info::alpha() const
{
	return m_impl->mAlpha;
}

std::optional<float> dssp::residue_info::kappa() const
{
	return m_impl->mKappa;
}

std::optional<float> dssp::residue_info::omega() const
{
	return m_impl->mOmega;
}

std::optional<float> dssp::residue_info::phi() const
{
	return m_impl->mPhi;
}

std::optional<float> dssp::residue_info::psi() const
{
	return m_impl->mPsi;
}

std::optional<float> dssp::residue_info::tco() const
{
	return m_impl->mTCO;
}

bool dssp::residue_info::is_pre_pro() const
{
	return m_impl->mType != kProline and m_impl->mNext != nullptr and m_impl->mNext->mType == kProline;
}

float dssp::residue_info::chiral_volume() const
{
	return m_impl->mChiralVolume;
}

const std::map<residue_type, std::vector<std::string>> kChiAtomsMap = {
	{ MapResidue("ASP"), { "CG", "OD1" } },
	{ MapResidue("ASN"), { "CG", "OD1" } },
	{ MapResidue("ARG"), { "CG", "CD", "NE", "CZ" } },
	{ MapResidue("HIS"), { "CG", "ND1" } },
	{ MapResidue("GLN"), { "CG", "CD", "OE1" } },
	{ MapResidue("GLU"), { "CG", "CD", "OE1" } },
	{ MapResidue("SER"), { "OG" } },
	{ MapResidue("THR"), { "OG1" } },
	{ MapResidue("LYS"), { "CG", "CD", "CE", "NZ" } },
	{ MapResidue("TYR"), { "CG", "CD1" } },
	{ MapResidue("PHE"), { "CG", "CD1" } },
	{ MapResidue("LEU"), { "CG", "CD1" } },
	{ MapResidue("TRP"), { "CG", "CD1" } },
	{ MapResidue("CYS"), { "SG" } },
	{ MapResidue("ILE"), { "CG1", "CD1" } },
	{ MapResidue("MET"), { "CG", "SD", "CE" } },
	{ MapResidue("MSE"), { "CG", "SE", "CE" } },
	{ MapResidue("PRO"), { "CG", "CD" } },
	{ MapResidue("VAL"), { "CG1" } }
};

std::size_t dssp::residue_info::nr_of_chis() const
{
	auto i = kChiAtomsMap.find(m_impl->mType);

	return i != kChiAtomsMap.end() ? i->second.size() : 0;
}

float dssp::residue_info::chi(std::size_t index) const
{
	float result = 0;

	auto type = m_impl->mType;

	auto i = kChiAtomsMap.find(type);
	if (i != kChiAtomsMap.end() and index < i->second.size())
	{
		std::vector<std::string> atoms{ "N", "CA", "CB" };

		atoms.insert(atoms.end(), i->second.begin(), i->second.end());

		// in case we have a positive chiral volume we need to swap atoms
		if (m_impl->mChiralVolume > 0)
		{
			if (type == kLeucine)
				atoms.back() = "CD2";
			if (type == kValine)
				atoms.back() = "CG2";
		}

		result = static_cast<float>(dihedral_angle(
			m_impl->get_atom(atoms[index + 0]),
			m_impl->get_atom(atoms[index + 1]),
			m_impl->get_atom(atoms[index + 2]),
			m_impl->get_atom(atoms[index + 3])));
	}

	return result;
}

std::tuple<float, float, float> dssp::residue_info::ca_location() const
{
	return { m_impl->mCAlpha.mX, m_impl->mCAlpha.mY, m_impl->mCAlpha.mZ };
}

chain_break_type dssp::residue_info::chain_break() const
{
	return m_impl->mChainBreak;
}

int dssp::residue_info::nr() const
{
	return m_impl->mNumber;
}

structure_type dssp::residue_info::type() const
{
	return m_impl->mSecondaryStructure;
}

int dssp::residue_info::ssBridgeNr() const
{
	return m_impl->mSSBridgeNr;
}

helix_position_type dssp::residue_info::helix(helix_type helixType) const
{
	return m_impl->GetHelixFlag(helixType);
}

bool dssp::residue_info::is_alpha_helix_end_before_start() const
{
	bool result = false;

	if (m_impl->mNext != nullptr)
		result = m_impl->GetHelixFlag(helix_type::alpha) == helix_position_type::End and m_impl->mNext->GetHelixFlag(helix_type::alpha) == helix_position_type::Start;

	return result;
}

bool dssp::residue_info::bend() const
{
	return m_impl->IsBend();
}

double dssp::residue_info::accessibility() const
{
	return m_impl->mAccessibility;
}

std::tuple<dssp::residue_info, int, bool> dssp::residue_info::bridge_partner(int i) const
{
	auto bp = m_impl->GetBetaPartner(i);

	residue_info ri(bp.m_residue);

	return std::make_tuple(std::move(ri), bp.ladder, bp.parallel);
}

int dssp::residue_info::sheet() const
{
	return m_impl->GetSheet();
}

int dssp::residue_info::strand() const
{
	return m_impl->GetStrand();
}

std::tuple<dssp::residue_info, double> dssp::residue_info::acceptor(int i) const
{
	auto &a = m_impl->mHBondAcceptor[i];
	return { residue_info(a.res), a.energy };
}

std::tuple<dssp::residue_info, double> dssp::residue_info::donor(int i) const
{
	auto &d = m_impl->mHBondDonor[i];
	return { residue_info(d.res), d.energy };
}

// --------------------------------------------------------------------

dssp::iterator::iterator(residue *res)
	: m_current(res)
{
}

dssp::iterator &dssp::iterator::operator++()
{
	++m_current.m_impl;
	return *this;
}

dssp::iterator &dssp::iterator::operator--()
{
	--m_current.m_impl;
	return *this;
}

// --------------------------------------------------------------------

dssp::dssp(const cif::mm::structure &s, int min_poly_proline_stretch_length, bool calculateSurfaceAccessibility)
	: dssp(s.get_datablock(), static_cast<int>(s.get_model_nr()), min_poly_proline_stretch_length, calculateSurfaceAccessibility)
{
}

dssp::dssp(const cif::datablock &db, int model_nr, int min_poly_proline_stretch, bool calculateSurfaceAccessibility)
	: m_impl(new DSSP_impl(db, model_nr, min_poly_proline_stretch))
{
	if (calculateSurfaceAccessibility)
	{
		std::thread t(std::bind(&DSSP_impl::calculateSurface, m_impl));
		m_impl->calculateSecondaryStructure();
		t.join();
	}
	else
		m_impl->calculateSecondaryStructure();
}

dssp::~dssp()
{
	delete m_impl;
}

dssp::iterator dssp::begin() const
{
	return iterator(m_impl->mResidues.empty() ? nullptr : m_impl->mResidues.data());
}

dssp::iterator dssp::end() const
{
	// careful now, MSVC is picky when it comes to dereferencing iterators that are at the end.
	residue *res = nullptr;
	if (not m_impl->mResidues.empty())
	{
		res = m_impl->mResidues.data();
		res += m_impl->mResidues.size();
	}

	return iterator(res);
}

dssp::residue_info dssp::operator[](const key_type &key) const
{
	auto i = std::find_if(begin(), end(),
		[key](const residue_info &res) { return res.asym_id() == std::get<0>(key) and res.seq_id() == std::get<1>(key); });
	
	if (i == end())
		throw std::out_of_range("Could not find residue with supplied key");
	
	return *i;
}

dssp::statistics dssp::get_statistics() const
{
	return m_impl->mStats;
}

std::string dssp::get_pdb_header_line(pdb_record_type pdb_record) const
{
	switch (pdb_record)
	{
		case pdb_record_type::HEADER:
			return m_impl->GetPDBHEADERLine();
		case pdb_record_type::COMPND:
			return m_impl->GetPDBCOMPNDLine();
		case pdb_record_type::SOURCE:
			return m_impl->GetPDBSOURCELine();
		case pdb_record_type::AUTHOR:
			return m_impl->GetPDBAUTHORLine();
		default:
			return {};
	}
}

// --------------------------------------------------------------------

void dssp::write_legacy_output(std::ostream& os) const
{
	writeDSSP(*this, os);
}

void dssp::annotate(cif::datablock &db, bool writeOther, bool writeDSSPCategories) const
{
	annotateDSSP(db, *this, writeOther, writeDSSPCategories);
}


