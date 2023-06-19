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

#include <stdexcept>

#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "dssp.hpp"
#include "dssp-io.hpp"

namespace fs = std::filesystem;

// --------------------------------------------------------------------

cif::file operator""_cf(const char *text, size_t length)
{
	struct membuf : public std::streambuf
	{
		membuf(char *text, size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(text), length);

	std::istream is(&buffer);
	return cif::file(is);
}

// --------------------------------------------------------------------

fs::path gTestDir = fs::current_path();

bool init_unit_test()
{
	cif::VERBOSE = 1;

	// not a test, just initialize test dir
	if (boost::unit_test::framework::master_test_suite().argc == 2)
	{
		gTestDir = boost::unit_test::framework::master_test_suite().argv[1];

		cif::add_data_directory(gTestDir / ".." / "rsrc");
	}

	// // do this now, avoids the need for installing
	// cif::add_file_resource("mmcif_pdbx.dic", gTestDir / ".." / "rsrc" / "mmcif_pdbx.dic");

	// // initialize CCD location
	// cif::add_file_resource("components.cif", gTestDir / ".." / "data" / "ccd-subset.cif");

	return true;
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ut_dssp)
{
	using namespace std::literals;

	cif::file f(gTestDir / "1cbs.cif.gz");
	BOOST_ASSERT(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	std::stringstream test;

	writeDSSP(dssp, test);

	std::ifstream reference(gTestDir / "1cbs.dssp");

	BOOST_CHECK(reference.is_open());

	std::string line_t, line_r;
	BOOST_CHECK(std::getline(test, line_t) and std::getline(reference, line_r));

	const char *kHeaderLineStart = "==== Secondary Structure Definition by the program DSSP, NKI version 4.3                           ====";
	BOOST_CHECK(line_t.compare(0, std::strlen(kHeaderLineStart), kHeaderLineStart) == 0);
	// BOOST_CHECK(line_r.compare(0, std::strlen(kHeaderLineStart), kHeaderLineStart) == 0);

	for (int line_nr = 2;; ++line_nr)
	{
		bool done_t = not std::getline(test, line_t);
		bool done_r = not std::getline(reference, line_r);

		BOOST_CHECK_EQUAL(done_r, done_t);
		if (done_r)
			break;

		if (line_t != line_r)
			std::cerr << line_nr << std::endl
					  << line_t << std::endl
					  << line_r << std::endl;

		if (line_t != line_r)
		{
			BOOST_CHECK(line_t == line_r);
			break;
		}
	}

	BOOST_CHECK(test.eof());
	BOOST_CHECK(reference.eof());
}

BOOST_AUTO_TEST_CASE(ut_mmcif_2)
{
	using namespace std::literals;
	using namespace cif::literals;

	cif::file f(gTestDir / "1cbs.cif.gz");
	BOOST_ASSERT(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	std::stringstream test;

	dssp.annotate(f.front(), true, false);
	test << f.front();

	cif::file rf(gTestDir / "1cbs-dssp.cif");

	// structure.datablock()["software"].erase("name"_key == "dssp");
	// rs.datablock()["software"].erase("name"_key == "dssp");

	// generate some output on different files:
	cif::VERBOSE = 2;

	// BOOST_CHECK(f.front() == rf.front());
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dssp_1)
{
	cif::file f(gTestDir / "1cbs.cif.gz");

	BOOST_ASSERT(f.is_valid());

	std::ifstream t(gTestDir / "1cbs-dssp-test.tsv");

	dssp dssp(f.front(), 1, 3, true);

	for (auto residue : dssp)
	{
		std::string line;
		getline(t, line);

		// std::cout << line << std::endl;

		auto fld = cif::split(line, "\t");

		BOOST_CHECK_EQUAL(fld.size(), 3);
		if (fld.size() != 3)
			continue;

		int seqID;
		std::from_chars(fld[0].data(), fld[0].data() + fld[0].length(), seqID);
		std::string asymID{ fld[1] };
		std::string secstr{ fld[2] };
		if (secstr == "_")
			secstr = " ";

		BOOST_CHECK_EQUAL(residue.asym_id(), asymID);
		BOOST_CHECK_EQUAL(residue.seq_id(), seqID);
		BOOST_CHECK_EQUAL((char)residue.type(), secstr.front());
	}
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dssp_2)
{
	cif::file f(gTestDir / "1cbs.cif.gz");

	BOOST_ASSERT(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	std::ifstream t(gTestDir / "1cbs-dssp-test.tsv");
	std::string line;

	while (getline(t, line))
	{
		auto fld = cif::split(line, "\t");

		BOOST_CHECK_EQUAL(fld.size(), 3);
		if (fld.size() != 3)
			continue;

		int seqID;
		std::from_chars(fld[0].data(), fld[0].data() + fld[0].length(), seqID);
		std::string asymID{ fld[1] };
		std::string secstr{ fld[2] };
		if (secstr == "_")
			secstr = " ";

		dssp::key_type key{ asymID, seqID };
		auto ri = dssp[key];
		
		BOOST_CHECK_EQUAL(ri.asym_id(), asymID);
		BOOST_CHECK_EQUAL(ri.seq_id(), seqID);
		BOOST_CHECK_EQUAL((char)ri.type(), secstr.front());
	}
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dssp_3)
{
	cif::file f(gTestDir / "1cbs.cif.gz");

	BOOST_ASSERT(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	dssp.annotate(f.front(), true, true);

	BOOST_TEST(f.is_valid());
}