/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2022 NKI/AVL, Netherlands Cancer Institute
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

#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include <charconv>
#include <stdexcept>

#include <cif++/cif.hpp>
#include <cif++/dssp/DSSP.hpp>

namespace tt = boost::test_tools;

std::filesystem::path gTestDir = std::filesystem::current_path(); // filled in first test

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

bool init_unit_test()
{
	cif::VERBOSE = 1;

	// not a test, just initialize test dir
	if (boost::unit_test::framework::master_test_suite().argc == 2)
		gTestDir = boost::unit_test::framework::master_test_suite().argv[1];

	// do this now, avoids the need for installing
	cif::add_file_resource("mmcif_pdbx.dic", gTestDir / ".." / "rsrc" / "mmcif_pdbx.dic");

	// initialize CCD location
	cif::add_file_resource("components.cif", gTestDir / ".." / "data" / "ccd-subset.cif");

	return true;
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dssp_1)
{
	cif::file f(gTestDir / "1cbs.cif");

	BOOST_ASSERT(f.is_valid());

	std::ifstream t(gTestDir / "1cbs-dssp-test.tsv");

	dssp::DSSP dssp(f.front(), 1, 3, true);

	for (auto residue : dssp)
	{
		std::string line;
		getline(t, line);

		std::cout << line << std::endl;

		auto f = cif::split(line, "\t");

		BOOST_CHECK_EQUAL(f.size(), 3);
		if (f.size() != 3)
			continue;

		int seqID;
		std::from_chars(f[0].begin(), f[0].end(), seqID);
		std::string asymID{ f[1] };
		std::string secstr{ f[2] };
		if (secstr == "_")
			secstr = " ";

		BOOST_CHECK_EQUAL(residue.asym_id(), asymID);
		BOOST_CHECK_EQUAL(residue.seq_id(), seqID);
		BOOST_CHECK_EQUAL((char)residue.type(), secstr.front());
	}
}