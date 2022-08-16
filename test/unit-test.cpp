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

#define BOOST_TEST_MODULE DSSP_Test
#include <boost/test/included/unit_test.hpp>
#include <boost/algorithm/string.hpp>

#include <stdexcept>

#include "DSSP.hpp"

#include "dssp_wrapper.hpp"

namespace ba = boost::algorithm;

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ut_dssp)
{
	using namespace std::literals;

	cif::file f("1cbs.cif.gz");
	BOOST_ASSERT(f.is_valid());

	dssp::DSSP dssp(f.front(), 1, 3, true);

	std::stringstream test;

	writeDSSP(dssp, test);

	std::ifstream reference("1cbs.dssp");

	BOOST_CHECK(reference.is_open());

	std::string line_t, line_r;
	BOOST_CHECK(std::getline(test, line_t) and std::getline(reference, line_r));

	const char* kHeaderLineStart = "==== Secondary Structure Definition by the program DSSP, NKI version 4.0                           ====";
	BOOST_CHECK(line_t.compare(0, std::strlen(kHeaderLineStart), kHeaderLineStart) == 0);
	BOOST_CHECK(line_r.compare(0, std::strlen(kHeaderLineStart), kHeaderLineStart) == 0);

	for (int line_nr = 2; ; ++line_nr)
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

		BOOST_CHECK(line_t == line_r);
	}

	BOOST_CHECK(test.eof());
	BOOST_CHECK(reference.eof());
}

BOOST_AUTO_TEST_CASE(ut_mmcif_2)
{
	using namespace std::literals;
	using namespace cif::literals;

	cif::file f("1cbs.cif.gz");
	BOOST_ASSERT(f.is_valid());

	dssp::DSSP dssp(f.front(), 1, 3, true);

	std::stringstream test;

	annotateDSSP(f.front(), dssp, true, test);

	cif::file rf("1cbs-dssp.cif");

	// structure.datablock()["software"].erase("name"_key == "dssp");
	// rs.datablock()["software"].erase("name"_key == "dssp");
	
	// generate some output on different files:
	cif::VERBOSE = 2;

	// BOOST_CHECK(f.front() == rf.front());
}
