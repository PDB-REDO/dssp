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

#define CATCH_CONFIG_RUNNER

#if CATCH22
# include <catch2/catch.hpp>
#else
# include <catch2/catch_all.hpp>
#endif

#include "../libdssp/src/dssp-io.hpp"
#include "../src/revision.hpp"
#include "dssp.hpp"

#include <cif++/dictionary_parser.hpp>

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

std::filesystem::path gTestDir = std::filesystem::current_path();

int main(int argc, char *argv[])
{
	Catch::Session session; // There must be exactly one instance

	// Build a new parser on top of Catch2's
#if CATCH22
	using namespace Catch::clara;
#else
	// Build a new parser on top of Catch2's
	using namespace Catch::Clara;
#endif

	auto cli = session.cli()                                // Get Catch2's command line parser
	           | Opt(gTestDir, "data-dir")                  // bind variable to a new option, with a hint string
	                 ["-D"]["--data-dir"]                   // the option names it will respond to
	           ("The directory containing the data files"); // description string for the help output

	// Now pass the new composite back to Catch2 so it uses that
	session.cli(cli);

	// Let Catch2 (using Clara) parse the command line
	int returnCode = session.applyCommandLine(argc, argv);
	if (returnCode != 0) // Indicates a command line error
		return returnCode;

	return session.run();
}

// --------------------------------------------------------------------

TEST_CASE("ut_dssp")
{
	using namespace std::literals;

	cif::file f(gTestDir / "1cbs.cif.gz");
	REQUIRE(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	std::stringstream test;

	writeDSSP(dssp, test);

	std::ifstream reference(gTestDir / "1cbs.dssp");

	CHECK(reference.is_open());

	std::string line_t, line_r;
	CHECK((std::getline(test, line_t) and std::getline(reference, line_r)));

	char kHeaderLineStart[] = "==== Secondary Structure Definition by the program DSSP, NKI version 4.4.5                         ====";
	memcpy(kHeaderLineStart + 69, kVersionNumber, strlen(kVersionNumber));

	CHECK(line_t.compare(0, std::strlen(kHeaderLineStart), kHeaderLineStart) == 0);
	// CHECK(line_r.compare(0, std::strlen(kHeaderLineStart), kHeaderLineStart) == 0);

	for (int line_nr = 2;; ++line_nr)
	{
		bool done_t = not std::getline(test, line_t);
		bool done_r = not std::getline(reference, line_r);

		CHECK(done_r == done_t);
		if (done_r)
			break;

		if (line_t != line_r)
			std::cerr << line_nr << std::endl
					  << line_t << std::endl
					  << line_r << std::endl;

		if (line_t != line_r)
		{
			CHECK(line_t == line_r);
			break;
		}
	}

	CHECK(test.eof());
	CHECK(reference.eof());
}

TEST_CASE("ut_mmcif_2")
{
	using namespace std::literals;
	using namespace cif::literals;

	cif::file f(gTestDir / "1cbs.cif.gz");
	REQUIRE(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	std::stringstream test;

	dssp.annotate(f.front(), true, false);

	cif::file rf(gTestDir / "1cbs-dssp.cif");

	auto &db1 = f.front();
	auto &db2 = rf.front();

	db1["software"].erase("name"_key == "dssp");
	db1.erase(find_if(db1.begin(), db1.end(), [](cif::category &cat)
		{ return cat.name() == "audit_conform"; }));

	db2["software"].erase("name"_key == "dssp");
	db2.erase(find_if(db2.begin(), db2.end(), [](cif::category &cat)
		{ return cat.name() == "audit_conform"; }));

	// generate some output on different files:
	// cif::VERBOSE = 2;

	CHECK(f.front() == rf.front());
}

// --------------------------------------------------------------------

TEST_CASE("dssp_1")
{
	cif::file f(gTestDir / "1cbs.cif.gz");

	REQUIRE(f.is_valid());

	std::ifstream t(gTestDir / "1cbs-dssp-test.tsv");

	dssp dssp(f.front(), 1, 3, true);

	for (auto residue : dssp)
	{
		std::string line;
		getline(t, line);

		// std::cout << line << std::endl;

		auto fld = cif::split(line, "\t");

		CHECK(fld.size() == 3);
		if (fld.size() != 3)
			continue;

		int seqID;
		std::from_chars(fld[0].data(), fld[0].data() + fld[0].length(), seqID);
		std::string asymID{ fld[1] };
		std::string secstr{ fld[2] };
		if (secstr == "_")
			secstr = " ";

		CHECK(residue.asym_id() == asymID);
		CHECK(residue.seq_id() == seqID);
		CHECK((char)residue.type() == secstr.front());
	}
}

// --------------------------------------------------------------------

TEST_CASE("dssp_2")
{
	cif::file f(gTestDir / "1cbs.cif.gz");

	REQUIRE(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	std::ifstream t(gTestDir / "1cbs-dssp-test.tsv");
	std::string line;

	while (getline(t, line))
	{
		auto fld = cif::split(line, "\t");

		CHECK(fld.size() == 3);
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

		CHECK(ri.asym_id() == asymID);
		CHECK(ri.seq_id() == seqID);
		CHECK((char)ri.type() == secstr.front());
	}
}

// --------------------------------------------------------------------

TEST_CASE("dssp_3")
{
	cif::file f(gTestDir / "1cbs.cif.gz");

	REQUIRE(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	dssp.annotate(f.front(), true, true);

	// CHECK(f.is_valid());
}