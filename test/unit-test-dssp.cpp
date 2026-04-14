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

#define CATCH_CONFIG_RUNNER

#include "../libdssp/src/dssp-io.hpp"
#include "../src/revision.hpp"
#include "dssp.hpp"

#include <catch2/catch_all.hpp>
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

std::filesystem::path gTestDir;

int main(int argc, char *argv[])
{
	gTestDir = std::filesystem::current_path();

	Catch::Session session; // There must be exactly one instance

	// Build a new parser on top of Catch2's
	using namespace Catch::Clara;

	std::filesystem::path rsrc_dir;

	auto cli = session.cli()                                // Get Catch2's command line parser
	           | Opt(gTestDir, "data-dir")                  // bind variable to a new option, with a hint string
	                 ["-D"]["--data-dir"]                   // the option names it will respond to
	           | Opt(rsrc_dir, "rsrc-dir")                  // bind variable to a new option, with a hint string
	                 ["-D"]["--rsrc-dir"]                   // the option names it will respond to
	           ("The directory containing the data files"); // description string for the help output

	// Now pass the new composite back to Catch2 so it uses that
	session.cli(cli);

	// Let Catch2 (using Clara) parse the command line
	int returnCode = session.applyCommandLine(argc, argv);
	if (returnCode != 0) // Indicates a command line error
		return returnCode;

	if (not rsrc_dir.empty() and std::filesystem::exists(rsrc_dir))
		cif::add_data_directory(rsrc_dir);

	cif::add_data_directory(gTestDir / "libdssp" / "mmcif_pdbx");

	cif::add_file_resource("components.cif", gTestDir / "minimal-components.cif");

	return session.run();
}

// --------------------------------------------------------------------

TEST_CASE("ut_dssp")
{
	using namespace std::literals;

	auto f = cif::pdb::read(gTestDir / "1cbs.cif.gz");
	REQUIRE(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	std::stringstream test;

	writeDSSP(dssp, test);

	std::ifstream reference(gTestDir / "1cbs.dssp");

	CHECK(reference.is_open());

	std::string line_t, line_r;
	CHECK((std::getline(test, line_t) and std::getline(reference, line_r)));

	char kHeaderLineStart[] = "==== Secondary Structure Definition by the program DSSP, NKI version 4.5.4                         ====";
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
		{
			if (cif::starts_with(line_t, "REFERENCE ") and cif::starts_with(line_r, "REFERENCE "))
				continue;

			// std::cerr << line_nr << '\n'
			// 		  << line_t << '\n'
			// 		  << line_r << '\n';
		}

		REQUIRE(line_t == line_r);
	}

	CHECK(test.eof());
	CHECK(reference.eof());
}

TEST_CASE("ut_mmcif_2")
{
	using namespace std::literals;
	using namespace cif::literals;

	auto f = cif::pdb::read(gTestDir / "1cbs.cif.gz");
	REQUIRE(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	dssp.annotate(f.front(), true, true);

	cif::file rf(gTestDir / "1cbs-dssp.cif");

	auto &db1 = f.front();
	auto &db2 = rf.front();

	db1["software"].erase("name"_key == "DSSP");
	std::erase_if(db1, [](cif::category &cat)
		{ return cat.name() == "audit_conform"; });

	db2["software"].erase("name"_key == "DSSP");
	std::erase_if(db2, [](cif::category &cat)
		{ return cat.name() == "audit_conform"; });

	if (not (f.front() == rf.front()))
	{
		// std::ofstream a(std::filesystem::temp_directory_path() / "dssp-test-a.cif");
		// a << f.front();
		// a.close();
		// std::ofstream b(std::filesystem::temp_directory_path() / "dssp-test-b.cif");
		// b << rf.front();
		// b.close();

		CHECK(false);
	}
}

// --------------------------------------------------------------------

TEST_CASE("dssp_1")
{
	auto f = cif::pdb::read(gTestDir / "1cbs.cif.gz");

	REQUIRE(f.is_valid());

	std::ifstream t(gTestDir / "1cbs-dssp-test.tsv");

	dssp dssp(f.front(), 1, 3, true);

	for (auto residue : dssp)
	{
		std::string line;
		getline(t, line);

		// std::cout << line << '\n';

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
	auto f = cif::pdb::read(gTestDir / "1cbs.cif.gz");

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
	auto f = cif::pdb::read(gTestDir / "1cbs.cif.gz");

	REQUIRE(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	dssp.annotate(f.front(), true, true);

	// CHECK(f.is_valid());
}

// --------------------------------------------------------------------

TEST_CASE("dssp_pdb")
{
	auto f = cif::pdb::read(gTestDir / "pdb1cbs.ent.gz");

	REQUIRE(f.is_valid());

	dssp dssp(f.front(), 1, 3, true);

	dssp.annotate(f.front(), true, true);
}