/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2022 Maarten L. Hekkelman
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

#include <filesystem>

#include <cfg.hpp>

namespace tt = boost::test_tools;
namespace utf = boost::unit_test;
namespace fs = std::filesystem;

fs::path gTestDir = fs::current_path();

// --------------------------------------------------------------------

bool init_unit_test()
{
	// not a test, just initialize test dir
	if (boost::unit_test::framework::master_test_suite().argc == 2)
		gTestDir = boost::unit_test::framework::master_test_suite().argv[1];

	return true;
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(t_1, * utf::tolerance(0.001))
{
	int argc = 3;
	const char *const argv[] = {
		"test", "--flag", nullptr
	};

	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option("flag", ""),
		cfg::make_option<int>("param_int", ""),
		cfg::make_option("param_int_2", 1, ""),
		cfg::make_option<float>("param_float", ""),
		cfg::make_option("param_float_2", 3.14f, ""));
	
	config.parse(argc, argv);

	BOOST_CHECK(config.has("flag"));
	BOOST_CHECK(not config.has("flag2"));

	BOOST_CHECK_EQUAL(config.get<int>("param_int_2"), 1);
	BOOST_CHECK_THROW(config.get<float>("param_int_2"), std::bad_any_cast);
	BOOST_CHECK_THROW(config.get<int>("param_int"), std::system_error);

	BOOST_TEST(config.get<float>("param_float_2") == 3.14);
	BOOST_CHECK_THROW(config.get<int>("param_float_2"), std::bad_any_cast);
	BOOST_CHECK_THROW(config.get<float>("param_float"), std::system_error);
}

BOOST_AUTO_TEST_CASE(t_2)
{
	int argc = 3;
	const char *const argv[] = {
		"test", "-vvvv", "--verbose", nullptr
	};

	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option("verbose,v", ""));
	
	config.parse(argc, argv);

	BOOST_CHECK_EQUAL(config.count("verbose"), 5);
}

BOOST_AUTO_TEST_CASE(t_3)
{
	int argc = 2;
	const char *const argv[] = {
		"test", "--param_int=42", nullptr
	};

	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<int>("param_int", ""));
	
	config.parse(argc, argv);

	BOOST_CHECK(config.has("param_int"));
	BOOST_CHECK_EQUAL(config.get<int>("param_int"), 42);
}

BOOST_AUTO_TEST_CASE(t_4)
{
	int argc = 3;
	const char *const argv[] = {
		"test", "--param_int", "42", nullptr
	};

	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<int>("param_int", ""));
	
	config.parse(argc, argv);

	BOOST_CHECK(config.has("param_int"));
	BOOST_CHECK_EQUAL(config.get<int>("param_int"), 42);
}

BOOST_AUTO_TEST_CASE(t_5)
{
	const char *const argv[] = {
		"test", "-i", "42", "-j43", nullptr
	};
	int argc = sizeof(argv) / sizeof(char*);
	
	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<int>("nr1,i", ""),
		cfg::make_option<int>("nr2,j", ""));
	
	config.parse(argc, argv);

	BOOST_CHECK(config.has("nr1"));
	BOOST_CHECK(config.has("nr2"));

	BOOST_CHECK_EQUAL(config.get<int>("nr1"), 42);
	BOOST_CHECK_EQUAL(config.get<int>("nr2"), 43);
}

BOOST_AUTO_TEST_CASE(t_6)
{
	const char *const argv[] = {
		"test", "-i", "42", "-j43", "foo", "bar", nullptr
	};
	int argc = sizeof(argv) / sizeof(char*);
	
	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<int>("nr1,i", ""),
		cfg::make_option<int>("nr2,j", ""));
	
	config.parse(argc, argv);

	BOOST_CHECK(config.has("nr1"));
	BOOST_CHECK(config.has("nr2"));

	BOOST_CHECK_EQUAL(config.get<int>("nr1"), 42);
	BOOST_CHECK_EQUAL(config.get<int>("nr2"), 43);

	BOOST_CHECK_EQUAL(config.operands().size(), 2);
	BOOST_CHECK_EQUAL(config.operands().front(), "foo");
	BOOST_CHECK_EQUAL(config.operands().back(), "bar");
}

BOOST_AUTO_TEST_CASE(t_7)
{
	const char *const argv[] = {
		"test", "--", "-i", "42", "-j43", "foo", "bar", nullptr
	};
	int argc = sizeof(argv) / sizeof(char*) - 1;
	
	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<int>("nr1,i", ""),
		cfg::make_option<int>("nr2,j", ""));
	
	config.parse(argc, argv);

	BOOST_CHECK(not config.has("nr1"));
	BOOST_CHECK(not config.has("nr2"));

	BOOST_CHECK_EQUAL(config.operands().size(), 5);

	auto compare = std::vector<std::string>{ argv[2], argv[3], argv[4], argv[5], argv[6] };
	BOOST_CHECK(config.operands() == compare);
}

BOOST_AUTO_TEST_CASE(t_8)
{
	const char *const argv[] = {
		"test", "-i", "foo", "-jbar", nullptr
	};
	int argc = sizeof(argv) / sizeof(char*) - 1;
	
	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<const char*>("i", ""),
		cfg::make_option<std::string_view>("j", ""),
		cfg::make_option("k", "baz", ""));
	
	config.parse(argc, argv);

	BOOST_CHECK(config.has("i"));
	BOOST_CHECK_EQUAL(config.get<std::string>("i"), "foo");
	BOOST_CHECK(config.has("j"));
	BOOST_CHECK_EQUAL(config.get<std::string>("j"), "bar");

	BOOST_CHECK(config.has("k"));
	BOOST_CHECK_EQUAL(config.get<std::string>("k"), "baz");
}

BOOST_AUTO_TEST_CASE(t_9)
{
	auto &config = cfg::config::instance();

	config.set_usage("usage: test [options]");

	config.init(
		"test [options]",
		cfg::make_option<const char*>("i", "First option"),
		cfg::make_option<std::string_view>("j", "This is the second option"),
		cfg::make_option("a-very-long-option-name,k", "baz", "And, you guessed it, this must be option three."));
	
// 	std::stringstream ss;

// 	int fd = open("/dev/null", O_RDWR);
// 	dup2(fd, STDOUT_FILENO);

// 	ss << config << std::endl;

// 	const char kExpected[] = R"(usage: test [options]
//   -i arg                                First option
//   -j arg                                This is the second option
//   -k [ --a-very-long-option-name ] arg (=baz)
//                                         And, you guessed it, this must be
//                                         option three.

// )";

// 	std::cerr << '>' << kExpected << '<' << std::endl;
// 	std::cerr << '>' << ss.str() << '<' << std::endl;

// 	BOOST_CHECK_EQUAL(ss.str(), kExpected);
}

BOOST_AUTO_TEST_CASE(t_10)
{
	std::string s1 = R"(SPDX-License-Identifier: BSD-2-Clause

Copyright (c) 2022 Maarten L. Hekkelman

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
)";

	cfg::word_wrapper ww(s1, 80);

	std::ostringstream os;

	for (auto line : ww)
		os << line << std::endl;
	
	BOOST_CHECK_EQUAL(os.str(), R"(SPDX-License-Identifier: BSD-2-Clause

Copyright (c) 2022 Maarten L. Hekkelman

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this 
list of conditions and the following disclaimer
2. Redistributions in binary form must reproduce the above copyright notice, 
this list of conditions and the following disclaimer in the documentation and/
or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR 
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

)");
}

BOOST_AUTO_TEST_CASE(t_11)
{
	const char *const argv[] = {
		"test", "-faap", "-fnoot", "-fmies", nullptr
	};
	int argc = sizeof(argv) / sizeof(char*) - 1;

	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<std::vector<std::string>>("file,f", ""));
	
	config.parse(argc, argv);

	BOOST_CHECK_EQUAL(config.count("file"), 3);
	
	std::vector<std::string> files = config.get<std::vector<std::string>>("file");
	BOOST_CHECK_EQUAL(files.size(), 3);
	BOOST_CHECK_EQUAL(files[0], "aap");
	BOOST_CHECK_EQUAL(files[1], "noot");
	BOOST_CHECK_EQUAL(files[2], "mies");
}

BOOST_AUTO_TEST_CASE(t_12)
{
	const char *const argv[] = {
		"test", "--aap", nullptr
	};
	int argc = sizeof(argv) / sizeof(char*) - 1;

	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<std::vector<std::string>>("file,f", ""));
	
	std::error_code ec;
	config.parse(argc, argv, ec);
	BOOST_CHECK(ec == cfg::config_error::unknown_option);

	config.set_ignore_unknown(true);
	ec = {};

	config.parse(argc, argv, ec);
	BOOST_CHECK(ec == std::errc());
}

// --------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(file_1, * utf::tolerance(0.001))
{
	const std::string_view config_file{ R"(
# This is a test configuration
aap=1
noot = 2
mies = 	
pi = 3.14
s = hello, world!
verbose
	)" };

	struct membuf : public std::streambuf
	{
		membuf(char * text, size_t length)
		{
			this->setg(text, text, text + length);
		}
	} buffer(const_cast<char *>(config_file.data()), config_file.length());

	std::istream is(&buffer);

	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<const char*>("aap", ""),
		cfg::make_option<int>("noot", ""),
		cfg::make_option<std::string>("mies", ""),
		cfg::make_option<float>("pi", ""),
		cfg::make_option<std::string>("s", ""),
		cfg::make_option("verbose,v", ""));
	
	std::error_code ec;

	config.parse_config_file(is, ec);

	BOOST_CHECK(ec == std::errc());

	BOOST_CHECK(config.has("aap"));
	BOOST_CHECK_EQUAL(config.get<std::string>("aap"), "1");

	BOOST_CHECK(config.has("noot"));
	BOOST_CHECK_EQUAL(config.get<int>("noot"), 2);

	BOOST_CHECK(config.has("pi"));
	BOOST_TEST(config.get<float>("pi") == 3.14);

	BOOST_CHECK(config.has("s"));
	BOOST_CHECK_EQUAL(config.get<std::string>("s"), "hello, world!");

	BOOST_CHECK(config.has("verbose"));
}

BOOST_AUTO_TEST_CASE(file_2)
{
	auto &config = cfg::config::instance();

	std::tuple<std::string_view,std::string_view,std::error_code> tests[] = {
		{ "aap !", "aap", make_error_code(cfg::config_error::invalid_config_file) },
		{ "aap=aap", "aap", {} },
		{ "aap", "aap", make_error_code(cfg::config_error::missing_argument_for_option) },
		{ "verbose=1", "verbose", make_error_code(cfg::config_error::option_does_not_accept_argument) },
				
	};

	for (const auto &[config_file, option, err] : tests)
	{
		struct membuf : public std::streambuf
		{
			membuf(char * text, size_t length)
			{
				this->setg(text, text, text + length);
			}
		} buffer(const_cast<char *>(config_file.data()), config_file.length());

		std::istream is(&buffer);

		std::error_code ec;
		config.init(
			"test [options]",
			cfg::make_option<const char*>("aap", ""),
			cfg::make_option<int>("noot", ""),
			cfg::make_option<float>("pi", ""),
			cfg::make_option<std::string>("s", ""),
			cfg::make_option("verbose,v", ""));
		
		config.parse_config_file(is, ec);

		BOOST_CHECK(ec == err);

		if (ec == std::errc())
			BOOST_CHECK(config.has(option));
	}
}

BOOST_AUTO_TEST_CASE(file_3)
{
	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<const char*>("aap", ""),
		cfg::make_option<int>("noot", ""),
		cfg::make_option<std::string>("config", ""));
	
	std::error_code ec;

	const char *const argv[] = {
		"test", "--aap=aap", "--noot=42", "--config=unit-test.conf", nullptr
	};
	int argc = sizeof(argv) / sizeof(char*) - 1;

	config.parse(argc, argv);

	config.parse_config_file("config", "bla-bla.conf", { gTestDir.string() }, ec);

	BOOST_CHECK(ec == std::errc());

	BOOST_CHECK(config.has("aap"));
	BOOST_CHECK_EQUAL(config.get<std::string>("aap"), "aap");

	BOOST_CHECK(config.has("noot"));
	BOOST_CHECK_EQUAL(config.get<int>("noot"), 42);
}

BOOST_AUTO_TEST_CASE(file_4)
{
	auto &config = cfg::config::instance();

	config.init(
		"test [options]",
		cfg::make_option<const char*>("aap", ""),
		cfg::make_option<int>("noot", ""),
		cfg::make_option<std::string>("config", ""));
	
	std::error_code ec;

	const char *const argv[] = {
		"test", "--aap=aap", nullptr
	};
	int argc = sizeof(argv) / sizeof(char*) - 1;

	config.parse(argc, argv);

	config.parse_config_file("config", "unit-test.conf", { gTestDir.string() }, ec);

	BOOST_CHECK(ec == std::errc());

	BOOST_CHECK(config.has("aap"));
	BOOST_CHECK_EQUAL(config.get<std::string>("aap"), "aap");

	BOOST_CHECK(config.has("noot"));
	BOOST_CHECK_EQUAL(config.get<int>("noot"), 3);
}