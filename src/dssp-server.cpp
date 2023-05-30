/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2023 NKI/AVL, Netherlands Cancer Institute
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

#include <mcfp.hpp>

#include <zeep/http/html-controller.hpp>
#include <zeep/http/rest-controller.hpp>
#include <zeep/http/daemon.hpp>

// --------------------------------------------------------------------


class dssp_html_controller : public zeep::http::html_controller
{
  public:
	dssp_html_controller(const std::string &prefix)
		: zeep::http::html_controller(prefix)
	{
		mount("{css,scripts,fonts,images,favicon}/", &dssp_html_controller::handle_file);
		mount("{favicon.ico,browserconfig.xml,manifest.json}", &dssp_html_controller::handle_file);
		mount("", &dssp_html_controller::index);
	}

	void index(const zeep::http::request& request, const zeep::http::scope& scope, zeep::http::reply& reply)
	{
		get_template_processor().create_reply_from_template("index", scope, reply);
	}
};

// --------------------------------------------------------------------

class dssp_rest_controller : public zeep::http::rest_controller
{
  public:
	dssp_rest_controller()
		: zeep::http::rest_controller("")
	{
		// map_post_request("dssp", &dssp_rest_controller::calculate, "data");
	}
};

// --------------------------------------------------------------------

int main(int argc, char *argv[])
{
	using namespace std::literals;
	namespace zh = zeep::http;

	cif::compound_factory::init(true);

	int result = 0;

	auto &config = mcfp::config::instance();

	config.init("dsspd [options] start|stop|status|reload",
		mcfp::make_option("help,h", "Display help message"),
		mcfp::make_option("version", "Print version"),
		mcfp::make_option("verbose,v", "verbose output"),

		mcfp::make_option<std::string>("address", "0.0.0.0", "External address"),
		mcfp::make_option<uint16_t>("port", 10351, "Port to listen to"),
		mcfp::make_option<std::string>("user,u", "www-data", "User to run the daemon"),
		mcfp::make_option<std::string>("context", "", "Root context for web server"),

		mcfp::make_option("no-daemon,F", "Do not fork into background"),
		
		mcfp::make_option<std::string>("config", "Config file to use"));

	std::error_code ec;

	config.parse(argc, argv, ec);
	if (ec != std::errc())
	{
		std::cerr << "Error parsing command line arguments: " << ec.message() << std::endl
				  << std::endl
				  << config << std::endl;
		exit(1);
	}

	config.parse_config_file("config", "dsspd.conf", { ".", "/etc" }, ec);
	if (ec != std::errc())
	{
		std::cerr << "Error parsing config file: " << ec.message() << std::endl;
		exit(1);
	}

	// --------------------------------------------------------------------

	if (config.has("version"))
	{
		write_version_string(std::cout, config.has("verbose"));
		exit(0);
	}

	if (config.has("help"))
	{
		std::cout << config << std::endl;
		exit(0);
	}
	
	if (config.operands().empty())
	{
		std::cerr << "Missing command, should be one of start, stop, status or reload" << std::endl;
		exit(1);
	}

	cif::VERBOSE = config.count("verbose");

	std::string user = config.get<std::string>("user");
	std::string address = config.get<std::string>("address");
	uint16_t port = config.get<uint16_t>("port");

	zh::daemon server([&]()
	{
		auto s = new zeep::http::server();

#ifndef NDEBUG
		s->set_template_processor(new zeep::http::file_based_html_template_processor("docroot"));
#else
		s->set_template_processor(new zeep::http::rsrc_based_html_template_processor());
#endif
		s->add_controller(new dssp_rest_controller());
		s->add_controller(new dssp_html_controller(config.get("context")));
		return s;
	}, kProjectName );

	std::string command = config.operands().front();

	if (command == "start")
	{
		std::cout << "starting server at http://" << address << ':' << port << '/' << std::endl;

		if (config.has("no-daemon"))
			result = server.run_foreground(address, port);
		else
			result = server.start(address, port, 2, 2, user);
	}
	else if (command == "stop")
		result = server.stop();
	else if (command == "status")
		result = server.status();
	else if (command == "reload")
		result = server.reload();
	else
	{
		std::cerr << "Invalid command" << std::endl;
		result = 1;
	}

	return result;
}