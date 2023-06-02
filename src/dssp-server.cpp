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

#include <cif++.hpp>

#include <gxrio.hpp>
#include <mcfp.hpp>

#include <zeep/http/daemon.hpp>
#include <zeep/http/html-controller.hpp>
#include <zeep/http/rest-controller.hpp>
#include <zeep/streambuf.hpp>

// --------------------------------------------------------------------

class dssp_html_controller : public zeep::http::html_controller
{
  public:
	dssp_html_controller()
		: zeep::http::html_controller()
	{
		mount("{css,scripts,fonts,images,favicon}/", &dssp_html_controller::handle_file);
		mount("{favicon.ico,browserconfig.xml,manifest.json}", &dssp_html_controller::handle_file);
		map_get("", "index");
		map_get("about", "about");
		map_get("download", "download");
	}
};

// --------------------------------------------------------------------

class dssp_rest_controller : public zeep::http::rest_controller
{
  public:
	dssp_rest_controller()
		: zeep::http::rest_controller("")
	{
		map_post_request("do", &dssp_rest_controller::work, "data", "format");
	}

	zeep::http::reply work(const zeep::http::file_param &coordinates, std::optional<std::string> format);
};

zeep::http::reply dssp_rest_controller::work(const zeep::http::file_param &coordinates, std::optional<std::string> format)
{
	zeep::char_streambuf sb(coordinates.data, coordinates.length);

	gxrio::istream in(&sb);

	cif::file f = cif::pdb::read(in);
	if (f.empty())
		throw std::runtime_error("Invalid input file, is it empty?");

	// --------------------------------------------------------------------

	short pp_stretch = 3; //minPPStretch.value_or(3);

	std::string fmt = format.value_or("mmcif");

	dssp dssp(f.front(), 1, pp_stretch, true);

	std::ostringstream os;

	if (fmt == "dssp")
		writeDSSP(dssp, os);
	else
		annotateDSSP(f.front(), dssp, true, true, os);

	// --------------------------------------------------------------------
	
	std::string name = f.front().name();
	if (fmt == "dssp")
		name += ".dssp";
	else
		name += ".cif";

	zeep::http::reply rep{ zeep::http::ok };
	rep.set_content(os.str(), "text/plain");
	rep.set_header("content-disposition", "attachement; filename = \"" + name + '"');

	return rep;
}

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
	if (ec)
	{
		std::cerr << "Error parsing command line arguments: " << ec.message() << std::endl
				  << std::endl
				  << config << std::endl;
		exit(1);
	}

	config.parse_config_file("config", "dsspd.conf", { ".", "/etc" }, ec);
	if (ec)
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

	zeep::http::daemon server([&, context = config.get("context")]()
		{
		auto s = new zeep::http::server();

#ifndef NDEBUG
		s->set_template_processor(new zeep::http::file_based_html_template_processor("docroot"));
#else
		s->set_template_processor(new zeep::http::rsrc_based_html_template_processor());
#endif
		s->add_controller(new dssp_rest_controller());
		s->add_controller(new dssp_html_controller());

		s->set_context_name(context);

		return s; },
		kProjectName);

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