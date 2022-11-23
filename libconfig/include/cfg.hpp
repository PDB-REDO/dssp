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

#pragma once

/// \file cfg.hpp
/// This header-only library contains code to parse argc/argv and store the
/// values contained therein into a singleton object.

#include <cassert>
#include <cstring>

#include <any>
#include <charconv>
#include <deque>
#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <type_traits>
#include <vector>

#include <cfg/text.hpp>
#include <cfg/utilities.hpp>

namespace cfg
{

// we use the new system_error stuff.

enum class config_error
{
	unknown_option = 1,
	option_does_not_accept_argument,
	missing_argument_for_option,
	option_not_specified,
	invalid_config_file
};

class config_category_impl : public std::error_category
{
  public:
	const char *name() const noexcept override
	{
		return "configuration";
	}

	std::string message(int ev) const override
	{
		switch (static_cast<config_error>(ev))
		{
			case config_error::unknown_option:
				return "unknown option";
			case config_error::option_does_not_accept_argument:
				return "option does not accept argument";
			case config_error::missing_argument_for_option:
				return "missing argument for option";
			case config_error::option_not_specified:
				return "option was not specified";
			case config_error::invalid_config_file:
				return "config file contains a syntax error";
			default:
				assert(false);
				return "unknown error code";
		}
	}

	bool equivalent(const std::error_code &code, int condition) const noexcept override
	{
		return false;
	}
};

inline std::error_category &config_category()
{
	static config_category_impl instance;
	return instance;
}

inline std::error_code make_error_code(config_error e)
{
	return std::error_code(static_cast<int>(e), config_category());
}

inline std::error_condition make_error_condition(config_error e)
{
	return std::error_condition(static_cast<int>(e), config_category());
}

// --------------------------------------------------------------------
// Some template wizardry to detect containers, needed to have special
// handling of options that can be repeated.

template <typename T>
using iterator_t = typename T::iterator;

template <typename T>
using value_type_t = typename T::value_type;

template <typename T>
using std_string_npos_t = decltype(T::npos);

template <typename T, typename = void>
struct is_container_type : std::false_type
{
};

template <typename T>
struct is_container_type<T,
	std::enable_if_t<
		std::experimental::is_detected_v<value_type_t, T> and
		std::experimental::is_detected_v<iterator_t, T> and
		not std::experimental::is_detected_v<std_string_npos_t, T>>> : std::true_type
{
};

template <typename T>
inline constexpr bool is_container_type_v = is_container_type<T>::value;

// --------------------------------------------------------------------
// The options classes

namespace detail
{
	// The option traits classes are used to convert from the string-based
	// command line argument to the type that should be stored.
	// In fact, here is where the command line arguments are checked for
	// proper formatting.
	template <typename T, typename = void>
	struct option_traits;

	template <typename T>
	struct option_traits<T, typename std::enable_if_t<std::is_arithmetic_v<T>>>
	{
		using value_type = T;

		static value_type set_value(std::string_view argument, std::error_code &ec)
		{
			value_type value{};
			auto r = charconv<value_type>::from_chars(argument.begin(), argument.end(), value);
			if (r.ec != std::errc())
				ec = std::make_error_code(r.ec);
			return value;
		}

		static std::string to_string(const T &value)
		{
			char b[32];
			auto r = charconv<value_type>::to_chars(b, b + sizeof(b), value);
			if (r.ec != std::errc())
				throw std::system_error(std::make_error_code(r.ec));
			return { b, r.ptr };
		}
	};

	template <>
	struct option_traits<std::filesystem::path>
	{
		using value_type = std::filesystem::path;

		static value_type set_value(std::string_view argument, std::error_code &ec)
		{
			return value_type{ argument };
		}

		static std::string to_string(const std::filesystem::path &value)
		{
			return value.string();
		}
	};

	template <typename T>
	struct option_traits<T, typename std::enable_if_t<not std::is_arithmetic_v<T> and std::is_assignable_v<std::string, T>>>
	{
		using value_type = std::string;

		static value_type set_value(std::string_view argument, std::error_code &ec)
		{
			return value_type{ argument };
		}

		static std::string to_string(const T &value)
		{
			return { value };
		}
	};

	// The Options. The reason to have this weird constructing of
	// polymorphic options based on templates is to have a very
	// simple interface. The disadvantage is that the options have
	// to be copied during the construction of the config object.

	struct option_base
	{
		std::string m_name;        ///< The long argument name
		std::string m_desc;        ///< The description of the argument
		char m_short_name;         ///< The single character name of the argument, can be zero
		bool m_is_flag = true,     ///< When true, this option does not allow arguments
			m_has_default = false, ///< When true, this option has a default value.
			m_multi = false,       ///< When true, this option allows mulitple values.
			m_hidden;              ///< When true, this option is hidden from the help text
		int m_seen = 0;            ///< How often the option was seen on the command line

		option_base(const option_base &rhs) = default;

		option_base(std::string_view name, std::string_view desc, bool hidden)
			: m_name(name)
			, m_desc(desc)
			, m_short_name(0)
			, m_hidden(hidden)
		{
			if (m_name.length() == 1)
				m_short_name = m_name.front();
			else if (m_name.length() > 2 and m_name[m_name.length() - 2] == ',')
			{
				m_short_name = m_name.back();
				m_name.erase(m_name.end() - 2, m_name.end());
			}
		}

		virtual ~option_base() = default;

		virtual void set_value(std::string_view value, std::error_code &ec)
		{
			assert(false);
		}

		virtual std::any get_value() const
		{
			return {};
		}

		virtual std::string get_default_value() const
		{
			return {};
		}

		uint32_t width() const
		{
			uint32_t result = m_name.length();
			if (result <= 1)
				result = 2;
			else if (m_short_name != 0)
				result += 7;
			if (not m_is_flag)
			{
				result += 4;
				if (m_has_default)
					result += 4 + get_default_value().length();
			}
			return result + 6;
		}

		void write(std::ostream &os, uint32_t width) const
		{
			uint32_t w2 = 2;
			os << "  ";
			if (m_short_name)
			{
				os << '-' << m_short_name;
				w2 += 2;
				if (m_name.length() > 1)
				{
					os << " [ --" << m_name << " ]";
					w2 += 7 + m_name.length();
				}
			}
			else
			{
				os << "--" << m_name;
				w2 += 2 + m_name.length();
			}

			if (not m_is_flag)
			{
				os << " arg";
				w2 += 4;

				if (m_has_default)
				{
					auto default_value = get_default_value();
					os << " (=" << default_value << ')';
					w2 += 4 + default_value.length();
				}
			}

			auto leading_spaces = width;
			if (w2 + 2 > width)
				os << std::endl;
			else
				leading_spaces = width - w2;

			word_wrapper ww(m_desc, get_terminal_width() - width);
			for (auto line : ww)
			{
				os << std::string(leading_spaces, ' ') << line << std::endl;
				leading_spaces = width;
			}
		}
	};

	template <typename T>
	struct option : public option_base
	{
		using traits_type = option_traits<T>;
		using value_type = typename option_traits<T>::value_type;

		std::optional<value_type> m_value;

		option(const option &rhs) = default;

		option(std::string_view name, std::string_view desc, bool hidden)
			: option_base(name, desc, hidden)
		{
			m_is_flag = false;
		}

		option(std::string_view name, const value_type &default_value, std::string_view desc, bool hidden)
			: option(name, desc, hidden)
		{
			m_has_default = true;
			m_value = default_value;
		}

		void set_value(std::string_view argument, std::error_code &ec) override
		{
			m_value = traits_type::set_value(argument, ec);
		}

		std::any get_value() const override
		{
			std::any result;
			if (m_value)
				result = *m_value;
			return result;
		}

		std::string get_default_value() const override
		{
			if constexpr (std::is_same_v<value_type, std::string>)
				return *m_value;
			else
				return traits_type::to_string(*m_value);
		}
	};

	template <typename T>
	struct multiple_option : public option_base
	{
		using value_type = typename T::value_type;
		using traits_type = option_traits<value_type>;

		std::vector<value_type> m_values;

		multiple_option(const multiple_option &rhs) = default;

		multiple_option(std::string_view name, std::string_view desc, bool hidden)
			: option_base(name, desc, hidden)
		{
			m_is_flag = false;
			m_multi = true;
		}

		void set_value(std::string_view argument, std::error_code &ec) override
		{
			m_values.emplace_back(traits_type::set_value(argument, ec));
		}

		std::any get_value() const override
		{
			return { m_values };
		}
	};

	template <>
	struct option<void> : public option_base
	{
		option(const option &rhs) = default;

		option(std::string_view name, std::string_view desc, bool hidden)
			: option_base(name, desc, hidden)
		{
		}
	};

} // namespace detail

// --------------------------------------------------------------------
/// \brief A singleton class. Use config::instance to create an instance

class config
{
  public:
	using option_base = detail::option_base;

	void set_usage(std::string_view usage)
	{
		m_usage = usage;
	}

	/// \brief Initialise a config instance with a \a usage message and a set of \a options
	template <typename... Options>
	void init(std::string_view usage, Options... options)
	{
		m_usage = usage;
		m_impl.reset(new config_impl(std::forward<Options>(options)...));
	}

	void set_ignore_unknown(bool ignore_unknown)
	{
		m_ignore_unknown = ignore_unknown;
	}

	static config &instance()
	{
		static std::unique_ptr<config> s_instance;
		if (not s_instance)
			s_instance.reset(new config);
		return *s_instance;
	}

	bool has(std::string_view name) const
	{
		auto opt = m_impl->get_option(name);
		return opt != nullptr and (opt->m_seen > 0 or opt->m_has_default);
	}

	int count(std::string_view name) const
	{
		auto opt = m_impl->get_option(name);
		return opt ? opt->m_seen : 0;
	}

	template <typename T>
	auto get(std::string_view name) const
	{
		auto opt = m_impl->get_option(name);
		if (opt == nullptr)
			throw std::system_error(make_error_code(config_error::unknown_option), std::string{ name });

		std::any value = opt->get_value();

		if (not value.has_value())
			throw std::system_error(make_error_code(config_error::option_not_specified), std::string{ name });

		return std::any_cast<T>(value);
	}

	const std::vector<std::string> &operands() const
	{
		return m_impl->m_operands;
	}

	friend std::ostream &operator<<(std::ostream &os, const config &conf)
	{
		uint32_t terminal_width = get_terminal_width();

		if (not conf.m_usage.empty())
			os << conf.m_usage << std::endl;

		uint32_t options_width = conf.m_impl->get_option_width();

		if (options_width > terminal_width / 2)
			options_width = terminal_width / 2;

		conf.m_impl->write(os, options_width);

		return os;
	}

	// --------------------------------------------------------------------

	void parse(int argc, const char *const argv[])
	{
		std::error_code ec;
		parse(argc, argv, ec);
		if (ec)
			throw std::system_error(ec);
	}

	void parse_config_file(std::string_view config_option, std::string_view config_file_name,
		std::initializer_list<std::string_view> search_dirs)
	{
		std::error_code ec;
		parse_config_file(config_option, config_file_name, search_dirs, ec);
		if (ec)
			throw std::system_error(ec);
	}

	void parse_config_file(std::string_view config_option, std::string_view config_file_name,
		std::initializer_list<std::string_view> search_dirs, std::error_code &ec)
	{
		std::string file_name{ config_file_name };
		if (has(config_option))
			file_name = get<std::string>(config_option);

		for (std::filesystem::path dir : search_dirs)
		{
			std::ifstream file(dir / file_name);

			if (not file.is_open())
				continue;

			parse_config_file(file, ec);
			break;
		}
	}

	void parse_config_file(const std::filesystem::path &file, std::error_code &ec)
	{
		std::ifstream is(file);
		if (is.is_open())
			parse_config_file(is, ec);
	}

	static bool is_name_char(int ch)
	{
		return std::isalnum(ch) or ch == '_' or ch == '-';
	}

	static constexpr bool is_eoln(int ch)
	{
		return ch == '\n' or ch == '\r' or ch == std::char_traits<char>::eof();
	}

	void parse_config_file(std::istream &is, std::error_code &ec)
	{
		auto &buffer = *is.rdbuf();

		enum class State
		{
			NAME_START,
			COMMENT,
			NAME,
			ASSIGN,
			VALUE_START,
			VALUE
		} state = State::NAME_START;

		std::string name, value;

		for (;;)
		{
			auto ch = buffer.sbumpc();

			switch (state)
			{
				case State::NAME_START:
					if (is_name_char(ch))
					{
						name = { static_cast<char>(ch) };
						value.clear();
						state = State::NAME;
					}
					else if (ch == '#')
						state = State::COMMENT;
					else if (ch != ' ' and ch != '\t' and not is_eoln(ch))
						ec = make_error_code(config_error::invalid_config_file);
					break;

				case State::COMMENT:
					if (is_eoln(ch))
						state = State::NAME_START;
					break;

				case State::NAME:
					if (is_name_char(ch))
						name.insert(name.end(), static_cast<char>(ch));
					else if	(is_eoln(ch))
					{
						auto opt = m_impl->get_option(name);

						if (opt == nullptr)
						{
							if (not m_ignore_unknown)
								ec = make_error_code(config_error::unknown_option);
						}
						else if (not opt->m_is_flag)
							ec = make_error_code(config_error::missing_argument_for_option);
						else
							++opt->m_seen;

						state = State::NAME_START;
					}
					else
					{
						buffer.sungetc();
						state = State::ASSIGN;
					}
					break;

				case State::ASSIGN:
					if (ch == '=')
						state = State::VALUE_START;
					else if (is_eoln(ch))
					{
						auto opt = m_impl->get_option(name);

						if (opt == nullptr)
						{
							if (not m_ignore_unknown)
								ec = make_error_code(config_error::unknown_option);
						}
						else if (not opt->m_is_flag)
							ec = make_error_code(config_error::missing_argument_for_option);
						else
							++opt->m_seen;

						state = State::NAME_START;
					}
					else if (ch != ' ' and ch != '\t')
						ec = make_error_code(config_error::invalid_config_file);
					break;

				case State::VALUE_START:
				case State::VALUE:
					if (is_eoln(ch))
					{
						auto opt = m_impl->get_option(name);

						if (opt == nullptr)
						{
							if (not m_ignore_unknown)
								ec = make_error_code(config_error::unknown_option);
						}
						else if (opt->m_is_flag)
							ec = make_error_code(config_error::option_does_not_accept_argument);
						else if (not value.empty() and (opt->m_seen == 0 or opt->m_multi))
						{
							opt->set_value(value, ec);
							++opt->m_seen;
						}

						state = State::NAME_START;
					}
					else if (state == State::VALUE)
						value.insert(value.end(), static_cast<char>(ch));
					else if (ch != ' ' and ch != '\t')
					{
						value = { static_cast<char>(ch) };
						state = State::VALUE;
					}
					break;
			}

			if (ec or ch == std::char_traits<char>::eof())
				break;
		}
	}

	void parse(int argc, const char *const argv[], std::error_code &ec)
	{
		using namespace std::literals;

		enum class State
		{
			options,
			operands
		} state = State::options;

		for (int i = 1; i < argc and not ec; ++i)
		{
			const char *arg = argv[i];

			if (arg == nullptr) // should not happen
				break;

			if (state == State::options)
			{
				if (*arg != '-') // according to POSIX this is the end of options, start operands
				                 // state = State::operands;
				{                // however, people nowadays expect to be able to mix operands and options
					m_impl->m_operands.emplace_back(arg);
					continue;
				}
				else if (arg[1] == '-' and arg[2] == 0)
				{
					state = State::operands;
					continue;
				}
			}

			if (state == State::operands)
			{
				m_impl->m_operands.emplace_back(arg);
				continue;
			}

			option_base *opt = nullptr;
			std::string_view opt_arg;

			assert(*arg == '-');
			++arg;

			if (*arg == '-') // double --, start of new argument
			{
				++arg;

				assert(*arg != 0); // this should not happen, as it was checked for before

				std::string_view s_arg(arg);
				std::string_view::size_type p = s_arg.find('=');

				if (p != std::string_view::npos)
				{
					opt_arg = s_arg.substr(p + 1);
					s_arg = s_arg.substr(0, p);
				}

				opt = m_impl->get_option(s_arg);
				if (opt == nullptr)
				{
					if (not m_ignore_unknown)
						ec = make_error_code(config_error::unknown_option);
					continue;
				}

				if (opt->m_is_flag)
				{
					if (not opt_arg.empty())
						ec = make_error_code(config_error::option_does_not_accept_argument);

					++opt->m_seen;
					continue;
				}

				++opt->m_seen;
			}
			else // single character options
			{
				bool expect_option_argument = false;

				while (*arg != 0 and not ec)
				{
					opt = m_impl->get_option(*arg++);

					if (opt == nullptr)
					{
						if (not m_ignore_unknown)
							ec = make_error_code(config_error::unknown_option);
						continue;
					}

					++opt->m_seen;
					if (opt->m_is_flag)
						continue;

					opt_arg = arg;
					expect_option_argument = true;
					break;
				}

				if (not expect_option_argument)
					continue;
			}

			if (opt_arg.empty() and i + 1 < argc) // So, the = character was not present, the next arg must be the option argument
			{
				++i;
				opt_arg = argv[i];
			}

			if (opt_arg.empty())
				ec = make_error_code(config_error::missing_argument_for_option);
			else
				opt->set_value(opt_arg, ec);
		}
	}

  private:
	config() = default;
	config(const config &) = delete;
	config &operator=(const config &) = delete;

	struct config_impl_base
	{
		virtual ~config_impl_base() = default;

		virtual option_base *get_option(std::string_view name) = 0;
		virtual option_base *get_option(char short_name) = 0;

		virtual uint32_t get_option_width() const = 0;
		virtual void write(std::ostream &os, uint32_t width) const = 0;

		std::vector<std::string> m_operands;
	};

	template <typename... Options>
	struct config_impl : public config_impl_base
	{
		static constexpr size_t N = sizeof...(Options);

		config_impl(Options... options)
			: m_options(std::forward<Options>(options)...)
		{
		}

		option_base *get_option(std::string_view name) override
		{
			return get_option<0>(name);
		}

		template<size_t Ix>
		option_base *get_option([[maybe_unused]] std::string_view name)
		{
			if constexpr (Ix == N)
				return nullptr;
			else
			{
				option_base &opt = std::get<Ix>(m_options);
				return (opt.m_name == name) ? &opt : get_option<Ix + 1>(name);
			}
		}

		option_base *get_option(char short_name) override
		{
			return get_option<0>(short_name);
		}

		template<size_t Ix>
		option_base *get_option([[maybe_unused]] char short_name)
		{
			if constexpr (Ix == N)
				return nullptr;
			else
			{
				option_base &opt = std::get<Ix>(m_options);
				return (opt.m_short_name == short_name) ? &opt : get_option<Ix + 1>(short_name);
			}
		}

		virtual uint32_t get_option_width() const override
		{
			return std::apply([](Options const& ...opts) {
				uint32_t width = 0;
				((width = std::max(width, opts.width())), ...);
				return width;
			}, m_options);
		}

		virtual void write(std::ostream &os, uint32_t width) const override
		{
			std::apply([&os,width](Options const& ...opts) {
				(opts.write(os, width), ...);
			}, m_options);
		}

		std::tuple<Options...> m_options;
	};

	std::unique_ptr<config_impl_base> m_impl;
	bool m_ignore_unknown = false;
	std::string m_usage;
};

// --------------------------------------------------------------------

template <typename T = void, std::enable_if_t<not is_container_type_v<T>, int> = 0>
auto make_option(std::string_view name, std::string_view description)
{
	return detail::option<T>(name, description, false);
}

template <typename T = void, std::enable_if_t<not is_container_type_v<T>, int> = 0>
auto make_hidden_option(std::string_view name, std::string_view description)
{
	return detail::option<T>(name, description, true);
}

template <typename T, std::enable_if_t<not is_container_type_v<T>, int> = 0>
auto make_option(std::string_view name, const T &v, std::string_view description)
{
	return detail::option<T>(name, v, description, false);
}

template <typename T, std::enable_if_t<not is_container_type_v<T>, int> = 0>
auto make_hidden_option(std::string_view name, const T &v, std::string_view description)
{
	return detail::option<T>(name, v, description, true);
}

template <typename T, std::enable_if_t<is_container_type_v<T>, int> = 0>
auto make_option(std::string_view name, std::string_view description)
{
	return detail::multiple_option<T>(name, description, false);
}

template <typename T, std::enable_if_t<is_container_type_v<T>, int> = 0>
auto make_hidden_option(std::string_view name, std::string_view description)
{
	return detail::option<T>(name, description, true);
}

} // namespace cfg

namespace std
{

template <>
struct is_error_condition_enum<cfg::config_error>
	: public true_type
{
};

} // namespace std