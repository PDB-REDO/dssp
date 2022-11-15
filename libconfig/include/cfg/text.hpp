/*-
 * SPDX-License-Identifier: BSD-2-Clause
 *
 * Copyright (c) 2022 Maarten L. hekkelman
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

#include <algorithm>
#include <charconv>
#include <cmath>
#include <experimental/type_traits>

namespace cfg
{

template <typename T>
struct my_charconv
{
	using value_type = T;

	static std::from_chars_result from_chars(const char *first, const char *last, value_type &value)
	{
		std::from_chars_result result{ first, {} };

		enum State
		{
			IntegerSign,
			Integer,
			Fraction,
			ExponentSign,
			Exponent
		} state = IntegerSign;
		int sign = 1;
		unsigned long long vi = 0;
		long double f = 1;
		int exponent_sign = 1;
		int exponent = 0;
		bool done = false;

		while (not done and result.ec == std::errc())
		{
			char ch = result.ptr != last ? *result.ptr : 0;
			++result.ptr;

			switch (state)
			{
				case IntegerSign:
					if (ch == '-')
					{
						sign = -1;
						state = Integer;
					}
					else if (ch == '+')
						state = Integer;
					else if (ch >= '0' and ch <= '9')
					{
						vi = ch - '0';
						state = Integer;
					}
					else if (ch == '.')
						state = Fraction;
					else
						result.ec = std::errc::invalid_argument;
					break;

				case Integer:
					if (ch >= '0' and ch <= '9')
						vi = 10 * vi + (ch - '0');
					else if (ch == 'e' or ch == 'E')
						state = ExponentSign;
					else if (ch == '.')
						state = Fraction;
					else
					{
						done = true;
						--result.ptr;
					}
					break;

				case Fraction:
					if (ch >= '0' and ch <= '9')
					{
						vi = 10 * vi + (ch - '0');
						f /= 10;
					}
					else if (ch == 'e' or ch == 'E')
						state = ExponentSign;
					else
					{
						done = true;
						--result.ptr;
					}
					break;

				case ExponentSign:
					if (ch == '-')
					{
						exponent_sign = -1;
						state = Exponent;
					}
					else if (ch == '+')
						state = Exponent;
					else if (ch >= '0' and ch <= '9')
					{
						exponent = ch - '0';
						state = Exponent;
					}
					else
						result.ec = std::errc::invalid_argument;
					break;

				case Exponent:
					if (ch >= '0' and ch <= '9')
						exponent = 10 * exponent + (ch - '0');
					else
					{
						done = true;
						--result.ptr;
					}
					break;
			}
		}

		if (result.ec == std::errc())
		{
			long double v = f * vi * sign;
			if (exponent != 0)
				v *= std::pow(10, exponent * exponent_sign);

			if (std::isnan(v))
				result.ec = std::errc::invalid_argument;
			else if (std::abs(v) > std::numeric_limits<value_type>::max())
				result.ec = std::errc::result_out_of_range;

			value = static_cast<value_type>(v);
		}

		return result;
	}

	template <typename Iterator, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
	static std::to_chars_result to_chars(Iterator first, Iterator last, const T &value)
	{
		int size = last - first;
		int r;

		if constexpr (std::is_same_v<T, long double>)
			r = snprintf(first, last - first, "%lg", value);
		else
			r = snprintf(first, last - first, "%g", value);

		std::to_chars_result result;
		if (r < 0 or r >= size)
			result = { first, std::errc::value_too_large };
		else
			result = { first + r, std::errc() };

		return result;
	}
};

template <typename T>
struct std_charconv
{
	static std::from_chars_result from_chars(const char *a, const char *b, T &d)
	{
		return std::from_chars(a, b, d);
	}

	template<typename Iterator>
	static std::to_chars_result to_chars(Iterator a, Iterator b, const T &value)
	{
		return std::to_chars(a, b, value);
	}
};

template <typename T>
using from_chars_function = decltype(std::from_chars(std::declval<const char *>(), std::declval<const char *>(), std::declval<T &>()));

template <typename T>
using charconv = typename std::conditional_t<std::experimental::is_detected_v<from_chars_function, T>, std_charconv<T>, my_charconv<T>>;

// --------------------------------------------------------------------
/// Simplified line breaking code taken from a decent text editor.
/// In this case, simplified means it only supports ASCII.
/// The algorithm uses dynamic programming to find the optimal
/// separation in lines.

class word_wrapper : public std::vector<std::string_view>
{
  public:
	word_wrapper(std::string_view text, size_t width)
		: m_width(width)
	{
		std::string_view::size_type line_start = 0, line_end = text.find('\n');
		
		for (;;)
		{
			auto line = text.substr(line_start, line_end - line_start);
			if (line.empty())
				this->push_back(line);
			else
			{
				auto lines = wrap_line(line);
				this->insert(this->end(), lines.begin(), lines.end());
			}

			if (line_end == std::string_view::npos)
				break;

			line_start = line_end + 1;
			line_end = text.find('\n', line_start);
		}
	}

  private:
	std::vector<std::string_view> wrap_line(std::string_view line)
	{
		std::vector<std::string_view> result;
		std::vector<size_t> offsets = { 0 };

		auto b = line.begin();
		while (b != line.end())
		{
			auto e = next_line_break(b, line.end());

			offsets.push_back(e - line.begin());

			b = e;
		}

		size_t count = offsets.size() - 1;

		std::vector<size_t> minima(count + 1, std::numeric_limits<size_t>::max());
		minima[0] = 0;
		std::vector<size_t> breaks(count + 1, 0);

		for (size_t i = 0; i < count; ++i)
		{
			size_t j = i + 1;
			while (j <= count)
			{
				size_t w = offsets[j] - offsets[i];

				if (w > m_width)
					break;

				while (w > 0 and std::isspace(line[offsets[i] + w - 1]))
					--w;

				size_t cost = minima[i];
				if (j < count) // last line may be shorter
					cost += (m_width - w) * (m_width - w);

				if (cost < minima[j])
				{
					minima[j] = cost;
					breaks[j] = i;
				}

				++j;
			}
		}

		size_t j = count;
		while (j > 0)
		{
			size_t i = breaks[j];
			result.push_back(line.substr(offsets[i], offsets[j] - offsets[i]));
			j = i;
		}

		reverse(result.begin(), result.end());

		return result;
	}

	std::string_view::const_iterator next_line_break(std::string_view::const_iterator text, std::string_view::const_iterator end)
	{
		if (text == end)
			return text;

		enum LineBreakClass
		{
			OP, // OpenPunctuation,
			CL, // ClosePunctuation,
			CP, // CloseParenthesis,
			QU, // Quotation,
			// GL, // NonBreaking,
			// NS, // Nonstarter,
			EX, // Exlamation,
			SY, // SymbolAllowingBreakAfter,
			IS, // InfixNumericSeparator,
			PR, // PrefixNumeric,
			PO, // PostfixNumeric,
			NU, // Numeric,
			AL, // Alphabetic,
			// ID, // Ideographic,
			// IN, // Inseperable,
			HY, // Hyphen,
			BA, // BreakAfter,
			// BB, // BreakBefor,
			// B2, // BreakOpportunityBeforeAndAfter,
			// ZW, // ZeroWidthSpace,
			CM, // CombiningMark,
			WJ, // WordJoiner,
			// H2, // HangulLVSyllable,
			// H3, // HangulLVTSyllable,
			// JL, // HangulLJamo,
			// JV, // HangulVJamo,
			// JT, // HangulTJamo,

			MB, // MandatoryBreak,
			CR, // CarriageReturn,
			LF, // LineFeed,
			// NL, // NextLine,
			// SU, // Surrogate,
			SP, // Space,
			// CB, // ContigentBreakOpportunity,
			// AM, // Ambiguous,
			// CC, // ComplexContext,
			UN, // Unknown
		};

		static const LineBreakClass kASCII_LineBreakTable[128] = {
			CM, CM, CM, CM, CM, CM, CM, CM,
			CM, BA, LF, MB, MB, CR, CM, CM,
			CM, CM, CM, CM, CM, CM, CM, CM,
			CM, CM, CM, CM, CM, CM, CM, CM,
			SP, EX, QU, AL, PR, PO, AL, QU,
			OP, CP, AL, PR, IS, HY, IS, SY,
			NU, NU, NU, NU, NU, NU, NU, NU,
			NU, NU, IS, IS, AL, AL, AL, EX,
			AL, AL, AL, AL, AL, AL, AL, AL,
			AL, AL, AL, AL, AL, AL, AL, AL,
			AL, AL, AL, AL, AL, AL, AL, AL,
			AL, AL, AL, OP, PR, CP, AL, AL,
			AL, AL, AL, AL, AL, AL, AL, AL,
			AL, AL, AL, AL, AL, AL, AL, AL,
			AL, AL, AL, AL, AL, AL, AL, AL,
			AL, AL, AL, OP, BA, CL, AL, CM
		};

		enum BreakAction
		{
			DBK = 0, // direct break 	(blank in table)
			IBK,     // indirect break	(% in table)
			PBK,     // prohibited break (^ in table)
			CIB,     // combining indirect break
			CPB      // combining prohibited break
		};

		// const breakAction brkTable[27][27] = {
		static const BreakAction brkTable[15][15] = {
			//         OP   CL   CP   QU   /*GL   NS  */ EX   SY   IS   PR   PO   NU   AL   /*ID   IN  */ HY   BA   /*BB   B2   ZW  */ CM   WJ   /*H2   H3   JL   JV   JT */
			/* OP */ { PBK, PBK, PBK, PBK, /*PBK, PBK,*/ PBK, PBK, PBK, PBK, PBK, PBK, PBK, /*PBK, PBK,*/ PBK, PBK, /*PBK, PBK, PBK,*/ CPB, PBK, /*PBK, PBK, PBK, PBK, PBK*/ },
			/* CL */ { DBK, PBK, PBK, IBK, /*IBK, PBK,*/ PBK, PBK, PBK, IBK, IBK, DBK, DBK, /*DBK, DBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* CP */ { DBK, PBK, PBK, IBK, /*IBK, PBK,*/ PBK, PBK, PBK, IBK, IBK, IBK, IBK, /*DBK, DBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* QU */ { PBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, IBK, IBK, IBK, IBK, /*IBK, IBK,*/ IBK, IBK, /*IBK, IBK, PBK,*/ CIB, PBK, /*IBK, IBK, IBK, IBK, IBK*/ },
			// /* GL */ { IBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, IBK, IBK, IBK, IBK, /*IBK, IBK,*/ IBK, IBK, /*IBK, IBK, PBK,*/ CIB, PBK, /*IBK, IBK, IBK, IBK, IBK*/ },
			// /* NS */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, DBK, DBK, /*DBK, DBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* EX */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, DBK, DBK, /*DBK, DBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* SY */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, IBK, DBK, /*DBK, DBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* IS */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, IBK, IBK, /*DBK, DBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* PR */ { IBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, IBK, IBK, /*IBK, DBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*IBK, IBK, IBK, IBK, IBK*/ },
			/* PO */ { IBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, IBK, IBK, /*DBK, DBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* NU */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, IBK, IBK, IBK, IBK, /*DBK, IBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* AL */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, IBK, IBK, /*DBK, IBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			// /* ID */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, IBK, DBK, DBK, /*DBK, IBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			// /* IN */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, DBK, DBK, /*DBK, IBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* HY */ { DBK, PBK, PBK, IBK, /*DBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, IBK, DBK, /*DBK, DBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* BA */ { DBK, PBK, PBK, IBK, /*DBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, DBK, DBK, /*DBK, DBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			// /* BB */ { IBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, IBK, IBK, IBK, IBK, /*IBK, IBK,*/ IBK, IBK, /*IBK, IBK, PBK,*/ CIB, PBK, /*IBK, IBK, IBK, IBK, IBK*/ },
			// /* B2 */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, DBK, DBK, /*DBK, DBK,*/ IBK, IBK, /*DBK, PBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			// /* ZW */ { DBK, DBK, DBK, DBK, /*DBK, DBK,*/ DBK, DBK, DBK, DBK, DBK, DBK, DBK, /*DBK, DBK,*/ DBK, DBK, /*DBK, DBK, PBK,*/ DBK, DBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* CM */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, DBK, IBK, IBK, /*DBK, IBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, DBK*/ },
			/* WJ */ { IBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, IBK, IBK, IBK, IBK, /*IBK, IBK,*/ IBK, IBK, /*IBK, IBK, PBK,*/ CIB, PBK, /*IBK, IBK, IBK, IBK, IBK*/ },
			// /* H2 */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, IBK, DBK, DBK, /*DBK, IBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, IBK, IBK*/ },
			// /* H3 */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, IBK, DBK, DBK, /*DBK, IBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, IBK*/ },
			// /* JL */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, IBK, DBK, DBK, /*DBK, IBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*IBK, IBK, IBK, IBK, DBK*/ },
			// /* JV */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, IBK, DBK, DBK, /*DBK, IBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, IBK, IBK*/ },
			// /* JT */ { DBK, PBK, PBK, IBK, /*IBK, IBK,*/ PBK, PBK, PBK, DBK, IBK, DBK, DBK, /*DBK, IBK,*/ IBK, IBK, /*DBK, DBK, PBK,*/ CIB, PBK, /*DBK, DBK, DBK, DBK, IBK*/ },
		};

		uint8_t ch = *text;

		LineBreakClass cls;

		if (ch == '\n')
			cls = MB;
		else if (ch < 128)
		{
			cls = kASCII_LineBreakTable[ch];
			if (cls > MB and cls != SP) // duh...
				cls = AL;
		}
		else
			cls = UN;

		if (cls == SP)
			cls = WJ;

		LineBreakClass ncls = cls;

		while (++text != end and cls != MB)
		{
			ch = *text;

			LineBreakClass lcls = ncls;

			if (ch == '\n')
			{
				++text;
				break;
			}

			ncls = kASCII_LineBreakTable[ch];

			if (ncls == SP)
				continue;

			BreakAction brk = brkTable[cls][ncls];

			if (brk == DBK or (brk == IBK and lcls == SP))
				break;

			cls = ncls;
		}

		return text;
	}

	size_t m_width;
};

} // namespace cfg
