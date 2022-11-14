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

#include <climits>
#include <cstdint>

#if __has_include(<sys/ioctl.h>)
#include <sys/ioctl.h>
#include <fcntl.h>
#include <unistd.h>
#elif defined(_MSC_VER)
#include <windows.h>
#endif

namespace cfg
{

#if defined(_MSC_VER)
/// @brief Get the width in columns of the current terminal
/// @return number of columns of the terminal
inline uint32_t get_terminal_width()
{
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    ::GetConsoleScreenBufferInfo(::GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
    return csbi.srWindow.Right - csbi.srWindow.Left + 1;
}

#elif __has_include(<sys/ioctl.h>)
/// @brief Get the width in columns of the current terminal
/// @return number of columns of the terminal
inline uint32_t get_terminal_width()
{
	uint32_t result = 80;

	if (::isatty(STDOUT_FILENO))
	{
		struct winsize w;
		::ioctl(0, TIOCGWINSZ, &w);
		result = w.ws_col;
	}
	return result;
}
#else
#warning "Could not find the terminal width, falling back to default"
inline uint32_t get_terminal_width()
{
	return 80;
}
#endif

} // namespace cfg