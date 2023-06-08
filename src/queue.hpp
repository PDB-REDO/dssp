/*-
 * SPDX-License-Identifier: BSD-2-Clause
 * 
 * Copyright (c) 2021 NKI/AVL, Netherlands Cancer Institute
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

#include <condition_variable>
#include <mutex>
#include <queue>

template <typename T, size_t N = 100>
class blocking_queue
{
  public:
	void push(T const &value)
	{
		std::unique_lock<std::mutex> lock(m_guard);

		while (m_queue.size() >= N)
			m_full_signal.wait(lock);

		m_queue.push(value);

		m_empty_signal.notify_one();
	}

	T pop()
	{
		std::unique_lock<std::mutex> lock(m_guard);
		while (m_queue.empty())
			m_empty_signal.wait(lock);

		auto value = m_queue.front();
		m_queue.pop();

		m_full_signal.notify_one();

		return value;
	}

	template<class Rep, class Period>
	std::tuple<bool,T> pop(const std::chrono::duration<Rep, Period>& wait_for)
	{
		std::unique_lock<std::mutex> lock(m_guard);
		while (m_queue.empty())
		{
		    auto now = std::chrono::system_clock::now();
			if (m_empty_signal.wait_until(lock, now + wait_for) == std::cv_status::timeout)
				return { true , T{} };
		}

		auto value = m_queue.front();
		m_queue.pop();

		m_full_signal.notify_one();

		return { false, value };
	}

	bool is_full() const
	{
		std::unique_lock<std::mutex> lock(m_guard);
		return m_queue.size() >= N;
	}

  private:
	std::queue<T> m_queue;
	mutable std::mutex m_guard;
	std::condition_variable m_empty_signal, m_full_signal;
};

template <typename T, size_t N = 10>
class non_blocking_queue
{
  public:
	bool push(T const &value)
	{
		bool result = false;

		std::unique_lock<std::mutex> lock(m_guard);

		if (m_queue.size() < N)
		{
			m_queue.push(value);
			m_empty_signal.notify_one();
			result = true;
		}
			
		return result;
	}

	T pop()
	{
		std::unique_lock<std::mutex> lock(m_guard);
		while (m_queue.empty())
			m_empty_signal.wait(lock);

		auto value = m_queue.front();
		m_queue.pop();

		return value;
	}

  private:
	std::queue<T> m_queue;
	mutable std::mutex m_guard;
	std::condition_variable m_empty_signal;
};