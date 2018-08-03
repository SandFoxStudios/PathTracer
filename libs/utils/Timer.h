#pragma once
#include <chrono>


class Timer
{
public:
	typedef std::chrono::microseconds us;
	typedef std::chrono::milliseconds ms;
	typedef std::chrono::nanoseconds ns;
public:
	Timer() {}
	~Timer() {}

	inline void start();

	template<typename Duration>
	int64_t elapsedTime();

	std::chrono::steady_clock::time_point begin;
};

template<typename Duration>
inline int64_t Timer::elapsedTime()
{
	return std::chrono::duration_cast<Duration>(std::chrono::steady_clock::now() - begin).count();
}
inline void Timer::start()
{
	begin = std::chrono::steady_clock::now();
}