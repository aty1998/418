/**
 * @file instrument.hpp
 * @brief Implements a key-based instrumentation timer.
 * @author Andrew Yang (atyang)
 */

#ifndef __INSTRUMENT_H__
#define __INSTRUMENT_H__

#include <map>
#include <vector>
#include <string>
#include <chrono>

using namespace std;

typedef chrono::steady_clock::time_point time_point;

/** Instrumentation timer class. */
class Timer
{
public:
    map<string, vector<long long>> times;   /** Time list for each key. */
    long long total;                        /** Total elapsed time. */

    Timer();

    void start(string key);
    void stop(string key);

private:
    map<string, time_point> starts;         /** Start times for each key. */
};

/** Instantiates an instrumentation timer. */
Timer::Timer() : total(0) { }

/** Starts the stopwatch for a specific key. */
void Timer::start(string key)
{
    starts[key] = chrono::steady_clock::now();
}

/** Stops the stopwatch for a specific key and records the result. */
void Timer::stop(string key)
{
    time_point end = chrono::steady_clock::now();
    auto it = starts.find(key);
    if (it != starts.end())
    {
        time_point start = it->second;
        if (times.count(key) == 0)
            times.insert(make_pair(key, vector<long long>()));
        auto diff = chrono::duration_cast<chrono::microseconds>(end - start);
        times[key].push_back(diff.count());
        total += diff.count();
    }
}

#endif /* __INSTRUMENT_H__ */