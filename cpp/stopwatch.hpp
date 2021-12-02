#pragma once
#include <chrono>

template<typename Clock = std::chrono::high_resolution_clock>
class Stopwatch {
    // Stopwatch for program timing. Supports resuming and averaging.
  public:
    using TimePoint = typename Clock::time_point;
    using Duration  = typename Clock::duration;
    using DefaultOutputRep  = double;
    using DefaultOutputUnit = std::ratio<1,1>; // seconds, no scaling

    Stopwatch() : running{true}, start_time{Clock::now()}, elapsed{0}, sum_elapsed{0}, num_elapsed{0} {}

    // Stopwatch operation
    void start() {
        start_time = Clock::now();
        running = true;
    }
    void stop() {
        if (!running) return;
        elapsed += Clock::now() - start_time;
        running = false;
    }
    template<typename Rep = DefaultOutputRep, typename Unit = DefaultOutputUnit>
    Rep save() {
        assert(!running);
        auto result = as<Rep, Unit>(elapsed);
        sum_elapsed += elapsed;
        num_elapsed += 1;
        elapsed = Duration{0};
        return result;
    }
    template<typename Rep = DefaultOutputRep, typename Unit = DefaultOutputUnit>
    Rep mark() {
        // Like save, but for a running stopwatch
        assert(running);
        auto now = Clock::now();
        elapsed += now - start_time;
        start_time = now;
        running = false;
        auto result = save<Rep, Unit>();
        running = true;
        return result;
    }

    // Read the stopwatch
    template<typename Rep = DefaultOutputRep, typename Unit = DefaultOutputUnit>
    Rep get_elapsed() const {
        auto result = elapsed;
        if (running) result += Clock::now() - start_time;
        return as<Rep, Unit>(result);
    }
    template<typename Rep = DefaultOutputRep, typename Unit = DefaultOutputUnit>
    Rep get_avg_elapsed() const {
        return as<Rep, Unit>(sum_elapsed)/static_cast<Rep>(num_elapsed);
    }

  private:
    bool         running;
    TimePoint    start_time;
    Duration     elapsed;
    Duration     sum_elapsed;
    std::size_t  num_elapsed;

    // Shorthand: convert Duration into a numeric value with desired units
    template<typename Rep, typename Unit>
    static Rep as(const Duration& d) {
        return std::chrono::duration_cast<std::chrono::duration<Rep, Unit>>(d).count();
    }
};
