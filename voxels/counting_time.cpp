#include "counting_time.h"

CountingTime* CountingTime::instanceCountingTime = nullptr;
std::chrono::steady_clock::time_point CountingTime::m_startTimePointer = std::chrono::steady_clock::now();
