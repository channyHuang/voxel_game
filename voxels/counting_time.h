#pragma once

#include <iostream>
#include <ctime>
#include <chrono>

class CountingTime {
public:
    static CountingTime *getInstance() {
        if (instanceCountingTime == nullptr) {
            instanceCountingTime = new CountingTime;
        }
        return instanceCountingTime;
    }

    CountingTime() {
        m_startTimePointer = std::chrono::steady_clock::now();
    }

    static uint64_t getSpendTime() {
        std::chrono::steady_clock::time_point nowTimePoint = std::chrono::steady_clock::now();
        std::chrono::milliseconds timeSpanMilli = std::chrono::duration_cast<std::chrono::milliseconds>(nowTimePoint - m_startTimePointer);
        return (uint64_t)timeSpanMilli.count();

        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        std::string sTime = std::string(std::ctime(&end_time));
        //outputVoronoiDiagram(nFaceSize, nVerticesSize, faceVertices, facecolors, "rough_" + sTime.substr(11, 2) + sTime.substr(14, 2) + sTime.substr(17, 2) + ".ply");
    }

    static uint64_t get_ticks_msec() {
        return getSpendTime();
    }

    static uint64_t get_ticks_usec() {
        return getSpendTime();
    }

private:
    static std::chrono::steady_clock::time_point m_startTimePointer;
    static CountingTime* instanceCountingTime;
};

