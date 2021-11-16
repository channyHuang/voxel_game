#pragma once

#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

class Vector2i
{
public:
    union {
        struct {
            int64_t x, y;
        };
        int64_t data[2];
    };

    Vector2i() : x(0), y(0) {}
    Vector2i(const float v) : x(v), y(v) {}
    Vector2i(const float _x, const float _y) : x(_x), y(_y) {}
    ~Vector2i() {}

    const int64_t &operator[](const int axis) const {
        return data[axis];
    }

    void operator=(const Vector2i &v) { x = v.x; y = v.y; }

    Vector2i operator+(const Vector2i &v) {
        return Vector2i (x + v.x, y + v.y);
    }

    Vector2i operator-(const Vector2i &v) {
        return Vector2i (x - v.x, y - v.y);
    }

    Vector2i operator/(const Vector2i &v) {
        if (v.x == 0 || v.y == 0) return *this;
        return Vector2i (x / v.x, y / v.y);
    }

    Vector2i operator*(const Vector2i &v) {
        return Vector2i (x * v.x, y * v.y);
    }

    Vector2i operator+(const int64_t v) {
        return *this + Vector2i(v);
    }

    Vector2i operator-(const int64_t v) {
        return *this - Vector2i(v);
    }

    Vector2i operator/(const int64_t v) {
        return *this / Vector2i(v);
    }

    Vector2i operator*(const int64_t v) {
        return *this * Vector2i(v);
    }

    int64_t volumn() {
        return x * y;
    }

    float len() {
        return std::sqrt(x * x + y * y);
    }

    void normalize() {
        float length = this->len();
        if (length == 0) length = 1.0f;
        x /= length;
        y /= length;
    }

    std::string toString() {
        std::string str = "(" + std::to_string(x) + ", " + std::to_string(y) + ")";
        return str;
    }

    friend Vector2i operator*(const float k, const Vector2i &p)
    {
        return Vector2i(k * p.x, k * p.y);
    }
};
