#pragma once

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>

class Vector2
{
public:
    union {
        struct {
            float x, y;
        };
        float data[2];
    };

    Vector2() : x(0), y(0) {}
    Vector2(const float v) : x(v), y(v) {}
    Vector2(const float _x, const float _y) : x(_x), y(_y) {}
    ~Vector2() {}

    const float &operator[](const int axis) const {
        return data[axis];
    }

    void operator=(const Vector2 &v) { x = v.x; y = v.y; }

    Vector2 operator+(const Vector2 &v) {
        return Vector2 (x + v.x, y + v.y);
    }

    Vector2 operator-(const Vector2 &v) {
        return Vector2 (x - v.x, y - v.y);
    }

    Vector2 operator/(const Vector2 &v) {
        if (v.x == 0 || v.y == 0) return *this;
        return Vector2 (x / v.x, y / v.y);
    }

    Vector2 operator*(const Vector2 &v) {
        return Vector2 (x * v.x, y * v.y);
    }

    Vector2 operator+(const float v) {
        return *this + Vector2(v);
    }

    Vector2 operator-(const float v) {
        return *this - Vector2(v);
    }

    Vector2 operator/(const float v) {
        return *this / Vector2(v);
    }

    Vector2 operator*(const float v) {
        return *this * Vector2(v);
    }

    float volumn() {
        return x * y;
    }

    float len() {
        return std::sqrt(x * x + y * y);
    }

    float dot(const Vector2 &v) {
        return (x * v.y + y * v.x);
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

    friend Vector2 operator*(const float k, const Vector2 &p)
    {
        return Vector2(k * p.x, k * p.y);
    }
};
