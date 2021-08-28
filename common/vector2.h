#ifndef VECTOR2_H
#define VECTOR2_H

#ifndef real
#define real float
#endif

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>

class Vector2
{
public:
    union {
        struct {
            real x, y;
        };
        real data[2];
    };

    Vector2() : x(0), y(0) {}
    Vector2(const real v) : x(v), y(v) {}
    Vector2(const real _x, const real _y) : x(_x), y(_y) {}
    ~Vector2() {}

    const real &operator[](const int axis) const {
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

    Vector2 operator+(const real v) {
        return *this + Vector2(v);
    }

    Vector2 operator-(const real v) {
        return *this - Vector2(v);
    }

    Vector2 operator/(const real v) {
        return *this / Vector2(v);
    }

    Vector2 operator*(const real v) {
        return *this * Vector2(v);
    }

    real volumn() {
        return x * y;
    }

    real len() {
        return std::sqrt(x * x + y * y);
    }

    void normalize() {
        real length = this->len();
        if (length == 0) length = 1.0f;
        x /= length;
        y /= length;
    }

    std::string toString() {
        std::string str = "(" + std::to_string(x) + ", " + std::to_string(y) + ")";
        return str;
    }

    friend Vector2 operator*(const real k, const Vector2 &p)
    {
        return Vector2(k * p.x, k * p.y);
    }
};

#endif // VECTOR2_H
