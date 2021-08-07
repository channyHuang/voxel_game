#ifndef Vector3I_H
#define Vector3I_H

#ifndef real
#define real float
#endif

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>

class Vector3i
{
public:
    union {
        struct {
            int64_t x, y, z;
        };
        int64_t data[3];
    };

    Vector3i() : x(0), y(0), z(0) {}
    Vector3i(const int64_t v) : x(v), y(v), z(v) {}
    Vector3i(const int64_t _x, const int64_t _y, const int64_t _z) : x(_x), y(_y), z(_z) {}
    ~Vector3i() {}

    void operator=(const Vector3i &v) { x = v.x; y = v.y; z = v.z; }

    Vector3i operator+(const Vector3i &v) {
        return Vector3i (x + v.x, y + v.y, z + v.z);
    }

    Vector3i operator-(const Vector3i &v) {
        return Vector3i (x - v.x, y - v.y, z - v.z);
    }

    Vector3i operator/(const Vector3i &v) {
        if (v.x == 0 || v.y == 0 || v.z == 0) return *this;
        return Vector3i (x / v.x, y / v.y, z / v.z);
    }

    Vector3i operator*(const Vector3i &v) {
        return Vector3i (x * v.x, y * v.y, z * v.z);
    }

    Vector3i operator+(const int64_t v) {
        return *this + Vector3i(v);
    }

    Vector3i operator-(const int64_t v) {
        return *this - Vector3i(v);
    }

    Vector3i operator/(const int64_t v) {
        return *this / Vector3i(v);
    }

    Vector3i operator*(const int64_t v) {
        return *this * Vector3i(v);
    }

    Vector3i dot(const Vector3i &v) {
        return (x * v.x + y * v.y + z * v.z);
    }

    Vector3i cross(const Vector3i &v) {
        return Vector3i(y * v.z - z * v.y, x * v.z - z * v.x, x * v.y - y * v.x);
    }

    int64_t volumn() {
        return x * y * z;
    }

    real len() {
        return std::sqrt(x * x + y * y + z * z);
    }

    void normalize() {
        real length = this->len();
        if (length == 0) length = 1.0f;
        x /= length;
        y /= length;
        z /= length;
    }

    std::string toString() {
        std::string str = "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
        return str;
    }
};


#endif // Vector3iI_H
