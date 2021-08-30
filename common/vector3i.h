#pragma once

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>

#include "math_funcs.h"

class Vector3i
{
public:
    union {
        struct {
            int64_t x, y, z;
        };
        int64_t data[3];
    };

    const int64_t &operator[](const int axis) const {
        return data[axis];
    }

    Vector3i() : x(0), y(0), z(0) {}
    Vector3i(const int64_t v) : x(v), y(v), z(v) {}
    Vector3i(const int64_t _x, const int64_t _y, const int64_t _z) : x(_x), y(_y), z(_z) {}
    Vector3i(const Vector3 &f) {
        x = std::floor(f.x);
        y = std::floor(f.y);
        z = std::floor(f.z);
    }

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

    Vector3i operator*(const Vector3i &v) const {
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

    float len() {
        return std::sqrt(x * x + y * y + z * z);
    }

    void normalize() {
        float length = this->len();
        if (length == 0) length = 1.0f;
        x /= length;
        y /= length;
        z /= length;
    }

    std::string toString() {
        std::string str = "(" + std::to_string(x) + ", " + std::to_string(y) + ", " + std::to_string(z) + ")";
        return str;
    }

    inline unsigned int get_zxy_index(const Vector3i area_size) const {
        return y + area_size.y * (x + area_size.x * z); // ZXY
    }

    Vector3 to_vec3() const {
        return Vector3(x, y, z);
    }

    static void sort_min_max(Vector3i &a, Vector3i &b) {
        Math::sort_min_max(a.x, b.x);
        Math::sort_min_max(a.y, b.y);
        Math::sort_min_max(a.z, b.z);
    }

    // Clamps between min and max, where max is excluded
    void clamp_to(const Vector3i min, const Vector3i max) {
        if (x < min.x) {
            x = min.x;
        }
        if (y < min.y) {
            y = min.y;
        }
        if (z < min.z) {
            z = min.z;
        }

        if (x >= max.x) {
            x = max.x - 1;
        }
        if (y >= max.y) {
            y = max.y - 1;
        }
        if (z >= max.z) {
            z = max.z - 1;
        }
    }
};

inline bool operator==(const Vector3i &a, const Vector3i &b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

inline bool operator!=(const Vector3i &a, const Vector3i &b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}

inline Vector3i operator+(const Vector3i &a, const Vector3i &b) {
    return Vector3i(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vector3i operator-(const Vector3i &a, const Vector3i &b) {
    return Vector3i(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vector3i operator>>(const Vector3i &a, int b) {
    return Vector3i(a.x >> b, a.y >> b, a.z >> b);
}

inline Vector3i operator<<(const Vector3i &a, int b) {
    return Vector3i(a.x << b, a.y << b, a.z << b);
}

inline Vector3i operator*(const int64_t p_scalar, const Vector3i &p_vec) {
    return p_vec * p_scalar;
}

