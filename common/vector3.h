#ifndef VECTOR3_H
#define VECTOR3_H

#ifndef real
#define real float
#endif

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>

class Vector3
{
public:
    union {
        struct {
            real x, y, z;
        };
        real data[3];
    };

    Vector3() : x(0), y(0), z(0) {}
    Vector3(const real v) : x(v), y(v), z(v) {}
    Vector3(const real _x, const real _y, const real _z) : x(_x), y(_y), z(_z) {}
    ~Vector3() {}

    const real &operator[](const int axis) const {
        return data[axis];
    }

    void operator=(const Vector3 &v) { x = v.x; y = v.y; z = v.z; }

    Vector3 operator+(const Vector3 &v) {
        return Vector3 (x + v.x, y + v.y, z + v.z);
    }

    Vector3 &operator+=(const Vector3 &v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    Vector3 operator-(const Vector3 &v) {
        return Vector3 (x - v.x, y - v.y, z - v.z);
    }

    Vector3 &operator -= (const Vector3 &v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    Vector3 operator/(const Vector3 &v) {
        if (v.x == 0 || v.y == 0 || v.z == 0) return *this;
        return Vector3 (x / v.x, y / v.y, z / v.z);
    }

    Vector3 operator*(const Vector3 &v) {
        return Vector3 (x * v.x, y * v.y, z * v.z);
    }

    Vector3 operator+(const real v) {
        return *this + Vector3(v);
    }

    Vector3 operator-(const real v) {
        return *this - Vector3(v);
    }

    Vector3 operator/(const real v) {
        return *this / Vector3(v);
    }

    Vector3 operator*(const real v) {
        return *this * Vector3(v);
    }

    real dot(const Vector3 &v) const {
        return (x * v.x + y * v.y + z * v.z);
    }

    Vector3 cross(const Vector3 &v) {
        return Vector3(y * v.z - z * v.y, x * v.z - z * v.x, x * v.y - y * v.x);
    }

    real volumn() {
        return x * y * z;
    }

    real len() {
        return std::sqrt(x * x + y * y + z * z);
    }

    real lenSqr() const {
        return (x * x + y * y + z * z);
    }

    real distanceTo(const Vector3 &v) {
        return (*this - v).len();
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

    friend Vector3 operator*(const real k, const Vector3 &p)
    {
        return Vector3(k * p.x, k * p.y, k * p.z);
    }
};

#endif // VECTOR3_H
