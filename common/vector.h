#ifndef VECTOR_H
#define VECTOR_H

#ifndef real
#define real float
#endif

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

template <typename T, int N>
class Vector
{
public:
    T data[N];

    Vector() {}
    Vector(const T &v) {
        for (unsigned int i = 0; i < N; i++) {
            data[i] = v;
        }
    }
    Vector(const std::vector<T> &v) {
        for (unsigned int i = 0; i < N; i++) {
            data[i] = v[i];
        }
    }
    Vector(const T &a, const T& b) {
        if (N >= 1) {
            data[0] = a;
        }
        if (N >= 2) {
            data[1] = b;
        }
    }
    ~Vector() {}

    void operator=(const Vector &v) {
        for (unsigned int i = 0; i < N; i++) {
            data[i] = v[i];
        }
    }

    Vector operator+(const Vector &v) {
        Vector<T, N> res;
        for (unsigned int i = 0; i < N; i++) {
            res[i] = data[i] + v[i];
        }
        return res;
    }

    Vector operator-(const Vector &v) {
        Vector<T, N> res;
        for (unsigned int i = 0; i < N; i++) {
            res[i] = data[i] - v[i];
        }
        return res;
    }

    Vector operator/(const Vector &v) {
        Vector<T, N> res;
        for (unsigned int i = 0; i < N; i++) {
            if (v[i] == 0) return *this;
            res[i] = data[i] / v[i];
        }
        return res;
    }

    Vector operator*(const Vector &v) {
        Vector<T, N> res;
        for (unsigned int i = 0; i < N; i++) {
            res[i] = data[i] * v[i];
        }
        return res;
    }

    Vector operator+(const T &v) {
        return *this + Vector(v);
    }

    Vector operator-(const T &v) {
        return *this - Vector(v);
    }

    Vector operator/(const T &v) {
        return *this / Vector(v);
    }

    Vector operator*(const T &v) {
        return *this * Vector(v);
    }

    Vector dot(const Vector &v) {
        Vector<T, N> res;
        for (unsigned int i = 0; i < N; i++) {
            res[i] = data[i] * v[i];
        }
        return res;
    }

    T volumn() {
        T res;
        for (unsigned int i = 0; i < N; i++) {
            res *= data[i];
        }
        return res;
    }

    real len() {
        real res;
        for (unsigned int i = 0; i < N; i++) {
            res += data[i] * data[i];
        }
        return std::sqrt(res);
    }

    void normalize() {
        real length = this->len();
        if (length == 0) length = 1.0f;
        for (unsigned int i = 0; i < N; i++) {
            data[i] /= length;
        }
    }

    std::string toString() {
        std::string str = "(";
        for (unsigned int i = 0; i < N; i++) {
            str += std::to_string(data[i]);
            if (i != N - 1) {
                str += ", ";
            }
        }
        return str;
    }
};

#endif // VECTOR_H
