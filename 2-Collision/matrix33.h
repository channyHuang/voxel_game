#ifndef MATRIX33_H
#define MATRIX33_H

#include "vector3.h"

class Matrix33 {
public:
    union {
        struct {
            Vector3 row[3];
        };
        float data[9];
    };

    inline Matrix33() {}
    inline Matrix33(float f00, float f01, float f02,
                    float f10, float f11, float f12,
                    float f20, float f21, float f22) {
        row[0].x = f00;
        row[0].y = f01;
        row[0].z = f02;
        row[1].x = f10;
        row[1].y = f11;
        row[1].z = f12;
        row[2].x = f20;
        row[2].y = f21;
        row[2].z = f22;
    }
    ~Matrix33() {};

    inline Vector3 operator[] (int r) const {
        return row[r];
    }

    void operator=(const Matrix33 &v) {
        for (int i = 0; i < 9; ++i) {
            data[i] = v.data[i];
        }
    }

    Matrix33 operator+(const Matrix33 &v) {
        Matrix33 res;
        for (int i = 0; i < 9; ++i) {
            data[i] += v.data[i];
        }
        return res;
    }

    inline Vector3 getRow(int r) const {
        return row[r];
    }

    inline Vector3 getCol(int c) const {
        return Vector3(row[0][c], row[1][c], row[2][c]);
    }

    inline Matrix33 Transpose() {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < i; ++j) {
                float tmp = row[i][j];
                row[i][j] = row[j][i];
                row[j][i] = tmp;
            }
        }
    }

    friend Matrix33 operator*(const float k, const Matrix33 &p)
    {
        Matrix33 res;
        for (int i = 0; i < 9; ++i) {
            res.data[i] = p.data[i] * k;
        }
        return res;
    }

    friend Matrix33 operator*(const Matrix33 &p, const Matrix33 &q)
    {
        Matrix33 res;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                res.row[i][j] = p.getRow(i).dot(q.getCol(j));
            }
        }
        return res;
    }
};

#endif // MATRIX33_H
