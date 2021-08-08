#ifndef MATH_FUNCS_H
#define MATH_FUNCS_H

class Math {
public:
	Math() {} // useless to instance

    template<class T>
    static T clamp(T v, T minn, T maxn) {
        if (v < minn) return minn;
        if (v > maxn) return maxn;
        return v;
    }

    template <typename T>
    static inline void sort_min_max(T &a, T &b) {
        if (a > b) {
            T temp = a;
            a = b;
            b = temp;
        }
    }
};

#endif
