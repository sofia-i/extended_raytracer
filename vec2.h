//
//  vec2.h
//  raytracer_ext
//
//  Created by Sofia Iannicelli on 3/6/23.
//

#ifndef vec2_h
#define vec2_h

#include <stdio.h>

#include <cmath>
#include <iostream>

template <typename T>
class vec2 {
public:
    T e[2];

    vec2() { }
    vec2(T e0, T e1) {
        e[0] = e0; e[1] = e1;
    }

    // allow access by x, y (acting as coordinates)
    T x() const { return e[0]; }
    T y() const { return e[1]; }

    vec2 operator-() const { return vec3(-e[0], -e[1]); }
    double operator[](int i) const { return e[i]; }
    // double& operator[](int i) { return e[i]; }

    // compound operators
    vec2& operator+=(const vec2<T> &v) {
        e[0] += v.e[0];
        e[1] += v.e[1];
        return *this;
    }

    vec2& operator*=(const vec2<T> &v) {
        e[0] *= v.e[0];
        e[1] *= v.e[1];
        return *this;
    }

    vec2& operator /=(const T val) {
        e[0] *= 1/val;
        e[1] *= 1/val;
        return *this;
    }

    double length() const {
        return sqrt(e[0]*e[0] + e[1]*e[1]);
    }
    
    bool isNormalized() const {
        double length = this->length();
        const auto relative_difference_factor = 0.0001;
        const auto greater_magnitude = std::max(1.0, length);
        
        return std::abs(1.0 - length) < relative_difference_factor * greater_magnitude;
    }
    
    bool checkInBounds(T min, T max) {
        return e[0] >= min && e[0] <= max &&
                e[1] >= min && e[1] <= max;
    }


};

// utility functions
template <typename T>
inline std::ostream& operator<<(std::ostream &out, const vec2<T> &v) {
    return out << v.e[0] << ' ' << v.e[1];
}

template <typename T>
inline vec2<T> operator+(const vec2<T> &u, const vec2<T> &v) {
    return vec2<T>(u.e[0] + v.e[0], u.e[1] + v.e[1]);
}

template <typename T>
inline vec2<T> operator-(const vec2<T> &u, const vec2<T> &v) {
    return vec2<T>(u.e[0] - v.e[0], u.e[1] - v.e[1]);
}

template <typename T>
inline vec2<T> operator*(const vec2<T> &u, const vec2<T> &v) {
    return vec2<T>(u.e[0] * v.e[0], u.e[1] * v.e[1]);
}

template <typename S, typename T>
inline vec2<T> operator*(const S k, const vec2,T> &v) {
    return vec2<T>(k * v.e[0], k * v.e[1]);
}

template <typename S, typename T>
inline vec2<T> operator*(const vec2<T> &v, const S k) {
    return vec2<T>(k * v.e[0], k * v.e[1]);
}

template <typename S, typename T>
inline vec2<T> operator/(const vec2<T> &v, const S k) {
    return vec2<T>((1/k) * v.e[0], (1/k) * v.e[1]);
}

template <typename T>
inline T dot(const vec2<T> &v, const vec2<T> &u) {
    return u.e[0] * v.e[0] + u.e[1] * v.e[1];
}

template <typename T>
inline vec2<T> getUnitVector(const vec2<T> &v) {
    return v / v.length();
}

inline vec2<int> toIntVec3(const vec2<double> &v) {
    return vec2<int>((int)(v.e[0]), (int)(v.e[1]));
}

inline vec2<int> clip(const vec3<int> &v, int low, int high) {
    int x = std::max(v.e[0], low);
    x = std::min(x, high);
    
    int y = std::max(v.e[1], low);
    y = std::min(y, high);
    
    return vec2<int>(x, y);
}

#endif /* vec2_h */
