#include <iostream>
#include <memory>
#include <cstdint>
#include <iomanip>
#include "array_types.hpp"

#include <functional>
#include <chrono>
#include <cmath>

using intptr_t = std::intptr_t;

class dfloat {
private:
    float a;
    float c;
public:
    dfloat(double x) {
        a = static_cast<float>(x);
        c = static_cast<float>(x - static_cast<double>(a));
    }
    dfloat(float x) {
        a = x;
        c = 0.0f;
    }
    dfloat(float new_a, float new_c) {
        a = new_a;
        c = new_c;
    }
    dfloat() = default;

    // реализация функций из научной работы для перегрузки операторов
    float RN(double t) const {
        return static_cast<float>(t);
    }
    dfloat Fast2Sum(float a, float b) const {
        float s = RN(a + b);
        float z = RN(s - a);
        float t = RN(b - z);
        dfloat ans(s, t);
        return ans;
    }
    dfloat toSum(float a, float b) const {
        float s = RN(a + b);
        float a_p = RN(s - b);
        float b_p = RN(s - a_p);
        float del_a = RN(a - a_p);
        float del_b = RN(b - b_p);
        float t = RN(del_a - del_b);
        dfloat ans(s, t);
        return ans;
    }
    dfloat toProd(float a, float b) const {
        float pi = RN(a * b);
        float rho = RN(a * b - pi);
        dfloat ans(pi, rho);
        return ans;
    }

    // перегрузка операторов
    dfloat& operator=(const dfloat& other) {
        if (this != &other) {
            a = other.a;
            c = other.c;
        }
        return *this;
    }

    dfloat& operator=(const double& x) {
        a = static_cast<float>(x);
        c = static_cast<float>(x - static_cast<double>(a));
        return *this;
    }

    dfloat operator+(const dfloat& other) const {
        dfloat s = toSum(a, other.a);
        float v = RN(c + other.c);
        float w = RN(s.c + v);
        dfloat z = Fast2Sum(s.a, w);
        return z;
    }

    dfloat operator-(const dfloat& other) const {
        dfloat x(-other.a, -other.c);
        return *this + x;
    }

    dfloat operator*(const dfloat& other) const {
        dfloat c_f = toProd(a, other.a);
        float tl1 = RN(a * other.c);
        float tl2 = RN(c * other.a);
        float cl2 = RN(tl1 + tl2);
        float cl3 = RN(c_f.c + cl2);
        dfloat z = Fast2Sum(c_f.a, cl3);
        return z;
    }

    dfloat operator/(const dfloat& other) const {
        float th = RN(a / other.a);

        // DWTimesFP1(yh, y, th);
        dfloat c_p = toProd(other.a, th);
        float cl2 = RN(other.c * th);
        dfloat t = Fast2Sum(c_p.a, cl2);
        float tl2 = RN(t.c + c_p.c);
        dfloat r = Fast2Sum(t.a, tl2);

        dfloat pi = toSum(a, -r.a);
        float del_h = RN(pi.c - r.c);
        float del_l = RN(del_h + c);
        float del = RN(pi.a + del_l);
        float tl = RN(del / other.a);
        dfloat z = Fast2Sum(th, tl);
        return z;
    }

    bool operator<(const dfloat& other) const {
        if (a < other.a) return true;
        else if (a == other.a) {
            return c < other.c;
        }
        return false;
    }
    bool operator>(const dfloat& other) const {
        if (a > other.a) return true;
        else if (a == other.a) {
            return c > other.c;
        }
        return false;
    }
    bool operator==(const dfloat& other) const {
        if (a == other.a && c == other.c) return true;
        return false;
    }

    double to_double() const {
        return static_cast<double>(a) + static_cast<double>(c);
    }
    void print() const {
        std::cout << "a = " << a << ", c = " << c << std::endl;
    }

    friend std::ostream& operator<<(std::ostream& os, const dfloat& df) {
        os << "(" << df.a << ", " << df.c << ")";
        return os;
    }
};

template <class T>
struct benchresult {
    T result;
    double btime;
};

template <class T, class input_type>
auto benchmark(std::function<T(input_type)> fn, input_type input, intptr_t nrepeat) {
    T result;
    auto start = std::chrono::steady_clock::now();
    for (intptr_t i = 0; i < nrepeat; i++) {
        result = fn(input);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration_s = end - start;
    double ms_per_run = duration_s.count() * 1000 / nrepeat;
    return benchresult<T> {result, ms_per_run};
}

template <class T>
T quadratic_form(matrix<T> A, vec<T> x){
    T ans = 0., s;
    intptr_t n = x.length();
    for (intptr_t i = 0; i < n; i++) {
        s = 0.;
        for (intptr_t j = 0; j < n; j++) {
            s = s + x(j) * A(i, j);
        }
        ans = ans + x(i) * s;
    }
    return ans;
}

int main(int argc, char* argv[])
{
    intptr_t n;

    std::cin >> n;
    matrix<double> a(n, n);
    vec<double> x(n);

    matrix<dfloat> a_df(n, n);
    vec<dfloat> x_df(n);

    for (intptr_t k = 0; k < n * n; k++) {
        std::cin >> a(k);
    }

    for (intptr_t k = 0; k < n; k++) {
        std::cin >> x(k);
    }

    for (intptr_t k = 0; k < n * n; k++) {
        a_df(k) = a(k);
    }

    for (intptr_t k = 0; k < n; k++) {
        x_df(k) = x(k);
    }


    std::function<double(int)> test_quadratic = [=](int idx) {return quadratic_form(a, x); };
    std::function<dfloat(int)> test_quadratic_df = [=](int idx) {return quadratic_form(a_df, x_df); };

    auto benchresult = benchmark(test_quadratic, 0, 1000);
    auto benchresult_df = benchmark(test_quadratic_df, 0, 1000);

    double err = std::abs(benchresult.result - benchresult_df.result.to_double());
    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 10)
        << "Timing: " << benchresult.btime << " ms\n"
        << "Answer = " << benchresult.result
        << std::endl
        << "Timing (dfloat): " << benchresult_df.btime << " ms\n"
        << "Answer (dfloat) = " << benchresult_df.result.to_double() << " " << benchresult_df.result
        << std::endl
        << "Error = " << err
        << std::endl << std::endl;
    return 0;
}