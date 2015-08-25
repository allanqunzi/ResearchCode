#include <array>

constexpr double PI = 3.14159265358979323846;
constexpr double twoPI = PI * 2.0;

template<typename T, std::size_t D>
using Vector = std::array<T, D>;

template< typename T, std::size_t D >
inline Vector<T, D> operator + ( const Vector<T, D>& l, const Vector<T, D>& r ){
    Vector<T, D> res;

    for (std::size_t i = 0; i < D; ++i)
    {
        res[i] = l[i] + r[i];
    }
    return res;
}

template< typename T, std::size_t D >
inline Vector<T, D> & operator += ( Vector<T, D>& l, const Vector<T, D>& r ){

    for (std::size_t i = 0; i < D; ++i)
    {
        l[i] += r[i];
    }
    return l;
}

template< typename T, std::size_t D >
inline T operator * ( const Vector<T, D>& l, const Vector<T, D>& r ){
    T res = T{};
    for (std::size_t i = 0; i < D; ++i)
    {
        res += l[i] * r[i];

    }
    return res;
}

template< typename T, std::size_t D >
inline Vector<T, D> operator - ( const Vector<T, D>& l, const Vector<T, D>& r ){
    Vector<T, D> res;

    for (std::size_t i = 0; i < D; ++i)
    {
        res[i] = l[i] - r[i];
    }
    return res;
}

template< typename T, std::size_t D >
inline Vector<T, D> & operator -= ( Vector<T, D>& l, const Vector<T, D>& r ){

    for (std::size_t i = 0; i < D; ++i)
    {
        l[i] -= r[i];
    }
    return l;
}