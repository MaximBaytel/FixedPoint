#ifndef RATIONAL_H
#define RATIONAL_H

#include <iostream>
#include <numeric>
#include <stdexcept>
#include <limits>
#include <string>
#include <sstream>


template <typename T>
class Rational
{
    static_assert(std::is_integral<T>::value);
    static_assert(std::is_unsigned<T>::value);
public:
    explicit Rational(T n = 0, T m = 1, bool isNegative = false): numenator(n), denominator(m), negative(isNegative)
    {
        if (!denominator)
        {
            throw std::invalid_argument("The deniminator can't be zero");
        }

        // if n is zero it isn't a problem, cause we'll get 0/1
        T gcd = std::gcd(numenator, denominator);

        numenator /= gcd;
        denominator /= gcd;

        // if it's zero the sign is irrelevant, so set it to the positive one
        if (!numenator)
        {
            negative = false;
        }

    }

    T GetNumenator() const
    {
        return numenator;
    }

    T GetDenominator() const
    {
        return denominator;
    }

    bool GetNegative() const
    {
        return negative;
    }

    operator std::string() const
    {
        std::stringstream out;

        if (negative)
        {
            out << '-';
        }

        if constexpr (std::is_same<T, std::uint8_t>())
        {
            out << (uint16_t)numenator << "/" << (uint16_t)denominator;
        }
        else
        {
            out << numenator << '/' << denominator;
        }


        return out.str();
    }



    Rational<T> operator - () const
    {
        // nothing is left to do if it's zero
        if (!numenator)
        {
            return *this;
        }

        return Rational(numenator, denominator, !negative);
    }

    Rational<T> GetReciprocal() const
    {
        if (!numenator)
        {
            throw std::invalid_argument("The deniminator can't be zero");
        }

        return Rational(denominator, numenator, negative);
    }

    Rational<T> operator + (const Rational<T>& r) const
    {
        T mult = ControlledMult(denominator, r.GetDenominator());
        T gcd = std::gcd(denominator, r.GetDenominator());
        T lcm = mult/gcd;

        T left  = ControlledMult(numenator, lcm / denominator);
        T right = ControlledMult(r.GetNumenator(), lcm / r.GetDenominator());

        if (negative != r.GetNegative())
        {
            if (left > right)
            {
                return Rational<T>(left - right ,lcm, negative);
            }
            else
            {
                return Rational<T>(right - left ,lcm, r.GetNegative());
            }
        }

        T max = std::max(left, right);
        T res = left + right;

        if (res < max)
        {
            throw std::overflow_error("The type T is overflowed here");
        }

        return Rational<T>(res ,lcm, negative);
    }

    Rational<T> operator - (const Rational<T>& r) const
    {
        return this->operator + (- r);
    }

    Rational<T> operator * (const Rational<T>& r) const
    {
        T gcd_left = std::gcd(numenator, r.GetDenominator());
        T numerator_left = numenator / gcd_left;
        T denominator_right = r.GetDenominator() / gcd_left;

        T gcd_right = std::gcd(r.GetNumenator(), denominator);
        T numerator_right = r.GetNumenator() / gcd_right;
        T denominator_left = denominator / gcd_right;

        return Rational<T>(ControlledMult(numerator_left, numerator_right), ControlledMult(denominator_right, denominator_left), negative != r.GetNegative());

    }


    Rational<T> operator / (const Rational<T>& r) const
    {
        return this->operator * (r.GetReciprocal());
    }

private:

    T ControlledMult(T l, T r) const
    {
        if (!l || !r)
        {
            return 0;
        }

        if (std::numeric_limits<T>::max() / l < r)
        {
            std::cerr << " l = " << l << " r = " << r << std::endl;
            throw std::overflow_error("The type T is overflowed here");
        }

        return l * r;
    }

    T numenator;
    T denominator;
    bool negative;
};

template <typename T>
std::ostream& operator << (std::ostream& out, const Rational<T>& r)
{
    out << std::string(r);
    return out;
}

using tiny_rational = Rational<std::uint8_t>;
using short_rational = Rational<std::uint16_t>;
using int_rational = Rational<std::uint32_t>;
using long_rational = Rational<std::uint64_t>;

#endif // RATIONAL_H
