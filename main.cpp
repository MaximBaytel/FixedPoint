#include <stdint.h>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <climits>
#include <iostream>
#include <sstream>

// FractionLength - how many bits are after a period
template <uint8_t FractionLength, typename Basetype = int32_t, typename HelperType = uint64_t>
class FixedPoint
{
public:
    static_assert(sizeof(HelperType) == 2 * sizeof(Basetype));
    static_assert(FractionLength > 0);
    static_assert(FractionLength < sizeof(Basetype) * CHAR_BIT);
    static_assert(std::is_integral<Basetype>::value);
    static_assert(std::is_integral<HelperType>::value);
    static_assert(std::is_signed<Basetype>::value);
    static_assert(std::is_unsigned<HelperType>::value);

    explicit FixedPoint(int decimal, unsigned int fraction = 0): val(0)
    {
        unsigned int decimal_abs = std::abs(decimal);
        if (decimal_abs > max_dec || decimal < min_dec || (decimal == min_dec && fraction) || (fraction >= full_fraction))
        {
            throw std::invalid_argument("It won't fit our so few bits");
        }

        val |= fraction;
        val |= (decimal_abs << FractionLength);
        val *= isign(decimal);
    }

    explicit FixedPoint(double v): val(0)
    {
        double decimal_double = 0.0;
        double fraction = std::modf(v, &decimal_double);

        if (decimal_double > max_dec || decimal_double < min_dec)
        {
            throw std::invalid_argument("It won't fit our so few bits");
        }

        Basetype decimal_part = static_cast<Basetype>(decimal_double);
        int8_t sign = isign(v);
        UnsignedBasetype decimal_part_abs = std::abs(decimal_part);
        double fraction_abs = std::abs(fraction);

        UnsignedBasetype count = static_cast<UnsignedBasetype>(std::scalbln(fraction_abs, FractionLength + 1));

        count += count & least_bit_mask;
        count >>= 1;

        if (count && decimal_part == min_dec_abs)
        {
            throw std::invalid_argument("It won't fit our so few bits");
        }

        if (count == full_fraction)
        {
            decimal_part_abs += 1;

            auto decimal_part_signed = sign * decimal_part;

            if (decimal_part_abs > max_dec || decimal_part_signed < min_dec)
            {
                throw std::invalid_argument("It won't fit our so few bits");
            }
        }
        else
        {
            val |= count;
        }

        val |= (decimal_part_abs << FractionLength);
        val *= sign;
    }


    FixedPoint operator + (const FixedPoint& r) const
    {
        return makeFx(val + r.val);
    }

    FixedPoint operator - (const FixedPoint& r) const
    {
        return makeFx(val - r.val);
    }

    FixedPoint operator - () const
    {
        return makeFx(-val);
    }

    FixedPoint operator * (const FixedPoint& r) const
    {
        auto total_sign = isign(val) * isign(r.val);
        HelperType temp = HelperType(std::abs(val)) * HelperType(std::abs(r.val));
        temp >>= (FractionLength - 1);
        uint8_t unit = temp & least_bit_mask;

        return makeFx(total_sign * ((temp >> 1) + unit));
    }

    FixedPoint operator / (const FixedPoint& r) const
    {
        auto total_sign = isign(val) * isign(r.val);

        HelperType left = std::abs(val);
        left <<= (FractionLength + 1);
        UnsignedBasetype temp = left / std::abs(r.val);
        uint8_t unit = temp & least_bit_mask;

        return makeFx(total_sign * ((temp >> 1) + unit ));
    }

    operator std::string() const
    {
        std::stringstream res;
        UnsignedBasetype temp = val;
        auto fraction = temp & fraction_mask;

        if (sign_mask & temp)
        {
            res << '-';

            temp = ~temp + 1;
            fraction = temp & fraction_mask;
        }

        res << (temp >> FractionLength);


        if (fraction)
        {
            uint64_t current_fraction = fraction;
            uint64_t current_unit = fraction_unit;

            while ((std::numeric_limits<uint64_t>::max() / current_fraction) < current_unit)
            {
                uint8_t least_bit = current_fraction & least_bit_mask;
                current_unit /= 5;
                current_fraction = (current_fraction + least_bit)/ 2;
            }


            res << '.' << current_fraction * current_unit;
        }

        return res.str();
    }

private:
    Basetype val;

    using UnsignedBasetype = std::make_unsigned_t<Basetype>;

    // fills and return given number of ones in the least significant bits of a byte
    static constexpr UnsignedBasetype mask(uint8_t num)
    {
        return num == 0 ? 0 : (mask(num - 1) << 1) | 1;
    }

    static constexpr uint64_t ipow(uint8_t num, unsigned int pow)
    {
        return (pow >= sizeof(unsigned int)*8) ? 0 :
                   pow == 0 ? 1 : num * ipow(num, pow-1);
    }

    template<typename T>
    static int8_t isign(T v)
    {
        return v >= 0? +1 : -1;
    }

    static constexpr uint8_t full_length = sizeof(Basetype) * CHAR_BIT;
    static constexpr UnsignedBasetype decimal_lenght = full_length - FractionLength;
    static constexpr UnsignedBasetype max_dec = (1 << (decimal_lenght - 1)) - 1;
    static constexpr Basetype  min_dec_abs = (1 << (decimal_lenght - 1));
    static constexpr Basetype  min_dec = -min_dec_abs;
    static constexpr HelperType five = 5;
    static constexpr uint64_t fraction_unit = ipow(five, FractionLength);

    static constexpr UnsignedBasetype full_fraction = ipow(2, FractionLength);
    static constexpr UnsignedBasetype sign_mask = 1 << (full_length - 1);
    static constexpr UnsignedBasetype fraction_mask = mask(FractionLength);
    static constexpr UnsignedBasetype decimal_mask = mask(full_length) & ~FractionLength;
    static constexpr uint8_t least_bit_mask = 0x01;


    explicit FixedPoint(): val(0)
    {
    }

    static FixedPoint makeFx(Basetype v)
    {
        FixedPoint fp;
        fp.val = v;
        return fp;
    }
};

template <uint8_t FractionLength>
std::ostream& operator << (std::ostream& out, const FixedPoint<FractionLength>& number)
{
    out << std::string(number);

    return out;
}

template <uint8_t FractionLength>
const FixedPoint<FractionLength> operator + (const FixedPoint<FractionLength>& l, int r)
{
    return l + FixedPoint<FractionLength>(r);
}

template <uint8_t FractionLength>
const FixedPoint<FractionLength> operator + (int l, const FixedPoint<FractionLength>& r)
{
    return r + FixedPoint<FractionLength>(l);
}

template <uint8_t FractionLength>
const FixedPoint<FractionLength> operator - (const FixedPoint<FractionLength>& l, int r)
{
    return l - FixedPoint<FractionLength>(r);
}

template <uint8_t FractionLength>
const FixedPoint<FractionLength> operator - (int l, const FixedPoint<FractionLength>& r)
{
    return FixedPoint<FractionLength>(l) - r;
}

template <uint8_t FractionLength>
const FixedPoint<FractionLength> operator * (const FixedPoint<FractionLength>& l, int r)
{
    return l * FixedPoint<FractionLength>(r);
}

template <uint8_t FractionLength>
const FixedPoint<FractionLength> operator * (int l, const FixedPoint<FractionLength>& r)
{
    return r * FixedPoint<FractionLength>(l);
}

template <uint8_t FractionLength>
const FixedPoint<FractionLength> operator / (const FixedPoint<FractionLength>& l, int r)
{
    return l / FixedPoint<FractionLength>(r);
}

template <uint8_t FractionLength>
const FixedPoint<FractionLength> operator / (int l, const FixedPoint<FractionLength>& r)
{
    return FixedPoint<FractionLength>(l) / FixedPoint<FractionLength>(r);
}

template <uint8_t FractionLength>
const FixedPoint<FractionLength> operator *= (FixedPoint<FractionLength>& l, const FixedPoint<FractionLength>& r)
{
    l = l * r;
    return l;
}

template <uint8_t FractionLength>
const FixedPoint<FractionLength> operator += (FixedPoint<FractionLength>& l, const FixedPoint<FractionLength>& r)
{
    l = l + r;
    return l;
}

// let's introduce aliases for convenience
using FixedPoint_2 = FixedPoint<2>;
using FixedPoint_4 = FixedPoint<4>;
using FixedPoint_10 = FixedPoint<10>;
using FixedPoint_12 = FixedPoint<12>;
using FixedPoint_14 = FixedPoint<14>;
using FixedPoint_16 = FixedPoint<16>;
using FixedPoint_18 = FixedPoint<18>;
using FixedPoint_19 = FixedPoint<18>;
using FixedPoint_20 = FixedPoint<20>;
using FixedPoint_22 = FixedPoint<20>;
using FixedPoint_24 = FixedPoint<24>;
using FixedPoint_28 = FixedPoint<28>;
using FixedPoint_30 = FixedPoint<30>;

using FixedPoint_not_compile = FixedPoint<4, int, unsigned long long>;

//using FixedPoint_15 = FixedPoint<9>;

template <typename T>
T sinus(T x, uint8_t n = 3)
{
    T t = x;
    T res = t;
    for ( int i=1; i<n; ++i)
    {
        T mult = -x*x/((2*i+1)*(2*i));
        t *= mult;
        res += t;
    }
    return res;
}

template <typename T>
T exponenta(T x, uint8_t n = 5)
{
    T res{1.0};
    T num {1.0};
    uint64_t den {1};

    for ( int i=1; i<n; ++i)
    {
        num *= x;
        den *= i;

        res += num/den;
    }

    return res;
}

template <typename T>
T poly(T x)
{
    return (2*x + 1)*x - 3;
}

template <typename T>
T poly1(T x)
{
    return 2*x;
}

template <typename T>
T tricky_poly(T y, T z)
{
    return 108 - (815 - 1500/z)/y;
}

int main()
{
    //FixedPoint_not_compile a{1.25};
    //FixedPoint_2 p1{3.25}; //it must be 1.5
    // FixedPoint_2 p2{1.25};  //it must be 4.25
    // //FixedPoint p6{0.25};  //it must be 4.25


     //std::cout << p1 << ' ' /*<< p2*/ << std::endl;
    // std::cout << "sum = " << p1 + p2 << std::endl;
    // std::cout << "p1 - p2 = " << p1 - p2 << std::endl;
    // std::cout << "p2 - p1 = " << p2 - p1 << std::endl;
    // std::cout << "p2 * p1 = " << p2 * p1 << std::endl;
    // std::cout << "p2 / p1 = " << p2 / p1 << std::endl;
    // std::cout << "p1 / p2 = " << p1 / p2 << std::endl;

    // std::cout << "-p1 = " << -p1 << std::endl;
    // std::cout << "-p2 = " << -p2 << std::endl;
    // std::cout << "2*p1 = " << 2*p1 << std::endl;


    // FixedPoint_16 p3{0.75}; //it must be 1.5

    // std::cout << p3 << ' ' << p3*p3 << ' ' << p3 * p3 * p3 << ' ' <<  p3 * p3 * p3 * p3 << std::endl;
    /*
    FixedPoint p3{0b00010000}; //it must be 1.5
    FixedPoint p4{0b00001000};  //it must be 4.25

    std::cout << "p3 / p4 = " << p3 / p4 << std::endl;
    std::cout << "p4 / p3 = " << p4 / p3 << std::endl;
    */

    //double x = 0.261799;

    //std::cout << "Original number: " << x << std::endl << FixedPoint_20{x} << std::endl;
    //std::cout << "Its fixed point versions: " << FixedPoint_2{x} << ' ' << FixedPoint_4{x} << ' ' << FixedPoint_16{x} << ' ' << FixedPoint_20{x} << ' ' << FixedPoint_24{x}  << ' ' << FixedPoint_30{x} << std::endl;

    // double x = 0.261799;
    // FixedPoint_16 fx{x};

    // std::cout << std::setprecision(20) << exponenta(x) << ' ' << exponenta(fx) << ' ' << std::exp(x) << std::endl;
    // std::cout << std::setprecision(20) << sinus(x) << ' ' << sinus(fx) << ' ' << std::sin(x) << std::endl;

    FixedPoint_20 x_prev{4};
    FixedPoint_20 x_curr{4.25};

    for (int i=2; i <= 100; ++i)
    {
        auto temp = tricky_poly(x_curr, x_prev);
        x_prev = x_curr;
        x_curr = temp;

        std::cout << "tricky poly res = " << x_curr << std::endl;

    }

    std::cout << "tricky poly res = " << x_curr << std::endl;

    //std::cout << std::setprecision(20) << poly(x) << ' ' << poly(fx) << std::endl;
     //std::cout << poly1(x) << ' ' << poly1(fx) << std::endl;


    return 0;
}
