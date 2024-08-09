#include <stdint.h>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <climits>
#include <iostream>
#include <sstream>

#include "decimal.h"

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
    //static_assert(std::is_signed<Basetype>::value);
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

    FixedPoint& operator << (uint8_t count)
    {
        val <<= count;
        return *this;
    }

    FixedPoint& operator >> (uint8_t count)
    {
        val >>= count;
        return *this;
    }


    FixedPoint operator * (const FixedPoint& r) const
    {
        // auto total_sign = isign(val) * isign(r.val);
        // HelperType temp = HelperType(std::abs(val)) * HelperType(std::abs(r.val));
        // temp >>= (FractionLength - 1);
        // uint8_t unit = temp & least_bit_mask;

        // return makeFx(total_sign * ((temp >> 1) + unit));

        //auto total_sign = isign(val) * isign(r.val);
        HelperType temp = HelperType(HelperType(val) * HelperType(r.val));

        return makeFx(temp >> FractionLength);
        // temp >>= (FractionLength - 1);
        // uint8_t unit = temp & least_bit_mask;

        // return makeFx(((temp >> 1) + unit));

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

        // if (sign_mask & temp)
        // {
        //     res << '-';

        //     temp = ~temp + 1;
        //     fraction = temp & fraction_mask;
        // }

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


            res << '.' << std::setfill('0') << std::setw(FractionLength) << current_fraction * current_unit;
        }

        return res.str();
    }

    Basetype raw() const
    {
        return val;
    }

    static FixedPoint makeFx(Basetype v)
    {
        FixedPoint fp;
        fp.val = v;
        return fp;
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

};

// template <uint8_t FractionLength>
// std::ostream& operator << (std::ostream& out, const FixedPoint<FractionLength>& number)
// {
//     out << std::string(number);

//     return out;
// }

template <uint8_t FractionLength,typename Basetype, typename HelperType>
std::ostream& operator << (std::ostream& out, const FixedPoint<FractionLength, Basetype, HelperType> & number)
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
const FixedPoint<FractionLength> operator - (unsigned int l, const FixedPoint<FractionLength>& r)
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
using FixedPoint_1 = FixedPoint<15, uint16_t, uint32_t>;
using FixedPoint_16_16 = FixedPoint<16, uint32_t, uint64_t>;

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

template<typename InputInteger, typename OutputInteger>
std::pair<OutputInteger, uint8_t> getDivisionMultiplier(InputInteger divisor)
{
    if (!divisor)
    {
        throw std::invalid_argument("Division by zero is impossible");
    }

    if (divisor == 1)
    {
        return {1,0};
    }

    constexpr uint8_t n = sizeof(InputInteger) * CHAR_BIT;

    const double log_d_temp = std::log2(static_cast<double>(divisor));
    const uint8_t log_d = std::ceil(log_d_temp);

    if (log_d == std::floor(log_d_temp))
    {
        return {1, log_d};
    }

    OutputInteger res = std::ceil(static_cast<double>(static_cast<OutputInteger>(1) << (log_d + n)) / double(divisor));

    return {res, n + log_d};
}


template<typename InputInteger>
std::tuple<InputInteger, uint8_t, bool> getDivisionMultiplier(InputInteger divisor)
{
    if (!divisor)
    {
        throw std::invalid_argument("Division by zero is impossible");
    }

    if (divisor == 1)
    {
        return {1,0, false};
    }

    constexpr uint8_t n = sizeof(InputInteger) * CHAR_BIT;

    const double log_d_temp = std::log2(static_cast<double>(divisor));
    const uint8_t log_d = std::ceil(log_d_temp);

    if (log_d == std::floor(log_d_temp))
    {
        return {1, log_d, false};
    }

    uint64_t temp_low = (1UL << (log_d + n));
    uint64_t temp_hight = (1UL << log_d) | (1UL << (log_d + n));

    temp_hight /= divisor;
    temp_low /= divisor;

    uint8_t additionla_shift = log_d;

    while (additionla_shift)
    {

        if (temp_low /2 >= temp_hight/2)
        {
            break;
        }

        temp_low /= 2;
        temp_hight /= 2;

        --additionla_shift;
    }

    return {temp_hight, n + additionla_shift, temp_hight > std::numeric_limits<uint8_t>::max()};
}

void print_reciprocal(uint8_t bit_number)
{
    if (!bit_number || bit_number > CHAR_BIT - 1)
    {
        throw std::invalid_argument("It's expected to fit ony byte");
    }

    uint8_t count = 1 << bit_number;

    std::cout << "uint16_t reciprocals_" << (uint16_t)count << "[] = {";

    for(uint8_t i = 0; i < count; ++i)
    {
        uint8_t denominator = 1 << bit_number;
        denominator |= i;

        double reciprocal = static_cast<double>(1.) / static_cast<double>(denominator);

        uint8_t first_byte = std::scalbln(reciprocal, CHAR_BIT + bit_number);

        if (!first_byte)
        {
            first_byte = 0xFF;
        }

        // align it accrodint to Q 1.15 format used later in the calculations
        // particularly, we should leave one bit for a whole part
        std::cout  << ' ' << std::uppercase << std::hex << "0x" << ((uint16_t)first_byte << 7) << ", ";

        if (!(count % 15))
        {
            std::cout << std::endl;
        }
    }

    std::cout << "};" << std::endl;
}

//uint8_t reciprocals_8[] = { 0xFF,  0xE3,  0xCC,  0xBA,  0xAA,  0x9D,  0x92,  0x88, };
//uint8_t reciprocals_8[] = { 0x7F,  0x71,  0x66,  0x5D,  0x55,  0x4E,  0x49,  0x44, };

uint16_t reciprocals_8[] = { 0x7F80,  0x7180,  0x6600,  0x5D00,  0x5500,  0x4E80,  0x4900,  0x4400, };

constexpr uint8_t uint16_size = sizeof(uint16_t) * CHAR_BIT;
constexpr uint8_t shift_from_int = (sizeof(int) - sizeof(uint16_t)) * CHAR_BIT;

template<uint8_t BitCount = 3,uint16_t Reciprocals[] = reciprocals_8>
uint16_t divide(uint16_t u, uint16_t v)
{
    // __builtin_clz returns the number of the first not zero bit counting from the left, and the argument is widened to 4 bytes int
    int shift_to_left =  __builtin_clz(v) - shift_from_int;

    std::cout << "degree = " << shift_to_left << std::endl;

    // the first not zero bit should be the most left in uint16_t
    v <<= shift_to_left;


    FixedPoint_1 d = FixedPoint_1::makeFx(v);
    std::cout << "array index " << (v >> (8 + BitCount + 1)) - 8 << std::endl;
    // to look it up in the table we should first move the significant for us part to the right and then to zero the most significant bit
    FixedPoint_1 x = FixedPoint_1::makeFx(reciprocals_8[(v >> (8 + BitCount + 1)) - 8]);
    FixedPoint_1 one = FixedPoint_1::makeFx(1 << 15); //just one written in Q 1.15

    // two steps of Newton
    std::cout << "d " << d << std::endl;
    std::cout << "x " << x << std::endl;
    x = x + x * (one - d * x);
    std::cout << "x " << x << " d*x " << d * x  << " (one - d * x) " << (one - d * x) <<std::endl;
    x = x + x * (one - d * x);
    std::cout << "x " << x << std::endl;

    std::cout << "reciprocal of " << (v >> shift_to_left)  << " is " << x << std::endl;

    // just to wrap x into Q 16.16 format though x must be Q 1.15 probably the first bit is always zero...just move it one bit left
    FixedPoint_16_16 x_16_16 = FixedPoint_16_16::makeFx(x.raw() << 1);

    std::cout << "reciprocal of " << (v >> shift_to_left)  << " is " << x_16_16 << std::endl;

    //that's just for some tests because you can't wrap 32 bits u to this format
    FixedPoint_16_16 u_16_16 = FixedPoint_16_16::makeFx(u << 16);

    // when we wrapped v to Q 1.15 we made a tricky thing:
    // as a physical value we moved it to the left, but as a logical value we moved it to the right

    auto q = (x_16_16 * u_16_16) >> (uint16_size - shift_to_left - 1);
    std::cout << "u/v = " << q  << std::endl;


    //std::cout << "in raw types" << std::endl;

    uint64_t almostResult = uint64_t(x.raw() << 1) * uint64_t(u);
    uint16_t result = almostResult >> (uint16_size + (uint16_size - shift_to_left - 1));
    //std::cout << result << std::endl;

    v >>= shift_to_left;

    uint32_t reminder = u - result * v;

    if (reminder >= v)
    {
        reminder -= v;
        ++result;

        if (reminder >= v)
        {
            reminder -= v;
            ++result;

            if (reminder >= v)
            {
                reminder -= v;
                ++result;
            }
        }
    }

    return result;
}
// void divide(uint32_t n, uint32_t v)
// {
//     if (!v)

// }

int main()
{
    //divide(3, 7);
    //divide(35, 17);
    //divide(779, 19);
    print_reciprocal(3);

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

    // FixedPoint_20 x_prev{4};
    // FixedPoint_20 x_curr{4.25};

    // for (int i=2; i <= 100; ++i)
    // {
    //     auto temp = tricky_poly(x_curr, x_prev);
    //     x_prev = x_curr;
    //     x_curr = temp;

    //     std::cout << "tricky poly res = " << x_curr << std::endl;

    // }

    // std::cout << "tricky poly res = " << x_curr << std::endl;

    //std::cout << std::setprecision(20) << poly(x) << ' ' << poly(fx) << std::endl;
     //std::cout << poly1(x) << ' ' << poly1(fx) << std::endl;

    // Decimal<5,2> d{"-031.25"};
    // std::cout << std::string(d) << std::endl;

    // auto p = getDisionMultiplier<uint8_t, uint16_t>(static_cast<uint8_t>(10));

    // std::cout << p.first << " " << (uint16_t)(p.second) << std::endl;

    // for(uint8_t divisor = 1; divisor > 0; divisor++)
    // {
    //     auto [multiplier1, shift1] = getDisionMultiplier<uint8_t, uint16_t>(divisor);

    //     auto [multiplier2, shift2, overflow] = getDisionMultiplier(divisor);

    //     for(uint8_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         uint32_t res1 = static_cast<uint32_t>(numenator * multiplier1) >> shift1;

    //         if (res1 != numenator / divisor)
    //         {
    //             std::cout << "panic: did something went wrong?" << std::endl;
    //         }

    //         uint32_t res2 = 0;

    //         if (!overflow)
    //         {
    //             res2 = (numenator * multiplier2) >> shift2;
    //         }
    //         else
    //         {
    //             uint16_t temp = (numenator * multiplier2);
    //             uint8_t temp2 = temp >> 8;
    //             uint16_t temp3 = (numenator - temp2) >> 1;
    //             res2 = (temp2 + temp3) >> (shift2 - 8 - 1);
    //         }

    //         if (res2 != numenator / divisor)
    //         {
    //             std::cout << "panic 2!" << std::endl;
    //         }
    //     }
    // }


    // for(uint8_t divisor = 1; divisor > 0; divisor++)
    // {
    //     auto [multiplier, shift] = getDivisionMultiplier<uint8_t, uint16_t>(divisor);

    //     for(uint8_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         uint32_t res = static_cast<uint32_t>(numenator * multiplier) >> shift;

    //         if (res != numenator / divisor)
    //         {
    //             std::cout << "panic: did something went wrong?" << std::endl;
    //         }
    //     }
    // }

    // auto [coeff, shift, _] = getDivisionMultiplier(static_cast<uint8_t>(10));
    // std::cout << (uint16_t)coeff << " " << (uint16_t)(shift) << std::endl;



    for(uint16_t divisor = 1; divisor > 0; divisor++)
    {
        for(uint16_t numenator = 1; numenator > 0; numenator++)
        {
            if (divide(numenator, divisor) != numenator / divisor)
            {
                std::cout << "panic: did something went wrong?" << std::endl;
            }
        }
    }


    return 0;
}
