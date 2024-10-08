#include <chrono>
#include <stdint.h>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <climits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <x86intrin.h>

#include "decimal.h"
#include "uint16_reciprocals.h"
#include "uint32_reciprocals.h"
#include "Rational.h"

// FractionLength - how many bits are after a period
template <uint8_t FractionLength, typename Basetype = int32_t, typename HelperType = uint64_t>
class FixedPoint
{
public:
    static_assert(sizeof(HelperType) == 2 * sizeof(Basetype));
    static_assert(FractionLength > 0);
    //static_assert(FractionLength < sizeof(Basetype) * CHAR_BIT);
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

using FixedPoint_0_16 = FixedPoint<16, uint16_t, uint32_t>;

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

    const uint8_t count = 1 << bit_number;

    //double max_error = 0.0;

    std::cout << "uint16_t reciprocals_" << (uint16_t)count << "[] = {";

    for(uint8_t i = 0; i < count; ++i)
    {
        uint8_t denominator = count;
        denominator |= i;

        double reciprocal = static_cast<double>(1.) / static_cast<double>(denominator);

        uint8_t first_byte = std::scalbln(reciprocal, CHAR_BIT + bit_number);

        if (!first_byte)
        {
            first_byte = 0xFF;
        }

        std::cout  << ' ' << std::uppercase << std::hex << "0x" << ((uint16_t)first_byte << 8) << ", ";

        if (!(count % 15))
        {
            std::cout << std::endl;
        }

        //uint16_t nextDenominator = ((denominator + 1) << 12);
        //std::cout << std::endl << (first_byte / 256.) << " " << (nextDenominator / 32768.) << std::endl;
        //max_error = std::max(max_error, std::abs(1. - (first_byte / 256.) * (nextDenominator / 32768.)));
        //std::cout << "error " << max_error << std::endl;
    }

    std::cout << "};" << std::endl;

    //std::cout << "max_error = " << max_error << std::endl;
}

void print_reciprocal_uint32t(uint8_t bit_number, const char* file_name)
{
    if (!file_name || ! std::strlen(file_name))
    {
        return;
    }


    if (!bit_number || bit_number > 2 * CHAR_BIT - 1)
    {
        throw std::invalid_argument("It's expected to fit two bytes");
    }

    const uint16_t count = 1 << bit_number;
    std::ofstream out(file_name);


    out << "uint32_t reciprocals_" << (uint16_t)count << "[] = {";

    for(uint16_t i = 0; i < count; ++i)
    {
        uint16_t denominator = count;
        denominator |= i;

        double reciprocal = static_cast<double>(1.) / static_cast<double>(denominator);

        uint16_t first_bytes = std::scalbln(reciprocal, 2 * CHAR_BIT + bit_number);

        if (!first_bytes)
        {
            first_bytes = 0xFFFF;
        }

        out  << ' ' << std::uppercase << std::hex << "0x" << ((uint32_t)first_bytes << 16) << ", ";

        if (!(i % 15))
        {
            out << std::endl;
        }
    }

    out << "};" << std::endl;

    out.close();
}


void print_all_reciprocals(const char* file_name)
{
    if (!file_name || ! std::strlen(file_name))
    {
        return;
    }

    std::ofstream out(file_name);

    uint16_t i = 0;

    out << "uint16_t all_reciprocals" << "[] = {";

    while (true)
    {
        if (!i)
        {
            out << "0x0000, ";
            i++;
            continue;
        }

        if (i == 1)
        {
            out << "0xFFFF, ";
            i++;
            continue;
        }

        double reciprocal = static_cast<double>(1.) / static_cast<double>(i);

        uint16_t two_bytes = std::scalbln(reciprocal, CHAR_BIT * sizeof(uint16_t));

        out << "0x" << std::uppercase << std::hex << std::setfill('0') << std::setw(4) << two_bytes << ", ";

        if (!(i % 15))
        {
            out << std::endl;
        }

        if (i == std::numeric_limits<uint16_t>::max())
        {
            break;
        }

        i++;
    }

    out << "};" << std::endl;

    out.close();
}

uint16_t reciprocals_8[] = { 0xFF00,  0xE300,  0xCC00,  0xBA00,  0xAA00,  0x9D00,  0x9200,  0x8800, };

uint16_t reciprocals_128[] = { 0xFF00,  0xFE00,  0xFC00,  0xFA00,  0xF800,  0xF600,  0xF400,  0xF200,  0xF000,  0xEF00,  0xED00,  0xEB00,  0xEA00,  0xE800,  0xE600,  0xE500,
                               0xE300,  0xE100,  0xE000,  0xDE00,  0xDD00,  0xDB00,  0xDA00,  0xD900,  0xD700,  0xD600,  0xD400,  0xD300,  0xD200,  0xD000,  0xCF00,  0xCE00,
                               0xCC00,  0xCB00,  0xCA00,  0xC900,  0xC700,  0xC600,  0xC500,  0xC400,  0xC300,  0xC100,  0xC000,  0xBF00,  0xBE00,  0xBD00,  0xBC00,  0xBB00,
                               0xBA00,  0xB900,  0xB800,  0xB700,  0xB600,  0xB500,  0xB400,  0xB300,  0xB200,  0xB100,  0xB000,  0xAF00,  0xAE00,  0xAD00,  0xAC00,  0xAB00,
                               0xAA00,  0xA900,  0xA800,  0xA800,  0xA700,  0xA600,  0xA500,  0xA400,  0xA300,  0xA300,  0xA200,  0xA100,  0xA000,  0x9F00,  0x9F00,  0x9E00,
                               0x9D00,  0x9C00,  0x9C00,  0x9B00,  0x9A00,  0x9900,  0x9900,  0x9800,  0x9700,  0x9700,  0x9600,  0x9500,  0x9400,  0x9400,  0x9300,  0x9200,
                               0x9200,  0x9100,  0x9000,  0x9000,  0x8F00,  0x8F00,  0x8E00,  0x8D00,  0x8D00,  0x8C00,  0x8C00,  0x8B00,  0x8A00,  0x8A00,  0x8900,  0x8900,
                               0x8800,  0x8700,  0x8700,  0x8600,  0x8600,  0x8500,  0x8500,  0x8400,  0x8400,  0x8300,  0x8300,  0x8200,  0x8200,  0x8100,  0x8100,  0x8000, };

constexpr uint8_t uint16_size = sizeof(uint16_t) * CHAR_BIT;
constexpr uint8_t shift_from_int = (sizeof(int) - sizeof(uint16_t)) * CHAR_BIT;

uint16_t multHigherHalf(uint16_t a, uint16_t b)
{
    uint32_t res = a * b;
    res >>= 16;
    return res;
}

uint32_t multHigherHalf32(uint32_t a, uint32_t b)
{
    uint64_t res = (uint64_t)a * (uint64_t)b;
    res >>= 32;
    return res;
}


template<uint8_t BitCount = 3,uint16_t Reciprocals[] = reciprocals_8>
uint16_t divide(uint16_t u, uint16_t v)
{
    // __builtin_clz returns the number of the first not zero bit counting from the left, and the argument is widened to 4 bytes int
    int shift_to_left =  __builtin_clz(v) - shift_from_int;

    // the first not zero bit should be the most left in uint16_t
    v <<= shift_to_left;
    // to look it up in the table we should first move the significant for us part to the right and then to zero the most significant bit
    uint16_t x = Reciprocals[(v >> (8 + BitCount + 1)) - 8];

    // two steps of Newton
    x = multHigherHalf(static_cast<uint16_t>(- multHigherHalf(v, x)), x) << 1;
    x = multHigherHalf(static_cast<uint16_t>(- multHigherHalf(v, x)), x) << 1;

    uint16_t q = multHigherHalf(x, u);
    q >>= uint16_size - shift_to_left - 1;

    if (q)
    {
        --q;
    }

    v >>= shift_to_left;

    uint32_t reminder = u - q * v;

    if (reminder >= v)
    {
        reminder -= v;
        ++q;

        if (reminder >= v)
        {
            reminder -= v;
            ++q;

            if (reminder >= v)
            {
                ++q;
            }
        }
    }

    return q;
}

template<uint8_t BitCount = 7,uint16_t Reciprocals[] = reciprocals_128>
uint16_t divide7(uint16_t u, uint16_t v)
{
    // __builtin_clz returns the number of the first not zero bit counting from the left, and the argument is widened to 4 bytes int
    int shift_to_left =  __builtin_clz(v) - shift_from_int;

    // the first not zero bit should be the most left in uint16_t
    v <<= shift_to_left;

    uint16_t x = Reciprocals[((v >> (8)) & 0x7F)];

    // a step of Newton
    x = multHigherHalf(static_cast<uint16_t>(- multHigherHalf(v, x)), x) << 1;

    uint16_t q = multHigherHalf(x, u);
    q >>= uint16_size - shift_to_left - 1;

    if (q)
    {
        --q;
    }

    v >>= shift_to_left;

    uint32_t reminder = u - q * v;

    if (reminder >= v)
    {
        reminder -= v;
        ++q;

        if (reminder >= v)
        {
            reminder -= v;
            ++q;

            if (reminder >= v)
            {
                ++q;
            }
        }
    }

    return q;
}


uint16_t new_divide(uint16_t u, uint16_t v)
{
    // __builtin_clz returns the number of the first not zero bit counting from the left, and the argument is widened to 4 bytes int
    int shift_to_left =  __builtin_clz(v) - 16;

    // the first not zero bit should be the most left in uint16_t (normalization)
    v <<= shift_to_left;

    // to look it up in the table we should first move the significant for us part to the right and then to zero the most significant bit
    uint16_t x = reciprocals_8[(v >> (8 + 3 + 1)) - 8];

    // two steps of Newton
    x = multHigherHalf(static_cast<uint16_t>((-(v * x)) >> 16), x) << 1;
    x = multHigherHalf(static_cast<uint16_t>((-(v * x)) >> 16), x) << 1;

    uint16_t q = multHigherHalf(x, u);
    q >>= 16 - shift_to_left - 1;

    //get normalization back
    v >>= shift_to_left;

    uint32_t reminder = u - q * v;

    // the reminder should be always less than v
    if (reminder >= v)
    {
        reminder -= v;
        ++q;

        if (reminder >= v)
        {
            ++q;

        }
    }

    return q;
}


uint16_t new_divide7(uint16_t u, uint16_t v)
{
    // __builtin_clz returns the number of the first not zero bit counting from the left, and the argument is widened to 4 bytes int
    int shift_to_left =  __builtin_clz(v) - 16;

    // the first not zero bit should be the most left in uint16_t
    v <<= shift_to_left;

    // to look it up in the table we should first move the significant for us part to the right and then to zero the most significant bit
    uint16_t x = reciprocals_128[((v >> (8)) & 0x7F)];

    // a step of Newton
    x = multHigherHalf(static_cast<uint16_t>((-(v * x)) >> 16), x) << 1;

    uint16_t q = multHigherHalf(x, u);
    q >>= 16 - shift_to_left - 1;

    v >>= shift_to_left;

    uint32_t reminder = u - q * v;

    if (reminder >= v)
    {
        reminder -= v;
        ++q;

        if (reminder >= v)
        {
            ++q;
        }
    }

    return q;
}

uint16_t divide_wo(uint16_t u, uint16_t v)
{
    uint16_t x = all_reciprocals[v];
    uint16_t q = multHigherHalf(x, u);

    uint32_t reminder = u - q * v;

    if (reminder >= v)
    {
        q++;
    }

    return q;
}

uint32_t divide32(uint32_t u, uint32_t v)
{
    // __builtin_clz returns the number of the first not zero bit counting from the left
    int shift_to_left =  __builtin_clz(v);

    // the first not zero bit should be the most left in uint32_t
    v <<= shift_to_left;

    // to look it up in the table we should first move the significant for us part to the right and then to zero the most significant bit
    uint32_t x = reciprocals_32768[(v >> 16) & 0x7FFF];

    // a step of Newton
    x = multHigherHalf32((-(static_cast<uint64_t>(v) * x) >> 32) , x) << 1;

    uint32_t q = multHigherHalf32(x, u);
    q >>= 32 - shift_to_left - 1;

    v >>= shift_to_left;

    uint32_t reminder = u - q * v;

    if (reminder >= v)
    {
        reminder -= v;
        ++q;

        if (reminder >= v)
        {
            ++q;
        }
    }

    return q;
}

template <typename T>
Rational<T> muller_func(Rational<T> x0, Rational<T> x1, uint8_t count)
{
    Rational<T> x_new;

    for (uint8_t i = 0; i < count; ++i)
    {
        x_new =  Rational<T>(108) - (Rational<T>(815) - Rational<T>(1500) / x0) / x1;

        x0 = x1;
        x1 = x_new;
    }

    return x_new;
}

int main()
{

    std::cout << muller_func(long_rational(4), long_rational(17, 4), 22) << std::endl;
    // tiny_rational zero;
    // tiny_rational one{1};

    // tiny_rational one_third{1, 3};
    // tiny_rational one_half{1, 2};
    // tiny_rational minus_one_half{1, 2, true};

    // std::cout << tiny_rational(3, 4) * tiny_rational(16, 9) << std::endl;



    // std::cout << zero << " " << (zero + zero) << " " << (zero * zero) << std::endl;

    // std::cout << one << " " << (one + one) << " " << (one * one) <<
    //           " " << (one + zero) << " " << (one * zero)
    //           << std::endl;

    // std::cout << one_third << " " << one_half << " " << (one_third + one_half) << " " <<
    //         (one_third - one_half) << " " << (one_third * one_half) << " " <<
    //          (one_third / one_half) << std::endl;

    // std::cout << one_half << " " << minus_one_half << " " << (one_half + minus_one_half)
    //           << " " << (one_half - minus_one_half) << " " << ( minus_one_half - one_half)
    //           << " " << (one_half * minus_one_half)  << " " << (one_half / minus_one_half) << std::endl;

     // std::cout << divide32(2147483649, 1) << std::endl;
     // return 0;
    //print_reciprocal_uint32t(15, "uint32_reciprocals.h");
    //print_all_reciprocals("uint16_reciprocals.h");
    //divide(3, 7);
    // std::cout << divide(36198, 53) << std::endl;
    // return -1;
    //std::cout << divide(36198, 53) << std::endl;
    //std::cout << divide(32769, 1) << std::endl;
    //print_reciprocal(3);
     //return 0;
    //std::cout << divide7(102, 3) << std::endl;
    //return 0;

    //std::cout << divide(32770, 2) << std::endl;

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



    // for(uint16_t divisor = 1; divisor > 0; divisor++)
    // {
    //     for(uint16_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         if (new_divide(numenator, divisor) != numenator / divisor)
    //         {
    //             std::cout << "panic: did something went wrong?" << std::endl;
    //             std::cout << "numenator = "  << numenator << " divisor = " << divisor  << " result = " << new_divide(numenator, divisor) << std::endl;
    //             return -1;
    //         }
    //     }
    // }

    // for(uint16_t divisor = 1; divisor > 0; divisor++)
    // {
    //     for(uint16_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         volatile uint16_t res = 0;
    //         res = new_divide(numenator, divisor);
    //         (void)(res);
    //     }
    // }

    // uint64_t avg = 0, max = 0, min = -1;

    // for(uint16_t divisor = 1; divisor > 0; divisor++)
    // {
    //     for(uint16_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         const uint64_t st = __rdtsc();
    //         volatile uint16_t res = numenator / divisor;
    //         const uint64_t et = __rdtsc() - st;
    //         (void)(res);

    //         avg += et;

    //         max = std::max(max, et);
    //         min = std::min(min, et);
    //     }
    // }

    // std::cout << "max = " << max << " min = " << min << " avg = " << (avg / ((uint64_t)std::numeric_limits<uint16_t>::max() * std::numeric_limits<uint16_t>::max())) << std::endl;



    // uint64_t avg = 0, max = 0, min = -1;

    // for(uint16_t divisor = 1; divisor > 0; divisor++)
    // {
    //     for(uint16_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         const uint64_t st = __rdtsc();
    //         volatile uint16_t res = new_divide7(numenator, divisor);
    //         const uint64_t et = __rdtsc() - st;
    //         (void)(res);

    //         avg += et;

    //         max = std::max(max, et);
    //         min = std::min(min, et);
    //     }
    // }

    // std::cout << "max = " << max << " min = " << min << " avg = " << (avg / ((uint64_t)std::numeric_limits<uint16_t>::max() * std::numeric_limits<uint16_t>::max())) << std::endl;

    // auto start = std::chrono::high_resolution_clock::now();

    // for(uint16_t divisor = 1; divisor > 0; divisor++)
    // {
    //     for(uint16_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         volatile uint16_t res = divide_wo(numenator, divisor);
    //         //volatile uint16_t res = numenator / divisor;
    //         (void)(res);
    //     }
    // }

    // auto stop = std::chrono::high_resolution_clock::now();

    // std::cout << "duration = " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;

    // for(uint16_t divisor = 1; divisor > 0; divisor++)
    // {
    //     for(uint16_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         if (new_divide7(numenator, divisor) != numenator / divisor)
    //         {
    //             std::cout << "panic: did something went wrong?" << std::endl;
    //             std::cout << "numenator = "  << numenator << " divisor = " << divisor  << " result = " << new_divide7(numenator, divisor) << std::endl;
    //             return -1;
    //         }
    //     }
    // }


    // auto start = std::chrono::high_resolution_clock::now();

    // for(uint16_t divisor = 1; divisor > 0; divisor++)
    // {
    //     for(uint16_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         volatile uint16_t res = divide_wo(numenator, divisor);
    //         (void)(res);
    //     }
    // }

    // auto stop = std::chrono::high_resolution_clock::now();

    // std::cout << "duration = " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;


    // auto start = std::chrono::high_resolution_clock::now();

    // for(uint16_t divisor = 1; divisor > 0; divisor++)
    // {
    //     for(uint16_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         volatile uint16_t res = numenator / divisor;
    //         (void)(res);
    //     }
    // }

    // auto stop = std::chrono::high_resolution_clock::now();

    // std::cout << "duration = " << std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count() << std::endl;


    // for(uint32_t divisor = 1; divisor > 0; divisor++)
    // {
    //     for(uint32_t numenator = 1; numenator > 0; numenator++)
    //     {
    //         if (divide32(numenator, divisor) != numenator / divisor)
    //         {
    //             std::cout << "panic: did something went wrong?" << std::endl;
    //             std::cout << "numenator = "  << numenator << " divisor = " << divisor  << " result = " << divide32(numenator, divisor) << std::endl;
    //             return -1;
    //         }
    //     }
    // }


    return 0;
}
