#ifndef DECIMAL_H
#define DECIMAL_H

#include <stdint.h>
#include <array>
#include <cstring>
#include <stdexcept>
#include <sstream>

template<uint8_t Lenght, uint8_t FractionLength>
class Decimal
{
    static_assert(FractionLength > 0);
    static_assert(Lenght > FractionLength);
public:
    explicit Decimal(const char* str)
    {
        val.fill(0);

        if (!str)
        {
            throw std::invalid_argument("nullptr isn't expected here");
        }

        auto len = std::strlen(str);

        if (!len ||len > Lenght + 2)
        {
            throw std::invalid_argument("Given string isn't suitable");
        }

        if (!std::strcmp(str, "."))
        {
           throw std::invalid_argument("It can'be just a dot");
        }

        int8_t sign = 1;

        const auto leading = str[0];

        if (leading == '-')
        {
            sign = -1;
            --len;
            str++;
        }
        else if (leading == '+')
        {
            --len;
            str++;
        }

        while (*str == '0')
        {
            ++str;
            --len;
        }

        int16_t dot_pos = -1;
        uint8_t start_index = 0;
        auto dot = std::strchr(str, '.');

        if (dot)
        {
            dot_pos = dot - str;
        }

        if (dot_pos >= 0)
        {
            start_index = DecimalLength - dot_pos;
        }
        else
        {
            start_index = DecimalLength - len;
        }

        for (uint8_t i = 0; i < len; ++i)
        {
            const auto ch = str[i];

            if (ch == '.')
            {
                continue;
            }

            val[start_index++] = ch - '0';
        }

        if (sign == -1)
        {
            val[0] |= SignMask;
        }
    }

    operator std::string() const
    {
        std::stringstream res;
        bool found_not_zero = false;

        for (uint8_t i = 0; i < val.size(); ++i)
        {
            uint8_t elem = val[i];

            if (!i && elem & SignMask)
            {
                res << '-';
                elem &= ResetSignMask;
            }

            if (i == DecimalLength)
            {
                res << '.';
                found_not_zero = true;
            }

            if (!elem && !found_not_zero)
            {
                continue;
            }

            found_not_zero = true;
            res << static_cast<char>((elem) + '0');
        }

        return res.str();
    }

private:
    std::array<uint8_t, Lenght> val;
    constexpr static uint8_t DecimalLength = Lenght - FractionLength;
    constexpr static uint8_t SignMask = 0x80;
    constexpr static uint8_t ResetSignMask = 0x7F;
};


template<uint8_t Lenght, uint8_t FractionLength>
std::ostream& operator << (std::ostream& out, const Decimal<Lenght,FractionLength>& number)
{
    out << std::string(number);

    return out;
}



#endif // DECIMAL_H
