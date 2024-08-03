#ifndef LONGTYPES_H
#define LONGTYPES_H

#include <vector>
#include <stdexcept>
#include <cstring>

template<typename BaseType>
class UnsignedLong
{
public:
    UnsignedLong(const char* str)
    {
        if (!str || !std::strlen(str))
        {
             throw std::invalid_argument("nullptr isn't expected here");
        }

        if ('0' == *str)
        {
            throw std::invalid_argument("leading zeroes aren't allowed");
        }

        auto len = std::strlen(str);

        val.reserve(len / decimal_binary);



    }
private:
    std::vector<BaseType> val;
    constexpr static double decimal_binary = 3.5;
};

#endif // LONGTYPES_H
