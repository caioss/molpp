#ifndef MATCHERS_HPP
#define MATCHERS_HPP

#include <gmock/gmock.h>
#include <cmath>

using ::std::get;
using namespace testing;

MATCHER_P(Prop, func, "")
{
    return (std::get<0>(arg).*func)() == std::get<1>(arg);
}

MATCHER_P2(PropFloat, func, epsilon, "")
{
    return std::fabs((std::get<0>(arg).*func)() - std::get<1>(arg)) < epsilon;
}

#endif // MATCHERS_HPP
