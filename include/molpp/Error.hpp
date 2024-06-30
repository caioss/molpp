#ifndef MOLPP_ERROR_HPP
#define MOLPP_ERROR_HPP

#include <stdexcept>

namespace mol
{

class Error : public std::runtime_error
{
public:
    Error(const std::string& what = "")
    : std::runtime_error(what)
    {}
};

} // namespace mol

#endif // MOLPP_ERROR_HPP
