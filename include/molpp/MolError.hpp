#ifndef MOLPP_MOLERROR_HPP
#define MOLPP_MOLERROR_HPP

#include <stdexcept>

namespace mol
{

class MolError : public std::runtime_error
{
public:
    MolError(const std::string& what = "")
    : std::runtime_error(what)
    {}
};

} // namespace mol

#endif // MOLPP_MOLERROR_HPP
