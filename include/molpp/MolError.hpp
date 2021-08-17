#ifndef MOLERROR_HPP
#define MOLERROR_HPP

#include <stdexcept>

namespace mol
{

class MolError : public std::runtime_error
{
public:
    MolError(const std::string &what = "")
    : std::runtime_error(what) {}
};

} // namespace mol

#endif // MOLARTISTERROR_HPP
