#ifndef ATOMBONDGUESSER_HPP
#define ATOMBONDGUESSER_HPP

#include <memory>

namespace mol {

class AtomSel;

namespace internal {

class AtomBondGuesser
{
public:
    void apply(std::shared_ptr<AtomSel> atoms) const;

private:
};

} // namespace internal
} // namespace mol

#endif // ATOMBONDGUESSER_HPP
