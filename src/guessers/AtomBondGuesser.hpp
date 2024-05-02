#ifndef MOLPP_GUESSERS_ATOMBONDGUESSER_HPP
#define MOLPP_GUESSERS_ATOMBONDGUESSER_HPP

#include <memory>

namespace mol
{

class AtomSel;

namespace internal
{

class AtomBondGuesser
{
public:
    void apply(AtomSel& atoms) const;

private:
};

} // namespace internal
} // namespace mol

#endif // MOLPP_GUESSERS_ATOMBONDGUESSER_HPP
