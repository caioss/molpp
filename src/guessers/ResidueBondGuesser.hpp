#ifndef RESIDUEBONDGUESSER_HPP
#define RESIDUEBONDGUESSER_HPP

#include <memory>

namespace mol {

class ResidueSel;

namespace internal {

class ResidueBondGuesser
{
public:
    void apply(std::shared_ptr<ResidueSel> residues) const;

private:
};

} // namespace internal
} // namespace mol

#endif // RESIDUEBONDGUESSER_HPP
