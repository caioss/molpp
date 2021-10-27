#ifndef RESIDUEBONDGUESSER_HPP
#define RESIDUEBONDGUESSER_HPP

namespace mol {

class ResidueSel;

namespace internal {

class ResidueBondGuesser
{
public:
    void apply(ResidueSel &residues) const;

private:
};

} // namespace internal
} // namespace mol

#endif // RESIDUEBONDGUESSER_HPP
