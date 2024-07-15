#ifndef MOLPP_GUESSERS_RESIDUEBONDGUESSER_HPP
#define MOLPP_GUESSERS_RESIDUEBONDGUESSER_HPP

namespace mol
{

class ResidueSel;

namespace internal
{

class ResidueBondGuesser
{
public:
    void apply(ResidueSel& residue_sel) const;

private:
};

} // namespace internal
} // namespace mol

#endif // MOLPP_GUESSERS_RESIDUEBONDGUESSER_HPP
