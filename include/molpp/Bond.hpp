#ifndef BOND_HPP
#define BOND_HPP

#include <cstdlib>

namespace mol {

class Bond
{
public:
    Bond() = delete;
    Bond(index_t const atom1, index_t const atom2);
    index_t atom1() const { return m_atom1; }
    index_t atom2() const { return m_atom2; }
    int order() const { return m_order; }
    void set_order(int const order) { m_order = order; }
    bool guessed() const { return m_guessed; }
    void set_guessed(float const is_guessed) { m_guessed = is_guessed; }
    bool guessed_order() const { return m_guessed_order; }
    void set_guessed_order(float const is_guessed) { m_guessed_order = is_guessed; }
    bool aromatic() const { return m_aromatic; }
    void set_aromatic(bool const is_aromatic) { m_aromatic = is_aromatic; }

private:
    bool m_guessed;
    bool m_guessed_order;
    bool m_aromatic;
    int m_order;
    index_t m_atom1;
    index_t m_atom2;
};

inline Bond::Bond(index_t const atom1, index_t const atom2)
: m_guessed { true },
  m_guessed_order { true },
  m_aromatic { false },
  m_order { 0 },
  m_atom1 { (atom1 < atom2) ? atom1 : atom2 },
  m_atom2 { (atom1 > atom2) ? atom1 : atom2 }
{}

} // namespace mol

#endif // BOND_HPP
