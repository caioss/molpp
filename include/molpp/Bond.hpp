#ifndef BOND_HPP
#define BOND_HPP

#include <cstdlib>

namespace mol {

class Bond
{
public:
    enum Order {Unknown = 0, Single = 1, Double = 2, Triple = 3, Aromatic = 5};

    Bond() = delete;
    Bond(size_t const atom1, size_t const atom2);
    size_t atom1() { return m_atom1; }
    size_t atom2() { return m_atom2; }
    Order order() const { return m_order; }
    void set_order(float const order);
    bool guessed() const { return m_guessed; }
    void set_guessed(float const is_guessed) { m_guessed = is_guessed; }
    bool guessed_order() const { return m_guessed_order; }
    void set_guessed_order(float const is_guessed) { m_guessed_order = is_guessed; }

private:
    bool m_guessed;
    bool m_guessed_order;
    Order m_order;
    size_t m_atom1;
    size_t m_atom2;
};

inline Bond::Bond(size_t const atom1, size_t const atom2)
: m_guessed { true },
  m_guessed_order { true },
  m_order { Bond::Unknown },
  m_atom1 { (atom1 < atom2) ? atom1 : atom2 },
  m_atom2 { (atom1 > atom2) ? atom1 : atom2 }
{}

inline void Bond::set_order(float const order)
{
    if (order == 1)
    {
        m_order = Bond::Single;
    }
    else if (order == 2)
    {
        m_order = Bond::Double;
    }
    else if (order == 3)
    {
        m_order = Bond::Triple;
    }
    else if (order > 1 && order < 2)
    {
        m_order = Bond::Aromatic;
    }
    else
    {
        m_order = Bond::Unknown;
    }
}

} // namespace mol

#endif // BOND_HPP
