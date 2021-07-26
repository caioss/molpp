#ifndef RESIDUE_HPP
#define RESIDUE_HPP

#include "MolppCore.hpp"
#include <memory>

namespace mol {

class Atom;
class Bond;
class AtomSel;

namespace internal {
    class AtomData;
}

class Residue
{
public:
    Residue(size_t const index, size_t const frame,std::shared_ptr<internal::AtomData> data)
    : m_index { index },
      m_frame { frame },
      m_data { data }
    {}

    bool operator==(Residue const &other) const
    {
        return m_data == other.m_data && m_index == other.m_index && m_frame == other.m_frame;
    }

    size_t index() const
    {
        return m_index;
    }

    size_t frame() const
    {
        return m_frame;
    }

    int resid() const;
    void set_resid(int const &resid);

    std::string resname() const;
    void set_resname(std::string const &resname);

    std::string segid() const;
    void set_segid(std::string const &segid);

    std::string chain() const;
    void set_chain(std::string const &chain);

    void add_atom(size_t index);
    void add_atom(Atom &atom);
    std::shared_ptr<AtomSel> atoms();
    size_t size() const;

private:
    size_t m_index;
    size_t m_frame;
    std::shared_ptr<internal::AtomData> m_data;
};

} // namespace mol

#endif // RESIDUE_HPP
