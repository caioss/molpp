#ifndef MOLSYSTEM_HPP
#define MOLSYSTEM_HPP

#include <molpp/AtomSel.hpp>
#include <molpp/MolppCore.hpp>
#include <molpp/AtomSelector.hpp>
#include <string>
#include <memory>
#include <vector>

namespace mol {

namespace internal {
class MolData;
}

class MolSystem
{
public:
    MolSystem(std::string const& topology);
    MolSystem(MolSystem&& other);
    ~MolSystem();
    void add_trajectory(std::string const& file_name, int begin=0, int end=-1, int step=1);
    AtomSel atoms(Frame const frame = std::nullopt) const;
    AtomSel select(std::vector<index_t> const &indices, Frame const frame = std::nullopt) const;
    AtomSel select(std::string const &selection, Frame const frame = std::nullopt) const;
    AtomSelector selector(std::string const& selection) const;
    void reset_bonds();
    void guess_bonds(Frame const frame);

private:
    std::unique_ptr<internal::MolData> m_data;
};

} // namespace mol

#endif // MOLSYSTEM_HPP
