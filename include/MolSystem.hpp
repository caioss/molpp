#ifndef MOLSYSTEM_HPP
#define MOLSYSTEM_HPP

#include "AtomSel.hpp"
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
    void add_trajectory(std::string const& file_name, int begin=0, int end=-1, int step=1);
    std::shared_ptr<AtomSel> all() const;
    std::shared_ptr<AtomSel> select(std::vector<size_t> const &indices) const;

private:
    std::shared_ptr<internal::MolData> m_data;
    std::shared_ptr<AtomSel> m_all;
};

} // namespace mol

#endif // MOLSYSTEM_HPP
