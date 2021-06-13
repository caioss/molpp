#ifndef MOLSYSTEM_HPP
#define MOLSYSTEM_HPP

#include <string>
#include <memory>

namespace mol {

namespace internal {
class AtomData;
}

class AtomSel;

class MolSystem
{
public:
    MolSystem(std::string const& topology);
    void add_trajectory(std::string const& file_name, int begin=0, int end=-1, int step=1);

private:
    std::shared_ptr<internal::AtomData> m_atoms;
};

} // namespace mol

#endif // MOLSYSTEM_HPP
