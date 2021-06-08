#ifndef MOLSYSTEM_HPP
#define MOLSYSTEM_HPP

#include <string>
#include <memory>

namespace mol {

namespace internal {
class AtomData;
}


class MolSystem
{
public:
    MolSystem(std::string const& topology);
    void add_trajectory(std::string const& filename);

private:
    std::shared_ptr<internal::AtomData> m_atoms;
};

} // namespace mol

#endif // MOLSYSTEM_HPP
