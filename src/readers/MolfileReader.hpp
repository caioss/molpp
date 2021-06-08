#ifndef MOLFILEREADER_HPP
#define MOLFILEREADER_HPP

#include "MolReader.hpp"
#include "molfile_plugin.h"
#include <string>
#include <unordered_set>

namespace mol::internal {

class MolfileReader : public MolReader
{
public:
    MolfileReader() = delete;
    MolfileReader(std::string const &name);
    ~MolfileReader();
    bool can_read(std::string const &filename) const override;
    bool has_topology() const override;
    bool has_trajectory() const override;
    bool has_bonds() const override;
    bool has_trajectory_metadata() const override;
    bool open(const std::string &file_name) override;
    void close() override;
    std::shared_ptr<AtomData> read_atoms() override;
    bool skip_timestep(std::shared_ptr<AtomData> atom_data);
    bool read_timestep(std::shared_ptr<AtomData> atom_data);

private:
    int m_num_atoms;
    void *m_handle;
    std::string m_name;
    molfile_plugin_t *m_plugin;
    std::unordered_set<std::string> m_extensions;
};

} // namespace mol::internal

#endif // MOLFILEREADER_HPP
