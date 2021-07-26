#ifndef MOLFILEREADER_HPP
#define MOLFILEREADER_HPP

#include "MolReader.hpp"
#include "molfile_plugin.h"
#include <string>

namespace mol::internal {

class MolfileReader : public MolReader
{
public:
    MolfileReader() = delete;
    MolfileReader(std::string const &file_ext);
    ~MolfileReader();
    static bool can_read(std::string const &file_ext);
    bool has_topology() const override;
    bool has_trajectory() const override;
    bool has_bonds() const override;
    bool has_trajectory_metadata() const override;
    Status open(const std::string &file_name) override;
    void close() override;
    std::shared_ptr<MolData> read_atoms() override;
    Status check_timestep_read(std::shared_ptr<MolData> atom_data) override;
    Status skip_timestep(std::shared_ptr<MolData> atom_data) override;
    Status read_timestep(std::shared_ptr<MolData> atom_data) override;

private:
    int m_num_atoms;
    void *m_handle;
    std::string m_name;
    molfile_plugin_t *m_plugin;
};

} // namespace mol::internal

#endif // MOLFILEREADER_HPP
