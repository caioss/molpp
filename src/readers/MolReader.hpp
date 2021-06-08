#ifndef MOLREADER_HPP
#define MOLREADER_HPP

#include <string>
#include <memory>
#include <vector>

namespace mol::internal
{

class AtomData;

class MolReader
{
public:
    static std::shared_ptr<MolReader> from_file(std::string const &file_name);
    virtual ~MolReader() {};
    std::shared_ptr<AtomData> read_topology(std::string const &file_name);
    bool read_trajectory(std::string const &file_name,
                         std::shared_ptr<AtomData> atom_data,
                         int begin=0, int end=-1, int step=1);
    virtual bool can_read(const std::string &file_name) const = 0;
    virtual bool has_topology() const = 0;
    virtual bool has_trajectory() const = 0;
    virtual bool has_bonds() const = 0;
    virtual bool has_trajectory_metadata() const = 0;
    virtual bool open(const std::string &file_name) = 0;
    virtual void close() = 0;
    virtual std::shared_ptr<AtomData> read_atoms() = 0;
    virtual bool skip_timestep(std::shared_ptr<AtomData> atom_data) = 0;
    virtual bool read_timestep(std::shared_ptr<AtomData> atom_data) = 0;

private:
    static std::vector<std::shared_ptr<MolReader>> m_readers;
};

} // namespace mol::internal

#endif // MOLREADER_HPP
