#ifndef MOLREADER_HPP
#define MOLREADER_HPP

#include <string>
#include <memory>
#include <vector>

namespace mol::internal
{

class MolData;

class MolReader
{
public:
    enum Status
    {
        SUCCESS,
        INVALID,
        WRONG_ATOMS,
        END,
        FAILED,
    };

    static std::shared_ptr<MolReader> from_file_ext(std::string const &file_ext);
    virtual ~MolReader() {};
    std::shared_ptr<MolData> read_topology(std::string const &file_name);
    Status read_trajectory(std::string const &file_name,
                           std::shared_ptr<MolData> atom_data,
                           int begin=0, int end=-1, int step=1);
    virtual bool has_topology() const = 0;
    virtual bool has_trajectory() const = 0;
    virtual bool has_bonds() const = 0;
    virtual bool has_trajectory_metadata() const = 0;
    virtual Status open(const std::string &file_name) = 0;
    virtual void close() = 0;
    virtual std::shared_ptr<MolData> read_atoms() = 0;
    virtual Status check_timestep_read(std::shared_ptr<MolData> atom_data) = 0;
    virtual Status skip_timestep(std::shared_ptr<MolData> atom_data) = 0;
    virtual Status read_timestep(std::shared_ptr<MolData> atom_data) = 0;

private:
};

} // namespace mol::internal

#endif // MOLREADER_HPP
