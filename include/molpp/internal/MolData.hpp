#ifndef MOLDATA_HPP
#define MOLDATA_HPP

#include <molpp/MolppCore.hpp>
#include <molpp/internal/AtomData.hpp>
#include <molpp/internal/BondData.hpp>
#include <molpp/internal/ResidueData.hpp>
#include <molpp/internal/PropertyContainer.hpp>

namespace mol::internal {

class MolData
{
public:
    MolData(size_t const num_atoms);

    AtomData &atoms() { return m_properties_old; }
    AtomData const &atoms() const { return m_properties_old; }
    BondData &bonds() { return m_bonds; }
    BondData const& bonds() const { return m_bonds; }
    ResidueData &residues() { return m_residues; }
    ResidueData const &residues() const { return m_residues; }

    template<IsAtomAggregate AggregateType>
    bool add_entity(size_t const size);

    template<IsAtomAggregate AggregateType>
    size_t entity_size() const;

    template<IsAtomAggregate AggregateType>
    void set_entity_size(size_t const size);

    size_t num_frames() const;
    Frame add_frame();
    void remove_frame(Frame const frame);

    template<IsAtomAggregate AggregateType, IsProperty PropertyType>
    PropertyType const* property_at(Frame const frame) const;

    template<IsAtomAggregate AggregateType, IsProperty PropertyType>
    PropertyType* property_at(Frame const frame);

    template<IsAtomAggregate AggregateType, IsProperty PropertyType>
    PropertyType* add_property(bool const is_time_based);

private:
    std::unordered_map<std::type_index, size_t> m_sizes;
    PropertyContainer m_properties;
    AtomData m_properties_old;
    BondData m_bonds;
    ResidueData m_residues;
};

template<IsAtomAggregate AggregateType>
inline bool MolData::add_entity(size_t const size)
{
    std::type_index const key = typeid(AggregateType);
    auto result = m_sizes.insert({key, size});
    return result.second;
}

template<IsAtomAggregate AggregateType>
inline size_t MolData::entity_size() const
{
    std::type_index const key = typeid(AggregateType);
    auto const iter = m_sizes.find(key);

    if (iter == m_sizes.end())
    {
        throw MolError("Entity not registered");
    }

    return iter->second;
}

template<IsAtomAggregate AggregateType>
inline void MolData::set_entity_size(size_t const size)
{
    std::type_index const key = typeid(AggregateType);
    auto const iter = m_sizes.find(key);

    if (iter == m_sizes.end())
    {
        throw MolError("Entity not registered");
    }

    iter->second = size;
    m_properties.resize<AggregateType>(size);
}

template<IsAtomAggregate AggregateType, IsProperty PropertyType>
inline PropertyType const* MolData::property_at(Frame const frame) const
{
    return m_properties.get<AggregateType, PropertyType>(frame);
}

template<IsAtomAggregate AggregateType, IsProperty PropertyType>
inline PropertyType* MolData::property_at(Frame const frame)
{
    return m_properties.get<AggregateType, PropertyType>(frame);
}

template<IsAtomAggregate AggregateType, IsProperty PropertyType>
inline PropertyType* MolData::add_property(bool const is_time_based)
{
    // entity_size() will check if Aggregate is already registered.
    size_t const size = entity_size<AggregateType>();
    return m_properties.add<AggregateType, PropertyType>(is_time_based, size);
}

} // namespace mol::internal

#endif // MOLDATA_HPP
