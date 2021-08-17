#include "core/MolData.hpp"
#include <molpp/internal/BaseSel.hpp>

using namespace mol::internal;

bool mol::internal::check_valid_frame(size_t const frame, std::shared_ptr<MolData> data)
{
    return frame < data->num_frames();
}
