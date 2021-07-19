#include "internal/BaseSel.hpp"
#include "core/AtomData.hpp"

using namespace mol::internal;

bool mol::internal::check_valid_frame(size_t const frame, std::shared_ptr<AtomData> data)
{
    return frame < data->num_frames();
}
