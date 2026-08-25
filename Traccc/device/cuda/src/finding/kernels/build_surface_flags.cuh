/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/finding/device/build_surface_flags.hpp"

// CUDA include(s).
#include <cuda_runtime.h>

namespace traccc::cuda {

/// Fill the per-surface flags, once per event
///
/// @param grid_size  The number of blocks to launch
/// @param block_size The number of threads per block
/// @param stream     The stream to launch the kernel on
/// @param det        View of the detector to read the surfaces from
/// @param payload    The payload for the kernel
///
template <typename detector_t>
void build_surface_flags(const dim3& grid_size, const dim3& block_size,
                         const cudaStream_t& stream,
                         typename detector_t::const_view_type det,
                         const device::build_surface_flags_payload& payload);

}  // namespace traccc::cuda
