/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "../../../utils/global_index.hpp"
#include "../build_surface_flags.cuh"

namespace traccc::cuda {
namespace kernels {

template <typename detector_t>
__global__ void build_surface_flags(
    const __grid_constant__ typename detector_t::const_view_type det,
    const __grid_constant__ device::build_surface_flags_payload payload) {
  device::build_surface_flags<detector_t>(details::global_index1(), det,
                                          payload);
}

}  // namespace kernels

template <typename detector_t>
void build_surface_flags(const dim3& grid_size, const dim3& block_size,
                         const cudaStream_t& stream,
                         typename detector_t::const_view_type det,
                         const device::build_surface_flags_payload& payload) {
  kernels::build_surface_flags<detector_t>
      <<<grid_size, block_size, 0, stream>>>(det, payload);
}
}  // namespace traccc::cuda
