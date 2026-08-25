/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/global_index.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/finding/device/surface_flags.hpp"
#include "traccc/fitting/kalman_filter/is_line_visitor.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>
#include <vecmem/containers/device_vector.hpp>

// Detray include(s).
#include <detray/definitions/indexing.hpp>
#include <detray/geometry/tracking_surface.hpp>

// System include(s).
#include <type_traits>

namespace traccc::device {

/// Payload for the @c device::build_surface_flags function
struct build_surface_flags_payload {
  /// View object to the vector of flags to fill, one per surface
  vecmem::data::vector_view<std::underlying_type_t<surface_flags>>
      surface_flags_view;
};

/// Compute the properties of every surface, once per event
///
/// The track finding kernel needs to know whether a surface is a line, and
/// whether it has material. Both properties are constant for a whole event, so
/// we compute them once here instead of once per candidate.
///
/// @param[in] thread_id The index of the surface to process
/// @param[in] det_data  View object to the tracking detector description
/// @param[in] payload   The global memory payload
///
template <typename detector_t>
TRACCC_HOST_DEVICE inline void build_surface_flags(
    const global_index_t thread_id,
    typename detector_t::const_view_type det_data,
    const build_surface_flags_payload& payload) {
  using flags_t = std::underlying_type_t<surface_flags>;

  vecmem::device_vector<flags_t> out_flags(payload.surface_flags_view);

  if (thread_id >= out_flags.size()) {
    return;
  }

  detector_t det(det_data);
  const detray::tracking_surface sf{det,
                                    static_cast<detray::dindex>(thread_id)};

  flags_t flags = 0;

  if (detail::is_line(sf)) {
    flags |= static_cast<flags_t>(surface_flags::e_is_line);
  }

  if (sf.has_material()) {
    flags |= static_cast<flags_t>(surface_flags::e_has_material);
  }

  out_flags.at(thread_id) = flags;
}

}  // namespace traccc::device
