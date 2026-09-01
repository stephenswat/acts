// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/detail/material_accessor.hpp"
#include "detray/material/material.hpp"
#include "detray/navigation/accelerators/search_window.hpp"
#include "detray/navigation/concepts.hpp"

namespace detray::detail {

/// A functor to retrieve the volume material parameters at a given position
struct get_material_params {
  template <typename mat_group_t, typename index_t, concepts::point point_t>
  DETRAY_HOST_DEVICE inline auto operator()(const mat_group_t &mat_group,
                                            const index_t &idx,
                                            const point_t &loc_p) const {
    using material_t = typename mat_group_t::value_type;

    if constexpr (concepts::volume_material<material_t>) {
      if constexpr (concepts::homogeneous_material<material_t>) {
        // Homogeneous volume material
        return &(detail::material_accessor::get(mat_group, idx, loc_p));
      } else {
        // Volume material maps
        return &(detail::material_accessor::get(mat_group, idx, loc_p)
                     .get_material());
      }
    } else {
      using scalar_t = typename material_t::scalar_type;
      // Cannot be reached for volumes
      return static_cast<const material<scalar_t> *>(nullptr);
    }
  }
};

/// A functor to access all surfaces registered in an acceleration structure
/// @tparam functor_t functor that performs a task per surface
template <typename functor_t>
struct apply_to_surfaces {
  /// Call operator that forwards the functor to all contained surfaces
  template <concepts::accelerator_collection accel_coll_t,
            typename accel_index_t, typename... Args>
  DETRAY_HOST_DEVICE inline void operator()(const accel_coll_t &coll,
                                            const accel_index_t index,
                                            Args &&...args) const {
    using accel_type = typename accel_coll_t::value_type;

    if constexpr (concepts::surface_accelerator<accel_type>) {
      // Run over the surfaces in a single acceleration data structure
      for (const auto &sf : coll[index].all()) {
        functor_t{}(sf, std::forward<Args>(args)...);
      }
    }
  }
};

/// A functor to find surfaces in the neighborhood of a track position
/// @tparam functor_t functor that performs a task per surface in the
///                   neighbourhood (e.g. intersection)
template <typename functor_t>
struct apply_to_neighbourhood {
  /// Apply the functor to one contiguous run of surfaces
  ///
  /// @param first the first surface descriptor of the run
  /// @param n the number of surface descriptors in the run
  /// @param vol_idx index of the volume that the surfaces belong to
  ///
  /// @note The signature does not depend on the acceleration structure that
  /// produced the run, so the compiler emits this code once for all of them,
  /// instead of inlining it into every branch of the accelerator dispatch.
  template <typename surface_t, typename... Args>
  DETRAY_HOST_DEVICE
// Not inlining this function keeps a single copy of the (very large) functor
// body and lets NVIDIA GPUs reconverge threads that use different acceleration
// data structures. The call is made once per run, not once per surface.
#if defined(__CUDACC__)
      DETRAY_NO_INLINE
#else
      inline
#endif
      static void process_run(const surface_t *first, const dindex n,
                              [[maybe_unused]] const dindex vol_idx,
                              Args &&...args) {
    for (dindex i = 0u; i < n; ++i) {
      const surface_t &sf = first[i];
      assert(sf.volume() == vol_idx);
      functor_t{}(sf, args...);
    }
  }

  /// Call operator that forwards the neighborhood search call in a volume
  /// to a surface finder data structure
  template <concepts::accelerator_collection accel_coll_t,
            typename accel_index_t, typename detector_t, typename track_t,
            typename window_size_t, typename... Args>
  DETRAY_HOST_DEVICE inline void operator()(
      const accel_coll_t &coll, const accel_index_t index,
      const detector_t &det, const typename detector_t::volume_type &volume,
      const track_t &track, const search_window<window_size_t, 2> &win_size,
      const typename detector_t::geometry_context &ctx, Args &&...args) const {
    using accel_type = typename accel_coll_t::value_type;

    decltype(auto) accel = coll[index];

    if constexpr (concepts::surface_accelerator<accel_type>) {
      // Run over the surfaces in a single acceleration data structure, one
      // contiguous run at a time
      for (const auto &entry : accel.runs(det, volume, track, win_size, ctx)) {
        const auto run = entry.run();
        process_run(run.first, run.size, volume.index(), args...);
      }
    } else if constexpr (concepts::volume_accelerator<accel_type>) {
      // Run over the daughter volumes in a single acceleration data
      // structure
      for (const dindex daughter_idx :
           accel.search(det, volume, track, win_size, ctx)) {
        functor_t{}(daughter_idx, std::forward<Args>(args)...);
      }
    }
  }
};

/// A functor to access the daughter volumes of a volume
template <typename functor_t>
struct apply_to_volumes {
  /// Call operator that forwards the functor call to all contained daughter
  /// volumes
  template <concepts::accelerator_collection accel_coll_t,
            typename accel_index_t, typename... Args>
  DETRAY_HOST_DEVICE inline void operator()(const accel_coll_t & /*accel_col*/,
                                            const accel_index_t /*accel_idx*/,
                                            Args &&.../*args*/) const {
    using accel_type = typename accel_coll_t::value_type;

    if constexpr (concepts::volume_accelerator<accel_type>) {
      // Run over all the daughter volumes
      // TODO: Implement e.g. BVH
    }
  }
};

}  // namespace detray::detail
