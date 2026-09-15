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
#include "detray/definitions/algorithms.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/intersection_config.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/ranges.hpp"
#include "detray/utils/type_registry.hpp"

namespace detray::detail {

struct intersection_initialize_get_radius {
  template <typename mask_group_t, typename mask_range_t>
  DETRAY_HOST_DEVICE inline auto operator()(
      const mask_group_t &mask_group, const mask_range_t &mask_range) const {
    using mask_t = typename mask_group_t::value_type;
    using algebra_t = typename mask_t::algebra_type;
    using scalar_t = dscalar<algebra_t>;
    if constexpr (concepts::cylindrical<mask_t>) {
      dindex mask_idx{detail::invalid_value<dindex>()};

      if constexpr (concepts::interval<mask_range_t>) {
        mask_idx = mask_range.lower();
      } else {
        mask_idx = mask_range;
      }

      assert(mask_idx < mask_group.size());

      return static_cast<scalar_t>(mask_group[mask_idx][mask_t::shape::e_r]);
    } else {
      return 0.f;
    }
  };
};

/// A functor to add all valid intersections between the trajectory and
/// surface
template <template <typename, typename, bool> class intersector_t,
          bool contains_pos_v>
struct intersection_initialize {
  /// Operator function to initialize intersections
  ///
  /// @tparam mask_group_t is the input mask group type found by variadic
  /// unrolling
  /// @tparam is_container_t is the intersection container type
  /// @tparam traj_t is the input trajectory type (e.g. ray or helix)
  /// @tparam surface_t is the input surface type
  /// @tparam transform_container_t is the input transform store type
  ///
  /// @param mask_group is the input mask group
  /// @param is_container is the intersection container to be filled
  /// @param traj is the input trajectory
  /// @param surface is the input surface
  /// @param contextual_transforms is the input transform container
  /// @param mask_tolerance is the tolerance for mask size
  /// @param overstep_tol negative cutoff for the path
  ///
  /// @return the number of valid intersections
  template <typename mask_group_t, typename mask_range_t,
            typename is_container_t, typename traj_t, typename surface_t,
            typename intersection_result_t, typename transform_t,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline auto operator()(
      const mask_group_t &mask_group, const mask_range_t &mask_range,
      is_container_t &is_container, const traj_t &traj,
      const surface_t &sf_desc, const intersection_result_t &intersections,
      const transform_t &ctf, const intersection::config &cfg,
      const scalar_t external_mask_tolerance = 0.f) const {
    using mask_t = typename mask_group_t::value_type;
    using shape_t = typename mask_t::shape;
    using algebra_t = typename mask_t::algebra_type;

    // Find the point of intersection with the underlying geometry
    constexpr intersector_t<shape_t, algebra_t, contains_pos_v> intersector{};

    constexpr std::uint8_t n_sol{decltype(intersector)::n_solutions};

    for (std::size_t i = 0u; i < n_sol; ++i) {
      if constexpr (concepts::subscriptable<intersection_result_t>) {
        if (!intersections[i].is_valid()) [[unlikely]] {
          continue;
        }
      } else {
        if (!intersections.is_valid()) [[unlikely]] {
          continue;
        }
      }

      // Resolve the masks that belong to the surface
      for (const auto &mask :
           detray::ranges::subrange(mask_group, mask_range)) {
        // Build the resulting intersection(s) from the intersection point
        if constexpr (concepts::subscriptable<is_container_t>) {
          if constexpr (concepts::subscriptable<intersection_result_t>) {
            resolve_mask(is_container[i], traj, intersections[i], sf_desc, mask,
                         ctf, cfg, external_mask_tolerance);
          } else {
            resolve_mask(is_container[i], traj, intersections, sf_desc, mask,
                         ctf, cfg, external_mask_tolerance);
          }
          if (is_container[i].is_probably_inside()) {
            break;
          }
        } else {
          if constexpr (concepts::subscriptable<intersection_result_t>) {
            resolve_mask(is_container, traj, intersections[i], sf_desc, mask,
                         ctf, cfg, external_mask_tolerance);
          } else {
            resolve_mask(is_container, traj, intersections, sf_desc, mask, ctf,
                         cfg, external_mask_tolerance);
          }
          if (is_container.is_probably_inside()) {
            break;
          }
        }
      }
    }
  }
};

template <template <typename, typename, bool> typename intersector_t,
          bool contains_pos_v>
struct select_intersector {
  template <typename mask_t>
  using type = intersector_t<typename mask_t::shape,
                             typename mask_t::algebra_type, contains_pos_v>;
};

template <template <typename, typename, bool> class intersector_constructor_t,
          bool contains_pos_v, typename intersector_t>
struct intersection_initialize_surface_per_mask {
  template <typename mask_group_t, typename mask_range_t,
            typename is_container_t, typename traj_t, typename surface_t,
            typename intersection_result_t, typename transform_t,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline auto operator()(
      const mask_group_t &mask_group, const mask_range_t &mask_range,
      is_container_t &is_container, const traj_t &traj,
      const surface_t &sf_desc, const intersection_result_t &intersections,
      const transform_t &ctf, const intersection::config &cfg,
      const scalar_t external_mask_tolerance = 0.f) {
    using mask_t = typename mask_group_t::value_type;
    using shape_t = typename mask_t::shape;
    using algebra_t = typename mask_t::algebra_type;

    using local_intersector_t =
        intersector_constructor_t<shape_t, algebra_t, contains_pos_v>;

    if constexpr (std::same_as<local_intersector_t, intersector_t>) {
      intersection_initialize<intersector_constructor_t, contains_pos_v>{}(
          mask_group, mask_range, is_container, traj, sf_desc, intersections,
          ctf, cfg, external_mask_tolerance);
    }
  }
};

template <template <typename, typename, bool> class intersector_constructor_t,
          bool contains_pos_v>
struct intersection_initialize_surface_per_intersector {
  template <typename intersector_t, typename mask_store_t, typename surface_t,
            typename is_container_t, typename traj_t, typename transform_t,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline auto operator()(
      const intersector_t &intersector, const mask_store_t &mask_store,
      const surface_t &sf_desc, is_container_t &is_container,
      const traj_t &traj, const surface_t &_sf_desc, const transform_t &ctf,
      const intersection::config &cfg,
      const scalar_t external_mask_tolerance = 0.f) {
    typename intersector_t::result_type result{};

    if constexpr (concepts::cylindrical_frame<
                      typename intersector_t::frame_type>) {
      const auto radius =
          mask_store.template visit<intersection_initialize_get_radius>(
              sf_desc.mask());
      result = intersector.point_of_intersection(traj, ctf, radius,
                                                 cfg.overstep_tolerance);
    } else {
      result =
          intersector.point_of_intersection(traj, ctf, cfg.overstep_tolerance);
    }

    constexpr std::uint8_t n_sol{intersector_t::n_solutions};

    if constexpr (n_sol > 1) {
      bool any_valid = false;
      for (std::size_t i = 0u; i < n_sol; ++i) {
        if (result[i].is_valid()) {
          any_valid |= true;
        }
      }
      if (!any_valid) [[unlikely]] {
        return;
      }
    } else {
      if (!result.is_valid()) [[unlikely]] {
        return;
      }
    }

    // Keep in mind that this function body is called once for every
    // intersector type, not for every mask. What we will do now is call a
    // different per-mask function object for every mask type.
    //
    // We need to be careful here, because we are already visiting every
    // intersector type, so we must ensure that the function object we call
    // here knows what the intersector is.
    mask_store.template visit<intersection_initialize_surface_per_mask<
        intersector_constructor_t, contains_pos_v, intersector_t>>(
        sf_desc.mask(), is_container, traj, _sf_desc, result, ctf, cfg,
        external_mask_tolerance);
  }
};

template <typename T>
struct max_intersections_for_intersectors {};

template <typename... Ts>
struct max_intersections_for_intersectors<types::list<Ts...>> {
  static constexpr auto value = std::max({Ts::n_solutions...});
};

template <typename intersection_t, typename... allocator_t>
DETRAY_HOST_DEVICE void insert_sorted(
    const intersection_t &sfi,
    std::vector<intersection_t, allocator_t...> &intersections) {
  auto itr_pos =
      detray::upper_bound(intersections.cbegin(), intersections.cend(), sfi);

  intersections.insert(itr_pos, sfi);
}

/// Specialization for the navigation state cache
template <typename nav_state_t>
DETRAY_HOST_DEVICE void insert_sorted(
    const typename nav_state_t::value_type &sfi, nav_state_t &intersections) {
  auto itr_pos{intersections.cbegin()};

  // For just two candidates int the cache, the navigation state keeps
  // the first as the previously visited candidate -> no sorting needed
  if constexpr (nav_state_t::capacity() > 2u) {
    itr_pos =
        detray::upper_bound(intersections.cbegin(), intersections.cend(), sfi);
  }

  intersections.insert(itr_pos, sfi);
}

/// Intersect a surface with a trajectory and add all valid intersections to
/// the intersection container
///
/// @tparam intersector_t is the intersector template to be used
/// @tparam mask_store_t is the mask store type
/// @tparam is_container_t is the intersection container type
/// @tparam traj_t is the input trajectory type (e.g. ray or helix)
/// @tparam surface_t is the input surface type
/// @tparam transform_container_t is the input transform store type
///
/// @param mask_store is the mask store that holds the surface masks
/// @param is_container is the intersection container to be filled
/// @param traj is the input trajectory
/// @param sf_desc is the input surface
/// @param contextual_transforms is the input transform container
/// @param ctx is the geometry context
/// @param cfg is the intersection configuration
/// @param external_mask_tolerance additional mask tol. given by the caller
template <template <typename, typename, bool> class intersector_constructor_t,
          typename mask_store_t, typename is_container_t, typename traj_t,
          typename surface_t, typename transform_container_t,
          concepts::scalar scalar_t>
DETRAY_HOST_DEVICE inline void intersection_initialize_surface(
    const mask_store_t &mask_store, is_container_t &is_container,
    const traj_t &traj, const surface_t &sf_desc,
    const transform_container_t &contextual_transforms,
    const typename transform_container_t::context_type &ctx,
    const intersection::config &cfg,
    const scalar_t external_mask_tolerance = 0.f) {
  using masks_t = typename mask_store_t::value_types;
  static constexpr bool contains_pos =
      is_container_t::value_type::contains_pos();
  using registry_t = types::mapped_registry<
      masks_t, select_intersector<intersector_constructor_t, contains_pos>>;

  const auto &ctf = contextual_transforms.at(sf_desc.transform(), ctx);

  static constexpr auto max_n_results =
      max_intersections_for_intersectors<typename registry_t::type_list>::value;

  using single_output_t = typename is_container_t::value_type;
  using output_t = std::conditional_t<(max_n_results == 1), single_output_t,
                                      single_output_t[max_n_results]>;

  output_t found_intersections{};

  // We could, naively, visit the mask store directly, but the intersection
  // initializer does a lot of non-mask-dependent work. Compiling it once per
  // mask generates a lot of unwanted code. Instead, we visit the intersectors
  // that match each of the masks so we can generate the intersection code
  // once for every intersector, rather than once for every mask.
  types::visit<registry_t, intersection_initialize_surface_per_intersector<
                               intersector_constructor_t, contains_pos>>(
      sf_desc.mask().id(), mask_store, sf_desc, found_intersections, traj,
      sf_desc, ctf, cfg, external_mask_tolerance);

  if constexpr (concepts::subscriptable<output_t>) {
    for (std::size_t i = 0u; i < max_n_results; ++i) {
      if (found_intersections[i].is_probably_inside()) {
        insert_sorted(found_intersections[i], is_container);
      }
    }
  } else {
    if (found_intersections.is_probably_inside()) {
      insert_sorted(found_intersections, is_container);
    }
  }
}

/// A functor to update the closest intersection between the trajectory and
/// surface
template <template <typename, typename, bool> class intersector_t>
struct intersection_update {
  /// Operator function to update the intersection
  ///
  /// @tparam mask_group_t is the input mask group type found by variadic
  /// unrolling
  /// @tparam traj_t is the input trajectory type (e.g. ray or helix)
  /// @tparam surface_t is the input surface type
  /// @tparam transform_container_t is the input transform store type
  ///
  /// @param mask_group is the input mask group
  /// @param mask_range is the range of masks in the group that belong to the
  ///                   surface
  /// @param traj is the input trajectory
  /// @param surface is the input surface
  /// @param contextual_transforms is the input transform container
  /// @param mask_tolerance is the tolerance for mask size
  /// @param overstep_tol negative cutoff for the path
  ///
  /// @return the intersection
  template <typename mask_group_t, typename mask_range_t, typename traj_t,
            typename intersection_t, typename transform_container_t,
            concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline bool operator()(
      const mask_group_t &mask_group, const mask_range_t &mask_range,
      const traj_t &traj, intersection_t &sfi,
      const transform_container_t &contextual_transforms,
      const typename transform_container_t::context_type &ctx,
      const intersection::config &cfg,
      const scalar_t external_mask_tolerance = 0.f) const {
    using mask_t = typename mask_group_t::value_type;
    using shape_t = typename mask_t::shape;
    using algebra_t = typename mask_t::algebra_type;

    // Find the point of intersection with the underlying geometry
    const auto &ctf = contextual_transforms.at(sfi.surface().transform(), ctx);

    constexpr intersector_t<shape_t, algebra_t, intersection_t::contains_pos()>
        intersector{};
    constexpr std::uint8_t n_sol{decltype(intersector)::n_solutions};

    typename decltype(intersector)::result_type result{};

    if constexpr (concepts::cylindrical<mask_t>) {
      dindex mask_idx{detail::invalid_value<dindex>()};
      if constexpr (concepts::interval<mask_range_t>) {
        mask_idx = mask_range.lower();
      } else {
        mask_idx = mask_range;
      }
      assert(mask_idx < mask_group.size());

      result = intersector.point_of_intersection(
          traj, ctf, mask_group[mask_idx], cfg.overstep_tolerance);
    } else {
      result =
          intersector.point_of_intersection(traj, ctf, cfg.overstep_tolerance);
    }

    // Check if any valid solutions were found
    if constexpr (n_sol > 1) {
      bool found_any{false};
      for (const auto &ip : result) {
        if (ip.is_valid()) {
          found_any = true;
        }
      }
      if (!found_any) [[unlikely]] {
        return false;
      }
    } else {
      if (!result.is_valid()) [[unlikely]] {
        return false;
      }
    }

    // Run over the masks that belong to the surface
    for (const auto &mask : detray::ranges::subrange(mask_group, mask_range)) {
      // Build the resulting intersecion(s) from the intersection point
      if constexpr (n_sol > 1) {
        resolve_mask(sfi, traj, result[0], sfi.surface(), mask, ctf, cfg,
                     external_mask_tolerance);
      } else {
        resolve_mask(sfi, traj, result, sfi.surface(), mask, ctf, cfg,
                     external_mask_tolerance);
      }

      if (sfi.is_probably_inside()) {
        return true;
      }
    }

    return false;
  }
};

}  // namespace detray::detail
