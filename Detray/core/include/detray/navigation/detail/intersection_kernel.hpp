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
#include "detray/navigation/intersection/intersector_base.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/ranges.hpp"
#include "detray/utils/type_registry.hpp"

namespace detray::detail {

/// @returns the @param i -th element of @param t, or @param t itself if it
/// holds a single element
template <typename T>
DETRAY_HOST_DEVICE constexpr decltype(auto) at_solution(T &t,
                                                        const std::size_t i) {
  if constexpr (concepts::subscriptable<T>) {
    return (t[i]);
  } else {
    return (t);
  }
}

template <template <typename, typename, bool> typename intersector_t,
          bool contains_pos_v>
struct select_intersector {
  template <typename mask_t>
  using type = intersector_t<typename mask_t::shape,
                             typename mask_t::algebra_type, contains_pos_v>;
};

/// Number of mask types in @tparam mask_list_t that map to @tparam intersector_t
/// @{
template <typename intersector_t,
          template <typename, typename, bool> class intersector_constructor_t,
          bool contains_pos_v, typename mask_list_t>
struct n_masks_for_intersector {};

template <typename intersector_t,
          template <typename, typename, bool> class intersector_constructor_t,
          bool contains_pos_v, typename... mask_ts>
struct n_masks_for_intersector<intersector_t, intersector_constructor_t,
                               contains_pos_v, types::list<mask_ts...>>
    : std::integral_constant<
          std::size_t,
          (static_cast<std::size_t>(
               std::same_as<intersector_constructor_t<
                                typename mask_ts::shape,
                                typename mask_ts::algebra_type, contains_pos_v>,
                            intersector_t>) +
           ... + 0u)> {};
/// @}

/// @returns true if any of the solutions in @param result is valid
template <typename result_t>
DETRAY_HOST_DEVICE constexpr bool any_solution_valid(const result_t &result) {
  if constexpr (concepts::subscriptable<result_t>) {
    for (const auto &ip : result) {
      if (ip.is_valid()) {
        return true;
      }
    }
    return false;
  } else {
    return result.is_valid();
  }
}

/// Check the intersection point @param ip against the masks of a surface and
/// stop at the first mask that contains it
template <bool contains_pos_v, typename mask_check_t, typename mask_group_t,
          typename mask_range_t, typename traj_t, typename intersection_point_t,
          typename transform_t, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE inline void check_masks_in_range(
    mask_check_t &check, const mask_group_t &mask_group,
    const mask_range_t &mask_range, const traj_t &traj,
    const intersection_point_t &ip, const transform_t &ctf,
    const mask_tolerance<scalar_t> &tol) {
  for (const auto &mask : detray::ranges::subrange(mask_group, mask_range)) {
    check = check_intersection_mask<contains_pos_v>(traj, ip, mask, ctf, tol);

    // Stop at the first mask that contains the intersection point (the
    // edge check is only performed if an edge tolerance is given)
    if (detray::detail::any_of(check.inside) ||
        detray::detail::any_of(check.with_edge)) {
      break;
    }
  }
}

/// Resolve the masks for every solution in @param result: the mask
/// independent parts run here, the mask check is done by @param check_masks
///
/// @param is_container the intersection(s) to be filled, one per solution
/// @param result the solution(s) of the intersection with the unbounded surface
/// @param sf_desc the surface descriptor
/// @param cfg the intersection configuration
/// @param external_mask_tolerance additional mask tol. given by the caller
/// @param check_masks callable that checks a point against the surface masks
template <std::size_t n_sol, typename is_container_t, typename result_t,
          typename surface_t, concepts::scalar scalar_t, typename check_fn_t>
DETRAY_HOST_DEVICE inline void resolve_solutions(
    is_container_t &is_container, const result_t &result,
    const surface_t &sf_desc, const intersection::config &cfg,
    const scalar_t external_mask_tolerance, check_fn_t &&check_masks) {
  for (std::size_t i = 0u; i < n_sol; ++i) {
    const auto &ip = at_solution(result, i);
    auto &is = at_solution(is_container, i);

    // Mask independent part: status and path check (an invalid solution
    // leaves the intersection marked as outside)
    if (!init_intersection(is, ip, cfg) || !ip.is_valid()) [[unlikely]] {
      continue;
    }

    const mask_tolerance<scalar_t> tol =
        mask_tolerances(sf_desc, ip, cfg, external_mask_tolerance);

    // Mask dependent part
    const auto check = check_masks(ip, tol);

    // Mask independent part: fill the intersection
    finalize_intersection(is, ip, check);
  }
}

/// A functor to check one intersection point against the masks of a surface.
/// Only the mask types that belong to @tparam intersector_t are checked
template <template <typename, typename, bool> class intersector_constructor_t,
          bool contains_pos_v, typename intersector_t>
struct intersect_surface_per_mask {
  template <typename mask_group_t, typename mask_range_t, typename mask_check_t,
            typename traj_t, typename intersection_point_t,
            typename transform_t, concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline void operator()(
      const mask_group_t &mask_group, const mask_range_t &mask_range,
      mask_check_t &check, const traj_t &traj, const intersection_point_t &ip,
      const transform_t &ctf, const mask_tolerance<scalar_t> &tol) const {
    using mask_t = typename mask_group_t::value_type;
    using shape_t = typename mask_t::shape;
    using algebra_t = typename mask_t::algebra_type;

    using local_intersector_t =
        intersector_constructor_t<shape_t, algebra_t, contains_pos_v>;

    if constexpr (std::same_as<local_intersector_t, intersector_t>) {
      check_masks_in_range<contains_pos_v>(check, mask_group, mask_range, traj,
                                           ip, ctf, tol);
    }
  }
};

/// A functor to intersect a surface whose intersector needs the mask (the
/// cylinder intersectors read the radius). The whole intersection runs in the
/// mask branch, which loses nothing as long as exactly one mask type maps to
/// @tparam intersector_t
template <template <typename, typename, bool> class intersector_constructor_t,
          bool contains_pos_v, typename intersector_t>
struct intersect_surface_cylinder_per_mask {
  template <typename mask_group_t, typename mask_range_t,
            typename is_container_t, typename traj_t, typename surface_t,
            typename transform_t, concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline void operator()(
      const mask_group_t &mask_group, const mask_range_t &mask_range,
      is_container_t &is_container, const traj_t &traj,
      const surface_t &sf_desc, const transform_t &ctf,
      const intersection::config &cfg,
      const scalar_t external_mask_tolerance) const {
    using mask_t = typename mask_group_t::value_type;
    using shape_t = typename mask_t::shape;
    using algebra_t = typename mask_t::algebra_type;

    using local_intersector_t =
        intersector_constructor_t<shape_t, algebra_t, contains_pos_v>;

    if constexpr (std::same_as<local_intersector_t, intersector_t>) {
      using nav_link_t = typename mask_t::links_type;
      using mask_check_t =
          mask_check_result<algebra_t, nav_link_t, contains_pos_v>;

      // The cylinder intersectors need the radius of the mask
      dindex mask_idx{detail::invalid_value<dindex>()};
      if constexpr (concepts::interval<mask_range_t>) {
        mask_idx = mask_range.lower();
      } else {
        mask_idx = mask_range;
      }
      assert(mask_idx < mask_group.size());

      const auto radius = static_cast<dscalar<algebra_t>>(
          mask_group[mask_idx][mask_t::shape::e_r]);

      constexpr intersector_t intersector{};
      const auto result = intersector.point_of_intersection(
          traj, ctf, radius, cfg.overstep_tolerance);

      if (!any_solution_valid(result)) [[unlikely]] {
        return;
      }

      resolve_solutions<intersector_t::n_solutions>(
          is_container, result, sf_desc, cfg, external_mask_tolerance,
          [&](const auto &ip, const mask_tolerance<scalar_t> &tol) {
            mask_check_t check{};
            check_masks_in_range<contains_pos_v>(check, mask_group, mask_range,
                                                 traj, ip, ctf, tol);
            return check;
          });
    }
  }
};

template <template <typename, typename, bool> class intersector_constructor_t,
          bool contains_pos_v>
struct intersect_surface_per_intersector {
  template <typename intersector_t, typename mask_store_t,
            typename mask_range_t, typename surface_t, typename is_container_t,
            typename traj_t, typename transform_t, concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE inline void operator()(
      [[maybe_unused]] const intersector_t &intersector,
      const mask_store_t &mask_store, const typename mask_store_t::ids mask_id,
      const mask_range_t mask_range, const surface_t &sf_desc,
      is_container_t &is_container, const traj_t &traj, const transform_t &ctf,
      const intersection::config &cfg,
      const scalar_t external_mask_tolerance = 0.f) const {
    using mask_list_t = typename mask_store_t::value_types::type_list;

    // Keep in mind that this function body is called once for every
    // intersector type, not for every mask. What we will do now is call a
    // different per-mask function object for every mask type.
    //
    // We need to be careful here, because we are already visiting every
    // intersector type, so we must ensure that the function object we call
    // here knows what the intersector is.
    if constexpr (concepts::cylindrical_frame<
                      typename intersector_t::frame_type>) {
      // The intersector needs the mask, so the intersection is done in the
      // mask branch. This is only free of duplication for a single mask type
      static_assert(
          n_masks_for_intersector<intersector_t, intersector_constructor_t,
                                  contains_pos_v, mask_list_t>::value == 1u,
          "Cylindrical intersectors must map to exactly one mask type");

      mask_store.template visit<intersect_surface_cylinder_per_mask<
          intersector_constructor_t, contains_pos_v, intersector_t>>(
          mask_id, mask_range, is_container, traj, sf_desc, ctf, cfg,
          external_mask_tolerance);
    } else {
      using algebra_t = typename intersector_t::algebra_type;
      using nav_link_t = typename types::front<mask_list_t>::links_type;
      using mask_check_t =
          mask_check_result<algebra_t, nav_link_t, contains_pos_v>;

      const auto result =
          intersector.point_of_intersection(traj, ctf, cfg.overstep_tolerance);

      if (!any_solution_valid(result)) [[unlikely]] {
        return;
      }

      resolve_solutions<intersector_t::n_solutions>(
          is_container, result, sf_desc, cfg, external_mask_tolerance,
          [&](const auto &ip, const mask_tolerance<scalar_t> &tol) {
            mask_check_t check{};
            mask_store.template visit<intersect_surface_per_mask<
                intersector_constructor_t, contains_pos_v, intersector_t>>(
                mask_id, mask_range, check, traj, ip, ctf, tol);
            return check;
          });
    }
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

/// Type registry of the intersectors that are needed for the mask types in
/// @tparam mask_store_t
template <template <typename, typename, bool> class intersector_constructor_t,
          typename mask_store_t, typename intersection_t>
using intersect_surface_registry_t =
    types::mapped_registry<typename mask_store_t::value_types,
                           select_intersector<intersector_constructor_t,
                                              intersection_t::contains_pos()>>;

/// Maximum number of intersections that a surface in @tparam mask_store_t can
/// produce
template <template <typename, typename, bool> class intersector_constructor_t,
          typename mask_store_t, typename intersection_t>
inline constexpr std::size_t intersect_surface_max_n_results =
    max_intersections_for_intersectors<typename intersect_surface_registry_t<
        intersector_constructor_t, mask_store_t,
        intersection_t>::type_list>::value;

/// One intersection per solution: a single intersection or an array of them
template <template <typename, typename, bool> class intersector_constructor_t,
          typename mask_store_t, typename intersection_t>
using intersect_surface_output_t = std::conditional_t<
    (intersect_surface_max_n_results<intersector_constructor_t, mask_store_t,
                                     intersection_t> == 1),
    intersection_t,
    intersection_t[intersect_surface_max_n_results<
        intersector_constructor_t, mask_store_t, intersection_t>]>;

/// Intersect a surface with a trajectory and write one intersection per
/// solution into @param found_intersections
///
/// @tparam intersector_constructor_t is the intersector template to be used
///
/// @param mask_store is the mask store that holds the surface masks
/// @param found_intersections the intersection(s) to be filled
/// @param traj is the input trajectory
/// @param sf_desc is the input surface
/// @param ctf is the transform of the surface
/// @param cfg is the intersection configuration
/// @param external_mask_tolerance additional mask tol. given by the caller
template <template <typename, typename, bool> class intersector_constructor_t,
          typename mask_store_t, typename output_t, typename traj_t,
          typename surface_t, typename transform_t, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE inline void intersect_surface(
    const mask_store_t &mask_store, output_t &found_intersections,
    const traj_t &traj, const surface_t sf_desc, const transform_t &ctf,
    const intersection::config &cfg, const scalar_t external_mask_tolerance) {
  using intersection_t = std::remove_extent_t<output_t>;
  using registry_t = intersect_surface_registry_t<intersector_constructor_t,
                                                  mask_store_t, intersection_t>;

  // Decode the mask link once and pass the id and the range down, so that
  // the link is not re-read from the surface descriptor in every visit
  const auto mask_link = sf_desc.mask();
  const auto mask_id = mask_link.id();
  const auto mask_range = mask_link.index();

  // We could, naively, visit the mask store directly, but the intersection
  // initializer does a lot of non-mask-dependent work. Compiling it once per
  // mask generates a lot of unwanted code. Instead, we visit the intersectors
  // that match each of the masks so we can generate the intersection code
  // once for every intersector, rather than once for every mask.
  types::visit<registry_t,
               intersect_surface_per_intersector<
                   intersector_constructor_t, intersection_t::contains_pos()>>(
      mask_id, mask_store, mask_id, mask_range, sf_desc, found_intersections,
      traj, ctf, cfg, external_mask_tolerance);
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
  using intersection_t = typename is_container_t::value_type;
  using output_t = intersect_surface_output_t<intersector_constructor_t,
                                              mask_store_t, intersection_t>;
  constexpr std::size_t max_n_results =
      intersect_surface_max_n_results<intersector_constructor_t, mask_store_t,
                                      intersection_t>;

  const auto &ctf = contextual_transforms.at(sf_desc.transform(), ctx);

  output_t found_intersections{};

  intersect_surface<intersector_constructor_t>(mask_store, found_intersections,
                                               traj, sf_desc, ctf, cfg,
                                               external_mask_tolerance);

  // The surface link is set once here, so that it is not carried through
  // the intersector dispatch
  if constexpr (concepts::subscriptable<output_t>) {
    for (std::size_t i = 0u; i < max_n_results; ++i) {
      if (found_intersections[i].is_probably_inside()) {
        found_intersections[i].set_surface(sf_desc);
        insert_sorted(found_intersections[i], is_container);
      }
    }
  } else {
    if (found_intersections.is_probably_inside()) {
      found_intersections.set_surface(sf_desc);
      insert_sorted(found_intersections, is_container);
    }
  }
}

/// Update the intersection @param sfi of a surface with a trajectory. Only the
/// closest solution is kept.
///
/// @tparam intersector_t is the intersector template to be used
/// @tparam mask_store_t is the mask store type
/// @tparam traj_t is the input trajectory type (e.g. ray or helix)
/// @tparam intersection_t is the intersection type
/// @tparam transform_container_t is the input transform store type
///
/// @param mask_store is the mask store that holds the surface masks
/// @param traj is the input trajectory
/// @param sfi is the intersection to be updated (holds the surface)
/// @param contextual_transforms is the input transform container
/// @param ctx is the geometry context
/// @param cfg is the intersection configuration
/// @param external_mask_tolerance additional mask tol. given by the caller
///
/// @returns true if the trajectory (probably) hits the surface
template <template <typename, typename, bool> class intersector_constructor_t,
          typename mask_store_t, typename traj_t, typename intersection_t,
          typename transform_container_t, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE inline bool intersection_update_surface(
    const mask_store_t &mask_store, const traj_t &traj, intersection_t &sfi,
    const transform_container_t &contextual_transforms,
    const typename transform_container_t::context_type &ctx,
    const intersection::config &cfg,
    const scalar_t external_mask_tolerance = 0.f) {
  using output_t = intersect_surface_output_t<intersector_constructor_t,
                                              mask_store_t, intersection_t>;

  const auto &sf_desc = sfi.surface();
  const auto &ctf = contextual_transforms.at(sf_desc.transform(), ctx);

  // Start from the current intersection, so that it is left untouched if the
  // surface is not intersected at all
  output_t found_intersections{};
  at_solution(found_intersections, 0u) = sfi;

  intersect_surface<intersector_constructor_t>(mask_store, found_intersections,
                                               traj, sf_desc, ctf, cfg,
                                               external_mask_tolerance);

  // Only the closest solution is kept
  sfi = at_solution(found_intersections, 0u);

  return sfi.is_probably_inside();
}

}  // namespace detray::detail
