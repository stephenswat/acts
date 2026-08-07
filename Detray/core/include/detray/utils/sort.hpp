// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/find_bound.hpp"
#include "detray/utils/ranges/detail/iterator_functions.hpp"

// System include(s).
#include <algorithm>
#include <functional>

namespace detray::detail {

template <std::random_access_iterator RandomIt>
DETRAY_HOST_DEVICE inline void insertion_sort(RandomIt first, RandomIt last) {
  for (auto it = first; it != last; it++) {
    // Searching the upper bound, i.e., first
    // element greater than *it from beginning
    auto const insertion_point = detray::detail::upper_bound(first, it, *it);

    // Shifting the unsorted part
    std::ranges::rotate(insertion_point, it, it + 1);
  }
}

// Function to sort the array
template <template <typename...> class vector_t, typename TYPE>
DETRAY_HOST_DEVICE inline void insertion_sort(vector_t<TYPE> &vec) {
  detray::detail::insertion_sort(vec.begin(), vec.end());
}

template <std::random_access_iterator RandomIt, class Comp = std::less<void>>
DETRAY_HOST_DEVICE inline void selection_sort(RandomIt first, RandomIt last,
                                              Comp &&comp = Comp()) {
  for (RandomIt i = first; i < (last - 1); ++i) {
    RandomIt k = i;

    for (RandomIt j = i + 1; j < last; ++j) {
      if (comp(*j, *k)) {
        k = j;
      }
    }

    if (k != i) {
      auto t = *i;
      *i = *k;
      *k = t;
    }
  }
}

// Function to sort the array
template <template <typename...> class vector_t, typename TYPE>
DETRAY_HOST_DEVICE inline void selection_sort(vector_t<TYPE> &vec) {
  detray::detail::selection_sort(vec.begin(), vec.end());
}

template <std::size_t N, std::size_t P, std::size_t K, std::size_t J,
          std::size_t I, std::random_access_iterator RandomIt, class Comp>
DETRAY_HOST_DEVICE inline void batcher_odd_even_merge_network_sort_i(
    RandomIt first, Comp&& comp) {
  static constexpr std::size_t I0 = I + J;
  static constexpr std::size_t I1 = I + J + K;

  if constexpr (I0 / (2 * P) == I1 / (2 * P)) {
    if (comp(first[I1], first[I0])) {
      const auto tmp = first[I0];
      first[I0] = first[I1];
      first[I1] = tmp;
    }
  }

  if constexpr ((I + 1) <= std::min(K - 1, N - J - K - 1)) {
    batcher_odd_even_merge_network_sort_i<N, P, K, J, I + 1>(first, comp);
  }
}

template <std::size_t N, std::size_t P, std::size_t K, std::size_t J,
          std::random_access_iterator RandomIt, class Comp>
DETRAY_HOST_DEVICE inline void batcher_odd_even_merge_network_sort_j(
    RandomIt first, Comp&& comp) {
  batcher_odd_even_merge_network_sort_i<N, P, K, J, 0>(first, comp);

  if constexpr ((J + 2 * K) <= (N - 1 - K)) {
    batcher_odd_even_merge_network_sort_j<N, P, K, J + 2 * K>(first, comp);
  }
}

template <std::size_t N, std::size_t P, std::size_t K,
          std::random_access_iterator RandomIt, class Comp>
DETRAY_HOST_DEVICE inline void batcher_odd_even_merge_network_sort_k(
    RandomIt first, Comp&& comp) {
  batcher_odd_even_merge_network_sort_j<N, P, K, K % P>(first, comp);

  if constexpr (K >= 2) {
    batcher_odd_even_merge_network_sort_k<N, P, K / 2>(first, comp);
  }
}

template <std::size_t N, std::size_t P, std::random_access_iterator RandomIt,
          class Comp>
DETRAY_HOST_DEVICE inline void batcher_odd_even_merge_network_sort_p(
    RandomIt first, Comp&& comp) {
  batcher_odd_even_merge_network_sort_k<N, P, P>(first, comp);

  if constexpr ((2 * P) < N) {
    batcher_odd_even_merge_network_sort_p<N, 2 * P>(first, comp);
  }
}

template <std::size_t N, std::random_access_iterator RandomIt,
          class Comp = std::less<void>>
DETRAY_HOST_DEVICE inline void batcher_odd_even_merge_network_sort(
    RandomIt first, Comp&& comp = Comp()) {
  if constexpr (N > 1) {
    batcher_odd_even_merge_network_sort_p<N, 1>(first, comp);
  }
}

}  // namespace detray::detail
