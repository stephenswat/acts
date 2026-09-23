// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <algorithm>
#include <array>
#include <random>

#include <detray/definitions/algorithms.hpp>
#include <detray/navigation/detail/candidate_order.hpp>
#include <gtest/gtest.h>

template <std::size_t N>
void check() {
  using cache_t = detray::navigation::detail::ordered_candidate_cache<int, N>;
  static_assert(std::random_access_iterator<typename cache_t::iterator>);
  static_assert(std::random_access_iterator<typename cache_t::const_iterator>);
  cache_t cache;
  std::array<int, N> reference;
  for (std::size_t i = 0; i < N; ++i)
    cache[i] = reference[i] = 1000000;
  std::mt19937 rng(42);
  for (int round = 0; round < 10000; ++round) {
    const int value = static_cast<int>(rng() % 1000);
    auto pos = std::upper_bound(reference.begin(), reference.end(), value) -
               reference.begin();
    if (pos < N) {
      for (std::size_t i = N - 1; i > static_cast<std::size_t>(pos); --i)
        reference[i] = reference[i - 1];
      reference[pos] = value;
      cache.insert_at(pos, value);
    }
    if (round % 19 == 0) {
      for (std::size_t i = 0; i < N; ++i)
        cache[i] = reference[i] = static_cast<int>(rng() % 1000);
      detray::sequential_sort(cache.begin(), cache.end());
      std::sort(reference.begin(), reference.end());
    }
    EXPECT_TRUE(std::equal(cache.begin(), cache.end(), reference.begin()));
    const cache_t& c = cache;
    typename cache_t::const_iterator first = cache.begin();
    EXPECT_EQ(c.end() - first, N);
  }
}
TEST(candidate_order, insert_and_sort) {
  check<2>();
  check<4>();
  check<8>();
  check<16>();
}
