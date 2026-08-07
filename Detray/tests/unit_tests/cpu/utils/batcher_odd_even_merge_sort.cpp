// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/definitions/algorithms.hpp"

// Google Test include(s).
#include <gtest/gtest.h>

namespace {
template <typename comp_t, typename it_t>
bool is_sorted(it_t i, std::size_t n) {
  for (std::size_t j = 0; j < n - 1; ++j) {
    if (!comp_t{}(i[j], i[j + 1])) {
      return false;
    }
  }

  return true;
}

template <typename it_t>
void fill_unsorted_array(it_t i, std::size_t n) {
  for (std::size_t j = 0; j < n; ++j) {
    if (j % 2 == 0) {
      i[j] = static_cast<std::decay_t<decltype(i[j])>>(1000 - j);
    } else {
      i[j] = static_cast<std::decay_t<decltype(i[j])>>(j);
    }
  }
}
}  // namespace

template <typename T>
class detray_utils_network_sort : public ::testing::Test {};

using TestTypes = ::testing::Types<
    std::tuple<std::integral_constant<std::size_t, 1>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 1>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 1>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 1>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 1>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 1>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 2>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 2>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 2>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 2>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 2>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 2>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 3>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 3>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 3>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 3>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 3>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 3>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 4>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 4>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 4>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 4>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 4>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 4>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 5>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 5>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 5>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 5>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 5>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 5>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 6>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 6>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 6>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 6>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 6>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 6>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 7>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 7>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 7>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 7>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 7>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 7>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 8>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 8>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 8>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 8>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 8>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 8>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 9>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 9>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 9>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 9>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 9>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 9>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 10>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 10>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 10>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 10>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 10>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 10>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 11>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 11>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 11>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 11>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 11>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 11>, int, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 12>, double, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 12>, double, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 12>, float, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 12>, float, std::greater<>>,
    std::tuple<std::integral_constant<std::size_t, 12>, int, std::less<>>,
    std::tuple<std::integral_constant<std::size_t, 12>, int, std::greater<>>>;
TYPED_TEST_SUITE(detray_utils_network_sort, TestTypes, );

TYPED_TEST(detray_utils_network_sort, batcher_odd_even_merge_sort) {
  using size_constant_t = std::tuple_element_t<0, TypeParam>;
  using scalar_t = std::tuple_element_t<1, TypeParam>;
  using comp_t = std::tuple_element_t<2, TypeParam>;

  std::array<scalar_t, size_constant_t::value> vec;

  fill_unsorted_array(vec.begin(), size_constant_t::value);

  if constexpr (size_constant_t::value > 2) {
    ASSERT_FALSE(is_sorted<comp_t>(vec, size_constant_t::value));
  }

  detray::batcher_odd_even_merge_network_sort<size_constant_t::value>(
      vec.begin(), comp_t{});

  ASSERT_TRUE(is_sorted<comp_t>(vec, size_constant_t::value));
}
