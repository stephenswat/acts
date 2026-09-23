// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "detray/definitions/detail/qualifiers.hpp"
#include <utility>

namespace detray::detail {
/// Call @param f once per index in [0, N), with the index as a compile-time
/// constant
/// @{
template <typename F, std::size_t... I>
DETRAY_HOST_DEVICE constexpr void constexpr_for(F &&f,
                                                std::index_sequence<I...>) {
  (f.template operator()<I>(), ...);
}

template <std::size_t N, typename F>
DETRAY_HOST_DEVICE constexpr void constexpr_for(F &&f) {
  constexpr_for(std::forward<F>(f), std::make_index_sequence<N>{});
}
/// @}
}  // namespace detray::detail
