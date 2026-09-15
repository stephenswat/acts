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
#include "detray/geometry/coordinates/cartesian2D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"

namespace detray::detail {

/// @brief The frame in which the unbounded surface is intersected.
///
/// A polar plane is intersected exactly like a cartesian plane, only the mask
/// check differs. Mapping the frame here lets the two shapes share one
/// intersector type.
/// @{
template <typename frame_t>
struct intersection_frame {
  using type = frame_t;
};

template <concepts::algebra algebra_t>
struct intersection_frame<polar2D<algebra_t>> {
  using type = cartesian2D<algebra_t>;
};

template <typename frame_t>
using intersection_frame_t = typename intersection_frame<frame_t>::type;
/// @}

}  // namespace detray::detail
