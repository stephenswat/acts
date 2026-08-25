/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cstdint>

namespace traccc::device {
// Set of flags that we can assign to surfaces.
enum class surface_flags : std::uint8_t {
  e_is_line = 0b01,
  e_has_material = 0b10
};
}  // namespace traccc::device
