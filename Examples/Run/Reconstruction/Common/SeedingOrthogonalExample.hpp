// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Geometry/GeometryIdentifier.hpp"

namespace ActsExamples {
class IBaseDetector;
}

/// The Propagation example
///
///
/// @param argc the number of argumetns of the call
/// @param argv the argument list
/// @param detector The detector descriptor instance
int runSeedingOrthogonalExample(int argc, char* argv[],
                      std::shared_ptr<ActsExamples::IBaseDetector> detector);
