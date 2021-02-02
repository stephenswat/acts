// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Memory/NaiveMemoryManager.hpp"

#include "../Utilities/ErrorCheck.cuh"

// CUB include(s).
#include <cuda_runtime.h>

// System include(s).
#include <stdexcept>

namespace Acts::Cuda {
  NaiveMemoryManager& NaiveMemoryManager::instance() {
    static NaiveMemoryManager mm;
    return mm;
  }

  void NaiveMemoryManager::setMemorySize(std::size_t sizeInBytes, int device) {
    max_capacity = sizeInBytes;
  }

  std::size_t NaiveMemoryManager::availableMemory(int device) const {
    return max_capacity - total_size;
  }

  void* NaiveMemoryManager::allocate(std::size_t sizeInBytes, int device) {
    if (device == -1) {
      ACTS_CUDA_ERROR_CHECK(cudaGetDevice(&device));
    }

    if (sizeInBytes > availableMemory()) {
      throw std::bad_alloc();
    }

    void * rv;

    cudaMalloc(&rv, sizeInBytes);

    allocations.push_back(rv);

    total_size += sizeInBytes;

    return rv;
  }

  void NaiveMemoryManager::reset(int device) {
    for (void * p : allocations) {
      cudaFree(p);
    }

    allocations.clear();

    total_size = 0;
  }
}