// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Memory/CubMemoryManager.hpp"

#include "../Utilities/ErrorCheck.cuh"

// CUB include(s).
#include <cub/util_allocator.cuh>

// System include(s).
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace Acts::Cuda {
  CubMemoryManager::~CubMemoryManager() {
    allocator->~CachingDeviceAllocator();
  }

  CubMemoryManager& CubMemoryManager::instance() {
    static CubMemoryManager mm;
    return mm;
  }

  void CubMemoryManager::setMemorySize(std::size_t sizeInBytes, int device) {
    allocator->SetMaxCachedBytes(sizeInBytes);
    max_capacity = sizeInBytes;
  }

  std::size_t CubMemoryManager::availableMemory(int device) const {
    return max_capacity - allocator->cached_bytes[device].live;
  }

  void* CubMemoryManager::allocate(std::size_t sizeInBytes, int device) {
    if (device == -1) {
      ACTS_CUDA_ERROR_CHECK(cudaGetDevice(&device));
    }

    if (sizeInBytes > availableMemory()) {
      throw std::bad_alloc();
    }

    void * rv;

    allocator->DeviceAllocate(device, &rv, sizeInBytes);

    allocations.push_back(rv);

    return rv;
  }

  void CubMemoryManager::reset(int device) {
    for (void * p : allocations) {
      allocator->DeviceFree(p);
    }

    allocations.clear();

  }

  CubMemoryManager::CubMemoryManager() {
    allocator = new cub::CachingDeviceAllocator(4, 5, 15, -1);
  }
}