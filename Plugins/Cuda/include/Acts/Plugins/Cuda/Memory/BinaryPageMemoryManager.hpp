// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <vector>
#include <memory>

namespace Acts::Cuda {
  class BinaryPageMemoryManager {
  public:
    ~BinaryPageMemoryManager();

    BinaryPageMemoryManager(const BinaryPageMemoryManager&) = delete;
    BinaryPageMemoryManager(BinaryPageMemoryManager&&) = delete;

    BinaryPageMemoryManager& operator=(const BinaryPageMemoryManager&) = delete;
    BinaryPageMemoryManager& operator=(BinaryPageMemoryManager&&) = delete;

    static BinaryPageMemoryManager& instance();

    void setMemorySize(std::size_t sizeInBytes, int device = -1);

    std::size_t availableMemory(int device = -1) const;

    void* allocate(std::size_t sizeInBytes, int device = -1);

    void reset(int device = -1);
  private:
    BinaryPageMemoryManager() = default;

    enum class PageState {
        OCCUPIED,
        VACANT,
        SPLIT
    };

    struct Page {
        bool root = false;
        PageState state;
        std::size_t size;
        void * addr;
        Page * sibling = nullptr;
        Page * left = nullptr;
        Page * right = nullptr;
        Page * parent = nullptr;

        ~Page();
        std::size_t maximum_contiguous();
        void deep_free();
        void shallow_free();
        void split();
        void unsplit();
    };

    std::vector<Page *> pages;
  };
}
