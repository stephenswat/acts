// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// CUDA plugin include(s).
#include "Acts/Plugins/Cuda/Memory/BinaryPageMemoryManager.hpp"

#include "../Utilities/ErrorCheck.cuh"

// CUB include(s).
#include <cub/util_allocator.cuh>

// System include(s).
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <stack>

namespace {
  std::size_t round_down(std::size_t size) {
    for (unsigned short i = 0; i <= 32; i++) {
      std::size_t s = 2UL << i;

      if (s * 2 > size) {
        return s;
      }
    }

    return 0;
  }

  std::size_t round_up(std::size_t size) {
    for (unsigned short i = 0; i <= 32; i++) {
      std::size_t s = 2UL << i;

      if (s >= size) {
        return s;
      }
    }

    return 0;
  }
}

namespace Acts::Cuda {
  BinaryPageMemoryManager::~BinaryPageMemoryManager() {
    for (Page * p : pages) {
      if (p->root) {
        cudaFree(p->addr);
      }

      delete p;
    }
  }

  BinaryPageMemoryManager& BinaryPageMemoryManager::instance() {
    static BinaryPageMemoryManager mm;
    return mm;
  }

  void BinaryPageMemoryManager::setMemorySize(std::size_t sizeInBytes, int device) {
    void * alloc;

    cudaMalloc(&alloc, sizeInBytes);

    std::size_t remaining = sizeInBytes;

    void * cur_ptr = alloc;

    while (remaining >= 512) {
      Page * p = new Page();
      p->root = (cur_ptr == alloc);
      p->state = PageState::VACANT;
      p->addr = cur_ptr;
      p->size = round_down(remaining);

      remaining -= p->size;
      cur_ptr = static_cast<void *>(static_cast<char *>(cur_ptr) + p->size);
      pages.push_back(p);
    }
  }

  std::size_t BinaryPageMemoryManager::availableMemory(int device) const {
    std::size_t max_size = 0;

    for (Page * p : pages) {
      max_size = std::max(p->maximum_contiguous(), max_size);
    }

    return max_size;
  }

  void* BinaryPageMemoryManager::allocate(std::size_t sizeInBytes, int device) {
    Page * cand = nullptr;
    std::size_t goal = round_up(sizeInBytes);

    std::stack<Page *> rem;

    for (Page * p : pages) {
      rem.push(p);
    }

    while (!rem.empty()) {
      Page * c = rem.top();
      rem.pop();

      if (c->state == PageState::VACANT) {
        if (c->size >= goal && (cand == nullptr || c->size < cand->size)) {
          cand = c;
        }

        if (cand != nullptr && cand->size == goal) {
          break;
        }
      } else if (c->state == PageState::SPLIT) {
        rem.push(c->left);
        rem.push(c->right);
      }
    }

    if (cand == nullptr) {
      throw std::bad_alloc();
    }

    while (cand->size > goal) {
      cand->split();
      cand = cand->left;
    }

    cand->state = PageState::OCCUPIED;

    return cand->addr;
  }

  void BinaryPageMemoryManager::reset(int device) {
    for (Page * p : pages) {
      p->shallow_free();
    }
  }

  BinaryPageMemoryManager::Page::~Page() {
    if (state == PageState::SPLIT) {
      delete left;
      delete right;
    }
  }

  std::size_t BinaryPageMemoryManager::Page::maximum_contiguous() {
    if (state == PageState::OCCUPIED) {
      return 0;
    } else if (state == PageState::VACANT) {
      return size;
    } else {
      return std::max(left->maximum_contiguous(), right->maximum_contiguous());
    }
  }

  void BinaryPageMemoryManager::Page::deep_free() {
    shallow_free();
    unsplit();
  }

  void BinaryPageMemoryManager::Page::shallow_free() {
    if (state == PageState::OCCUPIED) {
      state = PageState::VACANT;
    } else if (state == PageState::SPLIT) {
      left->shallow_free();
      right->shallow_free();
    }
  }

  void BinaryPageMemoryManager::Page::split() {
    if (state != PageState::VACANT) {
      throw std::runtime_error("Can't split!");
    }

    state = PageState::SPLIT;

    left = new Page();
    right = new Page();

    left->state = PageState::VACANT;
    right->state = PageState::VACANT;

    left->size = size / 2;
    right->size = size / 2;

    left->addr = addr;
    right->addr = static_cast<void *>(static_cast<char *>(left->addr) + left->size);

    left->sibling = right;
    right->sibling = left;

    left->parent = this;
    right->parent = this;
  }

  void BinaryPageMemoryManager::Page::unsplit() {
    if (state != PageState::SPLIT) {
      return;
    }

    if (left->state != PageState::VACANT || right->state != PageState::VACANT) {
      return;
    }

    delete left;
    delete right;

    state = PageState::VACANT;
  }
}