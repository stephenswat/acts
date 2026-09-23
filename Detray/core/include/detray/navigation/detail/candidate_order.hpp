// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <cassert>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <type_traits>

namespace detray::navigation::detail {

/// @brief Maps the position of a candidate in the sorted navigation cache onto
///        the slot it occupies in the candidate storage.
///
/// The navigation cache is kept sorted by distance to the track position.
/// Instead of moving the (comparatively large) candidates themselves whenever
/// a new candidate is inserted, only this permutation is reordered while the
/// candidate storage is left untouched.
///
/// For a cache of up to eight candidates the entire permutation is packed into
/// a single 64-bit word, which can be held in a register: reordering the cache
/// then costs a handful of ALU instructions, whereas moving the candidates
/// costs a loop of candidate-sized copies through (local) memory.
///
/// @tparam N the capacity of the navigation cache
template <std::size_t N, bool = (N <= 8u)>
class candidate_order;

/// Packed permutation: one byte per cache position in a single word
template <std::size_t N>
class candidate_order<N, true> {
  using word_t = std::uint64_t;

  static constexpr std::size_t bits_per_slot{8u};
  static constexpr word_t slot_mask{0xffull};

  /// @returns the permutation that maps every cache position onto itself
  static constexpr word_t identity() {
    word_t order{0u};
    for (std::size_t i = 0u; i < N; ++i) {
      order |= static_cast<word_t>(i) << (bits_per_slot * i);
    }
    return order;
  }

  static constexpr word_t k_identity{identity()};

  /// @returns a mask of the @param n least significant bits
  DETRAY_HOST_DEVICE
  static constexpr word_t bits_below(const std::size_t n) {
    return (n >= bits_per_slot * sizeof(word_t))
               ? ~word_t{0}
               : ((word_t{1} << n) - word_t{1});
  }

  /// The storage slot of every cache position, one byte each
  word_t m_order{k_identity};

 public:
  /// @returns the storage slot of the candidate at cache position @param pos
  DETRAY_HOST_DEVICE
  constexpr std::size_t operator[](const std::size_t pos) const {
    assert(pos < N);
    return static_cast<std::size_t>((m_order >> (bits_per_slot * pos)) &
                                    slot_mask);
  }

  /// Map every cache position back onto itself
  DETRAY_HOST_DEVICE
  constexpr void reset() { m_order = k_identity; }

  /// Exchange the candidates at cache positions @param i and @param j
  DETRAY_HOST_DEVICE
  constexpr void swap(const std::size_t i, const std::size_t j) {
    assert(i < N);
    assert(j < N);

    const word_t diff{
        ((m_order >> (bits_per_slot * i)) ^ (m_order >> (bits_per_slot * j))) &
        slot_mask};

    m_order ^= (diff << (bits_per_slot * i)) | (diff << (bits_per_slot * j));
  }

  /// Free up cache position @param pos by moving the candidates from that
  /// position onwards one position towards the back of the cache. The
  /// candidate in the last position drops out of the cache.
  ///
  /// @returns the storage slot that is up for reuse, i.e. the slot that the
  ///          new candidate at position @param pos has to be written to
  DETRAY_HOST_DEVICE
  constexpr std::size_t make_room_at(const std::size_t pos) {
    assert(pos < N);

    // The evicted candidate leaves its storage slot behind
    const word_t free_slot{(m_order >> (bits_per_slot * (N - 1u))) & slot_mask};

    // Positions before @param pos keep their candidate...
    const word_t kept{m_order & bits_below(bits_per_slot * pos)};
    // ...the ones from @param pos onwards move up by one position
    const word_t moved{(m_order << bits_per_slot) &
                       (bits_below(bits_per_slot * N) &
                        ~bits_below(bits_per_slot * (pos + 1u)))};

    m_order = kept | (free_slot << (bits_per_slot * pos)) | moved;

    return static_cast<std::size_t>(free_slot);
  }
};

/// Fallback for caches that hold more candidates than fit into a single word
template <std::size_t N>
class candidate_order<N, false> {
  static_assert(N <= 256u, "Navigation cache capacity is limited to 256");

  using slot_t = std::uint8_t;

  darray<slot_t, N> m_order{};

 public:
  DETRAY_HOST_DEVICE
  constexpr candidate_order() { reset(); }

  /// @returns the storage slot of the candidate at cache position @param pos
  DETRAY_HOST_DEVICE
  constexpr std::size_t operator[](const std::size_t pos) const {
    assert(pos < N);
    return static_cast<std::size_t>(m_order[pos]);
  }

  /// Map every cache position back onto itself
  DETRAY_HOST_DEVICE
  constexpr void reset() {
    for (std::size_t i = 0u; i < N; ++i) {
      m_order[i] = static_cast<slot_t>(i);
    }
  }

  /// Exchange the candidates at cache positions @param i and @param j
  DETRAY_HOST_DEVICE
  constexpr void swap(const std::size_t i, const std::size_t j) {
    assert(i < N);
    assert(j < N);

    const slot_t tmp{m_order[i]};
    m_order[i] = m_order[j];
    m_order[j] = tmp;
  }

  /// @see candidate_order<N, true>::make_room_at
  DETRAY_HOST_DEVICE
  constexpr std::size_t make_room_at(const std::size_t pos) {
    assert(pos < N);

    const slot_t free_slot{m_order[N - 1u]};

    for (std::size_t i = N - 1u; i > pos; --i) {
      m_order[i] = m_order[i - 1u];
    }
    m_order[pos] = free_slot;

    return static_cast<std::size_t>(free_slot);
  }
};

/// @brief Iterates a navigation candidate cache in the order that is given by
///        a @c candidate_order permutation.
///
/// @tparam candidate_t the (possibly const qualified) candidate type
/// @tparam N the capacity of the navigation cache
template <typename candidate_t, std::size_t N>
class ordered_iterator {
 public:
  using order_type = candidate_order<N>;

  using iterator_category = std::random_access_iterator_tag;
  using iterator_concept = std::random_access_iterator_tag;
  using value_type = std::remove_cv_t<candidate_t>;
  using difference_type = std::ptrdiff_t;
  using pointer = candidate_t *;
  using reference = candidate_t &;

  constexpr ordered_iterator() = default;

  DETRAY_HOST_DEVICE
  constexpr ordered_iterator(candidate_t *const candidates,
                             const order_type order, const difference_type pos)
      : m_candidates{candidates}, m_order{order}, m_pos{pos} {}

  /// Allow the conversion from a mutable to a const iterator
  template <typename other_t>
    requires(!std::is_const_v<other_t> && std::is_const_v<candidate_t> &&
             std::same_as<std::remove_cv_t<other_t>, value_type>)
  DETRAY_HOST_DEVICE constexpr ordered_iterator(
      const ordered_iterator<other_t, N> &other)
      : m_candidates{other.data()},
        m_order{other.order()},
        m_pos{other.position()} {}

  /// @{ Access the underlying data (needed for the conversion above)
  DETRAY_HOST_DEVICE
  constexpr candidate_t *data() const { return m_candidates; }
  DETRAY_HOST_DEVICE
  constexpr order_type order() const { return m_order; }
  DETRAY_HOST_DEVICE
  constexpr difference_type position() const { return m_pos; }
  /// @}

  /// @{ Dereference
  DETRAY_HOST_DEVICE
  constexpr reference operator*() const {
    assert(m_candidates != nullptr);
    return m_candidates[m_order[static_cast<std::size_t>(m_pos)]];
  }
  DETRAY_HOST_DEVICE
  constexpr pointer operator->() const { return &(*(*this)); }
  DETRAY_HOST_DEVICE
  constexpr reference operator[](const difference_type i) const {
    assert(m_candidates != nullptr);
    return m_candidates[m_order[static_cast<std::size_t>(m_pos + i)]];
  }
  /// @}

  /// @{ Advance
  DETRAY_HOST_DEVICE
  constexpr ordered_iterator &operator++() {
    ++m_pos;
    return *this;
  }
  DETRAY_HOST_DEVICE
  constexpr ordered_iterator operator++(int) {
    ordered_iterator tmp{*this};
    ++m_pos;
    return tmp;
  }
  DETRAY_HOST_DEVICE
  constexpr ordered_iterator &operator--() {
    --m_pos;
    return *this;
  }
  DETRAY_HOST_DEVICE
  constexpr ordered_iterator operator--(int) {
    ordered_iterator tmp{*this};
    --m_pos;
    return tmp;
  }
  DETRAY_HOST_DEVICE
  constexpr ordered_iterator &operator+=(const difference_type i) {
    m_pos += i;
    return *this;
  }
  DETRAY_HOST_DEVICE
  constexpr ordered_iterator &operator-=(const difference_type i) {
    m_pos -= i;
    return *this;
  }
  /// @}

  /// @{ Arithmetic
  DETRAY_HOST_DEVICE
  friend constexpr ordered_iterator operator+(ordered_iterator itr,
                                              const difference_type i) {
    itr += i;
    return itr;
  }
  DETRAY_HOST_DEVICE
  friend constexpr ordered_iterator operator+(const difference_type i,
                                              ordered_iterator itr) {
    itr += i;
    return itr;
  }
  DETRAY_HOST_DEVICE
  friend constexpr ordered_iterator operator-(ordered_iterator itr,
                                              const difference_type i) {
    itr -= i;
    return itr;
  }
  DETRAY_HOST_DEVICE
  friend constexpr difference_type operator-(const ordered_iterator &lhs,
                                             const ordered_iterator &rhs) {
    return lhs.m_pos - rhs.m_pos;
  }
  /// @}

  /// @{ Comparison: iterators into the same cache only differ in position
  DETRAY_HOST_DEVICE
  friend constexpr bool operator==(const ordered_iterator &lhs,
                                   const ordered_iterator &rhs) {
    return lhs.m_pos == rhs.m_pos;
  }
  DETRAY_HOST_DEVICE
  friend constexpr auto operator<=>(const ordered_iterator &lhs,
                                    const ordered_iterator &rhs) {
    return lhs.m_pos <=> rhs.m_pos;
  }
  /// @}

 private:
  /// The candidate storage, in the order in which the candidates were found
  candidate_t *m_candidates{nullptr};
  /// Maps the cache position onto the storage slot
  order_type m_order{};
  /// The current position in the cache
  difference_type m_pos{0};
};

/// Fixed-size candidate storage whose logical order is a compact permutation.
/// Insertion invalidates iterators: each iterator snapshots the current order.
template <typename candidate_t, std::size_t N>
class ordered_candidate_cache {
 public:
  using value_type = candidate_t;
  using iterator = ordered_iterator<candidate_t, N>;
  using const_iterator = ordered_iterator<const candidate_t, N>;

  DETRAY_HOST_DEVICE constexpr iterator begin() {
    return {m_storage.data(), m_order, 0};
  }
  DETRAY_HOST_DEVICE constexpr iterator end() {
    return {m_storage.data(), m_order, N};
  }
  DETRAY_HOST_DEVICE constexpr const_iterator begin() const { return cbegin(); }
  DETRAY_HOST_DEVICE constexpr const_iterator end() const { return cend(); }
  DETRAY_HOST_DEVICE constexpr const_iterator cbegin() const {
    return {m_storage.data(), m_order, 0};
  }
  DETRAY_HOST_DEVICE constexpr const_iterator cend() const {
    return {m_storage.data(), m_order, N};
  }
  DETRAY_HOST_DEVICE constexpr candidate_t &operator[](std::size_t i) {
    return m_storage[m_order[i]];
  }
  DETRAY_HOST_DEVICE constexpr const candidate_t &operator[](
      std::size_t i) const {
    return m_storage[m_order[i]];
  }
  static constexpr std::size_t size() { return N; }
  static constexpr bool empty() { return N == 0; }

  DETRAY_HOST_DEVICE constexpr void insert_at(std::size_t pos,
                                              const candidate_t &candidate) {
    const auto slot = m_order.make_room_at(pos);
    m_storage[slot] = candidate;
  }

 private:
  darray<candidate_t, N> m_storage;
  candidate_order<N> m_order;
};

}  // namespace detray::navigation::detail
