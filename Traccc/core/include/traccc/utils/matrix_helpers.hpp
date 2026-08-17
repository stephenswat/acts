/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"

// Detray include(s).
#include <detray/utils/known_substructure_matrix.hpp>

// System include(s)
#include <cassert>
#include <cstddef>
#include <utility>

namespace traccc {

/// @brief Compute a masked inverse of a 2x2 matrix.
///
/// For 1D measurements only the [0, 0] element of the input matrix is
/// meaningful, so we compute 1 / M[0,0] and zero out the rest. For 2D
/// measurements a regular 2x2 inverse is computed.
///
/// @param M   the 2x2 matrix to invert
/// @param dim the meaningful dimension of @c M (1 or 2)
/// @returns the masked inverse of @c M
template <detray::concepts::algebra algebra_t>
TRACCC_HOST_DEVICE inline detray::dmatrix<algebra_t, 2, 2> masked_inverse(
    const detray::dmatrix<algebra_t, 2, 2>& M, const unsigned int dim) {
  assert(dim == 1u || dim == 2u);

  detray::dmatrix<algebra_t, 2, 2> M_inv;
  if (dim == 1u) {
    assert(getter::element(M, 0u, 0u) != 0.f);
    getter::element(M_inv, 0u, 0u) = 1.f / getter::element(M, 0u, 0u);
    getter::element(M_inv, 0u, 1u) = 0.f;
    getter::element(M_inv, 1u, 0u) = 0.f;
    getter::element(M_inv, 1u, 1u) = 0.f;
  } else {
    M_inv = matrix::inverse(M);
  }
  return M_inv;
}

/// @brief The substructures of the matrices a Kalman update runs on
///
/// These are the three facts that make the products of an update cheap. The
/// library skips the multiplications by a structural zero, so naming the
/// structure here is what removes the work.
/// @{
/// A track covariance, which is symmetric.
template <std::size_t N>
using track_covariance_substructure =
    typename detray::ksm::make_symmetric_substructure<N>::canonical_type;

/// An observation model. A measurement subspace only ever selects
/// @c e_bound_loc0 and @c e_bound_loc1, so every further column is zero.
template <std::size_t D, std::size_t N>
using observation_model_substructure =
    typename detray::ksm::make_left_columns_substructure<D, N,
                                                         2u>::canonical_type;

/// A measurement covariance, which is diagonal:
/// @c edm::get_measurement_covariance writes its off-diagonal elements as
/// exact zeros.
template <std::size_t D>
using measurement_covariance_substructure =
    typename detray::ksm::make_diagonal_substructure<D>::canonical_type;
/// @}

/// @brief Read a dense square matrix into a matrix of symmetric substructure.
///
/// A symmetric substructure keeps a single scalar for the pair (i, j) and
/// (j, i), so only the upper triangle of @p m is read and the lower one is
/// discarded.
///
/// @c ksm::matrix::from_dense cannot serve here, because it asserts that the
/// two mirrored elements hold the same value. A transported track covariance
/// is symmetric in exact arithmetic only: the generated
/// @c transport_covariance_to_bound_impl builds every element of J*C*J^T as
/// its own expression, so two mirrored elements differ in the last bits.
///
/// @tparam ksm_matrix_t the destination type, whose substructure must be
///                      square and symmetric
/// @param M the dense matrix to read
///
/// @returns @c M with its lower triangle discarded
template <typename ksm_matrix_t, typename dense_t>
TRACCC_HOST_DEVICE inline ksm_matrix_t symmetric_from_dense(const dense_t& M) {
  using substructure_t = typename ksm_matrix_t::canonical_substructure_type;
  using scalar_t = typename ksm_matrix_t::scalar_type;

  static_assert(detray::ksm::detail::checks::is_symmetric<substructure_t>,
                "the destination substructure is not symmetric");

  constexpr std::size_t N = substructure_t::rows;
  static_assert(static_cast<std::size_t>(detray::traits::rows<dense_t>) == N);
  static_assert(static_cast<std::size_t>(detray::traits::columns<dense_t>) ==
                N);

  ksm_matrix_t rv{};

  const auto read = [&rv, &M]<std::size_t I, std::size_t J>() {
    if constexpr (I <= J) {
      rv.template at<I, J>() = static_cast<scalar_t>(getter::element<I, J>(M));
    }
  };

  [&read]<std::size_t... Ks>(std::index_sequence<Ks...>) {
    (read.template operator()<Ks / N, Ks % N>(), ...);
  }(std::make_index_sequence<N * N>{});

  return rv;
}

}  // namespace traccc
