/*
 * KEnRef.cpp
 *
 *  Created on: 2023/5/8
 *      Author: amr
 */
#define VERBOSE false
// #include <omp.h>
#include <limits>
#include <memory>
#include <cstdlib>
#include <algorithm>
// #include <Eigen/Dense>
#include "core/KEnRef.h"
//#include <iostream>//for testing only

#include <cassert>

#include "core/IoUtils.h"
//#include <utility>

namespace {
// Minimum number of rows below which the row-blocked kernels stay serial (thread-spawn cost would
// dominate). Env-configurable via KENREF_RELAX_PARALLEL_THRESHOLD (default 256), read once.
Eigen::Index relax_parallel_threshold() {
    static const Eigen::Index threshold = []() -> Eigen::Index {
        if (const char *env = std::getenv("KENREF_RELAX_PARALLEL_THRESHOLD")) {
            char *end = nullptr;
            const long v = std::strtol(env, &end, 10);
            if (end != env && *end == '\0' && v >= 0) return static_cast<Eigen::Index>(v);
        }
        return 256;
    }();
    return threshold;
}

// Pair-row block size for the vector r_array_to_d_array's collapse(2) parallelization. Env-tunable via
// KENREF_DARRAY_BLOCK (default 256) so profiling can sweep it; a block >= N gives one block per model
// (i.e. model-axis-only parallelism), a clean A/B against per-pair blocking. Read once.
Eigen::Index darray_block_size() {
    static const Eigen::Index block = []() -> Eigen::Index {
        if (const char *env = std::getenv("KENREF_DARRAY_BLOCK")) {
            char *end = nullptr;
            const long v = std::strtol(env, &end, 10);
            if (end != env && *end == '\0' && v > 0) return static_cast<Eigen::Index>(v);
        }
        return 256;
    }();
    return block;
}
} // namespace

template<typename KEnRef_Real>
KEnRef_Real SpecDenData<KEnRef_Real>::getSigma(const std::tuple<std::string, std::string>& atomPairs) const {
    // const size_t pairId = atomNamePairs_to_atomPairIndex.at(atomPairs);
    // return getSigma(static_cast<const size_t>(pairId % sigma->size()));
    return getSigma(this->atomNamePairs_to_atomPairIndex.at(atomPairs) % sigma->size());
}

template<typename KEnRef_Real>
void SpecDenDataBase<KEnRef_Real>::setAtomPairs(std::vector<std::tuple<std::string, std::string>> &atomPairs){
    this->atom_pairs = atomPairs;
    this->atomNamePairs_to_atomPairIndex.clear();
    for (int i = 0; i < atomPairs.size(); ++i) {
        atomNamePairs_to_atomPairIndex[atomPairs[i]] = i;
    }
}

template<typename KEnRef_Real>
[[nodiscard]] KEnRef_Real SpecDenData<KEnRef_Real>::getSigma(const size_t atomPairId) const {
    return sigma.value()[atomPairId];
}

template<typename KEnRef_Real>
KEnRef<KEnRef_Real>::KEnRef() = default;

template<typename KEnRef_Real>
KEnRef<KEnRef_Real>::~KEnRef() = default;

//KEnRef::KEnRef(const KEnRef &other) = default;
//KEnRef::KEnRef(KEnRef &&other)  noexcept {}
//KEnRef& KEnRef::operator=(const KEnRef &other) {}
//KEnRef& KEnRef::operator=(KEnRef &&other) {}

template<typename KEnRef_Real>
std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >
KEnRef<KEnRef_Real>::array_shift(const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &arrays_to_shift,
                                 uint nShift, int numOmpThreads) {
    if (nShift == 1)
        return arrays_to_shift;
    if (arrays_to_shift.empty())
        return {};

    int numInteractions = arrays_to_shift[0].rows() / nShift;
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > ret(arrays_to_shift.size() * nShift);
#pragma omp parallel for num_threads(numOmpThreads)
    for (int i = 0; i < arrays_to_shift.size() * nShift; ++i) {
        ret.at(i) = std::move(Eigen::MatrixX<KEnRef_Real>(numInteractions, arrays_to_shift[0].cols()));
    }

#pragma omp parallel for collapse(3) num_threads(numOmpThreads)
    for (int interactionId = 0; interactionId < numInteractions; interactionId++) {
        //generalization of pairId
        for (int modelNo = 0; modelNo < arrays_to_shift.size(); modelNo++) {
            for (int interactionComponent = 0; interactionComponent < nShift; interactionComponent++) {
                // std::cout << "To ret[" <<(modelNo*nShift+interactionComponent) <<"]("<<interactionId<<", Eigen::placeholders::all)\t";
                // std::cout << "from arrays_to_shift["<<modelNo<<"]("<<(interactionId + interactionComponent * numInteractions)<<", Eigen::placeholders::all)" << std::endl;
                ret[modelNo * nShift + interactionComponent](interactionId, Eigen::placeholders::all) =
                        arrays_to_shift[modelNo](interactionId + interactionComponent * numInteractions, Eigen::placeholders::all);
                //TODO is there a way other than copying with "="?
            }
        }
    }
    return ret;
}

template<typename KEnRef_Real>
std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> >
KEnRef<KEnRef_Real>::r_array_to_d_array(const CoordsMatrixType<KEnRef_Real> &Nxyz, bool unit, bool gradient,
                                        bool addEpsilon, int /*numOmpThreads*/,
                                        Eigen::Index startRow, Eigen::Index numRows) {
    //	std::cout << "r_array_to_d_array(Nxyz) called" << std::endl;
    const bool dist = !unit;
    // Process the contiguous row-range [startRow, startRow+N). Default (startRow=0, numRows<0) => whole
    // matrix, identical to the original. Restricting to a sub-batch is bit-for-bit identical because every
    // operation below is per-row; this is what lets the vector overload fill its result in row-blocks.
    const Eigen::Index N = (numRows < 0) ? (Nxyz.rows() - startRow) : numRows;
    assert(startRow >= 0 && N >= 0 && startRow + N <= Nxyz.rows());

    const auto x = Nxyz.col(0).segment(startRow, N).array();
    const auto y = Nxyz.col(1).segment(startRow, N).array();
    const auto z = Nxyz.col(2).segment(startRow, N).array();
    constexpr auto half = static_cast<KEnRef_Real>(0.5);
    constexpr auto two = static_cast<KEnRef_Real>(2.0);
    constexpr auto sqrt3 = static_cast<KEnRef_Real>(1.7320508075688772);
    constexpr KEnRef_Real half_sqrt3 = sqrt3 * half;
    constexpr KEnRef_Real epsilon_val = std::numeric_limits<KEnRef_Real>::epsilon();

    // dist: coeff = -5 (r^5/r^7), unit: coeff = -2 (r^2/r^4)
    const KEnRef_Real neg_coeff = unit ? KEnRef_Real(-2) : KEnRef_Real(-5);
    const KEnRef_Real neg_coeff_sqrt3 = neg_coeff * sqrt3;
    const KEnRef_Real half_neg_coeff_sqrt3 = neg_coeff_sqrt3 * half;

    // Precompute common terms efficiently
    const auto x2 = x.square();
    const auto y2 = y.square();
    const auto z2 = z.square();
    const auto xy = x * y;
    const auto xz = x * z;
    const auto yz = y * z;
    const auto x2_minus_y2 = x2 - y2;
    const auto r2 = (x2 + y2 + z2).eval();
    // dist: r^5 = r^2*r^2*r, unit: r^2; if/else avoids sqrt for unit path
    const KEnRef_Real eps_addend = addEpsilon ? epsilon_val : KEnRef_Real(0);
    typename decltype(r2)::PlainObject r_pow;
    if (dist)
        r_pow = r2 * r2 * r2.sqrt() + eps_addend;
    else
        r_pow = r2 + eps_addend;
    const auto inv_r_pow = r_pow.inverse().eval();

    // Compute ret1 (d_array)
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> ret1(N, 5);
    // Manual assignment to avoid temporary matrices
    const auto half_minus_x2_y2_plus_z2 = ((-x2 - y2) * half) + z2;

    ret1.col(0) = half_minus_x2_y2_plus_z2 * inv_r_pow;
    ret1.col(1) = half_sqrt3 * x2_minus_y2 * inv_r_pow;
    ret1.col(2) = sqrt3 * xz * inv_r_pow;
    ret1.col(3) = sqrt3 * yz * inv_r_pow;
    ret1.col(4) = sqrt3 * xy * inv_r_pow;

    if (!gradient) {
        return {std::move(ret1), Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15>{}};
    }

    // dist: r^7 = r^5*r^2 + eps, unit: r^4 = r^2*r^2 + eps; eager so it's computed once for all gradient cols
    const auto inv_r_pow_grad = (r_pow * r2 + eps_addend).inverse().eval();
    // Precompute common gradient terms
    const auto sqrt3_inv_r_pow = sqrt3 * inv_r_pow;
    const auto neg_coeff_inv_r_pow_grad = neg_coeff * inv_r_pow_grad;
    const auto neg_coeff_sqrt3_inv_r_pow_grad = neg_coeff_sqrt3 * inv_r_pow_grad;
    const auto half_neg_coeff_sqrt3_inv_r_pow_grad = half_neg_coeff_sqrt3 * inv_r_pow_grad;

    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> ret2(N, 15);

    // Column 0-2
    const auto half_minus_x2_y2_plus_z2_times_neg_coeff_inv_r_pow_grad = half_minus_x2_y2_plus_z2 * neg_coeff_inv_r_pow_grad;
    ret2.col(0) = x * half_minus_x2_y2_plus_z2_times_neg_coeff_inv_r_pow_grad - x * inv_r_pow;
    ret2.col(1) = y * half_minus_x2_y2_plus_z2_times_neg_coeff_inv_r_pow_grad - y * inv_r_pow;
    ret2.col(2) = z * half_minus_x2_y2_plus_z2_times_neg_coeff_inv_r_pow_grad + two * z * inv_r_pow;

    // Column 3-5
    const auto x2_minus_y2_times_half_neg_coeff_sqrt3_inv_r_pow_grad = x2_minus_y2 * half_neg_coeff_sqrt3_inv_r_pow_grad;
    ret2.col(3) = x * x2_minus_y2_times_half_neg_coeff_sqrt3_inv_r_pow_grad + x * sqrt3_inv_r_pow;
    ret2.col(4) = y * x2_minus_y2_times_half_neg_coeff_sqrt3_inv_r_pow_grad - y * sqrt3_inv_r_pow;
    ret2.col(5) = z * x2_minus_y2_times_half_neg_coeff_sqrt3_inv_r_pow_grad;

    // Column 6-8
    const auto xz_times_neg_coeff_sqrt3_inv_r_pow_grad = xz * neg_coeff_sqrt3_inv_r_pow_grad;
    ret2.col(6) = x * xz_times_neg_coeff_sqrt3_inv_r_pow_grad + z * sqrt3_inv_r_pow;
    ret2.col(7) = y * xz_times_neg_coeff_sqrt3_inv_r_pow_grad;
    ret2.col(8) = z * xz_times_neg_coeff_sqrt3_inv_r_pow_grad + x * sqrt3_inv_r_pow;

    // Column 9-11
    const auto yz_times_neg_coeff_sqrt3_inv_r_pow_grad = yz * neg_coeff_sqrt3_inv_r_pow_grad;
    ret2.col(9) = x * yz_times_neg_coeff_sqrt3_inv_r_pow_grad;
    ret2.col(10) = y * yz_times_neg_coeff_sqrt3_inv_r_pow_grad + z * sqrt3_inv_r_pow;
    ret2.col(11) = z * yz_times_neg_coeff_sqrt3_inv_r_pow_grad + y * sqrt3_inv_r_pow;

    // Column 12-14
    const auto xy_times_neg_coeff_sqrt3_inv_r_pow_grad = xy * neg_coeff_sqrt3_inv_r_pow_grad;
    ret2.col(12) = x * xy_times_neg_coeff_sqrt3_inv_r_pow_grad + y * sqrt3_inv_r_pow;
    ret2.col(13) = y * xy_times_neg_coeff_sqrt3_inv_r_pow_grad + x * sqrt3_inv_r_pow;
    ret2.col(14) = z * xy_times_neg_coeff_sqrt3_inv_r_pow_grad;

    return {std::move(ret1), std::move(ret2)};
}

//	std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> KEnRef::r_array_to_d_array(const Eigen::MatrixX3f& Nxyz, bool gradient){
template<typename KEnRef_Real>
std::tuple<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >, std::vector<Eigen::Matrix<KEnRef_Real,
    Eigen::Dynamic, 15> > >
KEnRef<KEnRef_Real>::r_array_to_d_array(const std::vector<CoordsMatrixType<KEnRef_Real> > &models_Nxyz, bool unit,
                                        bool gradient, bool addEpsilon, int numOmpThreads) {
    //	std::cout << "r_array_to_d_array(models_Nxyz) called" << std::endl;
    // Fresh-allocation convenience wrapper: delegates to the buffer-reuse overload below. For a hot
    // per-step caller (e.g. MD refinement), prefer r_array_to_d_array_into() with persistent buffers to
    // avoid re-allocating (and re-faulting) the [N×5]/[N×15] results every call.
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > ret1;
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > ret2;
    r_array_to_d_array_into(models_Nxyz, unit, gradient, addEpsilon, numOmpThreads, ret1, ret2);
    return {std::move(ret1), std::move(ret2)};
}

template<typename KEnRef_Real>
void
KEnRef<KEnRef_Real>::r_array_to_d_array_into(const std::vector<CoordsMatrixType<KEnRef_Real> > &models_Nxyz,
                                             bool unit, bool gradient, bool addEpsilon, int numOmpThreads,
                                             std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &ret1,
                                             std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > &ret2) {
    const long M = static_cast<long>(models_Nxyz.size());
    // Buffer reuse: resize the caller's vectors/matrices to the needed shape. Eigen's resize() is a no-op
    // when the size is unchanged, so a caller that reuses the same ret1/ret2 across calls (same M, N) does
    // NOT re-allocate or re-fault those pages after the first call — eliminating the per-call result
    // allocation / first-touch cost. A fresh (empty) ret1/ret2 simply allocates once here, as before.
    ret1.resize(M);
    ret2.resize(M);
    if (M == 0) return;

    // All models share the same pair count (same atom_id_pairs), so N is taken once and the (model,
    // pair-block) grid is rectangular — required for collapse(2).
    const Eigen::Index N = models_Nxyz[0].rows();
    for (long m = 0; m < M; ++m) {
        assert(models_Nxyz[m].rows() == N);
        ret1[m].resize(N, 5);
        ret2[m].resize(gradient ? N : Eigen::Index(0), 15);
    }

    // The single-model kernel is pure (unthreaded) Eigen, so this wrapper is the only place its work gets
    // multi-threaded. Both the model axis AND the pair axis are free (no reduction), so collapse(2) over
    // (model, pair-block) exposes M*nBlocks independent units in one flat parallel region — this lifts the
    // ceiling of model-only parallelism for the common few-models/many-pairs case. Each unit fills disjoint
    // middleRows of its own model's ret1/ret2, so it is race-free and order-preserving => bitwise-identical
    // to serial. Block size is env-tunable (KENREF_DARRAY_BLOCK, default 256); a block >= N gives one block
    // per model, i.e. model-axis-only parallelism (a clean profiling A/B). Size-gated on N.
    const Eigen::Index BLOCK = darray_block_size();
    const Eigen::Index PARALLEL_THRESHOLD = relax_parallel_threshold();
    const long nBlocks = static_cast<long>((N + BLOCK - 1) / BLOCK);
#pragma omp parallel for collapse(2) num_threads(numOmpThreads) \
    if((M * nBlocks) > 1 && N >= PARALLEL_THRESHOLD) schedule(static)
    for (long m = 0; m < M; ++m) {
        for (long b = 0; b < nBlocks; ++b) {
            const Eigen::Index r0 = static_cast<Eigen::Index>(b) * BLOCK;
            const Eigen::Index len = std::min<Eigen::Index>(BLOCK, N - r0);
            auto arrs = r_array_to_d_array(models_Nxyz[m], unit, gradient, addEpsilon, numOmpThreads, r0, len);
            ret1[m].middleRows(r0, len) = std::move(std::get<0>(arrs));
            if (gradient) ret2[m].middleRows(r0, len) = std::move(std::get<1>(arrs));
        }
    }
}

template<typename KEnRef_Real>
std::tuple<std::vector<Eigen::VectorX<KEnRef_Real> >, std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic,
    5> > > >
KEnRef<KEnRef_Real>::d_array_to_g_multiple_groupings(
    const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &d_array,
    //vector (models<pairId, tensor_elements>) with interaction tensors
    const std::vector<std::vector<std::vector<int> > > &groupings,
    //groupings of models to average interaction tensors (per dipole-dipole interaction pair), i.e., outer list for pairId and inner list for modelId
    const bool gradient, const int numOmpThreads) {
    //	std::cout << "d_array_to_g_multiple_groupings() called" << std::endl;
    std::vector<Eigen::VectorX<KEnRef_Real> > ret1;
    std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > > ret2;
    ret1.reserve(d_array.size());

    for (const auto &grouping: groupings) {
        const auto &[ret1_temp, ret2_temp] = d_array_to_g(d_array, grouping, gradient, numOmpThreads);
        ret1.emplace_back(ret1_temp);
        if (gradient)
            ret2.emplace_back(ret2_temp);
    }
    return {ret1, ret2};
}

template<typename KEnRef_Real>
std::tuple<std::vector<Eigen::VectorX<KEnRef_Real> >, std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic,
    5> > > >
KEnRef<KEnRef_Real>::d_array_to_g_matrix(
    const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &d_array,
    //vector (models<pairId, tensor_elements>) with interaction tensors
    const Eigen::MatrixX<int> &grouping_mat,
    //groupings of models to average interaction tensors (per dipole-dipole interaction pair), i.e., outer list for pairId and inner list for modelId
    const bool gradient, const int numOmpThreads) {
    //	std::cout << "d_array_to_g_matrix() called" << std::endl;
    std::vector<Eigen::VectorX<KEnRef_Real> > ret1;
    std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > > ret2;
    ret1.reserve(d_array.size());

    const auto &groupings = IoUtils::grouping_mat_to_subset_idx(grouping_mat);
    for (int row = 0; row < grouping_mat.rows(); ++row) {
        const auto& grouping = groupings[row];
        const auto &[ret1_temp, ret2_temp] = d_array_to_g(d_array, grouping, gradient, numOmpThreads);
        ret1.emplace_back(ret1_temp);
        if (gradient)
            ret2.emplace_back(ret2_temp);
    }
    return {ret1, ret2};
}

template<typename KEnRef_Real>
std::tuple<Eigen::VectorX<KEnRef_Real>, std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > >
KEnRef<KEnRef_Real>::d_array_to_g(
    const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &d_arrays,
    const std::vector<std::vector<int> > &grouping,
    const bool gradient, const int numOmpThreads) {

    const auto num_models = d_arrays.size();
    const auto num_interactionIds = d_arrays[0].rows();

    // std::cout << "num_models " 	<< num_models << " num_pairIds " << num_pairIds << "num_groups " << num_groups << std::endl;
    Eigen::VectorX<KEnRef_Real> ret1 = Eigen::VectorX<KEnRef_Real>::Zero(num_interactionIds);

    //Every element of ret2 (every d_matrix_grad) is a matrix(num_pairIds x 5)
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > ret2;
    if (gradient) {
        ret2.resize(num_models);
#pragma omp parallel for num_threads(numOmpThreads)
        for (int i = 0; i < num_models; ++i) {
            ret2[i] = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(num_interactionIds, 5); // AKA d_matrix_grad
        }
    }// else: ret2 remains empty by default

    // Precompute constants outside grouping loop
    const auto TWO = static_cast<KEnRef_Real>(2.0); //precompute constant

    //for every grouping block
    for (const auto &currentGrouping: grouping) {
        const auto currentGroupSize = currentGrouping.size();
        if (currentGroupSize == 0) continue; //skip empty groups
        const auto CURRENT_GROUP_SIZE_real = static_cast<KEnRef_Real>(currentGroupSize);
        const auto currentGroupSize_OVER_num_models_real = CURRENT_GROUP_SIZE_real / num_models;

        //to carry "average dipole interaction tensor" every group
        // Use Eigen::Matrix instead of auto for better optimization
        Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> d_matrix = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(num_interactionIds, 5);

        // Thread-local accumulation then a single critical combine.
        // NOTE: a built-in `reduction(+:d_matrix)` here is unsafe: d_matrix is Eigen::Matrix<.,Dynamic,5>,
        // whose default ctor is 0x5, so OpenMP's default-constructed private reduction copy is empty ->
        // size mismatch on `+=`. (A user-defined reduction with an explicit zero-sized initializer would
        // work, but the gain over this thread-local path is negligible.) Each local accumulator below is
        // explicitly sized + zeroed (never default-constructed) for exactly that reason.
        #pragma omp parallel num_threads(numOmpThreads)
        {
            Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> local_d_matrix = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(num_interactionIds, 5);
            #pragma omp for nowait
            for (int j = 0; j < currentGroupSize; ++j) {
                local_d_matrix += d_arrays[currentGrouping[j]];
            }
            #pragma omp critical
            d_matrix += local_d_matrix;
        }

        if (gradient) {
            const auto TWO_OVER_num_models_currentGroupSize = TWO / (num_models * CURRENT_GROUP_SIZE_real); //single division
            // calculate  d_matrix_grad. d_matrix_grad shape is pairIDs * interaction tensor elements
            // Compute gradient once and assign to all group members
            const auto d_matrix_grad = d_matrix * TWO_OVER_num_models_currentGroupSize;
            // All models of the same group equally share the same value
#pragma omp parallel for num_threads(numOmpThreads)
            for (int j = 0; j < currentGroupSize; j++) {
                // All elements are equal (i.e., all models get the same overall (average ?) value at the end.
                ret2[currentGrouping[j]] = d_matrix_grad;
            }
        }

        // Average the d_matrix
        d_matrix /= CURRENT_GROUP_SIZE_real;

        // calculate self dot product (norm squared) and accumulate group's contribution to mean g
        // Each iteration writes a DISTINCT ret1(j); the cross-grouping accumulation happens in the
        // serial outer loop, so no reduction is needed. (A `reduction(+:ret1)` here was also broken:
        // OpenMP default-constructs the private Eigen vector empty (size 0) -> out-of-bounds.)
        #pragma omp parallel for num_threads(numOmpThreads)
        for (int j = 0; j < num_interactionIds; j++) {
            //It is safe to use the vector directly
            ret1(j) += d_matrix.row(j).squaredNorm() * currentGroupSize_OVER_num_models_real;
        }
    }
    //	std::cout << "ret1" << std::endl << ret1 << std::endl;
    return {std::move(ret1), std::move(ret2)}; //move semantics
}

template<typename KEnRef_Real>
std::tuple<Eigen::MatrixX<KEnRef_Real>, std::optional<Eigen::MatrixX<KEnRef_Real> > >
KEnRef<KEnRef_Real>::power_scaled_loss_function(
    const Eigen::MatrixX<KEnRef_Real>& g, // current group norm squared values
    const Eigen::MatrixX<KEnRef_Real>& g0, // target group norm squared values.
    KEnRef_Real k, // force constant
    KEnRef_Real n, // correction power.
    const bool gradient,
    int numOmpThreads
) {
    //Use Eigen's array operations for small matrices, manual for large ones
    const auto size = g.size();
    Eigen::MatrixX<KEnRef_Real> ret1;
    std::optional<Eigen::MatrixX<KEnRef_Real>> ret2;

    // Optimized: Special fast paths for common n values
    const bool is_n_half = (n == static_cast<KEnRef_Real>(0.5));
    const bool is_n_quarter = (n == static_cast<KEnRef_Real>(0.25));

    if (size < 1000) { //Threshold for Eigen vs manual optimization
        // Small matrices: Let Eigen optimize
        const auto& g_arr = g.array();
        const auto& g0_arr = g0.array();

        if (is_n_half) {
            // Optimized: n = 0.5: sqrt(1 + |x|) - 1 * sign(x)
            auto term1 = ((1.0 + g_arr.abs()).sqrt() - 1.0) * g_arr.sign();
            auto term2 = ((1.0 + g0_arr.abs()).sqrt() - 1.0) * g0_arr.sign();
            auto common = term1 - term2;
            ret1 = (k * common.square()).matrix();

            if (gradient) {
                // Gradient for n=0.5: 0.5 / sqrt(1 + |x|)
                ret2 = (2.0 * k * common * (0.5 / (1.0 + g_arr.abs()).sqrt())).matrix();
            }
        }
        else if (is_n_quarter) {
            // Optimized: n = 0.25: (1 + |x|)^0.25 - 1 * sign(x)
            auto term1 = ((1.0 + g_arr.abs()).sqrt().sqrt() - 1.0) * g_arr.sign(); // sqrt(sqrt()) = ^0.25
            auto term2 = ((1.0 + g0_arr.abs()).sqrt().sqrt() - 1.0) * g0_arr.sign();
            auto common = term1 - term2;
            ret1 = (k * common.square()).matrix();

            if (gradient) {
                // Gradient for n=0.25: 0.25 / ((1 + |x|)^0.75)
                ret2 = (2.0 * k * common * (0.25 / ((1.0 + g_arr.abs()).sqrt().sqrt() * (1.0 + g_arr.abs()).sqrt()))).matrix();
            }
        }
        else {
            // General case
            auto abs_g = g_arr.abs();
            auto abs_g0 = g0_arr.abs();
            auto common = ((1.0 + abs_g).pow(n) - 1.0) * g_arr.sign()
                        - ((1.0 + abs_g0).pow(n) - 1.0) * g0_arr.sign();
            ret1 = (k * common.square()).matrix();

            if (gradient) {
                ret2 = (2.0 * k * common * (n * (1.0 + abs_g).pow(n - KEnRef_Real{1}))).matrix(); //KEnRef_Real{1} prevents deprecation warning
            }
        }
    } else {
        // Large matrices: Manual parallelization with optimized paths
        const auto rows = g.rows();
        const auto cols = g.cols();
        ret1.resize(rows, cols);
        //This value may become infinity if it excceds 3.402823466E38 in a single precision float
        if (gradient) {
            ret2 = Eigen::MatrixX<KEnRef_Real>(rows, cols);
        }
        // Optimized pow() calls for common n values
        #pragma omp parallel for num_threads(numOmpThreads) collapse(2) schedule(static)
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                const KEnRef_Real g_val = g(i, j);
                const KEnRef_Real g0_val = g0(i, j);
                const KEnRef_Real abs_g = std::abs(g_val);
                const KEnRef_Real abs_g0 = std::abs(g0_val);
                const KEnRef_Real sign_g = (g_val > 0) ? 1.0 : (g_val < 0) ? -1.0 : 0.0;
                const KEnRef_Real sign_g0 = (g0_val > 0) ? 1.0 : (g0_val < 0) ? -1.0 : 0.0;

                KEnRef_Real term1, term2;

                if (is_n_half) {
                    // Optimized: n = 0.5 using sqrt (much faster than pow)
                    term1 = (std::sqrt(1.0 + abs_g) - 1.0) * sign_g;
                    term2 = (std::sqrt(1.0 + abs_g0) - 1.0) * sign_g0;
                }
                else if (is_n_quarter) {
                    // Optimized: n = 0.25 using nested sqrt
                    term1 = (std::sqrt(std::sqrt(1.0 + abs_g)) - 1.0) * sign_g;
                    term2 = (std::sqrt(std::sqrt(1.0 + abs_g0)) - 1.0) * sign_g0;
                }
                else {
                    // General case
                    term1 = (std::pow(1.0 + abs_g, n) - 1.0) * sign_g;
                    term2 = (std::pow(1.0 + abs_g0, n) - 1.0) * sign_g0;
                }

                const KEnRef_Real common = term1 - term2;
                ret1(i, j) = k * common * common;

                if (gradient) {
                    KEnRef_Real grad_term;
                    if (is_n_half) {
                        // Optimized: 0.5 / sqrt(1 + |x|)
                        grad_term = 0.5 / std::sqrt(1.0 + abs_g);
                    }
                    else if (is_n_quarter) {
                        // Optimized: 0.25 / ((1 + |x|)^0.75)
                        grad_term = 0.25 / (std::sqrt(std::sqrt(1.0 + abs_g)) * std::sqrt(1.0 + abs_g));
                    }
                    else {
                        grad_term = n * std::pow(1.0 + abs_g, n - 1.0);
                    }
                    (*ret2)(i, j) = 2.0 * k * common * grad_term;
                }
            }
        }
    }
    return {std::move(ret1), std::move(ret2)};
}

template<typename KEnRef_Real>
std::tuple<Eigen::MatrixX<KEnRef_Real>, std::optional<Eigen::MatrixX<KEnRef_Real> > >
KEnRef<KEnRef_Real>::log_abs_diff_over_optimum_loss_function(
Eigen::MatrixX<KEnRef_Real> g, // current group norm squared values
Eigen::MatrixX<KEnRef_Real> g0, // target group norm squared values.
KEnRef_Real k, // force constant
const bool gradient,
int numOmpThreads
) {
    const auto &g_arr = g.array() /*+ std::numeric_limits<KEnRef_Real_t>::epsilon()*/;
    const auto &g0_arr = g0.array();
    //    std::cout << "g_arr " << g_arr << std::endl;
    Eigen::MatrixX<KEnRef_Real> ret1;
    Eigen::MatrixX<KEnRef_Real> ret2;

    const auto &common = g0_arr + Eigen::abs(g_arr - g0_arr);
    ret1 = (k * Eigen::log(common / g0_arr)).matrix();
    if (gradient) {
        ret2 = (k * Eigen::sign(g_arr - g0_arr) / common).matrix();
    }
    return {ret1, ret2};
}



template<typename KEnRef_Real>
Eigen::MatrixX<KEnRef_Real>
KEnRef<KEnRef_Real>::vectorOfVectors_to_Matrix(
    std::vector<Eigen::VectorX<KEnRef_Real> > g_vect, int numOmpThreads) {
    Eigen::MatrixX<KEnRef_Real> g_mat(g_vect[0].rows(), g_vect.size());
#pragma omp parallel for num_threads(numOmpThreads)
    for (int i = 0; i < g_vect.size(); ++i) {
        g_mat.col(i) = g_vect[i];
    }
    return g_mat;
}

template<typename KEnRef_Real>
std::vector<CoordsMatrixType<KEnRef_Real> >
KEnRef<KEnRef_Real>::coord_array_to_r_array(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    const std::vector<std::tuple<int, int> > &atomId_pairs, int numOmpThreads) {
    std::vector<CoordsMatrixType<KEnRef_Real> > ret(coord_array.size());
    for (int model_no = 0; model_no < coord_array.size(); ++model_no) {
        ret.at(model_no) = std::move(coord_array_to_r_array(coord_array[model_no], atomId_pairs, numOmpThreads));
    }
    return ret;
}

template<typename KEnRef_Real>
CoordsMatrixType<KEnRef_Real>
KEnRef<KEnRef_Real>::coord_array_to_r_array(
    const CoordsMatrixType<KEnRef_Real> &coord_array,
    const std::vector<std::tuple<int, int> > &atomId_pairs, int numOmpThreads) {
    CoordsMatrixType<KEnRef_Real> ret(atomId_pairs.size(), 3);
#pragma omp parallel for num_threads(numOmpThreads)
    for (int i = 0; i < atomId_pairs.size(); ++i) {
        auto [atom0, atom1] = atomId_pairs.at(i);
        ret.row(i) = coord_array.row(atom1) - coord_array.row(atom0);
    }
    return ret;
}

//Back-propagate energy derivative from r array to coordinates
template<typename KEnRef_Real>
std::vector<CoordsMatrixType<KEnRef_Real> >
KEnRef<KEnRef_Real>::coord_array_to_r_array_backprop(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    const std::vector<std::tuple<int, int> > &atomId_pairs,
    const std::vector<CoordsMatrixType<KEnRef_Real>> &d_energy_d_r_array,
    int numOmpThreads) {

    const auto num_interactions = atomId_pairs.size();
    const auto num_ensembleMembers = coord_array.size();
    const auto num_atoms = coord_array[0].rows();
    std::vector<CoordsMatrixType<KEnRef_Real> > gradients(num_ensembleMembers);
#pragma omp parallel for num_threads(numOmpThreads)
    for (int i = 0; i < num_ensembleMembers; i++) {
        gradients.at(i) = CoordsMatrixType<KEnRef_Real>::Zero(num_atoms, 3);
    }
    // propagate the internuclear vector derivatives back onto the atomic coordinates
    // An atom may appear in many pairs, so different (p) iterations can target the same gradients[m]
    // row concurrently. `#pragma omp atomic` is only valid on scalar lvalues, so each xyz component is
    // updated individually. This keeps the full (pairs x models) parallelism.
    // NOTE: atomic accumulation order is non-deterministic => results may differ at the ULP level
    // between runs / thread counts. For bitwise-reproducible output, disable parallelism: call with
    // numOmpThreads == 1 (or run with OMP_NUM_THREADS=1).
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
    for (int p = 0; p < num_interactions; ++p) {
        // seq_len(dim(d_energy_d_r_array)[1])
        for (int m = 0; m < num_ensembleMembers; m++) {
            const auto [atomId0, atomId1] = atomId_pairs[p];
            const auto &pair_grad = d_energy_d_r_array[m].row(p);
            for (int d = 0; d < 3; ++d) {
#pragma omp atomic
                gradients[m](atomId0, d) -= pair_grad(d);
#pragma omp atomic
                gradients[m](atomId1, d) += pair_grad(d);
            }
        }
    }
    //        std::cout << "gradients" << std::endl;
    //        for(int m = 0; m < num_models; m++){
    //            std::cout << "model " << m << " first 100 rows" << std::endl;
    //            std::cout << gradients[m].topRows(100) << std::endl;
    //        }
    return gradients;
}

//TODO this method is a duplicate of IoUtils::atomNamePairs_2_atomIdPairs(). Delete one of them
template<typename KEnRef_Real>
std::vector<std::tuple<int, int> >
KEnRef<KEnRef_Real>::atomNamePairs_2_atomIdPairs(
    const std::vector<std::tuple<std::string, std::string> > &atomName_pairs,
    const std::map<std::string, int> &atomNames_2_atomIds, int numOmpThreads) {
    auto atomId_pairs = std::vector<std::tuple<int, int> >(atomName_pairs.size());
    // Fill the vector using atomNames_2_atomIds
#pragma omp parallel for num_threads(numOmpThreads)
    for (int i = 0; i < atomName_pairs.size(); ++i) {
        auto [left, right] = atomName_pairs.at(i);
        // I use at() instead of operator[] to force an exception to be thrown
        atomId_pairs.at(i) = std::move(std::tuple{atomNames_2_atomIds.at(left), atomNames_2_atomIds.at(right)});
    }
    return atomId_pairs;
}

template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::optional<std::vector<CoordsMatrixType<KEnRef_Real> > > >
KEnRef<KEnRef_Real>::coord_array_to_g_energy(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
    const std::vector<std::tuple<std::string, std::string> > &atomName_pairs,
    // Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
    const std::vector<std::vector<std::vector<int> > > &grouping_list,
    // list of lists of integer vectors giving groupings of models to average interaction tensors
    const Eigen::MatrixX<KEnRef_Real> &g0,
    const std::map<std::string, int> &atomNames_2_atomIds,
    KEnRef_Real k, KEnRef_Real n, const bool gradient, const int numOmpThreads) {
    //	std::cout << "coord_array_to_energy(atomName_pairs_) called" << std::endl;
    const auto &atomId_pairs = atomNamePairs_2_atomIdPairs(atomName_pairs, atomNames_2_atomIds);
    return KEnRef::coord_array_to_g_energy(coord_array, atomId_pairs, grouping_list, g0, k, n, gradient, numOmpThreads);
}

template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::optional<std::vector<CoordsMatrixType<KEnRef_Real> > > >
KEnRef<KEnRef_Real>::coord_array_to_g_energy(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
    const std::vector<std::tuple<int, int> > &atomId_pairs,
    // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
    const std::vector<std::vector<std::vector<int> > > &grouping_list,
    // list of lists of integer vectors giving groupings of models to average interaction tensors
    const Eigen::MatrixX<KEnRef_Real> &g0, KEnRef_Real k, KEnRef_Real n, bool gradient,
    int numOmpThreads) {
#if VERBOSE
    Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "");
    std::cout << "========>\n";
    std::cout << "g0 0\tg0 1\n";
    std::cout << g0.format(fmt) << "\n" << std::endl;
#endif
    // calculate inter nuclear vectors
    const auto &r_arrays = coord_array_to_r_array(coord_array, atomId_pairs, numOmpThreads);

    // calculate dipole-dipole interaction tensors [and their derivatives]
    // NOTE: named locals (not a structured binding) — the OpenMP regions below reference these, and
    // clang does not support capturing a structured binding inside an OpenMP region.
    auto d_arrays_tuple = r_array_to_d_array(r_arrays, false, gradient, false, numOmpThreads);
    const auto &d_arrays = std::get<0>(d_arrays_tuple);
    const auto &d_arrays_grad = std::get<1>(d_arrays_tuple);

    // calculate norm squared for different groupings of dipole-dipole interaction tensors
    auto g_tuple = d_array_to_g_multiple_groupings(d_arrays, grouping_list, gradient, numOmpThreads);
    const auto &g_list = std::get<0>(g_tuple);
    const auto &g_list_grad = std::get<1>(g_tuple);

    const auto &g_matrix = vectorOfVectors_to_Matrix(g_list, numOmpThreads);
#if VERBOSE
    std::cout << "========>\n";
    std::cout << "g_matrix 0\tg_matrix 1\tg_matrix Z\n" << g_matrix.format(fmt) << "\n" << std::endl;
#endif
    // calculate energies from the norm squared values
    // (named locals, not a structured binding — see note above re: OpenMP capture)
    auto energy_tuple = power_scaled_loss_function(g_matrix, g0, k, n, gradient, numOmpThreads);
    const auto &energy_matrix = std::get<0>(energy_tuple);
    const auto &energy_matrix_grad = std::get<1>(energy_tuple);
    //	std::cout << "energy_matrix_grad" << std::endl << *energy_matrix_grad << std::endl;

    // return the sum of all the individual restraint energies
    KEnRef_Real sum = energy_matrix.sum();

    //Add derivatives using the chain rule: de/dr = de/dd  * dd/dr = de/dg * dg/dd * dd/dr
    if (gradient) {
        const auto num_pairs = atomId_pairs.size();
        const auto num_models = coord_array.size();
        const auto num_atoms = coord_array[0].rows();
        // First calculate de/dd = de/dg * dg/dd for all individual interaction tensor components
        std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > d_energy_d_d_vector(num_models);
        //<num_models(num_pairs, 5)>
#pragma omp parallel for num_threads(numOmpThreads)
        for (int i = 0; i < num_models; i++)
            d_energy_d_d_vector.at(i) = std::move(
                Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(static_cast<int>(num_pairs), 5));

        // de/dd[j] = sum over groupings i of ( energy_matrix_grad[:,i] (x) g_list_grad[i][j] ).
        // Parallelize the two FREE (non-reduction) axes - model j and pair-block - via collapse(2), keeping
        // the grouping reduction i as the serial inner loop. Distinct (j, pair-block) iterations write
        // disjoint row ranges of distinct d_energy_d_d_vector[j], so this is race-free with NO internal
        // accumulators. It scales whether num_models is small (pair-blocks supply the parallelism) or large
        // (the models do); per-block Eigen ops retain SIMD over the pairs dimension. The ragged grouping
        // membership is handled by the inner `j < g_list_grad[i].size()` guard (collapse stays rectangular).
        constexpr Eigen::Index DEDD_BLOCK = 128;
        const int nPairBlocks = static_cast<int>((static_cast<Eigen::Index>(num_pairs) + DEDD_BLOCK - 1) / DEDD_BLOCK);
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
        for (int j = 0; j < num_models; j++) {
            for (int b = 0; b < nPairBlocks; b++) {
                const Eigen::Index r0 = static_cast<Eigen::Index>(b) * DEDD_BLOCK;
                const Eigen::Index len = std::min<Eigen::Index>(DEDD_BLOCK, static_cast<Eigen::Index>(num_pairs) - r0);
                for (int i = 0; i < g_list.size(); i++) {
                    if (j < static_cast<int>(g_list_grad[i].size())) {
                        d_energy_d_d_vector[j].middleRows(r0, len).array() +=
                            energy_matrix_grad->col(i).segment(r0, len).rowwise().template replicate<5>().array() *
                            g_list_grad[i][j].middleRows(r0, len).array();
                    }
                }
            }
        }

        // Then calculate de/dr = de/dd * dd/dr for each xyz component of the internuclear vectors
        std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 3> > d_energy_d_r_array(num_models);
#pragma omp parallel for num_threads(numOmpThreads)
        for (int i = 0; i < num_models; i++)
            d_energy_d_r_array.at(i) = std::move(Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 3>{num_pairs, 3});
        // sum the individual interaction tensor component derivatives associated with x, y, and z
#pragma omp parallel for collapse(3) num_threads(numOmpThreads)
        for (int i = 0; i < num_models; i++) {
            for (int j = 0; j < 3; j++) {
                for (int p = 0; p < num_pairs; ++p) {
                    d_energy_d_r_array[i](p, j) =
                            (d_arrays_grad[i].row(p).array() *
                             d_energy_d_d_vector[i].row(p).template replicate<3, 1>().reshaped(1, 15).array())
                            (Eigen::seq(j, Eigen::fix<14>, Eigen::fix<3>)).sum();
                }
                //This line works, but it is slower than a loop one row at a time + OMP. You can delete it.
                // d_energy_d_r_array[i].col(j) = (d_arrays_grad[i].array() * d_energy_d_d_vector[i].replicate(Eigen::fix<3>, num_pairs).reshaped(num_pairs, Eigen::fix<15>).array()) (Eigen::placeholders::all, Eigen::seq(j, Eigen::fix<14>, Eigen::fix<3>)).rowwise().sum();
            }
            //std::cout << "d_energy_d_r_array[" << i<< "]\n" << d_energy_d_r_array[i] << std::endl << std::endl;
        }

        std::vector<CoordsMatrixType<KEnRef_Real> > gradients(num_models);
#pragma omp parallel for num_threads(numOmpThreads)
        for (int i = 0; i < num_models; i++) {
            gradients.at(i) = CoordsMatrixType<KEnRef_Real>::Zero(num_atoms, 3);
        }
        // propagate the internuclear vector derivatives back onto the atomic coordinates
        // Scalar atomics per xyz component (an atom is shared across pairs); keeps (pairs x models)
        // parallelism. Accumulation order is non-deterministic — see the note in
        // coord_array_to_r_array_backprop; use numOmpThreads == 1 for bitwise-reproducible output.
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
        for (int p = 0; p < num_pairs; ++p) {
            // seq_len(dim(d_energy_d_r_array)[1])
            for (int m = 0; m < num_models; m++) {
                const auto [atomId0, atomId1] = atomId_pairs[p];
                const auto &pair_grad = d_energy_d_r_array[m].row(p);
                for (int d = 0; d < 3; ++d) {
#pragma omp atomic
                    gradients[m](atomId0, d) -= pair_grad(d);
#pragma omp atomic
                    gradients[m](atomId1, d) += pair_grad(d);
                }
            }
        }

        //        std::cout << "gradients" << std::endl;
        //        for (int m = 0; m < num_models; m++){
        //            std::cout << "model " << m << " first 100 rows" << std::endl;
        //            std::cout << gradients[m].topRows(100) << std::endl;
        //        }
        return {sum, {gradients}};
    }
    return {sum, std::nullopt};
}


template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::vector<CoordsMatrixType<KEnRef_Real> > >
KEnRef<KEnRef_Real>::coord_array_to_g_energy_refactored(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
    const std::vector<std::tuple<std::string, std::string> > &atomName_pairs,
    // Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
    const std::vector<std::vector<std::vector<int> > > &grouping_list,
    // list of lists of integer vectors giving groupings of models to average interaction tensors
    const Eigen::MatrixX<KEnRef_Real> &g0,
    const std::map<std::string, int> &atomNames_2_atomIds,
    KEnRef_Real k, KEnRef_Real n, const bool gradient, const int numOmpThreads) {

    //	std::cout << "coord_array_to_energy(atomName_pairs_) called" << std::endl;
    const auto &atomId_pairs = atomNamePairs_2_atomIdPairs(atomName_pairs, atomNames_2_atomIds);
    return KEnRef::coord_array_to_g_energy_refactored(coord_array, atomId_pairs, grouping_list, g0, k, n, gradient, numOmpThreads);
}

template<typename KEnRef_Real>
NamedRowVector<KEnRef_Real>
KEnRef<KEnRef_Real>::calculateLambdaVector(const SpecDenDataBase<KEnRef_Real> &currentSpecDenData,
                                           const NamedRowVector<KEnRef_Real> &rates, int numOmpThreads) {
    // lambda_vector <- -colSums(rates[rownames(spec_den_data_list[[i]][["lambda_coef"]])]*spec_den_data_list[[i]][["lambda_coef"]])
    const auto &lambda_coef_matrix = currentSpecDenData.get_lambda_coef();
    const auto &rowNames = currentSpecDenData.get_lambda_coef().rowNames();

    // -colSums(rates[row_names,] * lambda_coef) is a vector-matrix product:
    //   lambda_vector = -( rates_row . lambda_coef )   where rates_row(j) = rates(rowNames[j]).
    // Building rates_row is a tiny serial gather (one entry per rate constant); the product itself is a
    // GEMV that Eigen runs (and parallelizes for large operands) race-free and deterministically, so it
    // keeps scaling if lambda_coef grows.
    Eigen::RowVectorX<KEnRef_Real> rates_row(static_cast<Eigen::Index>(rowNames.size()));
    for (size_t j = 0; j < rowNames.size(); ++j) {
        rates_row(static_cast<Eigen::Index>(j)) = rates(rowNames[j]);
    }
    NamedRowVector<KEnRef_Real> lambda_vector(-(rates_row * lambda_coef_matrix));
    lambda_vector.setColNames(lambda_coef_matrix.colNames());
    return lambda_vector;
}

template<typename KEnRef_Real>
std::tuple<NamedRowVector<KEnRef_Real>, std::optional<NamedMatrix<KEnRef_Real>>>
KEnRef<KEnRef_Real>::a_matrix_to_sigma_horizontal(
    const NamedMatrix<KEnRef_Real>& a_matrix,
    const NamedRowVector<KEnRef_Real>& lambda_prime_vec,
    KEnRef_Real proton_mhz,
    bool gradient){

    // K^2 = ((1.2566370614e-6 kilogram meters / (ampere^2 second^2))/(4*pi)*(1.054571726e-34 meter^2 kilograms / second)*(267.513e6 radian s^-1 T^-1)^2)^2
    constexpr auto K_sq = static_cast<KEnRef_Real>(5.69549944e-49);  //rad m^6 / s^2
    const auto omega0 = proton_mhz * static_cast<KEnRef_Real>(1e6) * static_cast<KEnRef_Real>(2 * M_PI);

    // Transpose a_matrix to match R's t(a_matrix)
    NamedMatrix<KEnRef_Real> a_matrix_T = a_matrix.transpose();

    // J0 <- -2 * colSums(a_matrix/lambda_prime_vec)
    Eigen::RowVectorX<KEnRef_Real> J0 = -2.0 * (a_matrix_T.array().colwise() / lambda_prime_vec.transpose().array()).colwise().sum();

    // Jomega < - -2 * colSums(a_matrix * lambda_prime_vec / (lambda_prime_vec ^ 2 + (omega0) ^ 2))
    // J2omega <- -2 * colSums(a_matrix * lambda_prime_vec / (lambda_prime_vec ^ 2 + (2*omega0)^2))
    Eigen::ArrayX<KEnRef_Real> denom = lambda_prime_vec.transpose().array().square() + std::pow(2.0 * omega0, 2);
    Eigen::RowVectorX<KEnRef_Real> J2omega = -2.0 * ((a_matrix_T.array().colwise() * lambda_prime_vec.transpose().array()).colwise() / denom).colwise().sum();

    // sigma = 0.1 * K_sq * (3 * J2omega - 0.5 * J0)
    NamedRowVector<KEnRef_Real> sigma = 0.1 * K_sq * (3.0 * J2omega - 0.5 * J0);
    if (a_matrix.hasRowNames())
        sigma.setNames(a_matrix_T.colNames());

    std::optional<NamedMatrix<KEnRef_Real>> gradient_matrix;

    if (gradient) {
        // Derivatives
        // d_J0_d_a_matrix <- -2/lambda_prime_vec
        Eigen::ArrayXX<KEnRef_Real> d_J0_d_a = -2.0 / lambda_prime_vec.transpose().array();
        // d_J2omega_d_a_matrix <- -2*lambda_prime_vec/(lambda_prime_vec^2 + (2*omega0)^2)
        Eigen::ArrayXX<KEnRef_Real> d_J2_d_a = -2.0 * lambda_prime_vec.transpose().array() / denom.array();
        // d_sigma_d_a_matrix <- 0.1*K_sq*(3*d_J2omega_d_a_matrix - 0.5*d_J0_d_a_matrix)
        Eigen::ArrayXX<KEnRef_Real> d_sigma_d_a = 0.1 * K_sq * (3.0 * d_J2_d_a - 0.5 * d_J0_d_a);

        // Transpose to match R's byrow=TRUE
        NamedMatrix<KEnRef_Real> grad(d_sigma_d_a.matrix().transpose().replicate(a_matrix.rows(), 1));

        if (a_matrix.hasRowNames())
            grad.setRowNames(a_matrix.rowNames());

        if (a_matrix.hasColNames())
            grad.setColNames(a_matrix.colNames());
        else if (lambda_prime_vec.hasColNames())
            grad.setColNames(lambda_prime_vec.colNames());

        gradient_matrix = grad;
    }
    return {sigma, gradient_matrix};
}

template<typename KEnRef_Real>
std::tuple<NamedVector<KEnRef_Real>, std::optional<NamedMatrix<KEnRef_Real>>>
KEnRef<KEnRef_Real>::a_matrix_to_sigma(
    const NamedMatrix<KEnRef_Real>& a_matrix,
    const NamedRowVector<KEnRef_Real>& lambda_prime_vec,
    KEnRef_Real proton_mhz,
    bool gradient){

    // K^2 = ((1.2566370614e-6 kilogram meters / (ampere^2 second^2))/(4*pi)*(1.054571726e-34 meter^2 kilograms / second)*(267.513e6 radian s^-1 T^-1)^2)^2
    constexpr auto K_sq = static_cast<KEnRef_Real>(5.69549944e-49);  //rad m^6 / s^2
    const auto omega0 = proton_mhz * static_cast<KEnRef_Real>(1e6) * static_cast<KEnRef_Real>(2 * M_PI);

    const auto& lambda_prime_vec_arr = lambda_prime_vec.array();
    // J0 <- -2 * colSums(a_matrix/lambda_prime_vec)
    const auto& J0 = -2.0 * (a_matrix.array().rowwise() / lambda_prime_vec_arr).rowwise().sum();

    // Jomega < - -2 * colSums(a_matrix * lambda_prime_vec / (lambda_prime_vec ^ 2 + (omega0) ^ 2))
    // J2omega <- -2 * colSums(a_matrix * lambda_prime_vec / (lambda_prime_vec ^ 2 + (2*omega0)^2))
    const auto& denom = lambda_prime_vec_arr.square() + std::pow(2.0 * omega0, 2);
    const auto& J2omega = -2.0 * ((a_matrix.array().rowwise() * lambda_prime_vec_arr).rowwise() / denom).rowwise().sum();

    // value <- 0.1*K_sq*(3*J2omega - 0.5*J0)
    NamedVector<KEnRef_Real> sigma = 0.1 * K_sq * (3.0 * J2omega - 0.5 * J0);
    if (a_matrix.hasRowNames())
        sigma.setNames(a_matrix.rowNames());

    std::optional<NamedMatrix<KEnRef_Real>> grad;

    if (gradient) {
        // Derivatives
        // d_J0_d_a_matrix <- -2/lambda_prime_vec
        const auto &d_J0_d_a = -2.0 / lambda_prime_vec_arr;
        // d_J2omega_d_a_matrix <- -2*lambda_prime_vec/(lambda_prime_vec^2 + (2*omega0)^2)
        const auto &d_J2_d_a = -2.0 * lambda_prime_vec_arr / denom.array();
        // d_sigma_d_a_matrix <- 0.1*K_sq*(3*d_J2omega_d_a_matrix - 0.5*d_J0_d_a_matrix)
        const auto &d_sigma_d_a = 0.1 * K_sq * (3.0 * d_J2_d_a - 0.5 * d_J0_d_a);

        grad = d_sigma_d_a.replicate(a_matrix.rows(), 1).matrix();

        if (a_matrix.hasRowNames())
            grad->setRowNames(a_matrix.rowNames());
        if (a_matrix.hasColNames())
            grad->setColNames(a_matrix.colNames());
        else if (lambda_prime_vec.hasColNames())
            grad->setColNames(lambda_prime_vec.colNames());
    }
    return {sigma, grad};
}



template<typename KEnRef_Real>
std::tuple<std::vector<NamedVector<KEnRef_Real>>, std::optional<std::vector<Eigen::MatrixX<KEnRef_Real> > > >
KEnRef<KEnRef_Real>::coord_array_to_sigma(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    // Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
    const NamedRowVector<KEnRef_Real> &rates,
    const std::vector<SpecDenData<KEnRef_Real> > &spec_den_data_list,
    KEnRef_Real proton_mhz,
    std::map<std::string, int> atomNames_2_atomIds, bool gradient, int numOmpThreads) {

    if (!rates.hasColNames()) {
        throw std::runtime_error("rates has no column names");
    }

    std::vector<NamedVector<KEnRef_Real>> sigma_vector_list(spec_den_data_list.size());
    std::optional<std::vector<Eigen::MatrixX<KEnRef_Real> > > sigma_grad_vector_list_optional = std::nullopt;
    if (gradient)
        sigma_grad_vector_list_optional = {std::vector<Eigen::MatrixX<KEnRef_Real>>(spec_den_data_list.size())};
    // std::vector<Eigen::MatrixX<KEnRef_Real>> sigma_vector_list(spec_den_data_list.size());

    int numModels = coord_array.size();

    for (int i = 0; i < spec_den_data_list.size(); ++i) {

        std::vector<CoordsMatrixType<KEnRef_Real> > coord_array_meter(numModels);
        for (int m = 0; m < numModels; ++m)
            coord_array_meter[m] = coord_array[m] * 1e-10;

        auto &currentSpecDenData = spec_den_data_list[i];
        const auto& cache = currentSpecDenData.get_atomIdPairs_to_sub0Atom_id_pairs_cache();
        const std::vector<std::tuple<int, int>> &atom_id_pairs =
            cache.has_value() ? cache.value() : atomNamePairs_2_atomIdPairs(currentSpecDenData.get_atom_pairs(), atomNames_2_atomIds);

        // calculate inter nuclear vectors
        const auto &r_arrays = coord_array_to_r_array( coord_array_meter, atom_id_pairs, numOmpThreads);

        // calculate dipole-dipole interaction tensors [and their derivatives]
        auto &&[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, false, gradient, false, numOmpThreads);
        // d_arrays_grad_list.at(i) = std::move(d_arrays_grad);

        const auto &grouping_matrix = IoUtils::subset_idx_to_grouping_mat(currentSpecDenData.get_multiple_grouping());
        //calculate the factor by which the number of models should be expanded
        uint nshift = grouping_matrix.cols() / numModels;
        const auto &d_arrays_shifted = array_shift(d_arrays, nshift);

        //calculate norm squared for different groupings of dipole-dipole interaction tensors
        //g_matrix_list[[i]] <- d_array_to_g_matrix(d_array_shifted, spec_den_data_list[[i]][["groupings"]], gradient=gradient)
        auto &&[g_list, g_matrix_grad] =
                d_array_to_g_matrix(d_arrays_shifted, grouping_matrix, gradient, numOmpThreads);
        const auto &g_matrix = vectorOfVectors_to_Matrix(g_list, numOmpThreads);
        // g_matrix_grad_list.at(i) = std::move(g_matrix_grad);

        //calculate matrix of a values
        const auto &a_matrix = g_matrix_to_a_matrix(g_matrix, currentSpecDenData.get_a_coef(), numOmpThreads);

        // calculate lambda eigenvalues
        const auto &lambda_vector = calculateLambdaVector(currentSpecDenData, rates, numOmpThreads);

        // update eigenvalues to account for molecular tumbling
        // NamedRowVector<KEnRef_Real> lambda_prime_vec = (lambda_vector.array() - rates(KEnRef<KEnRef_Real>::KC)).matrix();
        NamedRowVector<KEnRef_Real> lambda_prime_vec(lambda_vector);
        lambda_prime_vec.array() -= rates(KC);

        auto &&[sigma, sigma_grad_optional] = a_matrix_to_sigma(a_matrix, lambda_prime_vec, proton_mhz, gradient);
        sigma_vector_list.at(i) = std::move(sigma);

#if VERBOSE
        std::cout << "sigma\n" << sigma_vector_list[i] << std::endl;
        std::cout << "sigma_grad_optional\n";
        if (sigma_grad_optional.has_value()) {
            std::cout << sigma_grad_optional.value() << std::endl;
        }else {
            std::cout << "NOTHING" << std::endl;
        }
#endif

        if (sigma_grad_optional) {
            // const auto& sigma_matrix = vectorOfVectors_to_Matrix({sigma_grad_optional}, numOmpThreads);
            sigma_grad_vector_list_optional.value().at(i) = std::move(sigma_grad_optional.value());
        }
    }
    return {std::move(sigma_vector_list), std::move(sigma_grad_vector_list_optional)};
}


template<typename KEnRef_Real>
Eigen::MatrixX<KEnRef_Real>
KEnRef<KEnRef_Real>::g_matrix_to_a_matrix_backprop(
    const Eigen::MatrixX<KEnRef_Real> &a_coef,
    const Eigen::MatrixX<KEnRef_Real> &d_energy_d_a_matrix,
    int numOmpThreads) {

    // Backpropagation: d_energy/d_g(:,i) = sum_j d_energy/d_a(:,j) * a_coef(i,j) == d_energy_d_a_matrix * a_coef^T.
    // This is the transpose-GEMM mirror of the forward g_matrix_to_a_matrix (a_matrix = g_matrix * a_coef):
    // Eigen runs it parallel for large operands, SIMD-vectorized and deterministically, replacing the manual
    // per-column accumulation. The old `if (a_coef(i,j) != 0)` sparsity guard is unnecessary — exact-zero
    // coefficients contribute exactly zero to the product, so the result is unchanged.
    return d_energy_d_a_matrix * a_coef.transpose();
}

template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::optional<std::vector<CoordsMatrixType<KEnRef_Real> > > >
KEnRef<KEnRef_Real>::coord_array_to_sigma_energy(
    std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    const NamedRowVector<KEnRef_Real> &rates,
    const std::vector<SpecDenData<KEnRef_Real> > &spec_den_data_list,
    KEnRef_Real proton_mhz,
    KEnRef_Real k, KEnRef_Real n,
    const std::map<std::string, int>& atomNames_2_atomIds, bool gradient, int numOmpThreads /*, lossFunction lossFunc*/) {

    if (!rates.hasColNames()) {
        throw std::runtime_error("rates has no column names");
    }

    auto specDenDataSize = spec_den_data_list.size();
    std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > > d_arrays_grad_list(specDenDataSize);
    std::vector<std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > > > g_matrix_grad_list(specDenDataSize);
    std::vector<NamedVector<KEnRef_Real>> sigma_vector_list(specDenDataSize);
    std::vector<Eigen::MatrixX<KEnRef_Real> > sigma_grad_vector_list;
    if (gradient) {
        sigma_grad_vector_list = std::vector<Eigen::MatrixX<KEnRef_Real> >(specDenDataSize);
    }
    const int numModels = coord_array.size();
    for (int m = 0; m < numModels; ++m)
        coord_array[m] *= 1e-10;    // no need for coord_array_meter

    for (int i = 0; i < specDenDataSize; ++i) {
        auto &currentSpecDenData = spec_den_data_list[i];
        //Cache AtomIdPairs not's calculate them every time
        const auto& cache = currentSpecDenData.get_atomIdPairs_to_sub0Atom_id_pairs_cache();
        const std::vector<std::tuple<int, int>> &atom_id_pairs =
            cache.has_value() ? cache.value() : atomNamePairs_2_atomIdPairs(currentSpecDenData.get_atom_pairs(), atomNames_2_atomIds);
        // calculate inter nuclear vectors
        const auto &r_arrays =
            coord_array_to_r_array(coord_array, atom_id_pairs, numOmpThreads);

        // calculate dipole-dipole interaction tensors [and their derivatives]
        auto &&[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, false, gradient, false, numOmpThreads);
        d_arrays_grad_list.at(i) = std::move(d_arrays_grad);

        //TODO get rid of the grouping_matrix. find nShift someway else.
        const auto &grouping_matrix = IoUtils::subset_idx_to_grouping_mat(currentSpecDenData.get_multiple_grouping());
        //calculate the factor by which the number of models should be expanded
        uint nShift = grouping_matrix.cols() / numModels;
        const auto &d_arrays_shifted = array_shift(d_arrays, nShift);

        //calculate norm squared for different groupings of dipole-dipole interaction tensors
        //g_matrix_list[[i]] <- d_array_to_g_matrix(d_array_shifted, spec_den_data_list[[i]][["groupings"]], gradient=gradient)
        auto &&[g_list, g_matrix_grad] =
                d_array_to_g_matrix(d_arrays_shifted, grouping_matrix, gradient, numOmpThreads);
        const auto &g_matrix = vectorOfVectors_to_Matrix(g_list, numOmpThreads);
        g_matrix_grad_list.at(i) = std::move(g_matrix_grad);

        //calculate matrix of a values
        const auto &a_matrix = g_matrix_to_a_matrix(g_matrix, currentSpecDenData.get_a_coef(), numOmpThreads);

        // calculate lambda eigenvalues
        const auto &lambda_vector = calculateLambdaVector(currentSpecDenData, rates, numOmpThreads);

        // update eigenvalues to account for molecular tumbling
        // NamedRowVector<KEnRef_Real> lambda_prime_vec = (lambda_vector.array() - rates(KEnRef<KEnRef_Real>::KC)).matrix();
        NamedRowVector<KEnRef_Real> lambda_prime_vec(lambda_vector);
        lambda_prime_vec.array() -= rates(KC);

        // if (lossFunc == POWER_SCALED_LOSS_FUNCTION) //NO need anymore
        auto &&[sigma, sigma_grad_optional] = a_matrix_to_sigma(a_matrix, lambda_prime_vec, proton_mhz, gradient);
        sigma_vector_list.at(i) = std::move(sigma);
        if (gradient) {
            sigma_grad_vector_list.at(i) = std::move(sigma_grad_optional.value());
        }

#if VERBOSE
        std::cout << "sigma\n" << sigma_vector_list[i] << std::endl;
        std::cout << "sigma_grad_optional\n";
        if (sigma_grad_optional.has_value()) {
            std::cout << sigma_grad_optional.value() << std::endl;
        }else {
            std::cout << "NOTHING" << std::endl;
        }
#endif
    }

    //unlist sigma_vector_list
    int totalNumSigmas=0;
    for (int i = 0; i < sigma_vector_list.size(); ++i)
        totalNumSigmas += sigma_vector_list.at(i).size();
    Eigen::MatrixX<KEnRef_Real> sigma(totalNumSigmas, 1);

    std::vector<int> indexes(specDenDataSize);
    for (int i = 0, index = 0; i < specDenDataSize; ++i) {
        indexes.at(i) = index;
        const auto &sigmas = spec_den_data_list.at(i).sigmas();
        sigma(Eigen::seqN(index, sigmas.rows()), Eigen::placeholders::all) = sigma_vector_list[i];
        index += sigmas.rows();
    }

    //unlist spec_den_data_list sigmas
    int totalSigma0Values=0;
    for (int i = 0; i < specDenDataSize; ++i) {
        const auto &specDenData = spec_den_data_list.at(i);
        const auto& sigmas = specDenData.sigmas();

        const int numAtomPairs = specDenData.get_atom_pairs().size();
        //TODO this is not the proper way to calculate cols.
        // Groupings should be a matrix and cols should be its columns.
        const int cols = specDenData.get_multiple_grouping().at(0).at(0).size();
        uint nShift = cols / numModels;
        assert(sigmas.rows() == numAtomPairs / nShift);

        totalSigma0Values += sigmas.rows();
    }
    assert(totalNumSigmas == totalSigma0Values);

    Eigen::MatrixX<KEnRef_Real> sigma0(totalSigma0Values, 1);
    for (int i = 0, index = 0; i < specDenDataSize; ++i) {
        const auto &sigma0s = spec_den_data_list.at(i).sigmas();
        sigma0(Eigen::seqN(index, sigma0s.rows()), Eigen::placeholders::all) = sigma0s;
        index += sigma0s.rows();
    }

    const auto&[sigma_energy_vec, sigma_energy_grad_vec] = power_scaled_loss_function(sigma, sigma0, k, n, gradient, numOmpThreads);
    KEnRef_Real ret1 = sigma_energy_vec.sum();
    std::optional<std::vector<CoordsMatrixType<KEnRef_Real> >> d_energy_d_coord_array;

    if (gradient) {
        d_energy_d_coord_array = std::vector<CoordsMatrixType<KEnRef_Real> >(numModels);
        for (int i = 0; i < numModels; ++i) {
            d_energy_d_coord_array->at(i) = std::move(CoordsMatrixType<KEnRef_Real>::Zero(coord_array[0 /*i*/].rows(), 3));
        }

        for (int i = 0 ; i < specDenDataSize; ++i) {
            const auto &currentSpecDenData = spec_den_data_list[i];

            // extract particular segment of energy gradient
            const auto &d_energy_d_sigma = sigma_energy_grad_vec.value()(Eigen::seqN(indexes[i],currentSpecDenData.sigmas().rows()), Eigen::placeholders::all);

            const auto & d_energy_d_a_matrix = a_matrix_to_sigma_backprop(sigma_grad_vector_list.at(i), d_energy_d_sigma);
            // std::cout << "d_energy_d_a_matrix("<<d_energy_d_a_matrix.rows()<<", "<<d_energy_d_a_matrix.cols()<<")\n";
            // std::cout << d_energy_d_a_matrix << std::endl;

            const auto &d_energy_d_g_matrix = g_matrix_to_a_matrix_backprop(currentSpecDenData.get_a_coef(), d_energy_d_a_matrix);
            // std::cout << "d_energy_d_g_matrix("<<d_energy_d_g_matrix.rows()<<", "<<d_energy_d_g_matrix.cols()<<")\n";
            // std::cout << d_energy_d_g_matrix << std::endl;

            // back-propagate derivatives from g matrix to d array
            const auto &d_energy_d_d_array = d_array_to_g_matrix_backprop(g_matrix_grad_list.at(i), d_energy_d_g_matrix, 0);
            // std::cout << "d_energy_d_d_array Iteration #"<<i<<" length = "<< d_energy_d_d_array.size() << std::endl;
            // for (int j = 0; j < d_energy_d_d_array.size(); ++j) {
            //     std::cout << "d_energy_d_d_array["<< j <<"] ("<< d_energy_d_d_array.at(j).rows() <<", " << d_energy_d_d_array.at(j).cols() << ")"<< std::endl;
            //     std::cout << d_energy_d_d_array.at(j) << std::endl;
            // }

            // // //No need to shift back. We will handle the shifted arrays internally in the next function.
            // //shift tensor components from virtual models back to atom pairs
            // n_shift <- ncol(spec_den_data_list[[i]][["groupings"]])/dim(coord_array)[3]
            //if (n_shift > 1) {
            //	d_energy_d_d_array <- array_shift(d_energy_d_d_array, dnames=dimnames(d_array_list[[i]]))
            //}

            // back-propagate derivatives from d array to r array
            const auto &d_energy_d_r_array = r_array_to_d_array_backprop(d_arrays_grad_list[i], d_energy_d_d_array, 0);
            // std::cout << "d_energy_d_r_array  Iteration #" << i << " length = "<< d_energy_d_r_array.size() << std::endl;
            // for (int j = 0; j < d_energy_d_r_array.size(); ++j) {
            //     std::cout << "d_energy_d_r_array["<<j<< "] Dim =("<< d_energy_d_r_array.at(j).rows()<<", " <<d_energy_d_r_array.at(j).cols() << ")\n";
            //     std::cout << d_energy_d_r_array.at(j) << std::endl;
            // }

            const auto& cache = currentSpecDenData.get_atomIdPairs_to_sub0Atom_id_pairs_cache();
            const std::vector<std::tuple<int, int>> &atom_id_pairs =
                cache.has_value() ? cache.value() : atomNamePairs_2_atomIdPairs(currentSpecDenData.get_atom_pairs(), atomNames_2_atomIds);
            // accumulate back-propagated derivatives from r array into coord array gradient
            const auto & gradients =
                coord_array_to_r_array_backprop(coord_array, atom_id_pairs, d_energy_d_r_array,0);
            for (int j = 0; j < numModels; ++j)
                d_energy_d_coord_array->at(j) += gradients.at(j);
        }
        for (int i = 0; i < d_energy_d_coord_array->size(); ++i)
            d_energy_d_coord_array->at(i) *= 1e-10;
    }
    return {ret1, d_energy_d_coord_array};
}

template<typename KEnRef_Real>
Eigen::MatrixX<KEnRef_Real>
KEnRef<KEnRef_Real>::a_matrix_to_sigma_backprop(
    const Eigen::MatrixX<KEnRef_Real> &d_sigma_d_a_matrix,
    const Eigen::MatrixX<KEnRef_Real> &d_energy_d_sigma,
    int numOmpThreads) {
    return 	d_sigma_d_a_matrix.array() * d_energy_d_sigma.array().replicate(1,d_sigma_d_a_matrix.cols());
}



template<typename KEnRef_Real>
std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >
KEnRef<KEnRef_Real>::d_array_to_g_matrix_backprop(
    const std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > > &d_g_matrix_d_d_array,
    const Eigen::MatrixX<KEnRef_Real> &d_energy_d_g_matrix,
    /*const int num_interactions, const int num_ensembleMembers,*/ int numOmpThreads) {

    //num_pairs is a multiple of num_interactions. e.g. in case of 1-3 interaction, num_pairs = 3*num_interactions
    const int num_something = d_g_matrix_d_d_array.size(); //TODO what exactly is it? is it num_groupings?
    const int num_pairs = d_energy_d_g_matrix.rows(); // or d_g_matrix_d_d_array[0][0].rows() ??
    const int num_ensembleMembers = d_g_matrix_d_d_array[0].size();

    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > d_energy_d_d_vector(num_ensembleMembers);
#pragma omp parallel for num_threads(numOmpThreads)
    for (int i = 0; i < num_ensembleMembers; i++)
        d_energy_d_d_vector.at(i) = std::move(
            Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(static_cast<int>(num_pairs), 5));

    // d_energy_d_d_vector[j] accumulates across the grouping index i (the reduction axis). Parallelize the
    // two FREE axes - model j and pair-block - via collapse(2), keeping i serial inside. Distinct
    // (j, pair-block) iterations write disjoint row ranges of distinct d_energy_d_d_vector[j], so this is
    // race-free with NO internal accumulators, and scales whether num_ensembleMembers is small (pair-blocks
    // supply the parallelism) or large (the models do); per-block Eigen ops keep SIMD over the pairs axis.
    constexpr Eigen::Index DEDD_BLOCK = 128;
    const int nPairBlocks = static_cast<int>((static_cast<Eigen::Index>(num_pairs) + DEDD_BLOCK - 1) / DEDD_BLOCK);
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
    for (int j = 0; j < num_ensembleMembers; j++) {
        for (int b = 0; b < nPairBlocks; b++) {
            const Eigen::Index r0 = static_cast<Eigen::Index>(b) * DEDD_BLOCK;
            const Eigen::Index len = std::min<Eigen::Index>(DEDD_BLOCK, static_cast<Eigen::Index>(num_pairs) - r0);
            for (int i = 0; i < num_something; i++) {
                d_energy_d_d_vector[j].middleRows(r0, len).array() +=
                    d_energy_d_g_matrix.col(i).segment(r0, len).rowwise().template replicate<5>().array() *
                    d_g_matrix_d_d_array[i][j].middleRows(r0, len).array();
            }
        }
    }
    return d_energy_d_d_vector;
}

template<typename KEnRef_Real>
std::vector<CoordsMatrixType<KEnRef_Real> >
KEnRef<KEnRef_Real>::r_array_to_d_array_backprop(
    const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > &d_d_array_d_r_array,
    const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &d_energy_d_d_vector,
    int numOmpThreads) {

    const int num_ensembleMembers = d_energy_d_d_vector.size();
    const int num_models = d_d_array_d_r_array.size();
    const int num_pairs = d_d_array_d_r_array[0].rows();
    const int num_interactions = d_energy_d_d_vector[0].rows();
    const int shiftFactor = num_ensembleMembers / num_models;

    std::vector<CoordsMatrixType<KEnRef_Real> > d_energy_d_r_array(num_models);
#pragma omp parallel for num_threads(numOmpThreads)
    for (int i = 0; i < num_models; i++)
        d_energy_d_r_array.at(i) = std::move(Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 3>{num_pairs, 3});

    // sum the individual interaction tensor component derivatives associated with x, y, and z
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
    for (int m = 0; m < num_models; m++) {
        for (int p = 0; p < num_pairs; ++p) {
            const int modelShift = p / num_interactions;
            const int rowShift = p % num_interactions;

            const auto& dd_dr_row = d_d_array_d_r_array[m].row(p);
            const auto& de_dd_row = d_energy_d_d_vector[m * shiftFactor + modelShift].row(rowShift);

            for (int xyzIdx = 0; xyzIdx < 3; xyzIdx++) {
                // Directly multiply and sum the relevant 5 elements
                d_energy_d_r_array.at(m)(p, xyzIdx) =
                    (dd_dr_row(Eigen::seqN(xyzIdx, Eigen::fix<5>, Eigen::fix<3>)).array() * de_dd_row.array()).sum();
            }
        }
    }
    return d_energy_d_r_array;
}

template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::vector<CoordsMatrixType<KEnRef_Real> > >
KEnRef<KEnRef_Real>::coord_array_to_g_energy_refactored(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
    const std::vector<std::tuple<int, int> > &atomId_pairs,
    // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
    const std::vector<std::vector<std::vector<int> > > &grouping_list,
    // list of lists of integer vectors giving groupings of models to average interaction tensors
    const Eigen::MatrixX<KEnRef_Real> &g0, KEnRef_Real k, KEnRef_Real n, bool gradient,
    int numOmpThreads) {

    // calculate inter nuclear vectors
    const auto &r_arrays = coord_array_to_r_array(coord_array, atomId_pairs, numOmpThreads);

    // calculate dipole-dipole interaction tensors [and their derivatives]
    // Claude says that structured binding in an OMP region is dangerous. so Take care if you put the next line in an OMP region.
    const auto &[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, false, gradient, false, numOmpThreads);

    const auto &grouping_matrix = IoUtils::subset_idx_to_grouping_mat(grouping_list);
    // calculate norm squared for different groupings of dipole-dipole interaction tensors
    const auto &[g_list, g_list_grad] = d_array_to_g_matrix(d_arrays, grouping_matrix, gradient, numOmpThreads);

    const auto &g_matrix = vectorOfVectors_to_Matrix(g_list, numOmpThreads);
    // calculate energies from the norm squared values
    const auto &[energy_matrix, energy_matrix_grad] = power_scaled_loss_function(g_matrix, g0, k, n, gradient, numOmpThreads);

    // return the sum of all the individual restraint energies
    KEnRef_Real sum = energy_matrix.sum();

    //Add derivatives using the chain rule: de/dr = de/dd  * dd/dr = de/dg * dg/dd * dd/dr
    if (gradient) {
        const auto num_interactions = atomId_pairs.size();
        const auto num_ensembleMembers = coord_array.size();
        // First calculate de/dd = de/dg * dg/dd for all individual interaction tensor components
        const auto &d_energy_d_d_vector = d_array_to_g_matrix_backprop(g_list_grad, *energy_matrix_grad,
                                                                       /*num_interactions, num_ensembleMembers,*/
                                                                       numOmpThreads);

        // Then calculate de/dr = de/dd  * dd/dr for each xyz component of the internuclear vectors
        const auto &d_energy_d_r_array = r_array_to_d_array_backprop(d_arrays_grad, d_energy_d_d_vector,
                                                                     // num_interactions, num_ensembleMembers,
                                                                     numOmpThreads);

        const auto &gradients = coord_array_to_r_array_backprop(coord_array, atomId_pairs, d_energy_d_r_array,
                                                                numOmpThreads);
        return {sum, gradients};
    }
    return {sum, std::vector<CoordsMatrixType<KEnRef_Real> >{}};
}


template<typename KEnRef_Real>
Eigen::MatrixX<KEnRef_Real>
KEnRef<KEnRef_Real>::coord_array_to_g_matrix(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
    const std::vector<std::tuple<int, int> > &atomId_pairs,
    // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
    const std::vector<std::vector<std::vector<int> > > &grouping_list,
    // list of lists of integer vectors giving groupings of models to average interaction tensors
    int numOmpThreads) {
    //	std::cout << "coord_array_to_g_matrix() called" << std::endl;
    // calculate internuclear vectors
    const auto &r_arrays = coord_array_to_r_array(coord_array, atomId_pairs, numOmpThreads);

    // calculate dipole-dipole interaction tensors [and their derivatives]
    auto [d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, false, false, false, numOmpThreads);

    // calculate norm squared for different groupings of dipole-dipole interaction tensors
    //		g_list <- lapply(grouping_list, function(grouping) d_array_to_g(d_array, grouping, gradient=FALSE))
    auto [g_list, ignore] = d_array_to_g_multiple_groupings(d_arrays, grouping_list, numOmpThreads);

    return vectorOfVectors_to_Matrix(g_list, numOmpThreads);
}

/**
 * @brief Computes a_matrix = g_matrix * a_coef (only for non-zero a_coef entries)
 * @param g_matrix Input matrix (rows: observations, cols: groupings)
 * @param a_coef Coefficient matrix (rows: groupings, cols: output features)
 * @param numOmpThreads Number of OpenMP threads to use (if <=0, use default)
 * @return a_matrix Resulting matrix (rows: observations, cols: output features)
 */
template<typename KEnRef_Real>
NamedMatrix<KEnRef_Real>
KEnRef<KEnRef_Real>::g_matrix_to_a_matrix(
    const Eigen::MatrixX<KEnRef_Real> &g_matrix,
    const NamedMatrix<KEnRef_Real> &a_coef,
    int numOmpThreads) {
    // Initialize output matrix with zeros
    // NamedRowVector<KEnRef_Real> a_matrix = Eigen::MatrixX<KEnRef_Real>::Zero(g_matrix.rows(), a_coef.cols());
    // a_matrix(:,i) = sum_j g_matrix(:,j) * a_coef(j,i)  ==  g_matrix * a_coef : a plain matrix product.
    // Expressing it as a GEMM keeps the whole computation parallel (Eigen, internally threaded for large
    // operands) and deterministic, and removes the manual accumulation that raced under collapse(2).
    NamedMatrix<KEnRef_Real> a_matrix(g_matrix * a_coef);
    // a_matrix.setRowNames(g_matrix.rowNames()); //TODO check it
    a_matrix.setColNames(a_coef.colNames()); //TODO not the best way. But it is the best available (yet).
    return a_matrix;
}

//TODO Shall we move it (with its test of course) to GMXKenRefForceProvider?
template<typename KEnRef_Real>
void
KEnRef<KEnRef_Real>::saturate(CoordsMatrixType<KEnRef_Real> &derivatives_rectified, KEnRef_Real thresholdSquared,
                              int numOmpThreads) {
#pragma omp parallel for num_threads(numOmpThreads) default(shared)
    for (int i = 0; i < derivatives_rectified.rows(); i++) {
        auto scaleDown = derivatives_rectified.row(i).squaredNorm() / thresholdSquared;
        if (scaleDown > 1.0) {
            derivatives_rectified.row(i) /= std::sqrt(scaleDown);
        }
    }
}

template<typename KEnRef_Real>
Eigen::VectorX<KEnRef_Real>
KEnRef<KEnRef_Real>::s2OrderParams(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
    const std::vector<std::tuple<int, int> > &atomId_pairs,
    // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
    int numOmpThreads) {
    int numModels = coord_array.size();

    //    //calculate internuclear vectors
    //    const auto &r_arrays = coord_array_to_r_array(coord_array, atomId_pairs, numOmpThreads);
    //
    //    // calculate dipole-dipole interaction tensors [and their derivatives]
    //    const auto &[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, numOmpThreads);

    ////////////////////////////////////////////////////////////////////////////////

    //calculate array of internuclear vectors *group_r_array*
    // group_r_array <- ke::coord_array_to_r_array(aperm(group_mean_coord, c(2,1,3)), group_pairs)
    //const auto &group_r_array = KEnRef<KEnRef_Real>::coord_array_to_r_array(group_mean_coord, atomId_pairs, numOmpThreads);
    // As long as there are no "groups" (i.e. all are singltons), we can safely use *coord_array* for now.
    const auto &group_r_array = coord_array_to_r_array(coord_array, atomId_pairs, numOmpThreads);

    //calculate matrix of radii
    // group_r_mat <- sqrt(rowSums(group_r_array^2, dims=2))
    std::vector<Eigen::VectorX<KEnRef_Real> > group_r_mat;
    group_r_mat.reserve(numModels);
    for (int i = 0; i < numModels; ++i) {
        //        group_r_mat.at(i) = std::move(Eigen::VectorX<KEnRef_Real>(group_r_array[i].rows()));
        group_r_mat.emplace_back(Eigen::VectorX<KEnRef_Real>(group_r_array[i].rows()));
    }
    for (int m = 0; m < numModels; ++m) {
        for (int i = 0; i < group_r_array[m].rows(); ++i) {
            group_r_mat[m](i) = group_r_array[m].row(i).norm();
        }
    }

    //calculate order parameter using normalized array of internuclear vectors
    // group_rnorm_array <- group_r_array/as.vector(group_r_mat)
    // group_dnorm_array <- ke::r_array_to_d_array(group_rnorm_array)
    // group_s2 <- rowSums(colMeans(aperm(group_dnorm_array, c(2,1,3)))^2)
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> group_dnorm_array_tempMeans = Eigen::Matrix<KEnRef_Real,
        Eigen::Dynamic, 5>::Zero(group_r_mat[0].rows(), 5);
    for (int i = 0; i < numModels; ++i) {
        const auto &group_rnorm_array = group_r_array[i].array() / group_r_mat[i].rowwise().template replicate<3>().
                                        array();
        const auto &group_dnorm_array_1model = r_array_to_d_array(group_rnorm_array, false, false, false, numOmpThreads);
        group_dnorm_array_tempMeans += std::get<0>(group_dnorm_array_1model);
    }
    group_dnorm_array_tempMeans /= numModels;
    //    Eigen::VectorX<KEnRef_Real> group_s2 = group_dnorm_array_tempMeans.array().square().colwise().sum();
    Eigen::VectorX<KEnRef_Real> group_s2 = Eigen::VectorX<KEnRef_Real>::Zero(group_dnorm_array_tempMeans.rows());
    group_s2 = group_dnorm_array_tempMeans.array().square().rowwise().sum();

    //################################################

    //TODO Alternate approach directly from d_array
    // group_d_array <- ke::r_array_to_d_array(group_r_array)
    // group_dnorm_array_alt <- group_d_array/as.vector(sqrt(rowSums(group_d_array^2, dims=2)))
    // group_s2_alt <- rowSums(colMeans(aperm(group_dnorm_array_alt, c(2,1,3)))^2)
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > group_d_array =
            std::get<0>(r_array_to_d_array(group_r_array, false, false, false, numOmpThreads));

    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > group_dnorm_array_alt;
    group_dnorm_array_alt.reserve(numModels);
    for (int i = 0; i < numModels; ++i) {
        Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> temp = group_d_array[i].array().rowwise().norm().template
                replicate<1, 5>();
        group_dnorm_array_alt.emplace_back(group_d_array[i].array() / temp.array());
    }

    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> group_dnorm_array_alt_tempMeans = Eigen::Matrix<KEnRef_Real,
        Eigen::Dynamic, 5>::Zero(group_dnorm_array_alt[0].rows(), 5);
    for (int i = 0; i < numModels; ++i) {
        group_dnorm_array_alt_tempMeans.array() += group_dnorm_array_alt[i].array();
    }
    group_dnorm_array_alt_tempMeans.array() /= numModels;
    Eigen::VectorX<KEnRef_Real> group_s2_alt = group_dnorm_array_alt_tempMeans.array().square().colwise().sum();

    //    assert((group_s2 - group_s2_alt).array() < 1e7); //FIXME didn't work
    //#######################################################
    return group_s2;
}

template<typename KEnRef_Real>
std::tuple<Eigen::VectorX<KEnRef_Real>, std::optional<Eigen::MatrixX<KEnRef_Real> > >
KEnRef<KEnRef_Real>::a_matrix_to_relax(
    const Eigen::MatrixX<KEnRef_Real> &a_int_matrix,
    const Eigen::RowVectorX<KEnRef_Real> &lambda_int_vec,
    const Eigen::MatrixX<KEnRef_Real> &a_overall_matrix,
    const Eigen::RowVectorX<KEnRef_Real> &lambda_overall_vec,
    const SpecDenTermArray<KEnRef_Real> &spec_den_term_array,
    const bool gradient, const int numOmpThreads) {

    const Eigen::Index pairs = a_int_matrix.rows();
    const Eigen::Index n_int = a_int_matrix.cols();
    const Eigen::Index n_overall = a_overall_matrix.cols();

    assert(lambda_int_vec.size() == n_int);
    assert(lambda_overall_vec.size() == n_overall);
    assert(a_overall_matrix.rows() == pairs);
    assert(spec_den_term_array.coef.rows() == pairs);
    assert(spec_den_term_array.freq.rows() == spec_den_term_array.coef.rows()
        && spec_den_term_array.freq.cols() == spec_den_term_array.coef.cols());

    Eigen::VectorX<KEnRef_Real> value = Eigen::VectorX<KEnRef_Real>::Zero(pairs);
    std::optional<Eigen::MatrixX<KEnRef_Real> > grad;
    if (gradient) grad = Eigen::MatrixX<KEnRef_Real>::Zero(pairs, n_int);

    // freq^2 is loop-invariant across (i, j); precompute once (read-only, shared across threads).
    const Eigen::Array<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> freq_sq = spec_den_term_array.freq.array().square();

    // Parallelize over `pairs` in contiguous row-blocks. Each block keeps Eigen's SIMD over the
    // (contiguous, ColMajor) pairs dimension, and different blocks write disjoint rows of value/grad,
    // so there is no reduction and the result is bitwise-identical to the serial path. `numOmpThreads`
    // follows the project convention (0 => use as many threads as available). Small problems stay serial:
    // the minimum `pairs` to parallelize is env-configurable (KENREF_RELAX_PARALLEL_THRESHOLD, default 256).
    [[maybe_unused]] static const Eigen::Index PARALLEL_THRESHOLD = []() -> Eigen::Index {
        if (const char *env = std::getenv("KENREF_RELAX_PARALLEL_THRESHOLD")) {
            char *end = nullptr;
            const long v = std::strtol(env, &end, 10);
            if (end != env && *end == '\0' && v >= 0) return static_cast<Eigen::Index>(v);
        }
        return 256;
    }();
    constexpr Eigen::Index BLOCK = 128;   // rows per chunk: SIMD-friendly, cache-resident, load-balanced
    const long nChunks = static_cast<long>((pairs + BLOCK - 1) / BLOCK);

#pragma omp parallel for num_threads(numOmpThreads) if(pairs >= PARALLEL_THRESHOLD) schedule(static)
    for (long c = 0; c < nChunks; ++c) {
        const Eigen::Index r0 = static_cast<Eigen::Index>(c) * BLOCK;
        const Eigen::Index len = std::min<Eigen::Index>(BLOCK, pairs - r0);

        const auto coef_b = spec_den_term_array.coef.middleRows(r0, len).array();   // [len x terms]
        const auto freqsq_b = freq_sq.middleRows(r0, len);                          // [len x terms]
        Eigen::Array<KEnRef_Real, Eigen::Dynamic, 1> term(len);                     // reused across (i, j)

        for (Eigen::Index i = 0; i < n_overall; ++i) {
            const auto a_over_i = a_overall_matrix.col(i).segment(r0, len).array(); // [len]
            for (Eigen::Index j = 0; j < n_int; ++j) {
                const KEnRef_Real lambda_prime = lambda_int_vec(j) + lambda_overall_vec(i);
                // term = -a_overall[block,i] * lambda_prime * rowSum_t( coef / (lambda_prime^2 + freq^2) )
                term = (-a_over_i) * lambda_prime * (coef_b / (freqsq_b + lambda_prime * lambda_prime)).rowwise().sum();
                value.segment(r0, len).array() += a_int_matrix.col(j).segment(r0, len).array() * term;
                if (gradient) grad->col(j).segment(r0, len).array() += term;
            }
        }
    }
    return {std::move(value), std::move(grad)};
}

template<typename KEnRef_Real>
Eigen::MatrixX<KEnRef_Real>
KEnRef<KEnRef_Real>::a_matrix_to_relax_backprop(
    const Eigen::MatrixX<KEnRef_Real> &d_relax_d_a_matrix,
    const Eigen::VectorX<KEnRef_Real> &d_energy_d_relax,
    const int /*numOmpThreads*/) {
    assert(d_relax_d_a_matrix.rows() == d_energy_d_relax.rows());
    // R: d_relax_d_a_matrix * d_energy_d_relax — the length-`pairs` upstream gradient
    // multiplies every column elementwise (R recycles the vector down each column).
    return (d_relax_d_a_matrix.array().colwise() * d_energy_d_relax.array()).matrix();
}

template<typename KEnRef_Real>
std::tuple<Eigen::MatrixX<KEnRef_Real>, Eigen::RowVectorX<KEnRef_Real> >
KEnRef<KEnRef_Real>::dxyz_dunit_to_overall_modes(
    const Eigen::Matrix<KEnRef_Real, 3, 1> &dxyz_vec,
    const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &dunit_a_array_shifted,
    const KEnRef_Real tol, const int /*numOmpThreads*/) {

    assert(!dunit_a_array_shifted.empty());
    const KEnRef_Real dx = dxyz_vec(0);
    const KEnRef_Real dy = dxyz_vec(1);
    const KEnRef_Real dz = dxyz_vec(2);
    const KEnRef_Real sr3 = static_cast<KEnRef_Real>(1.73205080756888);   // sqrt(3), matches R literal

    const Eigen::Index n_pairs = dunit_a_array_shifted[0].rows();
    const size_t n_models = dunit_a_array_shifted.size();
    const KEnRef_Real inv_m = static_cast<KEnRef_Real>(1) / static_cast<KEnRef_Real>(n_models);

    using Arr = Eigen::Array<KEnRef_Real, Eigen::Dynamic, 1>;
    // R: abs(a - b) <= tol * max(1, abs(a), abs(b))
    auto isClose = [tol](const KEnRef_Real a, const KEnRef_Real b) {
        return std::abs(a - b) <= tol * std::max(std::max(static_cast<KEnRef_Real>(1), std::abs(a)), std::abs(b));
    };

    // ── Isotropic: Dx == Dy == Dz ────────────────────────────────────────────
    if (isClose(dx, dy) && isClose(dx, dz)) {
        Arr acc = Arr::Zero(n_pairs);
        for (size_t m = 0; m < n_models; ++m)
            acc += dunit_a_array_shifted[m].array().square().rowwise().sum();
        Eigen::MatrixX<KEnRef_Real> a_overall = (acc * inv_m).matrix();   // n_pairs x 1
        Eigen::RowVectorX<KEnRef_Real> lambda(1);
        lambda(0) = -6 * dx;
        return {std::move(a_overall), std::move(lambda)};
    }

    // ── Axially symmetric: Dx == Dy ──────────────────────────────────────────
    if (isClose(dx, dy)) {
        Arr c0 = Arr::Zero(n_pairs), c1 = Arr::Zero(n_pairs), c2 = Arr::Zero(n_pairs);
        for (size_t m = 0; m < n_models; ++m) {
            const auto &d = dunit_a_array_shifted[m];
            c0 += d.col(0).array().square();
            c1 += d.col(1).array().square() + d.col(4).array().square();
            c2 += d.col(2).array().square() + d.col(3).array().square();
        }
        Eigen::MatrixX<KEnRef_Real> a_overall(n_pairs, 3);
        a_overall.col(0) = (c0 * inv_m).matrix();
        a_overall.col(1) = (c1 * inv_m).matrix();
        a_overall.col(2) = (c2 * inv_m).matrix();
        Eigen::RowVectorX<KEnRef_Real> lambda(3);
        lambda(0) = -3 * (dx + dy);
        lambda(1) = -(2 * dx + 4 * dz);
        lambda(2) = -(5 * dx + dz);
        return {std::move(a_overall), std::move(lambda)};
    }

    // ── Fully anisotropic ────────────────────────────────────────────────────
    const KEnRef_Real block_a = 3 * (dx + dy);
    const KEnRef_Real block_b = sr3 * (dx - dy);
    const KEnRef_Real block_c = dx + dy + 4 * dz;
    const KEnRef_Real block_trace_half = static_cast<KEnRef_Real>(0.5) * (block_a + block_c);
    const KEnRef_Real block_diff_half = static_cast<KEnRef_Real>(0.5) * (block_a - block_c);
    const KEnRef_Real block_radius = std::sqrt(block_diff_half * block_diff_half + block_b * block_b);

    Eigen::RowVectorX<KEnRef_Real> lambda(5);
    lambda(0) = -(block_trace_half + block_radius);
    lambda(1) = -(block_trace_half - block_radius);
    lambda(2) = -(dx + 4 * dy + dz);
    lambda(3) = -(4 * dx + dy + dz);
    lambda(4) = -(dx + dy + 4 * dz);

    const KEnRef_Real theta = static_cast<KEnRef_Real>(0.5) * std::atan2(2 * block_b, block_a - block_c);
    const KEnRef_Real cos_t = std::cos(theta);
    const KEnRef_Real sin_t = std::sin(theta);

    Arr c0 = Arr::Zero(n_pairs), c1 = Arr::Zero(n_pairs), c2 = Arr::Zero(n_pairs),
        c3 = Arr::Zero(n_pairs), c4 = Arr::Zero(n_pairs);
    for (size_t m = 0; m < n_models; ++m) {
        const auto &d = dunit_a_array_shifted[m];
        const Arr mode1 = cos_t * d.col(0).array() + sin_t * d.col(1).array();
        const Arr mode2 = -sin_t * d.col(0).array() + cos_t * d.col(1).array();
        c0 += mode1.square();
        c1 += mode2.square();
        c2 += d.col(2).array().square();
        c3 += d.col(3).array().square();
        c4 += d.col(4).array().square();
    }
    Eigen::MatrixX<KEnRef_Real> a_overall(n_pairs, 5);
    a_overall.col(0) = (c0 * inv_m).matrix();
    a_overall.col(1) = (c1 * inv_m).matrix();
    a_overall.col(2) = (c2 * inv_m).matrix();
    a_overall.col(3) = (c3 * inv_m).matrix();
    a_overall.col(4) = (c4 * inv_m).matrix();
    return {std::move(a_overall), std::move(lambda)};
}

template<typename KEnRef_Real>
std::vector<NamedMatrix<KEnRef_Real> >
KEnRef<KEnRef_Real>::coord_array_to_relax(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    const NamedRowVector<KEnRef_Real> &rates,
    const std::vector<SpecDenRelaxData<KEnRef_Real> > &spec_den_relax_data_list,
    const std::map<std::string, int> &atomNames_2_atomIds,
    const int numOmpThreads) {

    if (!rates.hasColNames())
        throw std::runtime_error("rates has no column names");

    const int numModels = static_cast<int>(coord_array.size());
    Eigen::Matrix<KEnRef_Real, 3, 1> dxyz;
    dxyz << rates("Dx"), rates("Dy"), rates("Dz");

    // convert from Å to m once (shared by all substructures)
    std::vector<CoordsMatrixType<KEnRef_Real> > coord_array_meter(numModels);
    for (int m = 0; m < numModels; ++m)
        coord_array_meter[m] = coord_array[m] * static_cast<KEnRef_Real>(1e-10);

    std::vector<NamedMatrix<KEnRef_Real> > result;
    result.reserve(spec_den_relax_data_list.size());

    for (const auto &relaxData: spec_den_relax_data_list) {
        const bool unit = relaxData.is_unit();
        const auto &cache = relaxData.get_atomIdPairs_to_sub0Atom_id_pairs_cache();
        const std::vector<std::tuple<int, int> > &atom_id_pairs =
            cache.has_value() ? cache.value()
                              : atomNamePairs_2_atomIdPairs(relaxData.get_atom_pairs(), atomNames_2_atomIds);

        // internuclear vectors
        const auto &r_arrays = coord_array_to_r_array(coord_array_meter, atom_id_pairs, numOmpThreads);

        // dipole-dipole interaction tensors (R: dist = !unit)
        auto &&[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, unit, false, false, numOmpThreads);
        (void) d_arrays_grad;

        const auto &grouping_matrix = IoUtils::subset_idx_to_grouping_mat(relaxData.get_multiple_grouping());
        const uint nShift = grouping_matrix.cols() / numModels;
        const auto &d_arrays_shifted = array_shift(d_arrays, nShift);

        // internal-motion amplitudes and eigenvalues
        auto &&[g_list, g_grad] = d_array_to_g_matrix(d_arrays_shifted, grouping_matrix, false, numOmpThreads);
        (void) g_grad;
        const auto &g_matrix = vectorOfVectors_to_Matrix(g_list, numOmpThreads);
        const auto &a_int_matrix = g_matrix_to_a_matrix(g_matrix, relaxData.get_a_coef(), numOmpThreads);
        const auto &lambda_int_vec = calculateLambdaVector(relaxData, rates, numOmpThreads);

        // unit dipole-dipole tensors (reuse d_arrays if it is already the unit form)
        std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > d_unit_arrays;
        if (unit) {
            d_unit_arrays = d_arrays;
        } else {
            auto &&[du, du_grad] = r_array_to_d_array(r_arrays, true, false, false, numOmpThreads);
            (void) du_grad;
            d_unit_arrays = std::move(du);
        }
        const auto &d_unit_shifted = array_shift(d_unit_arrays, nShift);

        // overall rotational-diffusion tumbling modes
        auto &&[a_overall_matrix, lambda_overall_vec] = dxyz_dunit_to_overall_modes(dxyz, d_unit_shifted);

        // one relaxation rate per relax_data_list entry
        const auto &relax_list = relaxData.get_relax_data_list();
        const Eigen::Index n_pairs = a_int_matrix.rows();
        Eigen::MatrixX<KEnRef_Real> relax_eig(n_pairs, static_cast<Eigen::Index>(relax_list.size()));
        std::vector<std::string> colNames;
        colNames.reserve(relax_list.size());
        for (size_t r = 0; r < relax_list.size(); ++r) {
            auto &&[relax_vec, relax_grad] = a_matrix_to_relax(
                a_int_matrix, lambda_int_vec, a_overall_matrix, lambda_overall_vec,
                relax_list[r].spec_den_term_array, false, numOmpThreads);
            (void) relax_grad;
            relax_eig.col(static_cast<Eigen::Index>(r)) = relax_vec;
            colNames.push_back(relax_list[r].rate_name);
        }
        // 3 explicit args select NamedMatrixImpl's (matrix, colNames, rowNames) constructor over the
        // Eigen base constructors pulled in by `using Base::Base` (which have no matching 3-arg form).
        NamedMatrix<KEnRef_Real> relax_mat(relax_eig,
                                           std::optional<std::vector<std::string> >(colNames),
                                           std::nullopt);
        result.push_back(std::move(relax_mat));
    }
    return result;
}

template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::optional<std::vector<CoordsMatrixType<KEnRef_Real> > > >
KEnRef<KEnRef_Real>::coord_array_to_relax_energy(
    std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
    const NamedRowVector<KEnRef_Real> &rates,
    const std::vector<SpecDenRelaxData<KEnRef_Real> > &spec_den_relax_data_list,
    const KEnRef_Real k, const KEnRef_Real n,
    const std::map<std::string, int> &atomNames_2_atomIds, const bool gradient, const int numOmpThreads) {

    if (!rates.hasColNames())
        throw std::runtime_error("rates has no column names");

    const size_t specDenSize = spec_den_relax_data_list.size();
    const int numModels = static_cast<int>(coord_array.size());
    Eigen::Matrix<KEnRef_Real, 3, 1> dxyz;
    dxyz << rates("Dx"), rates("Dy"), rates("Dz");

    // scale coords from Å to m in place (matches coord_array_to_sigma_energy convention)
    for (int m = 0; m < numModels; ++m)
        coord_array[m] *= static_cast<KEnRef_Real>(1e-10);

    // captured for back-propagation
    std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > > d_arrays_grad_list(specDenSize);
    std::vector<std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > > > g_matrix_grad_list(specDenSize);
    std::vector<std::vector<Eigen::MatrixX<KEnRef_Real> > > relax_grad_list(specDenSize); // [i][rate] -> [pairs x n_int]

    // flattened predicted / target / force-constant vectors (order: substructure -> rate -> pair)
    std::vector<KEnRef_Real> relax_flat, relax0_flat, k_flat;

    for (size_t i = 0; i < specDenSize; ++i) {
        const auto &relaxData = spec_den_relax_data_list[i];
        const bool unit = relaxData.is_unit();
        const auto &relax_list = relaxData.get_relax_data_list();

        // R requires a target `value` for every relaxation rate
        for (const auto &entry: relax_list)
            if (!entry.value.has_value())
                throw std::runtime_error(
                    "coord_array_to_relax_energy() requires target `value` for all relaxation data; missing for rate "
                    + entry.rate_name);

        const auto &cache = relaxData.get_atomIdPairs_to_sub0Atom_id_pairs_cache();
        const std::vector<std::tuple<int, int> > &atom_id_pairs =
            cache.has_value() ? cache.value()
                              : atomNamePairs_2_atomIdPairs(relaxData.get_atom_pairs(), atomNames_2_atomIds);

        const auto &r_arrays = coord_array_to_r_array(coord_array, atom_id_pairs, numOmpThreads);

        auto &&[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, unit, gradient, false, numOmpThreads);
        d_arrays_grad_list.at(i) = std::move(d_arrays_grad);

        const auto &grouping_matrix = IoUtils::subset_idx_to_grouping_mat(relaxData.get_multiple_grouping());
        const uint nShift = grouping_matrix.cols() / numModels;
        const auto &d_arrays_shifted = array_shift(d_arrays, nShift);

        auto &&[g_list, g_matrix_grad] = d_array_to_g_matrix(d_arrays_shifted, grouping_matrix, gradient, numOmpThreads);
        const auto &g_matrix = vectorOfVectors_to_Matrix(g_list, numOmpThreads);
        g_matrix_grad_list.at(i) = std::move(g_matrix_grad);

        const auto &a_int_matrix = g_matrix_to_a_matrix(g_matrix, relaxData.get_a_coef(), numOmpThreads);
        const auto &lambda_int_vec = calculateLambdaVector(relaxData, rates, numOmpThreads);

        // unit dipole-dipole tensors -> overall tumbling modes
        std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > d_unit_arrays;
        if (unit) {
            d_unit_arrays = d_arrays;
        } else {
            auto &&[du, du_grad] = r_array_to_d_array(r_arrays, true, false, false, numOmpThreads);
            (void) du_grad;
            d_unit_arrays = std::move(du);
        }
        const auto &d_unit_shifted = array_shift(d_unit_arrays, nShift);
        auto &&[a_overall_matrix, lambda_overall_vec] = dxyz_dunit_to_overall_modes(dxyz, d_unit_shifted);

        relax_grad_list.at(i).resize(relax_list.size());
        for (size_t r = 0; r < relax_list.size(); ++r) {
            const auto &entry = relax_list[r];
            auto &&[relax_vec, relax_grad] = a_matrix_to_relax(
                a_int_matrix, lambda_int_vec, a_overall_matrix, lambda_overall_vec,
                entry.spec_den_term_array, gradient, numOmpThreads);

            const auto &value = entry.value.value();
            const Eigen::Index n_value = value.rows();
            assert(relax_vec.rows() == n_value);

            for (Eigen::Index p = 0; p < n_value; ++p) {
                relax_flat.push_back(relax_vec(p));
                relax0_flat.push_back(value(p, 0));
            }
            // per-rate k: length 1 broadcasts, else length n_value; default 1
            if (!entry.k.has_value()) {
                for (Eigen::Index p = 0; p < n_value; ++p) k_flat.push_back(static_cast<KEnRef_Real>(1));
            } else if (entry.k->rows() == 1) {
                const KEnRef_Real kv = (*entry.k)(0, 0);
                for (Eigen::Index p = 0; p < n_value; ++p) k_flat.push_back(kv);
            } else {
                assert(entry.k->rows() == n_value);
                for (Eigen::Index p = 0; p < n_value; ++p) k_flat.push_back((*entry.k)(p, 0));
            }

            if (gradient) relax_grad_list.at(i).at(r) = std::move(relax_grad.value());
        }
    }

    const Eigen::Index total = static_cast<Eigen::Index>(relax_flat.size());
    Eigen::MatrixX<KEnRef_Real> relax = Eigen::Map<Eigen::VectorX<KEnRef_Real> >(relax_flat.data(), total);
    Eigen::MatrixX<KEnRef_Real> relax0 = Eigen::Map<Eigen::VectorX<KEnRef_Real> >(relax0_flat.data(), total);
    Eigen::Array<KEnRef_Real, Eigen::Dynamic, 1> k_vec =
        Eigen::Map<Eigen::Array<KEnRef_Real, Eigen::Dynamic, 1> >(k_flat.data(), total);

    // The loss is linear in its force constant, so evaluate with k=1 and scale by the
    // per-rate k_vec (× the caller's scalar k) elementwise — matches R's vector-k loss.
    const auto &[base_energy, base_grad] = power_scaled_loss_function(
        relax, relax0, static_cast<KEnRef_Real>(1), n, gradient, numOmpThreads);

    const KEnRef_Real value = (k * k_vec * base_energy.array().col(0)).sum();

    std::optional<std::vector<CoordsMatrixType<KEnRef_Real> > > d_energy_d_coord_array;
    if (gradient) {
        // elementwise upstream gradient d(energy)/d(relax), scaled by k * k_vec
        const Eigen::VectorX<KEnRef_Real> d_energy_d_relax_all =
            (k * k_vec * base_grad.value().array().col(0)).matrix();

        d_energy_d_coord_array = std::vector<CoordsMatrixType<KEnRef_Real> >(numModels);
        for (int m = 0; m < numModels; ++m)
            d_energy_d_coord_array->at(m) = CoordsMatrixType<KEnRef_Real>::Zero(coord_array[0].rows(), 3);

        Eigen::Index offset = 0;
        for (size_t i = 0; i < specDenSize; ++i) {
            const auto &relaxData = spec_den_relax_data_list[i];
            const auto &relax_list = relaxData.get_relax_data_list();
            const Eigen::Index n_int = relaxData.get_a_coef().cols();
            // a_int_matrix has n_data_rows rows; recover it from any stored relax gradient
            const Eigen::Index n_pairs = relax_grad_list.at(i).at(0).rows();

            Eigen::MatrixX<KEnRef_Real> d_energy_d_a_int = Eigen::MatrixX<KEnRef_Real>::Zero(n_pairs, n_int);
            for (size_t r = 0; r < relax_list.size(); ++r) {
                const Eigen::VectorX<KEnRef_Real> d_energy_d_relax = d_energy_d_relax_all.segment(offset, n_pairs);
                offset += n_pairs;
                d_energy_d_a_int += a_matrix_to_relax_backprop(relax_grad_list.at(i).at(r), d_energy_d_relax);
            }

            const auto &d_energy_d_g_matrix = g_matrix_to_a_matrix_backprop(relaxData.get_a_coef(), d_energy_d_a_int);
            const auto &d_energy_d_d_array = d_array_to_g_matrix_backprop(g_matrix_grad_list.at(i), d_energy_d_g_matrix, 0);
            const auto &d_energy_d_r_array = r_array_to_d_array_backprop(d_arrays_grad_list[i], d_energy_d_d_array, 0);

            const auto &cache = relaxData.get_atomIdPairs_to_sub0Atom_id_pairs_cache();
            const std::vector<std::tuple<int, int> > &atom_id_pairs =
                cache.has_value() ? cache.value()
                                  : atomNamePairs_2_atomIdPairs(relaxData.get_atom_pairs(), atomNames_2_atomIds);
            const auto &gradients = coord_array_to_r_array_backprop(coord_array, atom_id_pairs, d_energy_d_r_array, 0);
            for (int m = 0; m < numModels; ++m)
                d_energy_d_coord_array->at(m) += gradients.at(m);
        }
        for (int m = 0; m < numModels; ++m)
            d_energy_d_coord_array->at(m) *= static_cast<KEnRef_Real>(1e-10);
    }
    return {value, std::move(d_energy_d_coord_array)};
}

template class SpecDenDataBase<float>;
template class SpecDenDataBase<double>;
template class SpecDenData<float>;
template class SpecDenData<double>;
template class SpecDenRelaxData<float>;
template class SpecDenRelaxData<double>;
template class KEnRef<float>;
template class KEnRef<double>;
#undef VERBOSE
