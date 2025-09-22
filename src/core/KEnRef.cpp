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
// #include <Eigen/Dense>
#include "core/KEnRef.h"
//#include <iostream>//for testing only

#include "core/IoUtils.h"
//#include <utility>

template<typename KEnRef_Real>
KEnRef_Real SpecDenData<KEnRef_Real>::getSigma(const std::tuple<std::string, std::string>& atomPairs) const {
    // const size_t pairId = atomNamePairs_to_atomPairIndex.at(atomPairs);
    // return getSigma(static_cast<const size_t>(pairId % sigma->size()));
    return getSigma(atomNamePairs_to_atomPairIndex.at(atomPairs) % sigma->size());
}

template<typename KEnRef_Real>
void SpecDenData<KEnRef_Real>::setAtomPairs(std::vector<std::tuple<std::string, std::string>> &atomPairs){
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
                // std::cout << "To ret[" <<(modelNo*nShift+interactionComponent) <<"]("<<interactionId<<", Eigen::all)\t";
                // std::cout << "from arrays_to_shift["<<modelNo<<"]("<<(interactionId + interactionComponent * numInteractions)<<", Eigen::all)" << std::endl;
                ret[modelNo * nShift + interactionComponent](interactionId, Eigen::all) =
                        arrays_to_shift[modelNo](interactionId + interactionComponent * numInteractions, Eigen::all);
                //TODO is there a way other than copying with "="?
            }
        }
    }
    return ret;
}

template<typename KEnRef_Real>
std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> >
KEnRef<KEnRef_Real>::r_array_to_d_array(const CoordsMatrixType<KEnRef_Real> &Nxyz, bool gradient, bool addEpsilon,
                                        int numOmpThreads) {
    //	std::cout << "r_array_to_d_array(Nxyz) called" << std::endl;
    auto N = Nxyz.rows();

    auto x = Nxyz.col(0).array();
    auto y = Nxyz.col(1).array();
    auto z = Nxyz.col(2).array();
    auto sqrt3 = static_cast<KEnRef_Real>(sqrt(3.0));

    enum ci {
        x2, y2, z2,
        //xyz,
        x2_y2_z2,
        xy, xz, yz, x2_minusy2,
        x2_y2_z2_p52, x2_y2_z2_p72,
        //sqrt3_x2_y2_z2_p52, sqrt3_x2_y2_z2_p72,
        half_minusx2_minusy2__z2, // NOLINT(*-reserved-identifier)
        sqrt3_over_x2_y2_z2_p52, neg5_over_x2_y2_z2_p72, neg5sqrt3_over_x2_y2_z2_p72, neg5sqrt3_over_2_x2_y2_z2_p72
    };

    //TODO use multiplying by inverse() instead of division
    Eigen::ArrayXX<KEnRef_Real> cache(N, 15); //must be Array (not Matrix) to allow item-wise operations
#define CACHE(a) (cache.col(a))
    CACHE(x2) = x.square();
    CACHE(y2) = y.square();
    CACHE(z2) = z.square();
    CACHE(xy) = x * y;
    CACHE(xz) = x * z;
    CACHE(yz) = y * z;
    //CACHE(xyz) 			= CACHE(xy) * z;
    CACHE(x2_y2_z2) = CACHE(x2) + CACHE(y2) + CACHE(z2); //x2 + y2 + z2;
    CACHE(x2_minusy2) = CACHE(x2) - CACHE(y2);
    // CACHE(x2_y2_z2_p52) = CACHE(x2_y2_z2).pow(5).sqrt() + std::numeric_limits<KEnRef_Real>::epsilon();
    CACHE(x2_y2_z2_p52) = CACHE(x2_y2_z2).pow(5).sqrt();
    if (addEpsilon)
        CACHE(x2_y2_z2_p52) += std::numeric_limits<KEnRef_Real>::epsilon();
    CACHE(half_minusx2_minusy2__z2) = ((-CACHE(x2) - CACHE(y2)) / 2) + CACHE(z2);
    // std::cout << "sqrt(2.38294e-20) " <<(sqrt(2.38294e-20)) << ",  sqrt(2.38294e-20))^5 " << pow(sqrt(2.38294e-20), 5)<< " ==> + epsilon("<< std::numeric_limits<KEnRef_Real>::epsilon() <<") " << (pow(sqrt(2.38294e-20), 5) + std::numeric_limits<KEnRef_Real>::epsilon() ) <<std::endl;
    // std::cout << "cache" << cache << std:: endl;

    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> ret1(N, 5);
    ret1.col(0) = CACHE(half_minusx2_minusy2__z2);
    ret1.col(1) = sqrt3 * CACHE(x2_minusy2) / 2;
    ret1.col(2) = sqrt3 * CACHE(xz);
    ret1.col(3) = sqrt3 * CACHE(yz);
    ret1.col(4) = sqrt3 * CACHE(xy);
    //	std::cout << "ret1 before division" << std::endl << ret1 << std::endl;
    //	std::cout << "x2_y2_z2" << std::endl << CACHE(x2_y2_z2) << std::endl;
    //	std::cout << "x2_y2_z2 power 5/2" << std::endl << CACHE(x2_y2_z2_p52).rowwise().template replicate<5>() << std::endl;

    //	std::cout << "ret1.array()" << ret1.array() << std::endl /*<< "CACHE(x2_y2_z2_p52)" << CACHE(x2_y2_z2_p52)*/ << std::endl << "CACHE(x2_y2_z2_p52).rowwise().replicate<5>()" << std::endl << CACHE(x2_y2_z2_p52).rowwise().replicate<5>() << std::endl;
    ret1.array() /= CACHE(x2_y2_z2_p52).rowwise().template replicate<5>()
            /* + std::numeric_limits<KEnRef_Real>::epsilon()*/;

    //    std:: cout << "ret1 after division" << std::endl << ret1 << std::endl;
    if (!gradient) {
        return {ret1, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15>{}};
    }
    //    std::cout << "CACHE after d ret1 " << cache << std::endl;

    CACHE(x2_y2_z2_p72) = CACHE(x2_y2_z2).pow(7).sqrt();
    if (addEpsilon)
        CACHE(x2_y2_z2_p72) += std::numeric_limits<KEnRef_Real_t>::epsilon();
    //	CACHE(sqrt3_x2_y2_z2_p52) = sqrt3 * CACHE(x2_y2_z2_p52);
    //	CACHE(sqrt3_x2_y2_z2_p72) = sqrt3 * CACHE(x2_y2_z2_p72);
    CACHE(sqrt3_over_x2_y2_z2_p52) = sqrt3 / CACHE(x2_y2_z2_p52);
    CACHE(neg5_over_x2_y2_z2_p72) = -5.0 / CACHE(x2_y2_z2_p72);
    CACHE(neg5sqrt3_over_x2_y2_z2_p72) = -5 * sqrt3 / CACHE(x2_y2_z2_p72);
    CACHE(neg5sqrt3_over_2_x2_y2_z2_p72) = CACHE(neg5sqrt3_over_x2_y2_z2_p72) / 2;
    //	CACHE(sqrt3_over_x2_y2_z2_p52)			=  sqrt3 / CACHE(x2_y2_z2).pow(5).sqrt();
    //	CACHE(_neg10_over_sqrt3_x2_y2_z2_p72) 	= -10 / (sqrt3 * CACHE(x2_y2_z2).pow(7).sqrt());
    //    std::cout << "CACHE after completion " << cache << std::endl;

    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> ret2 = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic,
        Eigen::Dynamic>(N, 15);

    ret2.col(0) = (x * CACHE(half_minusx2_minusy2__z2) * CACHE(neg5_over_x2_y2_z2_p72)) - (x / CACHE(x2_y2_z2_p52));
    ret2.col(1) = (y * CACHE(half_minusx2_minusy2__z2) * CACHE(neg5_over_x2_y2_z2_p72)) - (y / CACHE(x2_y2_z2_p52));
    ret2.col(2) = (z * CACHE(half_minusx2_minusy2__z2) * CACHE(neg5_over_x2_y2_z2_p72)) + (2 * z / CACHE(x2_y2_z2_p52));

    ret2.col(3) = (x * CACHE(x2_minusy2) * CACHE(neg5sqrt3_over_2_x2_y2_z2_p72)) + (x * CACHE(sqrt3_over_x2_y2_z2_p52));
    ret2.col(4) = (y * CACHE(x2_minusy2) * CACHE(neg5sqrt3_over_2_x2_y2_z2_p72)) - (y * CACHE(sqrt3_over_x2_y2_z2_p52));
    ret2.col(5) = (z * CACHE(x2_minusy2) * CACHE(neg5sqrt3_over_2_x2_y2_z2_p72));

    ret2.col(6) = (x * CACHE(xz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) + (z * CACHE(sqrt3_over_x2_y2_z2_p52));
    ret2.col(7) = (y * CACHE(xz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72));
    ret2.col(8) = (z * CACHE(xz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) + (x * CACHE(sqrt3_over_x2_y2_z2_p52));

    ret2.col(9) = (x * CACHE(yz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72));
    ret2.col(10) = (y * CACHE(yz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) + (z * CACHE(sqrt3_over_x2_y2_z2_p52));
    ret2.col(11) = (z * CACHE(yz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) + (y * CACHE(sqrt3_over_x2_y2_z2_p52));

    ret2.col(12) = (x * CACHE(xy) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) + (y * CACHE(sqrt3_over_x2_y2_z2_p52));
    ret2.col(13) = (y * CACHE(xy) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) + (x * CACHE(sqrt3_over_x2_y2_z2_p52));
    ret2.col(14) = (z * CACHE(xy) * CACHE(neg5sqrt3_over_x2_y2_z2_p72));

#undef CACHE
    //    std:: cout << "d ret2 " << std::endl << ret2 << std::endl;
    return {ret1, ret2};
}

//	std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> KEnRef::r_array_to_d_array(const Eigen::MatrixX3f& Nxyz, bool gradient){
template<typename KEnRef_Real>
std::tuple<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >, std::vector<Eigen::Matrix<KEnRef_Real,
    Eigen::Dynamic, 15> > >
KEnRef<KEnRef_Real>::r_array_to_d_array(const std::vector<CoordsMatrixType<KEnRef_Real> > &models_Nxyz, bool gradient,
                                        bool addEpsilon, int numOmpThreads) {
    //	std::cout << "r_array_to_d_array(models_Nxyz) called" << std::endl;
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > ret1;
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > ret2;
    ret1.reserve(models_Nxyz.size());
    ret2.reserve(models_Nxyz.size());
    for (const auto &Nxyz: models_Nxyz) {
        auto [arr1, arr2] = r_array_to_d_array(Nxyz, gradient, addEpsilon, numOmpThreads);
        ret1.emplace_back(arr1);
        ret2.emplace_back(arr2);
    }
    return {ret1, ret2};
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
    auto num_models = d_arrays.size();
    auto num_interactionIds = d_arrays[0].rows();

    // std::cout << "num_models " 	<< num_models << " num_pairIds " << num_pairIds << "num_groups " << num_groups << std::endl;
    Eigen::VectorX<KEnRef_Real> ret1 = Eigen::VectorX<KEnRef_Real>::Zero(num_interactionIds);

    //Every element of ret2 (every d_matrix_grad) is a matrix(num_pairIds x 5)
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > ret2;
    if (gradient) {
        ret2 = std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >(num_models);
#pragma omp parallel for num_threads(numOmpThreads)
        for (int i = 0; i < num_models; ++i) {
            ret2.at(i) = std::move(Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>(num_interactionIds, 5));
            // AKA d_matrix_grad
        }
    } else {
        ret2 = std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >(0);
    }

    //for every grouping block
    for (const auto &currentGrouping: grouping) {
        //create a new empty d_matrix (filled with 0) to carry "average dipole interaction tensor" every group
        Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> d_matrix = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(
            num_interactionIds, 5);

        const auto currentGroupSize = currentGrouping.size();
        const auto CURRENT_GROUP_SIZE_real = static_cast<KEnRef_Real_t>(currentGroupSize);
        const auto currentGroupSize_OVER_num_models_real = CURRENT_GROUP_SIZE_real / num_models;

        // sum the dipole interaction tensors within each group/block
#pragma omp parallel for num_threads(numOmpThreads) // schedule(static) //num_threads(gmx_omp_nthreads_get(ModuleMultiThread::Default))
        for (int j = 0; j < currentGroupSize; ++j) {
            //for every member of the grouping block
            //sum relevant models into relevant groups (e.g., model 1 & 2 into group 1, and models 3 & 4 into group 2)
            //#pragma omp atomic
            d_matrix += d_arrays[currentGrouping[j]];
        }

        if (gradient) {
            auto TWO_OVER_num_models_currentGroupSize =
                    static_cast<KEnRef_Real_t>(2.0) / num_models / CURRENT_GROUP_SIZE_real;
            // calculate  d_matrix_grad. d_matrix_grad shape is pairIDs * interaction tensor elements
            const auto &d_matrix_grad = d_matrix * TWO_OVER_num_models_currentGroupSize;
            // All models of the same group equally share the same value
#pragma omp parallel for num_threads(numOmpThreads)
            for (int j = 0; j < currentGroupSize; j++) {
                // All elements are equal (i.e., all models get the same overall (average ?) value at the end.
                ret2[currentGrouping[j]] = d_matrix_grad;
            }
        }

        // Divide d_matrix by currentGroupSize to get the average
        d_matrix /= CURRENT_GROUP_SIZE_real; //N.B. Dividing it line by line in OMP, was slower(!)
        //#pragma omp parallel for num_threads(numOmpThreads)
        //        for (int i = 0; i < d_matrix.rows(); i++) {
        //            d_matrix.row(i) /= CURRENT_GROUP_SIZE_real;
        //        }

        // calculate self dot product (norm squared) and accumulate group's contribution to mean g
#pragma omp parallel for num_threads(numOmpThreads)
        for (int j = 0; j < num_interactionIds; j++) {
            const auto contribution = d_matrix.row(j).squaredNorm() * currentGroupSize_OVER_num_models_real;
#pragma omp atomic
            ret1(j) += contribution;
        }
    }
    //	std::cout << "ret1" << std::endl << ret1 << std::endl;
    return {ret1, ret2};
}

template<typename KEnRef_Real>
std::tuple<Eigen::MatrixX<KEnRef_Real>, std::optional<Eigen::MatrixX<KEnRef_Real> > >
KEnRef<KEnRef_Real>::power_scaled_loss_function(
Eigen::MatrixX<KEnRef_Real> g, // current group norm squared values
Eigen::MatrixX<KEnRef_Real> g0, // target group norm squared values.
KEnRef_Real k, // force constant
KEnRef_Real n, // correction power.
const bool gradient,
int numOmpThreads
) {

    const auto &g_arr = g.array() /*+ std::numeric_limits<KEnRef_Real_t>::epsilon()*/;
    const auto &g0_arr = g0.array();
    //    std::cout << "g_arr " << g_arr << std::endl;
    Eigen::MatrixX<KEnRef_Real> ret1;
    Eigen::MatrixX<KEnRef_Real> ret2;

    const auto &common = (Eigen::pow(1 + Eigen::abs(g_arr), n) - 1) * Eigen::sign(g_arr) - (
             Eigen::pow(1 + Eigen::abs(g0_arr), n) - 1) * Eigen::sign(g0_arr);
    //          std::cout << "loss " << loss << std::endl;
    ret1 = k * common.square().matrix();
    //This value may become infinity if it excceds 3.402823466E38 in a single precision float
    if (gradient) {
        ret2 = (2.0 * k * common * (n * Eigen::pow(1 + Eigen::abs(g_arr), n - 1))).matrix();
    }
    return {ret1, ret2};
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
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
    for (int p = 0; p < num_interactions; ++p) {
        // seq_len(dim(d_energy_d_r_array)[1])
        for (int m = 0; m < num_ensembleMembers; m++) {
            const auto [atomId0, atomId1] = atomId_pairs[p];
            const auto &pair_grad = d_energy_d_r_array[m].row(p);
#pragma omp atomic
            gradients[m].row(atomId0) -= pair_grad;
#pragma omp atomic
            gradients[m].row(atomId1) += pair_grad;
        }
    }
    //        std::cout << "gradients" << std::endl;
    //        for(int m = 0; m < num_models; m++){
    //            std::cout << "model " << m << " first 100 rows" << std::endl;
    //            std::cout << gradients[m].topRows(100) << std::endl;
    //        }
    return gradients;
}

template<typename KEnRef_Real>
std::shared_ptr<std::vector<std::tuple<int, int> > >
KEnRef<KEnRef_Real>::atomNamePairs_2_atomIdPairs(
    const std::vector<std::tuple<std::string, std::string> > &atomName_pairs,
    const std::map<std::string, int> &atomNames_2_atomIds, int numOmpThreads) {
    auto atomId_pairs = std::make_shared<std::vector<std::tuple<int, int> > >(atomName_pairs.size());
    // Fill the vector using atomNames_2_atomIds
#pragma omp parallel for num_threads(numOmpThreads)
    for (int i = 0; i < atomName_pairs.size(); ++i) {
        auto [left, right] = atomName_pairs.at(i);
        // I use at() instead of operator[] to force an exception to be thrown
        atomId_pairs->at(i) = std::move(std::tuple{atomNames_2_atomIds.at(left), atomNames_2_atomIds.at(right)});
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
    return KEnRef::coord_array_to_g_energy(coord_array, *atomId_pairs, grouping_list, g0, k, n, gradient, numOmpThreads);
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
#if true
    Eigen::IOFormat fmt(Eigen::FullPrecision, 0, "\t", "\n", "", "");
#endif
#if VERBOSE
    std::cout << "========>\n";
    std::cout << "g0 0\tg0 1\n";
    std::cout << g0.format(fmt) << "\n" << std::endl;
#endif
    // calculate inter nuclear vectors
    const auto &r_arrays = coord_array_to_r_array(coord_array, atomId_pairs, numOmpThreads);

    // calculate dipole-dipole interaction tensors [and their derivatives]
    const auto &[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, gradient, false, numOmpThreads);

    // calculate norm squared for different groupings of dipole-dipole interaction tensors
    const auto &[g_list, g_list_grad] = d_array_to_g_multiple_groupings(d_arrays, grouping_list, gradient,
                                                                        numOmpThreads);

    const auto &g_matrix = vectorOfVectors_to_Matrix(g_list, numOmpThreads);
#if true
    std::cout << "========>\n";
    std::cout << "g_matrix 0\tg_matrix 1\tg_matrix Z\n" << g_matrix.format(fmt) << "\n" << std::endl;
#endif
    // calculate energies from the norm squared values
    const auto &[energy_matrix, energy_matrix_grad] = power_scaled_loss_function(g_matrix, g0, k, n, gradient, numOmpThreads);
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

        //#pragma omp parallel for collapse(3) num_threads(numOmpThreads)
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
        for (int i = 0; i < g_list.size(); i++) {
            //for each grouping
            for (int j = 0; j < g_list_grad[i].size(); j++) {
                d_energy_d_d_vector[j].array() += energy_matrix_grad->col(i).rowwise().template replicate<5>().array() *
                        g_list_grad[i][j].array();
//                // OMP line by line was slower (!)
//                for (int r = 0; r < d_energy_d_d_vector[j].rows(); ++r) {
//                    const auto &temp = (*energy_matrix_grad)(r, i) * g_list_grad[i][j].row(r);
//#pragma omp atomic
//                    d_energy_d_d_vector[j].row(r) += temp;
//                }
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
                // d_energy_d_r_array[i].col(j) = (d_arrays_grad[i].array() * d_energy_d_d_vector[i].replicate(Eigen::fix<3>, num_pairs).reshaped(num_pairs, Eigen::fix<15>).array()) (Eigen::all, Eigen::seq(j, Eigen::fix<14>, Eigen::fix<3>)).rowwise().sum();
            }
            //std::cout << "d_energy_d_r_array[" << i<< "]\n" << d_energy_d_r_array[i] << std::endl << std::endl;
        }

        std::vector<CoordsMatrixType<KEnRef_Real> > gradients(num_models);
#pragma omp parallel for num_threads(numOmpThreads)
        for (int i = 0; i < num_models; i++) {
            gradients.at(i) = CoordsMatrixType<KEnRef_Real>::Zero(num_atoms, 3);
        }
        // propagate the internuclear vector derivatives back onto the atomic coordinates
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
        for (int p = 0; p < num_pairs; ++p) {
            // seq_len(dim(d_energy_d_r_array)[1])
            for (int m = 0; m < num_models; m++) {
                const auto [atomId0, atomId1] = atomId_pairs[p];
                const auto &pair_grad = d_energy_d_r_array[m].row(p);
#pragma omp atomic
                gradients[m].row(atomId0) -= pair_grad;
#pragma omp atomic
                gradients[m].row(atomId1) += pair_grad;
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
    std::map<std::string, int> atomNames_2_atomIds,
    KEnRef_Real k, KEnRef_Real n, const bool gradient, const int numOmpThreads) {

    //	std::cout << "coord_array_to_energy(atomName_pairs_) called" << std::endl;
    const auto &atomId_pairs = atomNamePairs_2_atomIdPairs(atomName_pairs, atomNames_2_atomIds);
    return KEnRef::coord_array_to_g_energy_refactored(coord_array, *atomId_pairs, grouping_list, g0, k, n, gradient,
                                                      numOmpThreads);
}

template<typename KEnRef_Real>
NamedRowVector<KEnRef_Real>
KEnRef<KEnRef_Real>::calculateLambdaVector(const SpecDenData<KEnRef_Real> &currentSpecDenData,
                                           const NamedRowVector<KEnRef_Real> &rates, int numOmpThreads) {
    // lambda_vector <- -colSums(rates[rownames(spec_den_data_list[[i]][["lambda_coef"]])]*spec_den_data_list[[i]][["lambda_coef"]])
    const auto &lambda_coef_matrix = currentSpecDenData.get_lambda_coef();
    const auto &rowNames = currentSpecDenData.get_lambda_coef().rowNames();
    NamedRowVector<KEnRef_Real> lambda_vector(lambda_coef_matrix.cols());
    lambda_vector.setZero();
    lambda_vector.setColNames(lambda_coef_matrix.colNames());

    // Calculate -colSums(rates[row_names,] * lambda_coef)
#pragma omp parallel for num_threads(numOmpThreads)
    for (size_t j = 0; j < rowNames.size(); ++j) {
        // Get corresponding rate value
        auto rate = rates(rowNames[j]);
        // Multiply rate with lambda_coef row and accumulate
        lambda_vector -= rate * lambda_coef_matrix.row(j);
    }
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
        // calculate inter nuclear vectors
        const auto &r_arrays = coord_array_to_r_array( coord_array_meter, *atomNamePairs_2_atomIdPairs(currentSpecDenData.get_atom_pairs(), atomNames_2_atomIds), numOmpThreads);

        // calculate dipole-dipole interaction tensors [and their derivatives]
        auto &&[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, gradient, false, numOmpThreads);
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

    const int num_pairs = d_energy_d_a_matrix.rows();      // Number of atom pairs
    const int num_groupings = a_coef.rows();               // Number of groupings //TODO not sure about the name
    // const int num_eigenvalues = d_energy_d_a_matrix.cols(); // Number of eigenvalues

    // Initialize output matrix with zeros
    Eigen::MatrixX<KEnRef_Real> d_energy_d_g_matrix = Eigen::MatrixX<KEnRef_Real>::Zero(num_pairs, num_groupings);

    // Backpropagation: d_energy/d_g = sum_j (d_energy/d_a_j * da_j/d_g)
    for (int i = 0; i < num_groupings; ++i) {          // Loop over groupings (rows of a_coef)
        for (int j = 0; j < a_coef.cols(); ++j) {      // Loop over eigenvalues (columns of a_coef)
            if (a_coef(i, j) != KEnRef_Real(0)) {      // Only for non-zero coefficients
                d_energy_d_g_matrix.col(i) += d_energy_d_a_matrix.col(j) / a_coef(i, j);
            }
        }
    }
    return d_energy_d_g_matrix;
}

template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::optional<std::vector<CoordsMatrixType<KEnRef_Real> > > >
KEnRef<KEnRef_Real>::coord_array_to_sigma_energy(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array,
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

    for (int i = 0; i < specDenDataSize; ++i) {
        std::vector<CoordsMatrixType<KEnRef_Real> > coord_array_meter(numModels);
        for (int m = 0; m < numModels; ++m)
            coord_array_meter[m] = coord_array[m] * 1e-10;

        auto &currentSpecDenData = spec_den_data_list[i];
        std::shared_ptr<std::vector<std::tuple<int, int>>> atom_id_pairs = atomNamePairs_2_atomIdPairs(currentSpecDenData.get_atom_pairs(), atomNames_2_atomIds);
        // calculate inter nuclear vectors
        const auto &r_arrays =
            coord_array_to_r_array(coord_array_meter, *atom_id_pairs, numOmpThreads);

        // calculate dipole-dipole interaction tensors [and their derivatives]
        auto &&[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, gradient, false, numOmpThreads);
        d_arrays_grad_list.at(i) = std::move(d_arrays_grad);

        const auto &grouping_matrix = IoUtils::subset_idx_to_grouping_mat(currentSpecDenData.get_multiple_grouping());
        //calculate the factor by which the number of models should be expanded
        uint nshift = grouping_matrix.cols() / numModels;
        const auto &d_arrays_shifted = array_shift(d_arrays, nshift);

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
        sigma(Eigen::seqN(index, sigmas.rows()), Eigen::all) = sigma_vector_list[i];
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
        sigma0(Eigen::seqN(index, sigma0s.rows()), Eigen::all) = sigma0s;
        index += sigma0s.rows();
    }

    const auto&[sigma_energy_vec, sigma_energy_grad_vec] = power_scaled_loss_function(sigma, sigma0, k, n, gradient, numOmpThreads);
    KEnRef_Real ret1 = sigma_energy_vec.sum();
    std::optional<std::vector<CoordsMatrixType<KEnRef_Real> >> d_energy_d_coord_array;

    if (gradient) {
        d_energy_d_coord_array = std::vector<CoordsMatrixType<KEnRef_Real> >(coord_array.size());
        for (int i = 0; i < coord_array.size(); ++i) {
            d_energy_d_coord_array->at(i) = std::move(CoordsMatrixType<KEnRef_Real>::Zero(coord_array[0 /*i*/].rows(), 3));
        }

        for (int i = 0 ; i < specDenDataSize; ++i) {
            const auto &currentSpecDenData = spec_den_data_list[i];

            // extract particular segment of energy gradient
            const auto &d_energy_d_sigma = sigma_energy_grad_vec.value()(Eigen::seqN(indexes[i],currentSpecDenData.sigmas().rows()), Eigen::all);

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

            //TODO Cache the use the first call to this function and use it directly.
            std::shared_ptr<std::vector<std::tuple<int, int>>> atom_id_pairs = atomNamePairs_2_atomIdPairs(currentSpecDenData.get_atom_pairs(), atomNames_2_atomIds);
            // accumulate back-propagated derivatives from r array into coord array gradient
            const auto & gradients =
                coord_array_to_r_array_backprop(
                    coord_array, *atom_id_pairs, d_energy_d_r_array,0);
            for (int j = 0; j < numModels; ++j)
                d_energy_d_coord_array->at(j) += gradients.at(j);
        }
        for (int i = 0; i < d_energy_d_coord_array->size(); ++i)
            d_energy_d_coord_array->at(i) *= 1e-10; //TODO Shouldn't it be 1e+10?
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

    //#pragma omp parallel for collapse(3) num_threads(numOmpThreads)
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
    for (int i = 0; i < num_something; i++) {
        //for each grouping
        for (int j = 0; j < num_ensembleMembers; j++) {
            d_energy_d_d_vector[j].array() += d_energy_d_g_matrix.col(i).rowwise().template replicate<5>().array() *
                    d_g_matrix_d_d_array[i][j].array();
            //             // OMP line by line was slower (!)
            //             for (int r = 0; r < d_energy_d_d_vector[j].rows(); ++r) {
            //                 const auto &temp = d_energy_d_g_matrix(r, i) * d_g_matrix_d_d_array[i][j].row(r);
            // #pragma omp atomic
            //                 d_energy_d_d_vector[j].row(r) += temp;
            //             }
        }
    }
    return d_energy_d_d_vector;
}

template<typename KEnRef_Real>
std::vector<CoordsMatrixType<KEnRef_Real> >
KEnRef<KEnRef_Real>::r_array_to_d_array_backprop(
    const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > &d_d_array_d_r_array,
    const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &d_energy_d_d_vector,
    /*const int num_interactions, const int num_ensembleMembers, */int numOmpThreads) {

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
#pragma omp parallel for collapse(3) num_threads(numOmpThreads)
    for (int m = 0; m < num_models; m++) {
        //TODO evaluate the effect of swapping the inner 2 loops on performance
        for (int p = 0; p < num_pairs; ++p) {
            const int modelShift = p / num_interactions;
            const int rowShift = p % num_interactions;
            for (int xyzIdx = 0; xyzIdx < 3; xyzIdx++) {
                d_energy_d_r_array.at(m)(p, xyzIdx) =
                    (d_d_array_d_r_array[m].row(p).array() *
                        d_energy_d_d_vector[m * shiftFactor + modelShift].row(rowShift).template replicate<3, 1>().reshaped(1, 15).array())
                            (Eigen::seqN(xyzIdx, Eigen::fix<5>, Eigen::fix<3>)).sum();
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
    const auto &[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, gradient, false, numOmpThreads);

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
    auto [d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, false, false, numOmpThreads);

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
    NamedMatrix<KEnRef_Real> a_matrix = NamedMatrix<KEnRef_Real>::Zero(g_matrix.rows(), a_coef.cols());

    // Parallelize over output features (columns of a_matrix)
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
    for (int i = 0; i < a_coef.cols(); ++i) {
        // Iterate over groupings (columns of g_matrix / rows of a_coef)
        for (int j = 0; j < a_coef.rows(); ++j) {
            const auto coef = a_coef(j, i);
            if (coef != 0) {
                // SAFE: Column operations in Eigen are thread-safe
                a_matrix.col(i) += g_matrix.col(j) * coef;
            }
        }
    }
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
        const auto &group_dnorm_array_1model = r_array_to_d_array(group_rnorm_array, false, false, numOmpThreads);
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
            std::get<0>(r_array_to_d_array(group_r_array, false, false, numOmpThreads));

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

template class SpecDenData<float>;
template class SpecDenData<double>;
template class KEnRef<float>;
template class KEnRef<double>;
#undef VERBOSE
