/*
 * KEnRef.cpp
 *
 *  Created on: May 8, 2023
 *      Author: amr
 */

#include <omp.h>
#include <limits>
#include <memory>
#include <Eigen/Dense>
#include "core/KEnRef.h"
//#include <iostream>//FIXME for testing only
//#include <utility>

template<typename KEnRef_Real>
KEnRef<KEnRef_Real>::KEnRef() = default;

template<typename KEnRef_Real>
KEnRef<KEnRef_Real>::~KEnRef() = default;

//KEnRef::KEnRef(const KEnRef &other) = default;
//KEnRef::KEnRef(KEnRef &&other)  noexcept {}
//KEnRef& KEnRef::operator=(const KEnRef &other) {}
//KEnRef& KEnRef::operator=(KEnRef &&other) {}

template<typename KEnRef_Real>
std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> >
KEnRef<KEnRef_Real>::r_array_to_d_array(const CoordsMatrixType<KEnRef_Real> &Nxyz, bool gradient, int numOmpThreads) {
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
    CACHE(x2_y2_z2_p52) = CACHE(x2_y2_z2).pow(5).sqrt() + std::numeric_limits<KEnRef_Real_t>::epsilon();
    CACHE(half_minusx2_minusy2__z2) = ((-CACHE(x2) - CACHE(y2)) / 2) + CACHE(z2);
    //    std::cout << "cache" << cache << std:: endl;

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
    ret1.array() /= CACHE(x2_y2_z2_p52).rowwise().template replicate<5>()/* + std::numeric_limits<KEnRef_Real_t>::epsilon()*/;

    //    std:: cout << "ret1 after division" << std::endl << ret1 << std::endl;
    if (!gradient) {
        return {ret1, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15>{}};
    }
    //    std::cout << "CACHE after d ret1 " << cache << std::endl;

    CACHE(x2_y2_z2_p72) = CACHE(x2_y2_z2).pow(7).sqrt() + std::numeric_limits<KEnRef_Real_t>::epsilon();
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
                                        int numOmpThreads) {
    //	std::cout << "r_array_to_d_array(models_Nxyz) called" << std::endl;
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > ret1;
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15> > ret2;
    ret1.reserve(models_Nxyz.size());
    ret2.reserve(models_Nxyz.size());
    for (const auto &Nxyz: models_Nxyz) {
        auto [arr1, arr2] = r_array_to_d_array(Nxyz, gradient, numOmpThreads);
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
    //groupings of models to average interaction tensors (per dipole-dipole interaction pair), i.e. outer list for pairId and inner list for modelId
    bool gradient, int numOmpThreads) {
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
std::tuple<Eigen::VectorX<KEnRef_Real>, std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > >
KEnRef<KEnRef_Real>::d_array_to_g(
    const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > &d_arrays,
    const std::vector<std::vector<int> > &grouping,
    bool gradient, int numOmpThreads) {
    auto num_models = d_arrays.size();
    auto num_pairIds = d_arrays[0].rows();

    // std::cout << "num_models " 	<< num_models << " num_pairIds " << num_pairIds << "num_groups " << num_groups << std::endl;
    Eigen::VectorX<KEnRef_Real> ret1 = Eigen::VectorX<KEnRef_Real>::Zero(num_pairIds);

    //Every element of ret2 (every d_matrix_grad) is a matrix(num_pairIds x 5)
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> > ret2;
    if (gradient) {
        ret2 = std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >(num_models);
#pragma omp parallel for num_threads(numOmpThreads)
        for (int i = 0; i < num_models; ++i) {
            ret2.at(i) = std::move(Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>(num_pairIds, 5)); // AKA d_matrix_grad
        }
    } else {
        ret2 = std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> >(0);
    }

    //for every grouping block
    for (const auto & currentGrouping : grouping) {
        //create a new empty d_matrix (filled with 0) to carry "average dipole interaction tensor" every group
        Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> d_matrix = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(num_pairIds, 5);

        const auto currentGroupSize = currentGrouping.size();
        const auto CURRENT_GROUP_SIZE_real = static_cast<KEnRef_Real_t>(currentGroupSize);
        const auto currentGroupSize_OVER_num_models_real = CURRENT_GROUP_SIZE_real / num_models;

        // sum the dipole interaction tensors within each group/block
#pragma omp parallel for num_threads(numOmpThreads) // schedule(static) //num_threads(gmx_omp_nthreads_get(ModuleMultiThread::Default))
        for (int j = 0; j < currentGroupSize; ++j) {
            //for every member of the grouping block
            //sum relevant models into relevant groups (e.g. model 1 & 2 into group 1, and models 3 & 4 into group 2)
//#pragma omp atomic
            d_matrix += d_arrays[currentGrouping[j]];
        }

        if (gradient) {
            auto TWO_OVER_num_models_currentGroupSize =
                    static_cast<KEnRef_Real_t>(2.0) / num_models / CURRENT_GROUP_SIZE_real;
            // calculate  d_matrix_grad. d_matrix_grad shape is pairIDs * interaction tensor elements
            const auto& d_matrix_grad = d_matrix * TWO_OVER_num_models_currentGroupSize;
            // All models of the same group equally share the same value
#pragma omp parallel for num_threads(numOmpThreads)
            for (int j = 0; j < currentGroupSize; j++) {
                // All elements are equal (i.e. all models get the same overall (average ?) value at the end.
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
        for (int j = 0; j < num_pairIds; j++) {
            const auto contribution = d_matrix.row(j).squaredNorm() * currentGroupSize_OVER_num_models_real;
#pragma omp atomic
            ret1(j) += contribution;
        }
    }
    //	std::cout << "ret1" << std::endl << ret1 << std::endl;
    return {ret1, ret2};
}


// Calculate restraint energy from group norm squared values
// returns restraint energy (loss function) calculated using \eqn{k*(g-g0)^2} +/- gradient using \eqn{2*k(g-g0)}
template<typename KEnRef_Real>
std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic,
    Eigen::Dynamic> >
KEnRef<KEnRef_Real>::g_to_energy_uncorrected(
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g, // current group norm squared values
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g0, // target group norm squared values
    KEnRef_Real k, // force constant
    bool gradient, int numOmpThreads) {
    Eigen::Array<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g_minus_g0 = g.array() - g0.array();
    auto ret1 = k * g_minus_g0.square().matrix();
    if (gradient) {
        auto ret2 = 2.0 * k * g_minus_g0.matrix();
        return {ret1, ret2};
    } else {
        return {ret1, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>{}};
    }
}

// Calculate restraint energy from group norm squared values
// returns restraint energy (loss function) calculated using \eqn{k*(Sign[g]*Abs[g]^n - Sign[g0]*Abs[g0]^n)^2} +/- gradient using \eqn{2*k(g-g0)}
template<typename KEnRef_Real>
std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic,
    Eigen::Dynamic> >
KEnRef<KEnRef_Real>::g_to_energy(
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g, // current group norm squared values
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g0, // target group norm squared values. Cached first time, ignored next times.
    KEnRef_Real k, // force constant
    KEnRef_Real n, // correction power.
    bool gradient,
    int numOmpThreads, lossFunction lossFunc) {
    //	std::cout << "g   (" << g.rows() << " x " << g.cols() << ")" <<std::endl;
    //	std::cout << "g0  (" << g0.rows() << " x " << g0.cols() << ")" <<std::endl;

    auto g_arr = g.array() /*+ std::numeric_limits<KEnRef_Real_t>::epsilon()*/;
    auto g0_arr = g0.array();
    //    std::cout << "g_arr " << g_arr << std::endl;
    Eigen::Array<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> common; //TODO better keep it as an ArrayExpression (or simply auto) for better optimization
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> ret1;
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> ret2{}; //TODO don't initialize. use **nullOpt** instead

    switch (lossFunc) {
        case SQRT_ABS_POWER_N:
            common = (Eigen::pow(1 + Eigen::abs(g_arr), n) - 1) * Eigen::sign(g_arr) - (
                         Eigen::pow(1 + Eigen::abs(g0_arr), n) - 1) * Eigen::sign(g0_arr);
        //          std::cout << "loss " << loss << std::endl;
            ret1 = k * common.square().matrix(); //This value may become infinity if it excceds 3.402823466E38 in a single precision float
            if (gradient) {
                ret2 = (2.0 * k * common * (n * Eigen::pow(1 + Eigen::abs(g_arr), (n - 1)))).matrix();
            } /* else {
                ret2 = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>{};
            }*/
            break;
        case LOG_ABS_DIFFERENCE_OVER_NOE0:
            common = g0_arr + Eigen::abs(g_arr - g0_arr);
            ret1 = (k * Eigen::log(common / g0_arr)).matrix();
            if (gradient) {
                ret2 = (k * Eigen::sign(g_arr - g0_arr) / common).matrix();
            } /*else {
                ret2 = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>{};
            }*/
            break;
    }
//        std::cout << "g ret1 " << std::endl << ret1 << std::endl;
//        std::cout << "g ret2 " << std::endl << ret2 << std::endl;
    return {ret1, ret2};
}

template<typename KEnRef_Real>
Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>
KEnRef<KEnRef_Real>::vectorOfVectors_to_Matrix(
    std::vector<Eigen::VectorX<KEnRef_Real> > g_vect/*, int numOmpThreads*/) {
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g_mat(g_vect[0].rows(), g_vect.size());
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
    //	std::cout << "coord_array_to_r_array() called" << std::endl;
    std::vector<CoordsMatrixType<KEnRef_Real> > ret(coord_array.size());
//#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
    for (int model_no = 0; model_no < coord_array.size(); ++model_no) {
        ret.at(model_no) = {atomId_pairs.size(), 3};
#pragma omp parallel for num_threads(numOmpThreads)
        for (int i = 0; i < atomId_pairs.size(); ++i) {
            auto [atom0, atom1] = atomId_pairs.at(i);
            ret.at(model_no).row(i) = coord_array[model_no].row(atom1) - coord_array[model_no].row(atom0);
        }
    }
    return ret;
}


template<typename KEnRef_Real>
std::shared_ptr<std::vector<std::tuple<int, int> > >
KEnRef<KEnRef_Real>::atomNamePairs_2_atomIdPairs(
    const std::vector<std::tuple<std::string, std::string> > &atomName_pairs,
    std::map<std::string, int> &atomNames_2_atomIds) {
    auto atomId_pairs = std::make_shared<std::vector<std::tuple<int, int> > >(atomName_pairs.size());
    // Fill the vector using atomNames_2_atomIds
#pragma omp parallel for num_threads(numOmpThreads)
    for (int i = 0; i < atomName_pairs.size(); ++i) {
        auto [left, right] = atomName_pairs.at(i);
        // I use at() instead of operator[] to force an exception to be thrown
        atomId_pairs->at(i) = std::move(std::tuple<int, int>{
            atomNames_2_atomIds.at(left), atomNames_2_atomIds.at(right)
        });
    }
    return atomId_pairs;
}

template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::vector<CoordsMatrixType<KEnRef_Real> > >
KEnRef<KEnRef_Real>::coord_array_to_energy(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array, //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
    const std::vector<std::tuple<std::string, std::string> > &atomName_pairs, // Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
    const std::vector<std::vector<std::vector<int> > > &grouping_list, // list of lists of integer vectors giving groupings of models to average interaction tensors
    const Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> &g0,
    std::map<std::string, int> atomNames_2_atomIds,
    KEnRef_Real k, KEnRef_Real n, bool gradient, int numOmpThreads) {
    //	std::cout << "coord_array_to_energy(atomName_pairs_) called" << std::endl;
    const auto &atomId_pairs = atomNamePairs_2_atomIdPairs(atomName_pairs, atomNames_2_atomIds);
    return KEnRef::coord_array_to_energy(coord_array, *atomId_pairs, grouping_list, g0, k, n, gradient, numOmpThreads);
}

template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::vector<CoordsMatrixType<KEnRef_Real> > >
KEnRef<KEnRef_Real>::coord_array_to_energy(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array, //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
    const std::vector<std::tuple<int, int> > &atomId_pairs, // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
    const std::vector<std::vector<std::vector<int> > > &grouping_list, // list of lists of integer vectors giving groupings of models to average interaction tensors
    const Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> &g0, KEnRef_Real k, KEnRef_Real n, bool gradient,
    int numOmpThreads) {
    // calculate inter nuclear vectors
    const auto& r_arrays = coord_array_to_r_array(coord_array, atomId_pairs, numOmpThreads);

    // calculate dipole-dipole interaction tensors [and their derivatives]
    const auto &[d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, gradient, numOmpThreads);

    // calculate norm squared for different groupings of dipole-dipole interaction tensors
    const auto &[g_list, g_list_grad] = d_array_to_g_multiple_groupings(d_arrays, grouping_list, gradient, numOmpThreads);
    //	for(int i = 0; i < g_list_grad.size(); i++){
    //		auto g_list_grad_i = g_list_grad[i];
    //		for (int j = 0; j < g_list_grad_i.size(); j++){
    //			std::cout << "g_list_grad " << i+1 << " " << j+1 << std::endl;
    //			std::cout << g_list_grad_i[j] << std::endl;
    //		}
    //	}

    const auto& g_matrix = vectorOfVectors_to_Matrix(g_list/*, numOmpThreads*/);

    // calculate energies from the norm squared values
    //const auto &[energy_matrix, energy_matrix_grad] = g_to_energy_uncorrected(g_matrix, g0, k, gradient, numOmpThreads /*, KEnRef::SQRT_ABS_POWER_N*/);
    const auto &[energy_matrix, energy_matrix_grad] = g_to_energy(g_matrix, g0, k, n, gradient,
                                                              numOmpThreads /*, KEnRef::SQRT_ABS_POWER_N*/);
    //	std::cout << "energy_matrix" << std::endl << energy_matrix << std::endl;
    //	std::cout << "energy_matrix_grad" << std::endl << energy_matrix_grad << std::endl;

    // return the sum of all the individual restraint energies
    KEnRef_Real sum = energy_matrix.sum();
//    std::cout << "energy_matrix sum " << sum << std::endl;

    //Add derivatives using the chain rule: de/dr = de/dd  * dd/dr = de/dg * dg/dd * dd/dr
    if (gradient) {
        const auto num_pairs = atomId_pairs.size();
        const auto num_models = coord_array.size();
        const auto num_atoms = coord_array[0].rows();
        // First calculate de/dd = de/dg * dg/dd for all individual interaction tensor components
        std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> > d_energy_d_d_vector(num_models); //<num_models(num_pairs, 5)>
#pragma omp parallel for num_threads(numOmpThreads)
        for (int i = 0; i < num_models; i++)
            d_energy_d_d_vector.at(i) = std::move(
                Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>::Zero(static_cast<int>(num_pairs), 5));

//#pragma omp parallel for collapse(3) num_threads(numOmpThreads)
#pragma omp parallel for collapse(2) num_threads(numOmpThreads)
        for (int i = 0; i < g_list.size(); i++) { //for each grouping
            for (int j = 0; j < g_list_grad[i].size(); j++) {
                d_energy_d_d_vector[j].array() += (energy_matrix_grad.col(i).rowwise().template replicate<5>().array() * g_list_grad[i][j].array());
//                // OMP line by line was slower (!)
//                for (int r = 0; r < d_energy_d_d_vector[j].rows(); ++r) {
//                    const auto &temp = energy_matrix_grad(r, i) * g_list_grad[i][j].row(r);
//#pragma omp atomic
//                    d_energy_d_d_vector[j].row(r) += temp;
//                }
            }
        }

        // Then calculate de/dr = de/dd  * dd/dr for each xyz component of the internuclear vectors
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
                const auto& pair_grad = d_energy_d_r_array[m].row(p);
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
        return {sum, gradients};
    } else {
        return {sum, std::vector<CoordsMatrixType<KEnRef_Real> >{}};
    }
}


template<typename KEnRef_Real>
Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>
KEnRef<KEnRef_Real>::coord_array_to_g(
    const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array, //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
    const std::vector<std::tuple<int, int> > &atomId_pairs, // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
    const std::vector<std::vector<std::vector<int> > > &grouping_list, // list of lists of integer vectors giving groupings of models to average interaction tensors
    int numOmpThreads) {
    //	std::cout << "coord_array_to_g() called" << std::endl;
    // calculate internuclear vectors
    const auto &r_arrays = coord_array_to_r_array(coord_array, atomId_pairs, numOmpThreads);

    // calculate dipole-dipole interaction tensors [and their derivatives]
    auto [d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, numOmpThreads);

    // calculate norm squared for different groupings of dipole-dipole interaction tensors
    //		g_list <- lapply(grouping_list, function(grouping) d_array_to_g(d_array, grouping, gradient=FALSE))
    auto [g_list, ignore] = d_array_to_g_multiple_groupings(d_arrays, grouping_list, numOmpThreads);

    return vectorOfVectors_to_Matrix(g_list/*, numOmpThreads*/);
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
        const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array, //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
        const std::vector<std::tuple<int, int> > &atomId_pairs, // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
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
    const auto &group_r_array = KEnRef<KEnRef_Real>::coord_array_to_r_array(coord_array, atomId_pairs, numOmpThreads);

    //calculate matrix of radii
    // group_r_mat <- sqrt(rowSums(group_r_array^2, dims=2))
    std::vector<Eigen::VectorX<KEnRef_Real>> group_r_mat;
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
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> group_dnorm_array_tempMeans = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(group_r_mat[0].rows(), 5);
    for (int i = 0; i < numModels; ++i) {
        const auto& group_rnorm_array = group_r_array[i].array() / group_r_mat[i].rowwise().template replicate<3>().array();
        const auto &group_dnorm_array_1model = KEnRef<KEnRef_Real>::r_array_to_d_array(group_rnorm_array, false, numOmpThreads);
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
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>> group_d_array =
            std::get<0>(KEnRef<KEnRef_Real>::r_array_to_d_array(group_r_array, false, numOmpThreads));

    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>> group_dnorm_array_alt;
    group_dnorm_array_alt.reserve(numModels);
    for (int i = 0; i < numModels; ++i) {
        Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>temp = group_d_array[i].array().rowwise().norm().template replicate<1, 5>();
        group_dnorm_array_alt.emplace_back(group_d_array[i].array() / temp.array());
    }

    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> group_dnorm_array_alt_tempMeans = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(group_dnorm_array_alt[0].rows(), 5);
    for (int i = 0; i < numModels; ++i) {
        group_dnorm_array_alt_tempMeans.array() += group_dnorm_array_alt[i].array();
    }
    group_dnorm_array_alt_tempMeans.array() /= numModels;
    Eigen::VectorX<KEnRef_Real> group_s2_alt = group_dnorm_array_alt_tempMeans.array().square().colwise().sum();

//    assert((group_s2 - group_s2_alt).array() < 1e7); //FIXME didn't work
    //#######################################################
    return group_s2;
}

//template<typename KEnRef_Real>
//Eigen::VectorX<KEnRef_Real>
//s2OrderParamsPseudocode(
//        const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array, //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
//        const std::vector<std::string> &atomNames, //All atom names
//        const std::vector<std::tuple<int, int> > &atomId_pairs, // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
//        int numOmpThreads) {
//    //below are just some assumptions. We may need them later.
//    int numModels = 2;
//    int numHAtoms = 5000;
//    int numProtonGroups = 2500;
//    int numProtonGroupsAfterExclusions = 2300;
//
//    numModels = coord_array.size();
//
//    //TODO find *proton_groups*
//
//    //TODO find exclusions of proton_groups
//
//    //TODO # **group_list** is a list with each element giving the protons in the group, matching the input order
//    // implement it as vector<pair<string, vector<string>>>
//    // e.g.
//    // " CA  MET     1 " : [" HA  MET     1 "]
//    // " HB1 MET     1 " : [" HB1 MET     1 "]
//    // ...
//    // " CE  MET     1 " : [" HE1 MET     1 ", " HE2 MET     1 ", " HE3 MET     1 "]
//    std::vector<std::pair<std::string, std::vector<std::string>>> group_list{};
//
//    //TODO *group_mean_coord* includes the averages of the coordinates, in its model.
//    // each row represents the average coordinates in its corresponding group_list element.
//    // implement it as vector<CoordsMatrixType>.
//    std::vector<CoordsMatrixType<KEnRef_Real>> group_mean_coord;
//    group_mean_coord.reserve(numModels);
//    for (int m = 0; m < numModels; ++m) {
////        group_mean_coord.at(m) = CoordsMatrixType<KEnRef_Real>::Zero(group_list.size(), 3);
//        group_mean_coord.at(m) = CoordsMatrixType<KEnRef_Real>(group_list.size(), 3);
//    }
//
//    for (int m = 0; m < numModels; ++m) {
//        for (int groupIdx = 0; groupIdx < group_list.size(); ++groupIdx) {
//            auto& groupAtomNames = group_list[groupIdx].second;
//
//            std::string &atomName = groupAtomNames[0];
//            int atomId = 0; //TODO find Id of atomName
//            group_mean_coord[m].row(groupIdx) = coord_array[m].row(atomId);
//            if (groupAtomNames.size() > 1){
//                for (int i = 1; i < groupAtomNames.size(); ++i) {
//                    std::string &atomName = groupAtomNames[i];
//                    int atomId = 0; //TODO find Id of atomName
//                    group_mean_coord[m].row(groupIdx) += coord_array[m].row(atomId);
//                }
//                group_mean_coord[m].row(groupIdx) /= groupAtomNames.size();
//            }
//        }
//    }
//
//    //TODO calculate a distance matrix for every member of the ensemble (every Dimension separately) in *dist_array*
//    //TODO average *dist_array* over ensemble members to get average distances in the ensemble
//    //TODO determine set of atoms that are within the distance cutoff.
//    // Put them in *group_pairs*. This is equivalent to *atomId_pairs* in C++ code
//
//    //TODO calculate array of internuclear vectors *group_r_array*
//    // group_r_array <- ke::coord_array_to_r_array(aperm(group_mean_coord, c(2,1,3)), group_pairs)
//    const auto &group_r_array = KEnRef<KEnRef_Real>::coord_array_to_r_array(group_mean_coord, atomId_pairs, numOmpThreads);
//
//    //TODO calculate matrix of radii
//    // group_r_mat <- sqrt(rowSums(group_r_array^2, dims=2))
//    std::vector<Eigen::VectorX<KEnRef_Real>> group_r_mat;
//    group_r_mat.reserve(numModels);
//    for (int i = 0; i < numModels; ++i) {
//        group_r_mat.at(i) = std::move(Eigen::VectorX<KEnRef_Real>(group_r_array[i].rows()));
//    }
//    for (int i = 0; i < numModels; ++i) {
//        for (int j = 0; j < group_r_array[i].rows(); ++j) {
//            group_r_mat[i](j) = group_r_array[i].row(j).norm();
//        }
//    }
//
//    //TODO # calculate order parameter using normalized array of internuclear vectors
//    // group_rnorm_array <- group_r_array/as.vector(group_r_mat)
//    // group_dnorm_array <- ke::r_array_to_d_array(group_rnorm_array)
//    // group_s2 <- rowSums(colMeans(aperm(group_dnorm_array, c(2,1,3)))^2)
//    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> group_dnorm_array_tempMeans = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(group_r_mat[0].rows());
//    for (int i = 0; i < numModels; ++i) {
//        const auto& group_rnorm_array = group_r_array[i].array() / group_r_mat[i].array();
//        const auto &group_dnorm_array_1model = KEnRef<KEnRef_Real>::r_array_to_d_array(group_rnorm_array, false, numOmpThreads);
//        group_dnorm_array_tempMeans += std::get<0>(group_dnorm_array_1model);
//    }
//    group_dnorm_array_tempMeans /= numModels;
////    Eigen::VectorX<KEnRef_Real> group_s2 = group_dnorm_array_tempMeans.array().square().colwise().sum();
//    Eigen::VectorX<KEnRef_Real> group_s2 = Eigen::VectorX<KEnRef_Real>::Zero(group_dnorm_array_tempMeans.rows());
//    group_s2 = group_dnorm_array_tempMeans.array().square().colwise().sum();
//
//    //################################################
//
//    //TODO Alternate approach directly from d_array
//    // group_d_array <- ke::r_array_to_d_array(group_r_array)
//    // group_dnorm_array_alt <- group_d_array/as.vector(sqrt(rowSums(group_d_array^2, dims=2)))
//    // group_s2_alt <- rowSums(colMeans(aperm(group_dnorm_array_alt, c(2,1,3)))^2)
//    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>> group_d_array =
//            std::get<0>(KEnRef<KEnRef_Real>::r_array_to_d_array(group_r_array, false, numOmpThreads));
//
//    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>> group_dnorm_array_alt;
//    group_dnorm_array_alt.reserve(numModels);
//    for (int i = 0; i < numModels; ++i) {
//        group_dnorm_array_alt.at(i) = group_d_array[i].data() / group_d_array[i].data().rowwise().norm();
//    }
//
//    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> group_dnorm_array_alt_tempMeans = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>::Zero(group_dnorm_array_alt[0].rows());
//    for (int i = 0; i < numModels; ++i) {
//        group_dnorm_array_alt_tempMeans += group_dnorm_array_alt.data();
//    }
//    group_dnorm_array_alt_tempMeans.data() /= numModels;
//    Eigen::VectorX<KEnRef_Real> group_s2_alt = group_dnorm_array_alt_tempMeans.array().square().colwise().sum();
//
//    assert((group_s2.data() - group_s2_alt.data()) < 1e7);
//    //#######################################################
//    return group_s2;
//}

template class KEnRef<float>;
template class KEnRef<double>;
