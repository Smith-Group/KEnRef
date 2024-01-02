/*
 * KEnRef.cpp
 *
 *  Created on: May 8, 2023
 *      Author: amr
 */

#include <limits>
#include <memory>
#include <Eigen/Dense>
#include "KEnRef.h"
#include <iostream>//FIXME for testing only
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
std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15>>
KEnRef<KEnRef_Real>::r_array_to_d_array(const Eigen::MatrixX3<KEnRef_Real> &Nxyz, bool gradient) {
//	std::cout << "r_array_to_d_array(Nxyz) called" << std::endl;
    auto N = Nxyz.rows();

    auto x = Nxyz.col(0).array();
    auto y = Nxyz.col(1).array();
    auto z = Nxyz.col(2).array();
    auto sqrt3 = static_cast<KEnRef_Real>(sqrt(3.0));

    enum ci {
        x2, y2, z2,
//		xyz,
        x2_y2_z2,
        xy, xz, yz, x2_minusy2,
        x2_y2_z2_p52, x2_y2_z2_p72,
//		sqrt3_x2_y2_z2_p52, sqrt3_x2_y2_z2_p72,
        half_minusx2_minusy2__z2,
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
//	CACHE(xyz) 			= CACHE(xy) * z;
    CACHE(x2_y2_z2) = CACHE(x2) + CACHE(y2) + CACHE(z2); //x2 + y2 + z2;
    CACHE(x2_minusy2) = CACHE(x2) - CACHE(y2);
    CACHE(x2_y2_z2_p52) = CACHE(x2_y2_z2).pow(5).sqrt() + std::numeric_limits<KEnRef_Real_t>::epsilon();
    CACHE(half_minusx2_minusy2__z2) = ((-CACHE(x2) - CACHE(y2)) / 2) + CACHE(z2);
//    std::cout << "cache" << cache << std:: endl;

    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> ret1(N, 5);
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

	Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> ret2 = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>(N, 15);

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
std::tuple<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>>, std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15>>>
KEnRef<KEnRef_Real>::r_array_to_d_array(const std::vector<Eigen::MatrixX3<KEnRef_Real>> &models_Nxyz, bool gradient) {
//	std::cout << "r_array_to_d_array(models_Nxyz) called" << std::endl;
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>> ret1;
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15>> ret2;
    ret1.reserve(models_Nxyz.size());
    ret2.reserve(models_Nxyz.size());
    for (const auto &Nxyz: models_Nxyz) {
        auto[arr1, arr2] = r_array_to_d_array(Nxyz, gradient);
        ret1.emplace_back(arr1);
        ret2.emplace_back(arr2);
    }
    return {ret1, ret2};
}

template<typename KEnRef_Real>
std::tuple<std::vector<Eigen::VectorX<KEnRef_Real>>, std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>>>>
KEnRef<KEnRef_Real>::d_arrays_to_g(
        const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>> &d_array, //vector (models<pairId, tensor_elements>) with interaction tensors
        const std::vector<std::vector<std::vector<int>>> &groupings, //groupings of models to average interaction tensors (per dipole dipole interaction pair), i.e. outer list for pairId and inner list for modelId
        bool gradient) {
//	std::cout << "d_arrays_to_g() called" << std::endl;
    std::vector<Eigen::VectorX<KEnRef_Real>> ret1;
    std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>>> ret2;
    ret1.reserve(d_array.size());

    for (const auto & grouping : groupings) {
        const auto& [ret1_temp, ret2_temp] = d_array_to_g(d_array, grouping, gradient);
        ret1.emplace_back(ret1_temp);
        if (gradient)
            ret2.emplace_back(ret2_temp);
    }
    return {ret1, ret2};
}

template<typename KEnRef_Real>
std::tuple<Eigen::VectorX<KEnRef_Real>, std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>>>
KEnRef<KEnRef_Real>::d_array_to_g(
        const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>> &d_arrays, //vector (models<pairId, tensor_elements>) with interaction tensors
        const std::vector<std::vector<int>> &grouping, //groupings of models to average interaction tensors (per dipole-dipole interaction pair), i.e. outer list for pairId and inner list for modelId
        bool gradient) {

    auto num_models = d_arrays.size();
    auto num_pairIds = d_arrays[0].rows();
    auto num_groups = grouping.size();

//	std::cout << "num_models " 	<< num_models << std::endl;
//	std::cout << "num_pairIds "	<< num_pairIds << std::endl;
//	std::cout << "num_groups "	<< num_groups << std::endl;
    Eigen::VectorX<KEnRef_Real> ret1(num_pairIds);
    ret1.fill(0.0);

    //Every element of ret2 (every d_matrix_grad) is a matrix(num_pairIds x 5)
    std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>> ret2;
    if (gradient) {
        ret2.reserve(num_models);
        for (int i = 0; i < num_models; ++i) {
            ret2.emplace_back(num_pairIds, 5); // d_matrix_grad
        }
    } else {
        ret2 = std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>>(0);
    }

    for (int i = 0; i < num_groups; ++i) { //for every grouping block
        //create a new empty d_matrix (filled with 0) to carry "average dipole interaction tensor" every group
        std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>> d_matrix;// d_matrix is in the form num_pairIds<1 x 5> (previosly I thought it is num_pairIds<num_groups x 5> but it seems that I am wrong. in such case it would be just an Eigen::vector
        //TODO remove the vector altogether and keep only an Eigen::Matrix(num_pairIds, 5)
        d_matrix.reserve(num_pairIds);
        for (int temp = 0; temp < num_pairIds; ++temp) {
            Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> d_matrix_temp(num_pairIds, 5); // TODO 95% confirmed
            d_matrix_temp.fill(0);
            d_matrix.emplace_back(d_matrix_temp);
        }

        std::vector<int> currentGrouping = grouping[i];
        auto currentGroupSize = currentGrouping.size();

        // sum the dipole interaction tensors within each group/block
        for (int j = 0; j < currentGroupSize; ++j) { //for every member of the grouping block
            //sum relevant models into relevant groups (e.g. model 1 & 2 into group 1, and models 3 & 4 into group 2)
            //d_matrix += d_array[currentGrouping[j]]
            for (auto &d_matrix_temp: d_matrix) {
//				std::cout << "j = "<< j <<", d_matrix_temp        (" << d_matrix_temp.rows() << " x " << d_matrix_temp.cols() << ")" <<std::endl;
//				std::cout << "d_arrays[currentGrouping[j]](" << d_arrays[currentGrouping[j]].rows() << " x " << d_arrays[currentGrouping[j]].cols() << ")" <<std::endl;
                d_matrix_temp += d_arrays[currentGrouping[j]];
            }
        }

        if (gradient) {
            //# matrix is pairIDs * interaction tensor elements
            Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> d_matrix_grad_temp(num_pairIds, 5);
            // calculate  d_matrix_grad
            for (int j = 0; j < num_pairIds; j++) {
//				d_matrix_grad_temp.row(j) = d_matrix[j].row(i) * 2 / num_models / currentGroupSize; //just left the line in case I need it later
                d_matrix_grad_temp = d_matrix[j] * 2 / num_models / currentGroupSize;
            }
            for (int k = 0; k < currentGroupSize; k++) {
                // All models of the same group equally share the same value
                // All elements are equal (i.e. all models get the same overall (average ?) value at the end.
//				Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5> d_matrix_grad_lvalue = ret2[currentGrouping[k]]; //Matrix(num_pairIds x 5)
                ret2[currentGrouping[k]] = d_matrix_grad_temp;
            }
        }

        // Divide d_matrix by currentGroupSize to get the average
        for (int j = 0; j < num_pairIds; j++) {
            auto &d_matrix_temp = d_matrix[j];
            d_matrix_temp.array() /= currentGroupSize;
        }

        // calculate self dot product (norm squared) and accumulate group's contribution to mean g
        for (int j = 0; j < num_pairIds; j++) {
            auto &d_matrix_temp = d_matrix[j];
//			std::cout << "d_matrix_temp " << j << std::endl << d_matrix_temp << std::endl;
//            std::cout << "d_matrix_temp.row(j) " << d_matrix_temp.row(j) << " d_matrix_temp.row(j).squaredNorm() " << d_matrix_temp.row(j).squaredNorm() << " currentGroupSize " << currentGroupSize << " num_models " << num_models << std::endl;
            ret1(j) += d_matrix_temp.row(j).squaredNorm() * currentGroupSize / num_models;
        }
    }
//	std::cout << "ret1" << std::endl << ret1 << std::endl;
    return {ret1, ret2};
}


// Calculate restraint energy from group norm squared values
// returns restraint energy (loss function) calculated using \eqn{k*(g-g0)^2} +/- gradient using \eqn{2*k(g-g0)}
template<typename KEnRef_Real>
std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>>
KEnRef<KEnRef_Real>::g_to_energy_uncorrected(
        Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g,    // current group norm squared values
        Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g0,    // target group norm squared values
        KEnRef_Real k,    // force constant
        bool gradient) {
//	std::cout << "g   (" << g.rows() << " x " << g.cols() << ")" <<std::endl;
//	std::cout << "g0  (" << g0.rows() << " x " << g0.cols() << ")" <<std::endl;
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
std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>>
KEnRef<KEnRef_Real>::g_to_energy(
        Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g,    // current group norm squared values
        Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g0,    // target group norm squared values. Cached first time, ignored next times.
        KEnRef_Real k,    // force constant
        KEnRef_Real n,    // correction power. Cached first time, ignored next times.
        bool gradient) {
//	std::cout << "g   (" << g.rows() << " x " << g.cols() << ")" <<std::endl;
//	std::cout << "g0  (" << g0.rows() << " x " << g0.cols() << ")" <<std::endl;

    auto g_ = g.array() /*+ std::numeric_limits<KEnRef_Real_t>::epsilon()*/;
    auto g0_ = g0.array();
    std::cout << "g_ " << g_ << std::endl;
    Eigen::Array<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> loss =
            (Eigen::pow(1 + Eigen::abs(g_), n) - 1) * Eigen::sign(g_) - (Eigen::pow(1 + Eigen::abs(g0_), n) - 1) * Eigen::sign(g0_) ;
    std::cout << "loss " << loss << std::endl;
    auto ret1 = k * loss.square().matrix(); //This value may become infinity if it excceds 3.402823466E38 in a single precision float
    std::cout << "g ret1 " << ret1 << std::endl;


    if (gradient) {
        auto ret2 = 2.0 * k * loss * (n * Eigen::pow(1 + Eigen::abs(g_), (n-1)));
        std::cout << "g ret2 " << ret2 << std::endl;
        return {ret1, ret2.matrix()};
    } else {
        return {ret1, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>{}};
    }
}

template<typename KEnRef_Real>
Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>
KEnRef<KEnRef_Real>::vectorOfVectors_to_Matrix(std::vector<Eigen::VectorX<KEnRef_Real>> g_vect) {
    Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g_mat(g_vect[0].rows(), g_vect.size());
    for (int i = 0; i < g_vect.size(); ++i) {
        g_mat.col(i) = g_vect[i];
    }
    return g_mat;
}

template<typename KEnRef_Real>
std::vector<Eigen::MatrixX3<KEnRef_Real>>
KEnRef<KEnRef_Real>::coord_array_to_r_array( //TODO shall we change the return type to CoordsMatrixType ?
        const std::vector<Eigen::MatrixX3<KEnRef_Real>> &coord_array,
        const std::vector<std::tuple<int, int>> &atomId_pairs) {
//	std::cout << "coord_array_to_r_array() called" << std::endl;
    std::vector<Eigen::MatrixX3<KEnRef_Real>> ret;
    ret.reserve(coord_array.size());
    for (int model_no = 0; model_no < coord_array.size(); ++model_no) {
        Eigen::MatrixX3<KEnRef_Real> mat(atomId_pairs.size(), 3);
        for (int i = 0; i < mat.rows(); ++i) {
            auto [atom0, atom1] = atomId_pairs.at(i);
            mat.row(i) = coord_array[model_no].row(atom1) - coord_array[model_no].row(atom0);
        }
        ret.emplace_back(mat);
    }
    return ret;
}


template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::vector<CoordsMatrixType<KEnRef_Real>>>
KEnRef<KEnRef_Real>::coord_array_to_energy(
        const std::vector<Eigen::MatrixX3<KEnRef_Real>>& coord_array,    //Every vector item is an Nx3 Matrix representing atom coordinates of a model.
        const std::vector<std::tuple<std::string, std::string>>& atomName_pairs,    // Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
        const std::vector<std::vector<std::vector<int>>>& grouping_list,    // list of lists of integer vectors giving groupings of models to average interaction tensors
        const Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>& g0,
        KEnRef_Real k,
        std::map<std::string, int> atomNames_2_atomIds,
        KEnRef_Real n, bool gradient) {
//	std::cout << "coord_array_to_energy(atomName_pairs_) called" << std::endl;
    std::vector<std::tuple<int, int>> atomId_pairs{};
    // Fill the vector using atomNames_2_atomIds
    for (const auto &[left, right]: atomName_pairs) {
        // I use at() instead of operator[] to force an exception to be thrown
        atomId_pairs.emplace_back(atomNames_2_atomIds.at(left), atomNames_2_atomIds.at(right));
    }
//    std::tuple<KEnRef_Real, std::vector<CoordsMatrixType>> temp = KEnRef::coord_array_to_energy(coord_array, atomId_pairs, grouping_list, g0, k, gradient);
//    std::vector<CoordsMatrixType> &derivatives = std::get<1>(temp);
//    for (int i = 0; i < derivatives.size(); ++i) {
//        CoordsMatrixType &matrix = derivatives[i];
//        for (int j = 0; j < matrix.rows(); ++j)
//            for (int l = 0; l < 3; ++l)
//                matrix(j, l) = static_cast<KEnRef_Real>(i * 10000 + j * 10 + l);
//    }
//    return temp;
    return KEnRef::coord_array_to_energy(coord_array, atomId_pairs, grouping_list, g0, k, n, gradient);
}

template<typename KEnRef_Real>
std::tuple<KEnRef_Real, std::vector<CoordsMatrixType<KEnRef_Real>>>
KEnRef<KEnRef_Real>::coord_array_to_energy(
        std::vector<Eigen::MatrixX3<KEnRef_Real>> coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
		std::vector<std::tuple<int, int>> atomId_pairs, 	// Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
		const std::vector<std::vector<std::vector<int>>>& grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
		const Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> &g0, KEnRef_Real k, KEnRef_Real n, bool gradient)
{
	// calculate inter nuclear vectors
	auto r_arrays = coord_array_to_r_array(coord_array, atomId_pairs);

	// calculate dipole-dipole interaction tensors [and their derivates]
	auto [d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, gradient);

	// calculate norm squared for different groupings of dipole-dipole interaction tensors
	auto [g_list, g_list_grad] = d_arrays_to_g(d_arrays, grouping_list, gradient);
//	for(int i = 0; i < g_list_grad.size(); i++){
//		auto g_list_grad_i = g_list_grad[i];
//		for (int j = 0; j < g_list_grad_i.size(); j++){
//			std::cout << "g_list_grad " << i+1 << " " << j+1 << std::endl;
//			std::cout << g_list_grad_i[j] << std::endl;
//		}
//	}

	auto g_matrix = vectorOfVectors_to_Matrix(g_list);

	// calculate energies from the norm squared values
//    auto [energy_matrix, energy_matrix_grad] = g_to_energy_uncorrected(g_matrix, g0, k, gradient);
    auto [energy_matrix, energy_matrix_grad] = g_to_energy(g_matrix, g0, k, n, gradient);
//	std::cout << "energy_matrix" << std::endl << energy_matrix << std::endl;
//	std::cout << "energy_matrix_grad" << std::endl << energy_matrix_grad << std::endl;

	// return the sum of all the individual restraint energies
	KEnRef_Real sum = energy_matrix.sum();
//	std::cout << "energy_matrix sum" << std::endl << sum << std::endl;

	//Add derivatives using the chain rule: de/dr = de/dd  * dd/dr = de/dg * dg/dd * dd/dr
	if(gradient){
		auto num_pairs = atomId_pairs.size();
		auto num_models = coord_array.size();
		auto num_atoms = coord_array[0].rows();
		// First calculate de/dd = de/dg * dg/dd for all individual interaction tensor components
		std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>>d_energy_d_d_array; //<num_models(num_pairs, 5)>
		d_energy_d_d_array.reserve(num_models); //num_models
		for(int i = 0; i < num_models; i++){
			d_energy_d_d_array.emplace_back(Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>::Zero(static_cast<int>(num_pairs), 5));
		}

		for(int i = 0; i < g_list.size(); i++){//for each grouping
			auto e_matrix_grad_replicated = energy_matrix_grad.col(i).rowwise().template replicate<5>().array();
			auto g_list_grad_group_i = g_list_grad[i]; //<num_models>(num_pairs x 5)
//			std::cout << "d_energy_d_d_array" << " after iteration " << i << std::endl;
			for (int j = 0; j < g_list_grad_group_i.size(); j++) { //(num_pairs x 5)
				d_energy_d_d_array[j].array() += (e_matrix_grad_replicated * g_list_grad_group_i[j].array());
//				std::cout << d_energy_d_d_array[j] << std::endl;
			}
		}

		// Then calculate de/dr = de/dd  * dd/dr for each xyz component of the internuclear vectors
		std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15>> d_energy_d_r_array_all;
		d_energy_d_r_array_all.reserve(num_models);
		for(int i = 0; i < num_models; i++){
//			std::cout << "d_energy_d_d_array[i].replicate(3,1).reshaped(num_pairs, 15)" << std::endl << d_energy_d_d_array[i].replicate(3,1).reshaped(num_pairs, 15) << std::endl;
//			std::cout << "d_arrays_grad[i]" << std::endl << d_arrays_grad[i] << std::endl;
			d_energy_d_r_array_all.emplace_back(d_arrays_grad[i].array() * d_energy_d_d_array[i].replicate(3,1).reshaped(num_pairs, 15).array());
//			std::cout << "d_energy_d_d_array_all[" << i <<"]" << std::endl << d_energy_d_r_array_all[i] <<std::endl;
		}

		// sum the individual interaction tensor component derivatives associated with x, y, and z
		std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 3>> d_energy_d_r_array;
		d_energy_d_r_array.reserve(num_models);
		for(int i = 0; i< num_models; i++){
			d_energy_d_r_array.emplace_back(num_pairs, 3);
//			auto temp_array = d_energy_d_r_array[i];
			Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> temp;
//			std::cout << "d_energy_d_r_array" << std::endl;
			for(int j = 0; j< 3; j++){
				temp = d_energy_d_r_array_all[i](Eigen::all, Eigen::seq(j, 15, Eigen::fix<3>));
				d_energy_d_r_array[i].col(j) = temp.rowwise().sum();
			}
//			std::cout << d_energy_d_r_array[i] << std::endl;
		}

		std::vector<CoordsMatrixType<KEnRef_Real>> gradients;
		gradients.reserve(num_models);
		for(int i = 0; i < num_models; i++){
			gradients.emplace_back(num_atoms, 3);
			gradients[i].setZero(); //TOOD There should be a better way than this
		}
		// propagate the internuclear vector derivatives back onto the atomic coordinates
		for (int p = 0; p < num_pairs; ++p) { // seq_len(dim(d_energy_d_r_array)[1])
			for(int m = 0; m < num_models; m++){
				auto[atomId0, atomId1] = atomId_pairs[p];
				auto pair_grad = d_energy_d_r_array[m].row(p);
				gradients[m].row(atomId0) -= pair_grad;
				gradients[m].row(atomId1) += pair_grad;
			}
		}

//		std::cout << "gradients" << std::endl;
//		for(int m = 0; m < num_models; m++){
//			std::cout << "model " << m << std::endl;
//			std::cout << gradients[m] << std::endl;
//		}
		return {sum, gradients};
	}else{
		return {sum, std::vector<CoordsMatrixType<KEnRef_Real>>{}};
	}
}


template<typename KEnRef_Real>
Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>
KEnRef<KEnRef_Real>::coord_array_to_g(
		const std::vector<Eigen::MatrixX3<KEnRef_Real>>& coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
		const std::vector<std::tuple<int, int>>& atomId_pairs, 	// Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
		const std::vector<std::vector<std::vector<int>>>& grouping_list)	// list of lists of integer vectors giving groupings of models to average interaction tensors
{
//	std::cout << "coord_array_to_g() called" << std::endl;
	// calculate internuclear vectors
	auto r_arrays = coord_array_to_r_array(coord_array, atomId_pairs);

	// calculate dipole-dipole interaction tensors [and their derivatives]
	auto [d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays);

	// calculate norm squared for different groupings of dipole-dipole interaction tensors
	//		g_list <- lapply(grouping_list, function(grouping) d_array_to_g(d_array, grouping, gradient=FALSE))
	auto [g_list, ignore] = d_arrays_to_g(d_arrays, grouping_list);

	return vectorOfVectors_to_Matrix(g_list);
}

//TODO Shall we move it (with its test of course) to GMXKenRefForceProvider?
template<typename KEnRef_Real>
void
KEnRef<KEnRef_Real>::saturate(CoordsMatrixType<KEnRef_Real> &derivatives_rectified, int simulationIndex, KEnRef_Real energy) {
    constexpr KEnRef_Real threshold = 1000. * 1000.;
//    for(int i = 0; i < derivatives_rectified.rows(); i++){
//        auto rowBlock = derivatives_rectified.row(i);
    for(auto rowBlock: derivatives_rectified.rowwise()){
        KEnRef_Real scaleDown = rowBlock.squaredNorm() / threshold;
        if (scaleDown > 1){
            rowBlock /= std::sqrt(scaleDown);
//            std::cout << "Simulation # " << simulationIndex << " row: " << i << " scaleDown " << scaleDown <<
//                      (scaleDown > 1.0 ? " Scaled down " : " NOT USED ") << std::endl;
        }
    }
}

template class KEnRef<float>;
template class KEnRef<double>;