/*
 * KEnRef.h
 *
 *  Created on: May 8, 2023
 *      Author: amr
 */

#ifndef KENREF_H_
#define KENREF_H_

#include <memory>
#include <vector>
#include <map>
#include <tuple>

#include <Eigen/Dense>
// #include "../config/KEnRefConfig.h"

#ifdef DOUBLE
typedef double KEnRef_Real_t;
#else
typedef float KEnRef_Real_t;
#endif
template <typename KEnRef_Real>
using CoordsMatrixType = Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 3, Eigen::RowMajor>;
template<typename KEnRef_Real>
using CoordsMapType = Eigen::Map<CoordsMatrixType<KEnRef_Real>>;
template<typename KEnRef_Real>
using CoordsMapTypeConst = Eigen::Map<const CoordsMatrixType<KEnRef_Real>>;   // a read-only map

template<typename KEnRef_Real>
class KEnRef final {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	KEnRef();
    ~KEnRef();
//	KEnRef(const KEnRef &other);
//  KEnRef(KEnRef &&other) noexcept ;
//	KEnRef& operator=(const KEnRef &other);
//	KEnRef& operator=(KEnRef &&other);

    enum lossFunction{SQRT_ABS_POWER_N, LOG_ABS_DIFFERENCE_OVER_NOE0};

    //First element of the returned tuple is (numPairs, <d_1,d_2,d_3,d_4,d_5,>)
    //Second element of the returned tuple is (numPairs, <  x1d1,x2d1,x3d1,
    //                                                      x1d2,x2d2,x3d2,
    //                                                      x1d3,x2d3,x3d3,
    //                                                      x1d4,x2d4,x3d4,
    //                                                      x1d5,x2d5,x3d5>)
	static std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15>>
	r_array_to_d_array(const CoordsMatrixType<KEnRef_Real> &Nxyz, bool gradient=false, int numOmpThreads = 0);

	//return tuple where item0 is dipole-dipole interaction tensors (model<pairs, 5_tensor_elements>)
	//item1 is derivatives (It is a vector of 2D Matrix (models<pairId, (5_tensor_elements * XYZ)>).
	static std::tuple<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>>, std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 15>>>
	r_array_to_d_array(
			const std::vector<CoordsMatrixType<KEnRef_Real>>& models_Nxyz,	//model<pairID, XYZ>
			bool gradient=false, int numOmpThreads = 0
			);


    /** Calculate group norm squared from dipole-dipole interaction tensors, and optionally their gradients in the 5 tensor dimensions.
	 *
	 * @param  d_arrays vector (models<pairId, interaction tensor_elements>) with interaction tensors
	 * @param  grouping groupings of models to average interaction tensors (per dipole-dipole interaction pair), i.e. outer list for pairId and inner list for modelId
	 * @param  gradient whether to calculate & return the derivates
	 * @param  numOmpThreads number of OpenMP threads
	 * @return tuple containing
	 *			1) vector of norm squared for each atom pair,
	 *			2) derivates in a vector of Matrix (models<pairId, derivates of interaction tensor_elements>) or empty vector
	 */
	static std::tuple<Eigen::VectorX<KEnRef_Real>, std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>>>
	d_array_to_g(
			const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>> &d_arrays,
			const std::vector<std::vector<int>> &grouping,
			bool gradient = false, int numOmpThreads = 0
			);

	// Calculate group norm squared from dipole-dipole interaction tensors
	static std::tuple<std::vector<Eigen::VectorX<KEnRef_Real>>, std::vector<std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>>>>
	d_array_to_g_multiple_groupings(
			const std::vector<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 5>>& d_array, //vector (models<pairId, tensor_elements>) with interaction tensors
			const std::vector<std::vector<std::vector<int>>>& groupings, //groupings of models to average interaction tensors
			bool gradient=false, int numOmpThreads = 0
			);

	//Calculate internuclear vectors from atomic coordinates
	//!\param coord_array vector of Nx3 Matrix (models<atoms, xyz>) with atomic coordinates
	//!\param atomId_pairs matrix with each row having the indices of an atom pair (first dimension in `coord_array` Matrices)
	//!\param numOmpThreads number of available MpenMP threads
	//!\return a vector of Matrces (models<pairIndex, xyz>) with internuclear vectors.
	//Atom pair names in R code follow the format `resSeq:Atom-resSeq:Atom` (not implemented here).
	static std::vector<CoordsMatrixType<KEnRef_Real>>
	coord_array_to_r_array(
            const std::vector<CoordsMatrixType<KEnRef_Real>> &coord_array,
            const std::vector<std::tuple<int, int>> &atomId_pairs,  int numOmpThreads = 0
			);

    static
    std::shared_ptr<std::vector<std::tuple<int, int>>>
    atomNamePairs_2_atomIdPairs(const std::vector<std::tuple<std::string, std::string>> &atomName_pairs, std::map<std::string, int> &atomNames_2_atomIds);

    static std::tuple<KEnRef_Real, std::vector<CoordsMatrixType<KEnRef_Real>>>
	coord_array_to_energy(
			const std::vector<CoordsMatrixType<KEnRef_Real>>& coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
			const std::vector<std::tuple<std::string, std::string>>& atomName_pairs, 	// Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
			const std::vector<std::vector<std::vector<int>>>& grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
			const Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>& g0, //target group norm squared values
			std::map<std::string, int> atomNames_2_atomIds,
            KEnRef_Real k = 1.0, //force constant
            KEnRef_Real n = 0.25,
			bool gradient=false, int numOmpThreads = 0
			);

	static std::tuple<KEnRef_Real, std::vector<CoordsMatrixType<KEnRef_Real>>>
	coord_array_to_energy(
            const std::vector<CoordsMatrixType<KEnRef_Real>>& coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
			const std::vector<std::tuple<int, int>>& atomId_pairs, 	// Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
			const std::vector<std::vector<std::vector<int>>>& grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
			const Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> &g0, //target group norm squared values
			KEnRef_Real k = 1.0, //force constant
            KEnRef_Real n = 0.25,
			bool gradient=false, int numOmpThreads = 0,
            bool printStatistics = false
			);

	static Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>
	coord_array_to_g(
			const std::vector<CoordsMatrixType<KEnRef_Real>>& coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
			const std::vector<std::tuple<int, int>>& atomId_pairs, 	// Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
			const std::vector<std::vector<std::vector<int>>>& grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
            int numOmpThreads = 0
            );

	// Calculate restraint energy from group norm squared values
	// returns restraint energy calculated using \eqn{k*(g-g0)^2}
	static std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>>
	g_to_energy_uncorrected(
			Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g_list,	// current group norm squared values
			Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g0,	// target group norm squared values
			KEnRef_Real k = 1.0,	// force constant
            bool gradient=false, // whether to calculate the derivative
            int numOmpThreads = 0
            );

    // Calculate restraint energy from group norm squared values
    // returns restraint energy calculated using \eqn{k*(g^n -g0^n)^2}
    static std::tuple<Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>>
    g_to_energy(
            Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g_list,	// current group norm squared values
            Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic> g0,	// target group norm squared values
            KEnRef_Real k = 1.0,	// force constant
            KEnRef_Real n = 0.25,    // correction power
            bool gradient = false  // whether to calculate the derivative
            ,
            int numOmpThreads = 0, lossFunction lossFunc = KEnRef::SQRT_ABS_POWER_N
    );
	//Collects list/vector of norm squared of all groups in a single matrix (num_pairIds, num_models (or num of grouping vectors?))
	static Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, Eigen::Dynamic>
	vectorOfVectors_to_Matrix(std::vector<Eigen::VectorX<KEnRef_Real>> g_vect/*, int numOmpThreads = 0*/);
    static void saturate(CoordsMatrixType<KEnRef_Real> &derivatives_rectified, KEnRef_Real thresholdSquared,
                         int numOmpThreads = 0);

    static Eigen::VectorX<KEnRef_Real>
    s2OrderParams(
            const std::vector<CoordsMatrixType<KEnRef_Real> > &coord_array, //Every vector item is a Nx3 Matrix representing atom coordinates of a model.
            const std::vector<std::tuple<int, int> > &atomId_pairs, // Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
            int numOmpThreads);
};

#endif /* KENREF_H_ */
