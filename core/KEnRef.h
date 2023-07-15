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
#include "../config/KEnRefConfig.h"

typedef float KEnRef_Real;
typedef Eigen::Matrix<KEnRef_Real, Eigen::Dynamic, 3, Eigen::RowMajor> CoordsMatrixType;
typedef Eigen::Map<CoordsMatrixType> CoordsMapType;
typedef Eigen::Map<const CoordsMatrixType> CoordsMapTypeConst;   // a read-only map

class KEnRef {
public:
	KEnRef();
	virtual ~KEnRef();
//	KEnRef(const KEnRef &other);
//  KEnRef(KEnRef &&other) noexcept ;
//	KEnRef& operator=(const KEnRef &other);
//	KEnRef& operator=(KEnRef &&other);

	static std::tuple<Eigen::Matrix<float, Eigen::Dynamic, 5>, Eigen::Matrix<float, Eigen::Dynamic, 15>>
	r_array_to_d_array(const Eigen::MatrixX3<float>& Nxyz, bool gradient=false);

	//return tuple where item0 is dipole-dipole interaction tensors (model<pairs, 5_tensor_elements>)
	//item1 is derivatives (It is a vector of 2D Matrix (models<pairId, (5_tensor_elements * XYZ)>).
	static std::tuple<std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>, std::vector<Eigen::Matrix<float, Eigen::Dynamic, 15>>>
	r_array_to_d_array(
			const std::vector<Eigen::MatrixX3<float>>& r_array,	//model<pairID, XYZ>
			bool gradient=false
			);

	// Calculate group norm squared from dipole-dipole interaction tensors, and optionally their gradients in the 5 tensor dimensions.
	static std::tuple<Eigen::VectorX<float>, std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>>
	d_array_to_g(
			const std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>& d_arrays, //vector (models<pairId, tensor_elements>) with interaction tensors
			const std::vector<std::vector<int>>& grouping, //groupings of models to average interaction tensors
			bool gradient=false
			);

	// Calculate group norm squared from dipole-dipole interaction tensors
	static std::tuple<std::vector<Eigen::VectorX<float>>, std::vector<std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>>>
	d_arrays_to_g(
			const std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>& d_array, //vector (models<pairId, tensor_elements>) with interaction tensors
			const std::vector<std::vector<std::vector<int>>>& groupings, //groupings of models to average interaction tensors
			bool gradient=false
			);

	//Calculate internuclear vectors from atomic coordinates
	//!\param coord_array vector of Nx3 Matrix (models<atoms, xyz>) with atomic coordinates
	//!\param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array` Matrices)
	//!\return a vector of Matrces (models<pairIndex, xyz>) with internuclear vectors.
	//Atom pair names in R code follow the format `resSeq:Atom-resSeq:Atom` (not implemented here).
	static std::vector<Eigen::MatrixX3<float>>
	coord_array_to_r_array(
            const std::vector<Eigen::MatrixX3<float>> &coord_array,
            const std::vector<std::tuple<int, int>> &atomId_pairs
			);

	static std::tuple<float, std::vector<Eigen::MatrixX3<float>>>
	coord_array_to_energy(
			const std::vector<Eigen::MatrixX3<float>>& coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
			const std::vector<std::tuple<std::string, std::string>>& atomName_pairs, 	// Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
			const std::vector<std::vector<std::vector<int>>>& grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
			const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>& g0, //target group norm squared values
			float k, //force constant
			std::map<std::string, int> atomNames_2_atomIds,
			bool gradient=false
			);

	static std::tuple<float, std::vector<Eigen::MatrixX3<float>>>
	coord_array_to_energy(
            std::vector<Eigen::MatrixX3<float>> coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
			std::vector<std::tuple<int, int>> atomId_pairs, 	// Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
			const std::vector<std::vector<std::vector<int>>>& grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
			const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> &g0, //target group norm squared values
			float k, //force constant
			bool gradient=false
			);

	static Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>
	coord_array_to_g(
			const std::vector<Eigen::MatrixX3<float>>& coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
			const std::vector<std::tuple<int, int>>& atomId_pairs, 	// Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
			const std::vector<std::vector<std::vector<int>>>& grouping_list	// list of lists of integer vectors giving groupings of models to average interaction tensors
			);

	// Calculate restraint energy from group norm squared values
	// returns restraint energy calculated using \eqn{k*(g-g0)^2}
	static std::tuple<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>
	g_to_energy(
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> g_list,	// current group norm squared values
			Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> g0,	// target group norm squared values
			float k = 1.0,	// force constant
			bool gradient=false // whether to calculate the derivative
			);

	//Collects list/vector of norm squared of all groups in a single matrix (num_pairIds, num_models (or num of grouping vectors?))
	static Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>
	vectorOfVectors_to_Matrix(std::vector<Eigen::VectorX<float>> g_vect);

};

#endif /* KENREF_H_ */
