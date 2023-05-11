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

class KEnRef {
public:
	KEnRef();
	virtual ~KEnRef();
	KEnRef(const KEnRef &other);
	KEnRef(KEnRef &&other);
//	KEnRef& operator=(const KEnRef &other);
//	KEnRef& operator=(KEnRef &&other);

	static std::tuple<Eigen::MatrixXf, Eigen::MatrixXf>
	r_array_to_d_array(const Eigen::MatrixX3f& Nxyz, bool gradient=false);

	//return tuple where item0 is dipole-dipole interaction tensors (model<pairs, 5_tensor_elements>)
	//item1 is derivatives (It is a vector of 2D Matrix (models<pairId, (5_tensor_elements * XYZ)>).
	static std::tuple<std::vector<Eigen::MatrixXf>, std::vector<Eigen::MatrixXf>>
	r_array_to_d_array(
			std::vector<Eigen::MatrixX3<float>> r_array,	//model<pairID, XYZ>
			bool gradient
			);

	//Calculate internuclear vectors from atomic coordinates
	//!\param coord_array vector of Nx3 Matrix (models<atoms, xyz>) with atomic coordinates
	//!\param atom_pairs matrix with each row having the names or indices of an atom pair (first dimension in `coord_array` Matrices)
	//!\return a vector of Matrces (models<pairIndex, xyz>) with internuclear vectors.
	//Atom pair names in R code follow the format `resSeq:Atom-resSeq:Atom` (not implemented here).
	static std::vector<Eigen::MatrixX3<float>>
	coord_array_to_r_array(
			std::vector<Eigen::MatrixX3<float>> coord_array,
			std::vector<std::tuple<int, int>> atomId_pairs
			);

	static std::tuple<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>
	coord_array_to_energy(
			std::vector<Eigen::MatrixX3<float>> coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
			std::vector<std::tuple<std::string, std::string>> atomName_pairs, 	// Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
			std::vector<std::vector<int>> grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
			float g0, //target group norm squared values
			float k, //force constant
			std::map<std::string, int> atomNames_2_atomIds,
			bool gradient=false
			);

	static std::tuple<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>
	coord_array_to_energy(
			std::vector<Eigen::MatrixX3<float>> coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
			std::vector<std::tuple<int, int>> atomId_pairs, 	// Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
			std::vector<std::vector<int>> grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
			float g0, //target group norm squared values
			float k, //force constant
			bool gradient=false
			);
};

#endif /* KENREF_H_ */
