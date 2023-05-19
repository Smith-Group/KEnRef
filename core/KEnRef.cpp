/*
 * KEnRef.cpp
 *
 *  Created on: May 8, 2023
 *      Author: amr
 */

#include <memory>
#include <Eigen/Dense>
#include "KEnRef.h"
#include <iostream>//FIXME for testing only

KEnRef::KEnRef() {
	// TODO Auto-generated constructor stub
}
KEnRef::~KEnRef() {
	// TODO Auto-generated destructor stub
}
KEnRef::KEnRef(const KEnRef &other) {
	// TODO Auto-generated constructor stub
}
KEnRef::KEnRef(KEnRef &&other) {
	// TODO Auto-generated constructor stub
}
//KEnRef& KEnRef::operator=(const KEnRef &other) {}
//KEnRef& KEnRef::operator=(KEnRef &&other) {}

std::tuple<Eigen::Matrix<float, Eigen::Dynamic, 5>, Eigen::Matrix<float, Eigen::Dynamic, 15>> KEnRef::r_array_to_d_array(const Eigen::MatrixX3f& Nxyz, bool gradient){
	int N = Nxyz.rows();

	auto x 			= Nxyz.col(0).array();
	auto y 			= Nxyz.col(1).array();
	auto z 			= Nxyz.col(2).array();
	float sqrt3		= sqrt(3.0);

	enum ci{
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
	Eigen::ArrayXXf cache(N, 15); //must be Array (not Matrix) to allow item-wise operations
#define CACHE(a) (cache.col(a))
	CACHE(x2) = x.square();
	CACHE(y2) = y.square();
	CACHE(z2) = z.square();
	CACHE(xy) 			= x * y;
	CACHE(xz) 			= x * z;
	CACHE(yz) 			= y * z;
//	CACHE(xyz) 			= CACHE(xy) * z;
	CACHE(x2_y2_z2) 	= CACHE(x2) + CACHE(y2) + CACHE(z2); //x2 + y2 + z2;
	CACHE(x2_minusy2) 	= CACHE(x2) - CACHE(y2);
	CACHE(x2_y2_z2_p52) = CACHE(x2_y2_z2).pow(5).sqrt();
	CACHE(half_minusx2_minusy2__z2) = ((-CACHE(x2) - CACHE(y2)) / 2) + CACHE(z2);


	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> ret1(N, 5);
	ret1.col(0) = CACHE(half_minusx2_minusy2__z2);
	ret1.col(1) = sqrt3 * CACHE(x2_minusy2) / 2;
	ret1.col(2) = sqrt3 * CACHE(xz);
	ret1.col(3) = sqrt3 * CACHE(yz);
	ret1.col(4) = sqrt3 * CACHE(xy);
//	std::cout << "ret1 before division" << std::endl << ret1 << std::endl;
//	std::cout << "x2_y2_z2" << std::endl << CACHE(x2_y2_z2) << std::endl;
//	std::cout << "x2_y2_z2 power 5/2" << std::endl << CACHE(x2_y2_z2_p52).rowwise().replicate<5>() << std::endl;

	ret1.array() /= CACHE(x2_y2_z2_p52).rowwise().replicate<5>(); //TODO double check this line

	if(!gradient){
		return {ret1, Eigen::Matrix<float, Eigen::Dynamic, 15>{}};
	}

	CACHE(x2_y2_z2_p72) = CACHE(x2_y2_z2).pow(7).sqrt();
//	CACHE(sqrt3_x2_y2_z2_p52) = sqrt3 * CACHE(x2_y2_z2_p52);
//	CACHE(sqrt3_x2_y2_z2_p72) = sqrt3 * CACHE(x2_y2_z2_p72);
	CACHE(sqrt3_over_x2_y2_z2_p52)			= sqrt3      / CACHE(x2_y2_z2_p52);
	CACHE(neg5_over_x2_y2_z2_p72) 			= -5.0       / CACHE(x2_y2_z2_p72);
	CACHE(neg5sqrt3_over_x2_y2_z2_p72)		= -5 * sqrt3 / CACHE(x2_y2_z2_p72);
	CACHE(neg5sqrt3_over_2_x2_y2_z2_p72)	= CACHE(neg5sqrt3_over_x2_y2_z2_p72) / 2;
//	CACHE(sqrt3_over_x2_y2_z2_p52)			=  sqrt3 / CACHE(x2_y2_z2).pow(5).sqrt();
//	CACHE(_neg10_over_sqrt3_x2_y2_z2_p72) 	= -10 / (sqrt3 * CACHE(x2_y2_z2).pow(7).sqrt());

	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> ret2 = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(N, 15);

	ret2.col(0)  = (x * CACHE(half_minusx2_minusy2__z2) * CACHE(neg5_over_x2_y2_z2_p72))	- (    x / CACHE(x2_y2_z2_p52));
	ret2.col(1)  = (y * CACHE(half_minusx2_minusy2__z2) * CACHE(neg5_over_x2_y2_z2_p72))	- (    y / CACHE(x2_y2_z2_p52));
	ret2.col(2)  = (z * CACHE(half_minusx2_minusy2__z2) * CACHE(neg5_over_x2_y2_z2_p72))	+ (2 * z / CACHE(x2_y2_z2_p52));

	ret2.col(3)  = (x * CACHE(x2_minusy2) * CACHE(neg5sqrt3_over_2_x2_y2_z2_p72))	+ (x * CACHE(sqrt3_over_x2_y2_z2_p52));
	ret2.col(4)  = (y * CACHE(x2_minusy2) * CACHE(neg5sqrt3_over_2_x2_y2_z2_p72)) 	- (y * CACHE(sqrt3_over_x2_y2_z2_p52));
	ret2.col(5)  = (z * CACHE(x2_minusy2) * CACHE(neg5sqrt3_over_2_x2_y2_z2_p72));

	ret2.col(6)  = (x * CACHE(xz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) 		+ (z * CACHE(sqrt3_over_x2_y2_z2_p52));
	ret2.col(7)  = (y * CACHE(xz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72));
	ret2.col(8)  = (z * CACHE(xz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) 		+ (x * CACHE(sqrt3_over_x2_y2_z2_p52));

	ret2.col(9)  = (x * CACHE(yz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72));
	ret2.col(10) = (y * CACHE(yz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) 		+ (z * CACHE(sqrt3_over_x2_y2_z2_p52));
	ret2.col(11) = (z * CACHE(yz) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) 		+ (y * CACHE(sqrt3_over_x2_y2_z2_p52));

	ret2.col(12) = (x * CACHE(xy) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) 		+ (y * CACHE(sqrt3_over_x2_y2_z2_p52));
	ret2.col(13) = (y * CACHE(xy) * CACHE(neg5sqrt3_over_x2_y2_z2_p72)) 		+ (x * CACHE(sqrt3_over_x2_y2_z2_p52));
	ret2.col(14) = (z * CACHE(xy) * CACHE(neg5sqrt3_over_x2_y2_z2_p72));

#undef CACHE

	return {ret1, ret2};
}

//	std::tuple<Eigen::MatrixXf, Eigen::MatrixXf> KEnRef::r_array_to_d_array(const Eigen::MatrixX3f& Nxyz, bool gradient){
std::tuple<std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>, std::vector<Eigen::Matrix<float, Eigen::Dynamic, 15>>>
KEnRef::r_array_to_d_array(const std::vector<Eigen::MatrixX3f>& models_Nxyz, bool gradient){
	std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>> ret1;
	std::vector<Eigen::Matrix<float, Eigen::Dynamic, 15>> ret2;
	ret1.reserve(models_Nxyz.size());
	ret2.reserve(models_Nxyz.size());
	for (auto Nxyz: models_Nxyz) {
		auto [arr1, arr2] = r_array_to_d_array(Nxyz, gradient);
		ret1.emplace_back(arr1);
		ret2.emplace_back(arr2);
	}
	return {ret1, ret2};
}

std::vector<Eigen::MatrixX3<float>>
KEnRef::coord_array_to_r_array(
		std::vector<Eigen::MatrixX3<float>> coord_array,
		std::vector<std::tuple<int, int>> atomId_pairs)
{
	std::vector<Eigen::MatrixX3<float>> ret;
	ret.reserve(coord_array.size());
	for (int model_no = 0; model_no < coord_array.size(); ++model_no) {
		Eigen::MatrixX3<float> mat(atomId_pairs.size(), 3);
		for (int i = 0; i < mat.rows(); ++i) {
			auto [atom0, atom1] = atomId_pairs.at(i);
			mat.row(i) = coord_array[model_no].row(atom1) - coord_array[model_no].row(atom0);
		}
		ret.emplace_back(mat);
	}
	return ret;
}


std::tuple<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>
KEnRef::coord_array_to_energy(
		std::vector<Eigen::MatrixX3<float>> coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
		std::vector<std::tuple<std::string, std::string>> atomName_pairs, 	// Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
		std::vector<std::vector<int>> grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
		float g0, float k,
		std::map<std::string, int> atomNames_2_atomIds,
		bool gradient)
{
	std::vector<std::tuple<int, int>> atomId_pairs{};
	// Fill the vector using atomNames_2_atomIds
    for (const auto& [key, value] : atomName_pairs){
    	// I use at() instead of operator[] to force an exception to be thrown
    	atomId_pairs.emplace_back(std::tuple<int, int>{atomNames_2_atomIds.at(key), atomNames_2_atomIds.at(value)});
    }
	return KEnRef::coord_array_to_energy(coord_array, atomId_pairs, grouping_list, g0, k, gradient);
}


std::tuple<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>
KEnRef::coord_array_to_energy(
		std::vector<Eigen::MatrixX3<float>> coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
		std::vector<std::tuple<int, int>> atomId_pairs, 	// Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
		std::vector<std::vector<int>> grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
		float g0, float k, bool gradient)
{
	// calculate internuclear vectors
	auto r_arrays = coord_array_to_r_array(coord_array, atomId_pairs);

	// calculate dipole-dipole interaction tensors [and their derivates]
	auto [d_arrays, d_arrays_derivatives] = r_array_to_d_array(r_arrays, gradient);

	// calculate norm squared for different groupings of dipole-dipole interaction tensors
	//g_list <- lapply(grouping_list, function(grouping) d_array_to_g(d_array, grouping, gradient=gradient))

	// convert the list above into a matrix (pairs, groupings)
	//g_matrix <- simplify2array(g_list)

	// calculate energies from the norm squared values
	//energy_matrix <- g_to_energy(g_matrix, g0, k, gradient=gradient)

	// return the sum of all the individual restraint energies
	//value <- sum(energy_matrix)

	//add derivates

	return {Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(0,0), Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>(0,0)};
}
