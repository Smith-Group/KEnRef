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

	ret1.array() /= CACHE(x2_y2_z2_p52).rowwise().replicate<5>();

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

std::tuple<std::vector<Eigen::VectorX<float>>, std::vector<std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>>>
KEnRef::d_arrays_to_g(
		const std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>& d_array, //vector (models<pairId, tensor_elements>) with interaction tensors
		const std::vector<std::vector<std::vector<int>>>& groupings, //groupings of models to average interaction tensors (per dipole dipole interaction pair), i.e. outer list for pairId and inner list for modelId
		bool gradient)
{
	std::vector<Eigen::VectorX<float>> ret1;
	std::vector<std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>> ret2;
	ret1.reserve(d_array.size());

	for(int i = 0; i < groupings.size(); i++){
		auto grouping = groupings[i];
		auto [ ret1_temp, ret2_temp] = d_array_to_g(d_array, grouping, gradient);
		ret1.emplace_back(ret1_temp);
		if(gradient)
			ret2.emplace_back(ret2_temp);
	}
	return {ret1, ret2};
}

std::tuple<Eigen::VectorX<float>, std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>>
KEnRef::d_array_to_g(
		const std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>& d_arrays, //vector (models<pairId, tensor_elements>) with interaction tensors
		const std::vector<std::vector<int>>& grouping, //groupings of models to average interaction tensors (per dipole dipole interaction pair), i.e. outer list for pairId and inner list for modelId
		bool gradient){

	int num_models = d_arrays.size();
	int num_pairIds = d_arrays[0].rows();
	int num_groups = grouping.size();

//	std::cout << "num_models " 	<< num_models << std::endl;
//	std::cout << "num_pairIds "	<< num_pairIds << std::endl;
//	std::cout << "num_groups "	<< num_groups << std::endl;
	Eigen::VectorX<float> ret1(num_pairIds);
	ret1.fill(0.0);

	//Every element of ret2 (every d_matrix_grad) is a matrix(num_pairIds x 5)
	std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>> ret2;
	if(gradient){
		ret2.reserve(num_models);
		for (int i = 0; i < num_models; ++i) {
			ret2.emplace_back(Eigen::Matrix<float, Eigen::Dynamic, 5>(num_pairIds, 5)); // d_matrix_grad
		}
	}else{
		ret2 = std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>>(0);
	}

	for (int i = 0; i < num_groups; ++i) { //for every grouping block
		//create a new empty d_matrix (filled with 0) to carry "average dipole interaction tensor" every group
		std::vector<Eigen::Matrix<float, Eigen::Dynamic, 5>> d_matrix;// d_matrix is in the form num_pairIds<1 x 5> (previosly I thought it is num_pairIds<num_groups x 5> but it seems that I am wrong. in such case it would be just an Eigen::vector
		//TODO remove the vector altogether and keep only an Eigen::Matrix(num_pairIds, 5)
		d_matrix.reserve(num_pairIds);
		for (int temp = 0; temp < num_pairIds; ++temp) {
			Eigen::Matrix<float, Eigen::Dynamic, 5> d_matrix_temp(num_pairIds, 5); // TODO 95% confirmed
			d_matrix_temp.fill(0);
			d_matrix.emplace_back(d_matrix_temp);
		}

		std::vector<int> currentGrouping = grouping[i];
		int currentGroupSize = currentGrouping.size();

		// sum the dipole interaction tensors within each group/block
		for (int j = 0; j < currentGroupSize; ++j) { //for every member of the grouping block
			//sum relevant models into relevant groups (e.g. model 1 & 2 into group 1, and models 3 & 4 into group 2)
			//d_matrix += d_array[currentGrouping[j]]
			for(auto& d_matrix_temp : d_matrix){
//				std::cout << "d_matrix_temp               (" << d_matrix_temp.rows() << " x " << d_matrix_temp.cols() << ")" <<std::endl;
//				std::cout << "d_arrays[currentGrouping[j]](" << d_arrays[currentGrouping[j]].rows() << " x " << d_arrays[currentGrouping[j]].cols() << ")" <<std::endl;
				d_matrix_temp += d_arrays[currentGrouping[j]];
			}
		}

		if(gradient){
			//# matrix is pairIDs * interaction tensor elements
			Eigen::Matrix<float, Eigen::Dynamic, 5> d_matrix_grad_temp(num_pairIds, 5);
			// calculate  d_matrix_grad
			for (int j = 0; j < num_pairIds; j++) {
//				d_matrix_grad_temp.row(j) = d_matrix[j].row(i) * 2 / num_models / currentGroupSize; //just left the line in case I need it later
				d_matrix_grad_temp = d_matrix[j] * 2 / num_models / currentGroupSize;
			}
			for(int k = 0; k < currentGroupSize; k++){
				// All models of the same group equally share the same value
				// All elements are equal (i.e. all models get the same overall (average ?) value at the end.
//				Eigen::Matrix<float, Eigen::Dynamic, 5> d_matrix_grad_lvalue = ret2[currentGrouping[k]]; //Matrix(num_pairIds x 5)
				ret2[currentGrouping[k]] = d_matrix_grad_temp;
			}
		}

		// Divide d_matrix by currentGroupSize to get the average
		for(int j = 0; j < num_pairIds; j++){
			auto& d_matrix_temp = d_matrix[j];
			d_matrix_temp.array() /= currentGroupSize;
		}

		// calculate self dot product (norm squared) and accumulate group's contribution to mean g
		for(int j = 0; j < num_pairIds; j++){
			auto& d_matrix_temp = d_matrix[j];
//			std::cout << "d_matrix_temp " << j << std::endl << d_matrix_temp << std::endl;
			ret1(j) += d_matrix_temp.row(j).squaredNorm() * currentGroupSize / num_models;
		}
	}
//	std::cout << "ret1" << std::endl << ret1 << std::endl;
	return {ret1, ret2};
}



// Calculate restraint energy from group norm squared values
// returns restraint energy calculated using \eqn{k*(g-g0)^2} +/- gradient using \eqn{2*k(g-g0)}
std::tuple<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>>
KEnRef::g_to_energy(
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> g,	// current group norm squared values
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> g0,	// target group norm squared values
		float k,	// force constant
		bool gradient)
{
	std::cout << "g   (" << g.rows() << " x " << g.cols() << ")" <<std::endl;
	std::cout << "g0  (" << g0.rows() << " x " << g0.cols() << ")" <<std::endl;
	Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic> g_minus_g0 = g.array() - g0.array();
	auto ret1 = k * g_minus_g0.square().matrix();
    if (gradient) {
		auto ret2 = 2.0 * k *  g_minus_g0.matrix();
		return {ret1, ret2};
    }else{
    	return {ret1, Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>{}};
    }
}

Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>
KEnRef::vectorOfVectors_to_Matrix(std::vector<Eigen::VectorX<float>> g_vect){
	Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> g_mat(g_vect[0].rows(), g_vect.size());
	for (int i = 0; i < g_vect.size(); ++i) {
		g_mat.col(i) = g_vect[i];
	}
	return g_mat;
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


std::tuple<float, std::vector<Eigen::MatrixX3<float>>>
KEnRef::coord_array_to_energy(
		std::vector<Eigen::MatrixX3<float>> coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
		std::vector<std::tuple<std::string, std::string>> atomName_pairs, 	// Matrix with each row having the names of an atom pair (related to first dimension in `coord_array` matrices)
		std::vector<std::vector<std::vector<int>>> grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> g0,
		float k,
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


std::tuple<float, std::vector<Eigen::MatrixX3<float>>>
KEnRef::coord_array_to_energy(
		std::vector<Eigen::MatrixX3<float>> coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
		std::vector<std::tuple<int, int>> atomId_pairs, 	// Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
		std::vector<std::vector<std::vector<int>>> grouping_list,	// list of lists of integer vectors giving groupings of models to average interaction tensors
		Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> g0, float k, bool gradient)
{
	// calculate inter nuclear vectors
	auto r_arrays = coord_array_to_r_array(coord_array, atomId_pairs);

	// calculate dipole-dipole interaction tensors [and their derivates]
	auto [d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, gradient);

	// calculate norm squared for different groupings of dipole-dipole interaction tensors
	auto [g_list, g_list_grad] = d_arrays_to_g(d_arrays, grouping_list, gradient);

	auto g_matrix = vectorOfVectors_to_Matrix(g_list);

	// calculate energies from the norm squared values
	auto [energy_matrix, energy_matrix_grad] = g_to_energy(g_matrix, g0, k, gradient);
	std::cout << "energy_matrix" << std::endl << energy_matrix << std::endl;

	// return the sum of all the individual restraint energies
	float sum = energy_matrix.sum();
	std::cout << "sum" << std::endl << sum << std::endl;

	//add derivates
	if(gradient){
		// TODO calculate de/dd = de/dg * dg/dd for all individual interaction tensor components





		// TODO sum the contributions from the different norm squared values

		//TODO

		std::vector<Eigen::MatrixX3<float>> gradients;
		gradients.reserve(coord_array.size());
		for (int i = 0; i < coord_array.size(); ++i) { // seq_len(dim(d_energy_d_r_array)[1])

		//TODO
		}



		return {sum, gradient};
	}else{
		return {sum, std::vector<Eigen::MatrixX3<float>>{}};
	}
}


Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>
KEnRef::coord_array_to_g(
		std::vector<Eigen::MatrixX3<float>> coord_array,	//Every vector item is an Nx3 Matrix representing atom coordinates of a model.
		std::vector<std::tuple<int, int>> atomId_pairs, 	// Matrix with each row having the indices of an atom pair (first dimension in `coord_array` matrices)
		std::vector<std::vector<std::vector<int>>> grouping_list)	// list of lists of integer vectors giving groupings of models to average interaction tensors
{
	// calculate internuclear vectors
	auto r_arrays = coord_array_to_r_array(coord_array, atomId_pairs);

	// calculate dipole-dipole interaction tensors [and their derivatives]
	auto [d_arrays, d_arrays_grad] = r_array_to_d_array(r_arrays, false); //TODO AMR Why this does not have default value for gradient

	// calculate norm squared for different groupings of dipole-dipole interaction tensors
	//		g_list <- lapply(grouping_list, function(grouping) d_array_to_g(d_array, grouping, gradient=FALSE))
	auto [g_list, ignore] = d_arrays_to_g(d_arrays, grouping_list);

	return vectorOfVectors_to_Matrix(g_list);
}
