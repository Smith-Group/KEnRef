/*
 * ensembleutils.h
 *
 *      Author: amr
 */

#ifndef ENSEMBLEUTILS_H_
#define ENSEMBLEUTILS_H_

//#include <vector>
#include <optional>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include "../config/KEnRefConfig.h"
#include <Eigen/Dense>
#include <Eigen/CXX11/Tensor>

#include <typeinfo>
#include "KEnRef.h"

//using Eigen::VectorXf;

class EnsembleUtils{

public:
	static void get_ensemble_data(Eigen::Matrix3Xf coords_array, std::optional<std::map<std::string, float>> ensemble_params);

};

void EnsembleUtils::get_ensemble_data(Eigen::Matrix3Xf coords_array, std::optional<std::map<std::string, float>> ensemble_params){

//	array<float, 0> coords_array_vale;//[Dim, atoms, models]
	auto ensemble_params_value = (ensemble_params ? ensemble_params.value() : (std::map<std::string, float>) {{"proton_mhz", 800}, {"temperature", 298.15}, /*{"tau_c", 4e-9}*/});

	Eigen::Matrix3f matrix3f{{1, 2, 3.}, {4, 5.5, 6}, {7, 8, 9}};
//	std::cout << Eigen::matrix3f << std::endl;
//	std::cout << typeid(Eigen::matrix3f).name() << std::endl;
	std::cout << typeid(Eigen::Matrix4f).name() << std::endl;
	std::cout << typeid(Eigen::Matrix3f).name() << std::endl;

	Eigen::Tensor<double, 3> t(3, 100, 20);
	t.setZero();
	std::cout << typeid(t).name() << std::endl;
	auto chip = t.chip(0, 2);
	std::cout << typeid(chip).name() << "\t" << chip.dim() << "\t" << chip.NumDimensions << std::endl;
	std::cout << chip;
	std::cout << std::endl;
//	Matrix<double, 3, 100> test_mat;
//	test_mat << chip;
//	cout << test_mat << endl;


	Eigen::Vector3d v(1,2,3);
	Eigen::Vector3d w(0,1,2);

	std::cout << "Dot product: " << v.dot(w) << std::endl;
	double dp = v.adjoint()*w; // automatic conversion of the inner product to a scalar
	std::cout << "Dot product via a matrix product: " << dp << std::endl;
	dp = v.dot(w); // automatic conversion
	std::cout << "Dot product into double: " << dp << std::endl;
	std::cout << "Cross product:\n" << v.cross(w) << std::endl;


	Eigen::MatrixX3<float> d {{0, 0, 0}, {1, 1, 1}, {1, 2, 3}, {1, 4, 9}};
//	d << 1 << 2 << 3;
	std::cout << "input matrix:" << std::endl << d << std::endl;
	auto [d_array, d_array_derivatives] = KEnRef::r_array_to_d_array(d);
	std::cout << "returned value:" << std::endl << d_array << std::endl;
}


#endif /* ENSEMBLEUTILS_H_ */
