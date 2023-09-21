/*
 * gmxwrapper.cpp
 *
 *  Created on: Nov 29, 2022
 *      Author: amr
 */
#include <iostream>
#include <filesystem>
#include <string>
#include <map>
//#include <optional>
//#include "gromacs/mdtypes/iforceprovider.h"
//#include "../gmxinterface/gmxkenrefinitializer.h"

#include "gmxpre.h"

#include "mdrun/mdrun_main.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"

#include "../gmxinterface/gmxwrapper.h"
#include "../gmxinterface/gmxkenrefrunner.h"
#include "../gmxinterface/KEnRefMDModule.h"
#include "../core/IoUtils.h"
// Remember to put all local imports (especially those which reference Eigen) after
// all Gromacs imports; to avoid compiler error "reference to 'real' is ambiguous"
#include "../core/ensembleutils.h"


GMXWrapper::GMXWrapper() = default;

GMXWrapper::~GMXWrapper() = default;

GMXWrapper::GMXWrapper(const GMXWrapper &other) = default;

GMXWrapper::GMXWrapper(GMXWrapper &&other)  noexcept {}

//GMXWrapper& GMXWrapper::operator=(const GMXWrapper &other) {}
//GMXWrapper& GMXWrapper::operator=(GMXWrapper &&other) {}

#if GROMACS_WRAPPER && ! _TEST_

//! Initializer for a module that defaults to nice level zero. //copied from mdrun_main.cpp
void initSettingsNoNice(gmx::CommandLineModuleSettings* settings)
{
	settings->setDefaultNiceLevel(0);
}

//! Post-On Selef-Test
void POST(){
//	GmxKEnRefRunner gmxRunner;

	Eigen::Matrix3f matrix3f{{1, 2, 3}, {4, 5.5, 6}, {7, 8, 9}};
	std::map<std::string, float> test_map;
	EnsembleUtils::get_ensemble_data(test_map);



	std::vector<std::vector<std::vector<int>>> toy_grouping_list {{{0, 1, 2, 3}}, {{0, 1}, {2, 3}}, {{0}, {1}, {2}, {3}}};
	Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 3> toy_r_mat(4, 3);
	toy_r_mat <<
			0.848351683690084, -0.529433112659379, 0,
			0.966177888683851, 0.257876496444355, 0,
			0.966177888683851, -0.257876496444355, 0,
			0.848351683690084, 0.529433112659379, 0;

	std::cout << "toy_r_mat" << std::endl << toy_r_mat << std::endl;

	auto [toy_d_array, toy_d_array_grad] = KEnRef<KEnRef_Real_t>::r_array_to_d_array(toy_r_mat, true);
	std::cout << "toy_d_array" << std::endl << toy_d_array << std::endl;
	std::cout << "toy_d_array_grad" << std::endl << toy_d_array_grad<< std::endl;

	std::vector<Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 5>> toy_d_array_vec;
	toy_d_array_vec.reserve((toy_d_array.rows()));
//	std::cout << "toy_d_array_vec size\t" << toy_d_array_vec.size() << std::endl;
	for(int i = 0; i < toy_d_array.rows(); i++){
		toy_d_array_vec.emplace_back(toy_d_array.row(i));
	};

	for (int gg = 0; gg < toy_grouping_list.size(); gg++) {
		std::cout << "Calculating g" << gg+1 << std::endl;
		auto [toy_g_array, toy_g_array_grad] = KEnRef<KEnRef_Real_t>::d_array_to_g(toy_d_array_vec, toy_grouping_list[gg], true);
		std::cout << "toy_g_array" << std::endl << toy_g_array << std::endl;
		std::cout << "toy_g_array_grad" << std::endl;
		for(const auto& matrix: toy_g_array_grad){
			std::cout << matrix << /*std::endl <<*/ std::endl;
		}
		std::cout << "----------" << std::endl;
	}

	/////////////////////////////////////////////////////////////////////////

	Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 3> model1(5, 3);
	model1 <<
			32.708, 53.484, 20.701,
			32.284, 52.123, 22.636,
			31.277, 51.654, 21.284,
			31.852, 49.646, 22.312,
			32.854, 49.716, 20.812;
	Eigen::Matrix<KEnRef_Real_t, Eigen::Dynamic, 3> model2(5, 3);
	model2 <<
			32.733, 52.960, 22.152,
			33.130, 51.220, 23.736,
			31.878, 50.694, 22.613,
			33.471, 49.251, 21.415,
			34.819, 49.854, 22.481;
	std::vector<Eigen::Matrix<KEnRef_Real_t , Eigen::Dynamic, 3>> eros3_sub_coord{model1, model2};
	std::vector<std::tuple<int, int>> eros3_sub_atom_idPairs{{0, 1}, {0, 2}, {0, 3}, {0, 4}};
	std::vector<std::vector<std::vector<int>>> eros3_grouping_list {{{0}, {1}}, {{0, 1}}};

	auto eros3_sub_r_array = KEnRef<KEnRef_Real_t>::coord_array_to_r_array(eros3_sub_coord, eros3_sub_atom_idPairs);
	std::cout << "eros3_sub_r_array" << std::endl;
	for (int i = 0; i < eros3_sub_r_array.size(); i++) {
		auto r_array = eros3_sub_r_array[i];
		std::cout << "Model " << i+1 << std::endl << r_array << std::endl;
	}
	auto [eros3_sub_d_array, eros3_sub_d_array_grad] = KEnRef<KEnRef_Real_t>::r_array_to_d_array(eros3_sub_r_array, true);
	std::cout << "eros3_sub_d_array" << std::endl;
	for(int i = 0; i < eros3_sub_d_array.size(); i++){
		auto matrix = eros3_sub_d_array[i];
		std::cout << "eros3_sub_d_array " << i+1 << std::endl << matrix << std::endl;
	}
	std::cout << "eros3_sub_d_array_grad" << std::endl;
	for(int i = 0; i < eros3_sub_d_array_grad.size(); i++){
		auto matrix = eros3_sub_d_array_grad[i];
		std::cout << "eros3_sub_d_array_grad " << i+1 << std::endl << matrix << std::endl;
	}

	std::vector<Eigen::VectorX<KEnRef_Real_t>> g_list;
	for (int gg = 0; gg < eros3_grouping_list.size(); gg++) {
		std::cout << "Calculating eros3_sub g" << gg+1 << std::endl;
		auto [eros3_sub_g_list, eros3_sub_g_list_grad] = KEnRef<KEnRef_Real_t>::d_array_to_g(eros3_sub_d_array, eros3_grouping_list[gg], true); //TODO d_arrays_to_g ???
		std::cout << "eros3_sub_g_list" << std::endl << eros3_sub_g_list << std::endl;
		std::cout << "eros3_sub_g_list_grad" << std::endl;
		for(const auto& matrix: eros3_sub_g_list_grad){
			std::cout << matrix << /*std::endl <<*/ std::endl;
			std::cout << "----------" << std::endl;
		}
		std::cout << "=============" << std::endl;
		g_list.emplace_back(eros3_sub_g_list);
	}

	auto g_matrix = KEnRef<KEnRef_Real_t>::vectorOfVectors_to_Matrix(g_list);
	std::cout << "g_matrix" << std::endl << g_matrix << std::endl;


	auto eros3_sub_g = KEnRef<KEnRef_Real_t>::coord_array_to_g(eros3_sub_coord, eros3_sub_atom_idPairs, eros3_grouping_list);
	std::cout << "eros3_sub_g" << std::endl << eros3_sub_g << std::endl;

	std::vector<Eigen::MatrixX3<KEnRef_Real_t>> m1_twice = {eros3_sub_coord[0], eros3_sub_coord[0]};
	//get g values using M1 twice
	auto eros3_sub_1_g = KEnRef<KEnRef_Real_t>::coord_array_to_g(m1_twice, eros3_sub_atom_idPairs, eros3_grouping_list);
	std::cout << "eros3_sub_1_g" << std::endl << eros3_sub_1_g << std::endl;


	auto[eros3_sub_energy, eros3_sub_energy_grad] = KEnRef<KEnRef_Real_t>::coord_array_to_energy(eros3_sub_coord, eros3_sub_atom_idPairs, eros3_grouping_list, eros3_sub_1_g, 1.0, true);


	//	auto [g_list, eros3_sub_g_list_grad] = KEnRef::d_array_to_g(eros3_sub_d_array, eros3_grouping_list, true);
	//	KEnRef::g_to_energy(g_matrix, eros3_sub_1_g, 1.0, true);

	//	exit(0);

}


int main(int argc, char* argv[]) {// This should be the main function that loads GMXWrapper (which should replace GROMACS main class)

	POST();

	//		KEnRefMDModule kEnRefMDModule;
	//		return 0;

	//next line is copied from mdrun_main.cpp
	return gmx::CommandLineModuleManager::runAsMainCMainWithSettings(
			argc, argv, &gmx::/*KEnRef_*/gmx_mdrun, &initSettingsNoNice);

}

#endif
