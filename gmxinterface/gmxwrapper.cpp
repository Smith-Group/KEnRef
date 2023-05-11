/*
 * gmxwrapper.cpp
 *
 *  Created on: Nov 29, 2022
 *      Author: amr
 */
#include  <iostream>
#include <filesystem>
#include <string>
#include <map>
#include <optional>

#include "gmxpre.h"

#include "mdrun/mdrun_main.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/mdtypes/iforceprovider.h"

#include "../gmxinterface/gmxwrapper.h"
#include "../gmxinterface/gmxkenrefinitializer.h"
#include "../gmxinterface/gmxkenrefrunner.h"
#include "../gmxinterface/KEnRefMDModule.h"
#include "../core/IoUtils.h"
// Remember to put all local imports (especially those which reference Eigen) after
// all Gromacs imports; to avoid compiler error "reference to 'real' is ambiguous"
#include "../core/ensembleutils.h"


GMXWrapper::GMXWrapper() {
	// TODO Auto-generated constructor stub

}

GMXWrapper::~GMXWrapper() {
	// TODO Auto-generated destructor stub
}

GMXWrapper::GMXWrapper(const GMXWrapper &other) {
	// TODO Auto-generated constructor stub

}

GMXWrapper::GMXWrapper(GMXWrapper &&other) {
	// TODO Auto-generated constructor stub

}

//GMXWrapper& GMXWrapper::operator=(const GMXWrapper &other) {
//	// TODO Auto-generated method stub
//
//}
//
//GMXWrapper& GMXWrapper::operator=(GMXWrapper &&other) {
//	// TODO Auto-generated method stub
//
//}

#if GROMACS_WRAPPER && ! _TEST_

//! Initializer for a module that defaults to nice level zero. //copied from mdrun_main.cpp
void initSettingsNoNice(gmx::CommandLineModuleSettings* settings)
{
	settings->setDefaultNiceLevel(0);
}

int main(int argc, char* argv[]) {// This should be the main function that loads GMXWrapper (which should replace GROMACS main class)

	std::filesystem::path p1;
	p1 = "./";
	std::cout << p1 << std::endl;

	GmxKEnRefRunner gmxRunner;

	Eigen::Matrix3f matrix3f{{1, 2, 3}, {4, 5.5, 6}, {7, 8, 9}};
	std::map<std::string, float> test_map;
	EnsembleUtils::get_ensemble_data(matrix3f, test_map);


	//		KEnRefMDModule kEnRefMDModule;
	//		return 0;

	//next line is copied from mdrun_main.cpp
	return gmx::CommandLineModuleManager::runAsMainCMainWithSettings(
			argc, argv, &gmx::/*KEnRef_*/gmx_mdrun, &initSettingsNoNice);

}
#endif
