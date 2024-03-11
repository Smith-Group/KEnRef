/*
 * gmxwrapper.cpp
 *
 *  Created on: Nov 29, 2022
 *      Author: amr
 */
#include <string>
#include <map>
//#include <optional>
//#include "gromacs/mdtypes/iforceprovider.h"
//#include "../gmxinterface/gmxkenrefinitializer.h"
#include "gmxpre.h"
#include "mdrun/mdrun_main.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gmxinterface/gmxwrapper.h"
#include "gmxinterface/gmxkenrefrunner.h"
#include "gmxinterface/KEnRefMDModule.h"
#include "core/IoUtils.h"
// Remember to put all local imports (especially those which reference Eigen) after
// all Gromacs imports; to avoid compiler error "reference to 'real' is ambiguous"
#include "core/ensembleutils.h"


GMXWrapper::GMXWrapper() = default;
GMXWrapper::~GMXWrapper() = default;
GMXWrapper::GMXWrapper(const GMXWrapper &other) = default;
GMXWrapper::GMXWrapper(GMXWrapper &&other)  noexcept  = default;
GMXWrapper& GMXWrapper::operator=(const GMXWrapper &other)  = default;
GMXWrapper& GMXWrapper::operator=(GMXWrapper &&other) noexcept = default;

#if GROMACS_WRAPPER && ! _TEST_

//! Initializer for a module that defaults to nice level zero. //copied from mdrun_main.cpp
void initSettingsNoNice(gmx::CommandLineModuleSettings* settings)
{
	settings->setDefaultNiceLevel(0);
}

//! Post-On Selef-Test
void POST(){
//	GmxKEnRefRunner gmxRunner;
    std::map<std::string, float> test_map;
    EnsembleUtils::get_ensemble_data(test_map);

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
