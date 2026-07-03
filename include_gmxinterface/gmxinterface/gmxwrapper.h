/*
 * gmxwrapper.h
 *
 *  Created on: Nov 29, 2022
 *      Author: amr
 */

#ifndef GMXWRAPPER_H_
#define GMXWRAPPER_H_

#include "config/KEnRefConfig.h"
#include "core/KEnRef.h"

class GMXWrapper { //This class should replace GROMACS main class
public:
    GMXWrapper();
    virtual ~GMXWrapper();
    GMXWrapper(const GMXWrapper &other);
    GMXWrapper(GMXWrapper &&other) noexcept;
    GMXWrapper &operator=(const GMXWrapper &other);
    GMXWrapper &operator=(GMXWrapper &&other) noexcept;
};

class Settings {
public:
    Settings() = delete;

    inline static bool debug = false;
    inline static std::string alt_out_path;
    inline static KEnRef_Real_t k = 1;
    inline static KEnRef_Real_t n = 0.25;
    inline static KEnRef_Real_t max_force = 9999;
    // inline static int max_atom_to_read = -1;
    inline static std::string GUIDE_C_ALPHA = "guideC-alpha";
    inline static std::string indexFileName; // = "../../res/cleanstart/KEnRefAtomIndex.ndx";
    inline static std::string refFileName; //"ref.pdb"; // = "../../res/cleanstart/6v5d_step0_for_atomname_mapping.pdb";
    inline static std::string atomName_mapping_fileName; // = "../../res/cleanstart/6v5d_step0_for_atomname_mapping.pdb";
    inline static std::string experimentalDataFileName; // = "../../res/cleanstart/singleton_data_step0_model01.csv";
    inline static std::string experimentalDataFolder; // = "../../res/cleanstart/";
    inline static KEnRef_Real_t proton_mhz = 700.0;
    inline static std::string selected_energy_model; // registry model name (SIGMA / PLATEAUS / RELAX / ...); validated against ModelRegistry
};

#endif /* GMXWRAPPER_H_ */
