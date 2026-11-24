/*
 * KEnRefSettingsOptions.h
 *
 * KEnRef's settings, shared by bin/KEnRef and the GROMACS plugin.
 */

#ifndef KENREFSETTINGSOPTIONS_H_
#define KENREFSETTINGSOPTIONS_H_

#include <CLI11/CLI11.hpp>

namespace kenref {

//! Declares KEnRef's settings on \p app, each also read from its KENREF_* environment variable
void addSettingsOptions(CLI::App *app);

/*! \brief Reads the settings from the KENREF_* environment variables and KENREF_PARAMS.
 *
 * For the GROMACS plugin, which has no KEnRef command line. Throws std::runtime_error when they
 * are incomplete or wrong.
 */
void readSettingsFromEnvironment();

//! Uses the atom-name mapping file as the reference file when no reference file was given
void useAtomNameMappingAsDefaultReference();

} // namespace kenref

#endif /* KENREFSETTINGSOPTIONS_H_ */
