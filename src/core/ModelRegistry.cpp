/*
 * ModelRegistry.cpp
 *
 *  The single definition of the model registry's storage.
 *
 *  table() cannot be defined inline in the header. A function-local static inside a class-template
 *  member function has vague linkage: every translation unit that touches the registry emits its own
 *  COMDAT copy, and merging them is left to the toolchain. We cannot rely on that here, because the
 *  library and its consumer are routinely built by DIFFERENT compilers -- PLUMED compiles KEnRefBias
 *  in its own tree, with its own flags, against a prebuilt libkenref_core, and the GROMACS force
 *  provider does the same. When the copies fail to merge, bootstrapModels() registers into one table
 *  while has()/create() consult another, so every model looks unregistered at run time and nothing in
 *  the build says a word. That is precisely what happened on PLUMED's Intel (icpx) CI job, where the
 *  archive is built by g++ but the bias is compiled by icpx.
 *
 *  Defining it out of line here, plus the `extern template` declaration in the header, gives exactly
 *  one registry per process no matter which compilers were involved.
 */

#include "core/ModelRegistry.h"

namespace kenref {

template<typename Real>
std::map<std::string, typename ModelRegistry<Real>::Entry>& ModelRegistry<Real>::table() {
    static std::map<std::string, Entry> t;
    return t;
}

// Emit every member for the one instantiation the build actually uses. KEnRef_Real_t follows
// KENREF_DOUBLE, so this is double in a default build and float when that option is off.
template class ModelRegistry<KEnRef_Real_t>;

} // namespace kenref
