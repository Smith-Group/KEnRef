/*
 * ModelBootstrap.cpp
 *
 *  Registers every compile-time-enabled energy model into the ModelRegistry. This is a hand-written
 *  STUB for now; Phase E replaces it with a CMake-GENERATED version that, for each enabled model,
 *  forward-declares and calls its register*Model() under `#if KENREF_ENABLE_<MODEL>`:
 *
 *      #if KENREF_ENABLE_SIGMA
 *      void registerSigmaModel();
 *      #endif
 *      ...
 *      void bootstrapModels() {
 *      #if KENREF_ENABLE_SIGMA
 *          registerSigmaModel();
 *      #endif
 *          ...
 *      }
 *
 *  KEnRefDriver references bootstrapModels(), so the linker keeps every referenced register*Model()
 *  TU without --whole-archive (works for .a and .so). Models are added in Phase B.
 */

#include "core/ModelRegistry.h"

// Until Phase E wires the CMake `KENREF_ENABLE_*` options + generation, default every model ON so
// the hand-written bootstrap registers them. Phase E's generated TU will -D these explicitly.
#ifndef KENREF_ENABLE_SIGMA
#define KENREF_ENABLE_SIGMA 1
#endif
#ifndef KENREF_ENABLE_PLATEAUS
#define KENREF_ENABLE_PLATEAUS 1
#endif

namespace kenref {

#if KENREF_ENABLE_SIGMA
void registerSigmaModel();
#endif
#if KENREF_ENABLE_PLATEAUS
void registerPlateausModel();
#endif

void bootstrapModels() {
#if KENREF_ENABLE_SIGMA
    registerSigmaModel();
#endif
#if KENREF_ENABLE_PLATEAUS
    registerPlateausModel();
#endif
    // RelaxModel lands in Phase F (KENREF_ENABLE_RELAX).
}

} // namespace kenref
