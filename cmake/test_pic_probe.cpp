// Link probe for the position-independence check.
//
// This TU exists only to give CMake something to compile for the MODULE libraries that verify the
// static core archives can be linked INTO A SHARED OBJECT -- see the kenref_pic_probe_* targets in
// cmake/KEnRefCore.cmake. The body is deliberately trivial: the check is performed by the linker, not
// at run time, and --whole-archive means it does not depend on which archive members this file pulls.
#include "core/KEnRef.h"

extern "C" void kenref_pic_probe() {}
