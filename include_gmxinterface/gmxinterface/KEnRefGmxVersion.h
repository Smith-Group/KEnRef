/*
 * KEnRefGmxVersion.h
 *
 * Reads a GROMACS release out of the version string libgromacs reports at run time, in the encoding of
 * the compile-time GMX_VERSION macro, so the two can be compared.
 *
 * WHY THE TWO MUST MATCH. KEnRef compiles against GROMACS's INTERNAL headers (the MDModules notifiers,
 * SimulationContext, gmx_multisim_t), which carry no ABI promise between releases. The notifiers are the
 * sharpest case: subscribe() is inline template code, so the layout of MDModulesNotifiers is baked into
 * KEnRef when it is compiled -- one std::vector<std::function<...>> per payload type, in list order. The
 * payload list changes between releases (18 types in 2025.4, 19 in 2026.x, 720 B against 816 B), and a
 * module compiled against one list and run inside another does not crash: its subscriptions land on
 * unrelated slots, and haveSubscribers<const gmx_multisim_t*>() -- which decides whether the replicas are
 * kept in lockstep -- silently reads false. Measured in a reduced repro, not inferred.
 *
 * WHY A RELEASE CAN BE SWAPPED UNDER A BUILT BINARY. Releases of one year share a SONAME (2025.4 and 2025.5
 * are both libgromacs_mpi.so.10), and KEnRef links with RUNPATH, which LD_LIBRARY_PATH overrides. Loading
 * another GROMACS module, or reinstalling GROMACS in place, is enough. Across years the SONAME differs and
 * the loader already refuses; this check is for the case it lets through.
 *
 * The configure-time half (cmake/Dependencies.cmake, kenref_check_gromacs_versions) catches the other way
 * to mix them: a GROMACS source checkout moved to another release after GROMACS was built from it.
 *
 * Pure, so `google_tests/testGmxVersion.cpp` tests it without GROMACS or MPI.
 */

#ifndef KENREFGMXVERSION_H_
#define KENREFGMXVERSION_H_

#include <string_view>

namespace kenref {

/*! \brief The release in \p s as GMX_VERSION encodes it (year * 10000 + patch), or -1.
 *
 * \p s is what gmx_version() returns: "2026.3", "2026.3-dev-20260625-1210900145" (a build from git), or
 * "2027-dev-20261005-abcdef" (a development branch, patch 0). Anything after the release is ignored.
 * Returns -1 when \p s does not start with a four-digit year, so a format GROMACS has not used yet is
 * reported as unreadable rather than read as some other release. */
inline long gmxVersionNumber(std::string_view s) {
	auto digits = [&s](std::size_t from) {
		std::size_t to = from;
		while (to < s.size() && s[to] >= '0' && s[to] <= '9')
			++to;
		return to;
	};
	auto value = [&s](std::size_t from, std::size_t to) {
		long v = 0;
		for (std::size_t i = from; i < to; ++i)
			v = v * 10 + (s[i] - '0');
		return v;
	};
	const std::size_t yearEnd = digits(0);
	if (yearEnd != 4)
		return -1;
	long patch = 0;
	if (yearEnd < s.size() && s[yearEnd] == '.') {
		const std::size_t patchEnd = digits(yearEnd + 1);
		if (patchEnd == yearEnd + 1)
			return -1;
		patch = value(yearEnd + 1, patchEnd);
	}
	return value(0, yearEnd) * 10000 + patch;
}

} // namespace kenref

#endif // KENREFGMXVERSION_H_
