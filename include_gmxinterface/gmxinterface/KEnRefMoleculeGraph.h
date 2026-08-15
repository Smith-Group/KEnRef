/*
 * KEnRefMoleculeGraph.h
 *
 *  Bond connectivity for the molecules KEnRef restrains, and the make-whole that uses it.
 *
 *  Why this exists: under domain decomposition GROMACS hands out coordinates wrapped into the box, so
 *  a molecule straddling a periodic boundary arrives BROKEN. Every purely distance-based repair fails
 *  for a protein whose radius approaches half the box -- measured on res/sigma/md_single, the restrained
 *  set reaches 2.171 nm from its centroid while half the smallest box vector is 2.165 nm, so "nearest
 *  image to an anchor" picks the wrong image and silently distorts the structure.
 *
 *  The fix is the one both reference implementations use: walk the molecule's BOND graph. Bonds are
 *  ~0.1 nm, so the minimum image of a bond is never ambiguous no matter how large the molecule is.
 *  GROMACS does this in gmx::WholeMoleculeTransform (t_graph / mk_graph, restricted to a single PP
 *  rank); PLUMED does it in ActionAtomistic::makeWhole() via a PLMD::Tree spanning tree over MOLINFO
 *  connectivity. This is the same construction, in GLOBAL atom indexing so it works under DD once the
 *  molecule has been gathered.
 *
 *  The graph covers the WHOLE molecule, not the restrained subset: the guide/sub selections are sparse
 *  (the guide set jumps 116 -> 233), so they have no usable connectivity of their own.
 */

#ifndef KENREFMOLECULEGRAPH_H_
#define KENREFMOLECULEGRAPH_H_

#include <unordered_map>
#include <utility>
#include <vector>

#include <gromacs/math/vectypes.h>

#include "core/KEnRef.h"

struct gmx_mtop_t;

/*! \brief Bond graph + spanning tree over the molecules containing a set of atoms of interest.
 *
 * Built once at setup from the global topology; used every step to make the gathered coordinates
 * whole. All indices it exposes are ZERO-based GLOBAL atom indices, matching gmx::LocalAtomSet. */
class KEnRefMoleculeGraph {
public:
	/*! \brief Build the graph for every molecule that contains at least one atom of interest.
	 *
	 * \param[in] mtop              the global topology
	 * \param[in] atomsOfInterest0  zero-based global indices whose molecules must be made whole
	 *
	 * Whole molecules are pulled in, not just the atoms asked for, because connectivity is what makes
	 * the repair unambiguous. Molecules that contain none of the atoms of interest are ignored, so a
	 * protein in a box of water costs the protein, not the box. */
	void build(const gmx_mtop_t &mtop, const std::vector<int> &atomsOfInterest0);

	/*! \brief Build directly from an explicit atom list and edge list, for tests.
	 *
	 * Same spanning-forest construction as build(), minus the topology parsing, so the repair
	 * algorithm can be pinned on hand-made data with no tpr and no GROMACS topology involved.
	 *
	 * \param[in] atoms0  zero-based global indices, ascending
	 * \param[in] edges   connections as pairs of GLOBAL indices, both of which must be in \p atoms0 */
	void buildFromEdgesForTesting(const std::vector<int> &atoms0,
	                              const std::vector<std::pair<int, int>> &edges);

	//! Zero-based global indices of every atom in the covered molecules, ascending. The set to gather.
	[[nodiscard]] const std::vector<int> &atoms() const { return atoms_; }

	//! Row of \p global0 within atoms(), or -1 when that atom is not covered.
	[[nodiscard]] int rowOf(int global0) const {
		const auto it = globalToRow_.find(global0);
		return it == globalToRow_.end() ? -1 : it->second;
	}

	//! Whether build() found any connectivity at all.
	[[nodiscard]] bool isBuilt() const { return !atoms_.empty(); }

	//! Number of separately-connected fragments covered (each is re-imaged independently).
	[[nodiscard]] std::size_t numFragments() const { return fragmentRoots_.size(); }

	/*! \brief Make gathered coordinates whole, in place.
	 *
	 * \param[in,out] x    one row per entry of atoms(), in that order, in the same units as \p box
	 * \param[in]     box  GROMACS lower-triangular box
	 *
	 * Each atom is placed at the periodic image nearest its PARENT in the spanning tree, walking
	 * parents-before-children. Every step of that walk crosses one bond, so the image is unambiguous.
	 *
	 * A no-op, bit for bit, on input that is already whole: each bond difference is then already the
	 * minimum image, no shift is computed and no coordinate is written. That is what makes it safe to
	 * apply unconditionally rather than only under domain decomposition. */
	void makeWhole(CoordsMatrixType<KEnRef_Real_t> &x, const matrix box) const;

	/*! \brief Longest bond found by the last makeWhole(), for diagnostics.
	 *
	 * A bond much longer than a chemical bond means the graph and the coordinates disagree, which is
	 * the one way this construction can go wrong. */
	[[nodiscard]] KEnRef_Real_t checkBondLengths(const CoordsMatrixType<KEnRef_Real_t> &x) const;

private:
	/*! \brief Turn an edge list into the BFS spanning forest, given atoms_/globalToRow_ are set.
	 *
	 * Shared by build() and buildFromEdgesForTesting() so the tested construction is the real one. */
	void buildSpanningForest(const std::vector<std::pair<int, int>> &edges);

	//! Global atom indices covered, ascending.
	std::vector<int> atoms_;
	//! global atom index -> row in atoms_.
	std::unordered_map<int, int> globalToRow_;
	//! Spanning-tree edges as (parent row, child row), ordered parents-before-children.
	std::vector<std::pair<int, int>> treeEdges_;
	//! Row of the first atom of each connected fragment; these stay where GROMACS put them.
	std::vector<int> fragmentRoots_;
};

#endif /* KENREFMOLECULEGRAPH_H_ */
