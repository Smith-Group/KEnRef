/*
 * KEnRefMoleculeGraph.cpp
 *
 *  See the header for why a bond graph is required rather than a distance heuristic.
 */

#include "gmxinterface/KEnRefMoleculeGraph.h"

#include <algorithm>
#include <cmath>
#include <deque>
#include <unordered_set>

#include <gromacs/topology/topology.h>
#include <gromacs/topology/idef.h>
#include <gromacs/topology/ifunc.h>
#include <gromacs/utility/fatalerror.h>

namespace {

/*! \brief Whether an interaction type connects atoms rigidly enough to walk along.
 *
 * The same selection GROMACS uses to build its own molecule graph (pbcutil/mshift.cpp): everything
 * flagged IF_CHEMBOND, plus constraints and SETTLE. The latter two matter because a constrained or
 * rigid-water bond is a real connection that carries no IF_CHEMBOND flag, and leaving it out would
 * split a molecule into fragments that then get re-imaged independently. */
inline bool kenrefIsConnectingInteraction(int ftype) {
    if (ftype == F_CONSTR || ftype == F_CONSTRNC || ftype == F_SETTLE)
        return true;
    return (interaction_function[ftype].flags & IF_CHEMBOND) != 0;
}

/*! \brief Append the atom pairs of one interaction list, shifted into global indexing.
 *
 * SETTLE is the exception that has to be special-cased: its parameter list holds ONE entry per rigid
 * water with three atoms (O, H, H) rather than a sequence of pairs, so the generic "consecutive atoms
 * are bonded" reading would connect one water's hydrogen to the next water's oxygen and fuse every
 * water in the system into a single fragment. */
void kenrefAppendBonds(std::vector<std::pair<int, int>> &edges, int ftype, const InteractionList &il,
                       int atomOffset) {
    const int nral = NRAL(ftype);
    if (nral < 2)
        return;
    for (int i = 0; i < il.size(); i += 1 + nral) {
        const int *a = il.iatoms.data() + i + 1;
        if (ftype == F_SETTLE) {
            edges.emplace_back(atomOffset + a[0], atomOffset + a[1]);
            edges.emplace_back(atomOffset + a[0], atomOffset + a[2]);
        } else {
            for (int k = 1; k < nral; ++k)
                edges.emplace_back(atomOffset + a[k - 1], atomOffset + a[k]);
        }
    }
}

} // namespace

void KEnRefMoleculeGraph::build(const gmx_mtop_t &mtop, const std::vector<int> &atomsOfInterest0) {
    atoms_.clear();
    globalToRow_.clear();
    treeEdges_.clear();
    fragmentRoots_.clear();
    if (atomsOfInterest0.empty())
        return;

    const std::unordered_set<int> wanted(atomsOfInterest0.begin(), atomsOfInterest0.end());

    /* Walk the molecule blocks in global atom order. Each block repeats one molecule type `nmol`
     * times, and the atoms of every copy are contiguous, so the global index of a copy's local atom
     * is simply its offset plus the local index. Only copies that actually contain an atom of
     * interest are kept -- that is what stops a protein-in-water system from dragging in the water. */
    std::vector<std::pair<int, int>> edges;
    int atomOffset = 0;
    for (const gmx_molblock_t &molb : mtop.molblock) {
        const gmx_moltype_t &molt = mtop.moltype[molb.type];
        const int natomsPerMol = molt.atoms.nr;
        for (int mol = 0; mol < molb.nmol; ++mol) {
            const int base = atomOffset + mol * natomsPerMol;

            bool molIsWanted = false;
            for (int a = 0; a < natomsPerMol && !molIsWanted; ++a)
                molIsWanted = wanted.count(base + a) != 0;
            if (!molIsWanted)
                continue;

            for (int a = 0; a < natomsPerMol; ++a)
                atoms_.push_back(base + a);
            for (int ftype = 0; ftype < F_NRE; ++ftype)
                if (kenrefIsConnectingInteraction(ftype))
                    kenrefAppendBonds(edges, ftype, molt.ilist[ftype], base);
        }
        atomOffset += molb.nmol * natomsPerMol;
    }

    if (atoms_.empty())
        return;
    std::sort(atoms_.begin(), atoms_.end());
    for (std::size_t r = 0; r < atoms_.size(); ++r)
        globalToRow_[atoms_[r]] = static_cast<int>(r);

    // Every atom of interest must be covered; otherwise the caller's set spans a molecule we missed.
    for (const int g : atomsOfInterest0) {
        if (globalToRow_.find(g) == globalToRow_.end()) {
            gmx_fatal(FARGS,
                      "KEnRef could not find global atom %d in any molecule of the topology while "
                      "building the connectivity used to repair periodic images. The restrained atom "
                      "selection does not match the topology this simulation was built from.",
                      g);
        }
    }

    buildSpanningForest(edges);
}

void KEnRefMoleculeGraph::buildFromEdgesForTesting(const std::vector<int> &atoms0,
                                                   const std::vector<std::pair<int, int>> &edges) {
    atoms_ = atoms0;
    globalToRow_.clear();
    treeEdges_.clear();
    fragmentRoots_.clear();
    std::sort(atoms_.begin(), atoms_.end());
    for (std::size_t r = 0; r < atoms_.size(); ++r)
        globalToRow_[atoms_[r]] = static_cast<int>(r);
    buildSpanningForest(edges);
}

void KEnRefMoleculeGraph::buildSpanningForest(const std::vector<std::pair<int, int>> &edges) {
    std::vector<std::vector<int>> adjacency(atoms_.size());
    for (const auto &[ga, gb] : edges) {
        const auto ia = globalToRow_.find(ga);
        const auto ib = globalToRow_.find(gb);
        if (ia == globalToRow_.end() || ib == globalToRow_.end())
            continue;               // an edge of a molecule copy we did not keep
        adjacency[ia->second].push_back(ib->second);
        adjacency[ib->second].push_back(ia->second);
    }

    /* Breadth-first spanning forest. Emitting edges in BFS discovery order guarantees a parent is
     * always placed before its children, which is exactly the precondition makeWhole() needs. */
    std::vector<char> visited(atoms_.size(), 0);
    treeEdges_.reserve(atoms_.size());
    for (std::size_t start = 0; start < atoms_.size(); ++start) {
        if (visited[start])
            continue;
        visited[start] = 1;
        fragmentRoots_.push_back(static_cast<int>(start));
        std::deque<int> queue{static_cast<int>(start)};
        while (!queue.empty()) {
            const int parent = queue.front();
            queue.pop_front();
            for (const int child : adjacency[parent]) {
                if (visited[child])
                    continue;
                visited[child] = 1;
                treeEdges_.emplace_back(parent, child);
                queue.push_back(child);
            }
        }
    }
}

void KEnRefMoleculeGraph::makeWhole(CoordsMatrixType<KEnRef_Real_t> &x, const matrix box) const {
    for (const auto &[parent, child] : treeEdges_) {
        KEnRef_Real_t d[3];
        for (int k = 0; k < 3; ++k)
            d[k] = x(child, k) - x(parent, k);

        bool shifted = false;
        for (int k = 2; k >= 0; --k) {           // highest dimension first: triclinic reduction
            const auto boxkk = static_cast<KEnRef_Real_t>(box[k][k]);
            if (boxkk <= KEnRef_Real_t(0))
                continue;
            const KEnRef_Real_t half = boxkk / KEnRef_Real_t(2);
            while (d[k] > half) {
                for (int j = 0; j <= k; ++j)
                    d[j] -= static_cast<KEnRef_Real_t>(box[k][j]);
                shifted = true;
            }
            while (d[k] < -half) {
                for (int j = 0; j <= k; ++j)
                    d[j] += static_cast<KEnRef_Real_t>(box[k][j]);
                shifted = true;
            }
        }
        // Leave already-correct coordinates untouched so whole input stays bit-for-bit unchanged.
        if (shifted)
            for (int k = 0; k < 3; ++k)
                x(child, k) = x(parent, k) + d[k];
    }
}

KEnRef_Real_t KEnRefMoleculeGraph::checkBondLengths(const CoordsMatrixType<KEnRef_Real_t> &x) const {
    KEnRef_Real_t longest = 0;
    for (const auto &[parent, child] : treeEdges_) {
        const KEnRef_Real_t dx = x(child, 0) - x(parent, 0);
        const KEnRef_Real_t dy = x(child, 1) - x(parent, 1);
        const KEnRef_Real_t dz = x(child, 2) - x(parent, 2);
        longest = std::max(longest, std::sqrt(dx * dx + dy * dy + dz * dz));
    }
    return longest;
}
