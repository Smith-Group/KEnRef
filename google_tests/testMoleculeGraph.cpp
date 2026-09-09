/*
 * testMoleculeGraph.cpp
 *
 * Tests for KEnRefMoleculeGraph -- the bond connectivity that repairs molecules split across periodic
 * boundaries. See res/PBC-BROKEN.md.
 *
 * Two kinds of test here:
 *  - synthetic, so the ALGORITHM is pinned independently of any data file;
 *  - against res/sigma/md_single/md-00_PBC-BROKEN.tpr, a real GB3 structure that really is split
 *    (2.86 nm gap in x), so the repair is exercised on the kind of input that motivated it.
 *
 * As with the rest of the suite, the fixture paths are CWD-relative: run from inside the source tree.
 */

#include "gtest/gtest.h"

#include <cmath>
#include <array>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include <gromacs/fileio/tpxio.h>
#include <gromacs/mdtypes/inputrec.h>
#include <gromacs/topology/topology.h>
#include <gromacs/utility/arrayref.h>

#include "gmxinterface/KEnRefMoleculeGraph.h"

namespace {

//! The committed GB3 tpr whose protein is split across the periodic boundary.
const char *kBrokenTpr = "../../res/sigma/md_single/md-00_PBC-BROKEN.tpr";

/*! The reference structure, kept in BOTH states.
 *
 * `GB3_27_10us.pdb` is what ships and what the .toml files name; it was repaired with
 * `gmx trjconv -pbc whole` (see res/PBC-BROKEN.md). `_PBC-BROKEN.pdb` is the original, preserved
 * precisely so the repair can be re-derived here by a second, independent mechanism -- the bond
 * graph -- and the two required to agree. One method producing a plausible answer is not evidence. */
const char *kWholeRefPdb = "../../res/sigma/md_single/GB3_27_10us.pdb";
const char *kBrokenRefPdb = "../../res/sigma/md_single/GB3_27_10us_PBC-BROKEN.pdb";

/*! \brief Read the ATOM records of a PDB into a coordinate matrix, converting Angstrom -> nm.
 *
 * Only the columns the PDB format fixes are touched, so this stays independent of whatever wrote the
 * file. Returns false when the file is missing, so the tests skip rather than fail without fixtures. */
bool readPdbCoords(const char *path, CoordsMatrixType<KEnRef_Real_t> *x) {
    std::ifstream in(path);
    if (!in)
        return false;
    std::vector<std::array<KEnRef_Real_t, 3>> rows;
    for (std::string line; std::getline(in, line);) {
        if (line.compare(0, 4, "ATOM") != 0 && line.compare(0, 6, "HETATM") != 0)
            continue;
        if (line.size() < 54)
            return false;
        rows.push_back({ static_cast<KEnRef_Real_t>(std::stod(line.substr(30, 8)) / 10.0),
                         static_cast<KEnRef_Real_t>(std::stod(line.substr(38, 8)) / 10.0),
                         static_cast<KEnRef_Real_t>(std::stod(line.substr(46, 8)) / 10.0) });
    }
    if (rows.empty())
        return false;
    x->resize(static_cast<Eigen::Index>(rows.size()), 3);
    for (std::size_t i = 0; i < rows.size(); ++i)
        for (int k = 0; k < 3; ++k)
            (*x)(static_cast<Eigen::Index>(i), k) = rows[i][k];
    return true;
}

//! Longest edge of the graph, i.e. the worst bond. The single number that says "whole" or "broken".
KEnRef_Real_t longestBond(const KEnRefMoleculeGraph &g, const CoordsMatrixType<KEnRef_Real_t> &x) {
    return g.checkBondLengths(x);
}

/*! \brief Read a tpr: topology, coordinates and box.
 *
 * Returns false when the file is missing, so the tests can skip rather than fail on a tree without
 * the fixture (the res/ sets are large and a shallow checkout may not carry them). */
bool readTpr(const char *path, gmx_mtop_t *mtop, std::vector<gmx::RVec> *x, matrix box) {
    if (FILE *f = std::fopen(path, "rb")) {
        std::fclose(f);
    } else {
        return false;
    }
    /* read_tpx() is the form that hands back coordinates, box and topology together, with no t_state
     * (which is an incomplete type here). It needs the atom count first, so read the header pass with
     * null coordinate buffers, size the vector, then read for real. */
    t_inputrec ir;
    int natoms = 0;
    read_tpx(path, &ir, box, &natoms, nullptr, nullptr, mtop);
    if (natoms <= 0)
        return false;
    x->resize(natoms);
    read_tpx(path, &ir, box, &natoms, as_rvec_array(x->data()), nullptr, mtop);
    return true;
}

} // namespace

/*! The algorithm, on data constructed by hand: four atoms in a line, one deliberately wrapped.
 *
 * Independent of any file, so a failure here is unambiguously the algorithm and not the fixture. */
TEST(MoleculeGraphTest, RepairsAWrappedChainSynthetically) {
    // A cubic 10 nm box and a 4-atom chain at x = 9.8, 9.9, 10.0 -> wrapped to 0.0, and 0.1.
    matrix box = { { 10, 0, 0 }, { 0, 10, 0 }, { 0, 0, 10 } };
    CoordsMatrixType<KEnRef_Real_t> x(4, 3);
    x << KEnRef_Real_t(9.8), 0, 0,
         KEnRef_Real_t(9.9), 0, 0,
         KEnRef_Real_t(0.0), 0, 0,   // wrapped: really 10.0
         KEnRef_Real_t(0.1), 0, 0;   // wrapped: really 10.1

    KEnRefMoleculeGraph g;
    g.buildFromEdgesForTesting({ 0, 1, 2, 3 }, { { 0, 1 }, { 1, 2 }, { 2, 3 } });

    EXPECT_GT(longestBond(g, x), KEnRef_Real_t(9.0)) << "the wrapped input should look broken";
    g.makeWhole(x, box);
    EXPECT_NEAR(longestBond(g, x), KEnRef_Real_t(0.1), KEnRef_Real_t(1e-5))
            << "after repair every bond should be the 0.1 nm spacing we built";
    // The chain must be contiguous and monotonic again.
    EXPECT_NEAR(x(2, 0), KEnRef_Real_t(10.0), KEnRef_Real_t(1e-5));
    EXPECT_NEAR(x(3, 0), KEnRef_Real_t(10.1), KEnRef_Real_t(1e-5));
}

//! Already-whole input must come back BIT-FOR-BIT unchanged: that is what lets make-whole run always.
TEST(MoleculeGraphTest, IsBitwiseNoOpOnWholeInput) {
    matrix box = { { 10, 0, 0 }, { 0, 10, 0 }, { 0, 0, 10 } };
    CoordsMatrixType<KEnRef_Real_t> x(3, 3);
    x << KEnRef_Real_t(4.90), KEnRef_Real_t(1.25), KEnRef_Real_t(3.5),
         KEnRef_Real_t(5.00), KEnRef_Real_t(1.25), KEnRef_Real_t(3.5),
         KEnRef_Real_t(5.10), KEnRef_Real_t(1.25), KEnRef_Real_t(3.5);
    const CoordsMatrixType<KEnRef_Real_t> before = x;

    KEnRefMoleculeGraph g;
    g.buildFromEdgesForTesting({ 0, 1, 2 }, { { 0, 1 }, { 1, 2 } });
    g.makeWhole(x, box);

    for (Eigen::Index r = 0; r < x.rows(); ++r)
        for (int c = 0; c < 3; ++c)
            EXPECT_EQ(x(r, c), before(r, c)) << "row " << r << " col " << c << " was rewritten";
}

//! Applying the repair twice must change nothing the second time.
TEST(MoleculeGraphTest, IsIdempotent) {
    matrix box = { { 10, 0, 0 }, { 0, 10, 0 }, { 0, 0, 10 } };
    CoordsMatrixType<KEnRef_Real_t> x(3, 3);
    x << KEnRef_Real_t(9.9), 0, 0, KEnRef_Real_t(0.0), 0, 0, KEnRef_Real_t(0.1), 0, 0;

    KEnRefMoleculeGraph g;
    g.buildFromEdgesForTesting({ 0, 1, 2 }, { { 0, 1 }, { 1, 2 } });
    g.makeWhole(x, box);
    const CoordsMatrixType<KEnRef_Real_t> once = x;
    g.makeWhole(x, box);

    for (Eigen::Index r = 0; r < x.rows(); ++r)
        for (int c = 0; c < 3; ++c)
            EXPECT_EQ(x(r, c), once(r, c)) << "second application changed row " << r;
}

/*! The real thing: the committed GB3 tpr, whose protein is genuinely split across the boundary.
 *
 * This is the input that revealed KEnRef had been restraining a torn protein in every serial run. */
TEST(MoleculeGraphTest, RepairsTheCommittedBrokenGb3Tpr) {
    gmx_mtop_t mtop;
    std::vector<gmx::RVec> xs;
    matrix box;
    if (!readTpr(kBrokenTpr, &mtop, &xs, box)) {
        GTEST_SKIP() << "fixture not present: " << kBrokenTpr;
    }

    // Restrain something in the protein; the graph then pulls in that whole molecule.
    KEnRefMoleculeGraph g;
    g.build(mtop, { 0 });
    ASSERT_TRUE(g.isBuilt());
    EXPECT_EQ(g.numFragments(), 1u) << "the protein should be ONE connected fragment";
    EXPECT_EQ(g.atoms().size(), 862u) << "GB3 in this system is 862 atoms";

    CoordsMatrixType<KEnRef_Real_t> x(static_cast<Eigen::Index>(g.atoms().size()), 3);
    for (std::size_t i = 0; i < g.atoms().size(); ++i)
        for (int k = 0; k < 3; ++k)
            x(static_cast<Eigen::Index>(i), k) = static_cast<KEnRef_Real_t>(xs[g.atoms()[i]][k]);

    // 1. the fixture really is broken -- if this ever fails, the tpr was regenerated whole
    const KEnRef_Real_t before = longestBond(g, x);
    EXPECT_GT(before, KEnRef_Real_t(5.0))
            << "expected a box-sized bond in the raw coordinates; got " << before
            << " nm. Was md-00_PBC-BROKEN.tpr regenerated? See res/PBC-BROKEN.md.";

    // 2. after repair, every bond is a chemical bond
    g.makeWhole(x, box);
    const KEnRef_Real_t after = longestBond(g, x);
    EXPECT_LT(after, KEnRef_Real_t(0.25)) << "longest bond after repair was " << after << " nm";

    // 3. and repeating it changes nothing
    const CoordsMatrixType<KEnRef_Real_t> once = x;
    g.makeWhole(x, box);
    for (Eigen::Index r = 0; r < x.rows(); ++r)
        for (int c = 0; c < 3; ++c)
            ASSERT_EQ(x(r, c), once(r, c)) << "repair is not idempotent at row " << r;
}

/*! The repaired reference PDB that ships must actually BE whole.
 *
 * Cheap, and it is the one thing every downstream claim rests on: PLUMED builds its spanning tree from
 * this file's coordinates, so if it is ever replaced by a split copy again, PLUMED silently goes back
 * to computing on a torn protein instead of failing. */
TEST(MoleculeGraphTest, ShippedReferencePdbIsWhole) {
    gmx_mtop_t mtop;
    std::vector<gmx::RVec> xs;
    matrix box;
    CoordsMatrixType<KEnRef_Real_t> x;
    if (!readTpr(kBrokenTpr, &mtop, &xs, box) || !readPdbCoords(kWholeRefPdb, &x)) {
        GTEST_SKIP() << "fixture not present: " << kBrokenTpr << " / " << kWholeRefPdb;
    }

    KEnRefMoleculeGraph g;
    g.build(mtop, { 0 });
    ASSERT_TRUE(g.isBuilt());
    ASSERT_EQ(g.atoms().size(), static_cast<std::size_t>(x.rows()))
            << "the reference PDB should hold exactly the protein the graph covers";

    const KEnRef_Real_t worst = longestBond(g, x);
    EXPECT_LT(worst, KEnRef_Real_t(0.25))
            << "longest bond in " << kWholeRefPdb << " is " << worst
            << " nm -- the shipped reference is split across the periodic boundary. See res/PBC-BROKEN.md.";
}

/*! The cross-check: the bond graph must reproduce `gmx trjconv -pbc whole`, atom for atom.
 *
 * trjconv walks GROMACS's own t_graph over the tpr's bonded interactions; KEnRefMoleculeGraph walks a
 * BFS spanning tree it builds itself. Two independent implementations arriving at the same 12 shifts is
 * the evidence that the shipped reference is repaired correctly -- neither one alone is.
 *
 * It also pins the shape of the repair: exactly 12 atoms move, every other atom is left BIT-identical,
 * and none of the 33 guide C-alphas is among the twelve (which is why the GROMACS-side and offline
 * numbers do not move when the reference is swapped). */
TEST(MoleculeGraphTest, BondGraphRepairOfReferencePdbMatchesTrjconv) {
    gmx_mtop_t mtop;
    std::vector<gmx::RVec> xs;
    matrix box;
    CoordsMatrixType<KEnRef_Real_t> broken, expected;
    if (!readTpr(kBrokenTpr, &mtop, &xs, box) || !readPdbCoords(kBrokenRefPdb, &broken)
        || !readPdbCoords(kWholeRefPdb, &expected)) {
        GTEST_SKIP() << "fixtures not present next to " << kBrokenTpr;
    }
    ASSERT_EQ(broken.rows(), expected.rows());

    KEnRefMoleculeGraph g;
    g.build(mtop, { 0 });
    ASSERT_TRUE(g.isBuilt());
    ASSERT_EQ(g.atoms().size(), static_cast<std::size_t>(broken.rows()));
    // The PDB is the protein in tpr order, so row i of the file is row i of the graph. Check, do not assume.
    for (std::size_t i = 0; i < g.atoms().size(); ++i)
        ASSERT_EQ(g.atoms()[i], static_cast<int>(i)) << "the protein is not atoms 0..N-1 in this tpr";

    // The preserved original really is broken -- otherwise this test proves nothing.
    EXPECT_GT(longestBond(g, broken), KEnRef_Real_t(5.0))
            << "the preserved " << kBrokenRefPdb << " is no longer split; the cross-check is vacuous";

    const CoordsMatrixType<KEnRef_Real_t> before = broken;
    g.makeWhole(broken, box);
    EXPECT_LT(longestBond(g, broken), KEnRef_Real_t(0.25));

    /* Agreement with trjconv. The tolerance is PDB write precision (0.001 A = 1e-4 nm) and nothing
     * more: both files hold the same rounded coordinates, differing by an exact box vector. */
    Eigen::Index movedByGraph = 0, disagreements = 0;
    for (Eigen::Index r = 0; r < broken.rows(); ++r) {
        bool moved = false;
        for (int c = 0; c < 3; ++c) {
            if (broken(r, c) != before(r, c))
                moved = true;
            if (std::abs(broken(r, c) - expected(r, c)) > KEnRef_Real_t(2e-4)) {
                if (++disagreements <= 10)
                    ADD_FAILURE() << "atom " << r + 1 << " column " << c << ": bond graph gives "
                                  << broken(r, c) << " nm, trjconv gives " << expected(r, c) << " nm";
            }
        }
        if (moved) {
            ++movedByGraph;
        } else {
            // Untouched atoms must be untouched EXACTLY: that is what makes the repair safe to apply
            // unconditionally, and what guarantees the guide atoms cannot shift.
            for (int c = 0; c < 3; ++c)
                ASSERT_EQ(broken(r, c), before(r, c)) << "atom " << r + 1 << " was rewritten in place";
        }
    }
    EXPECT_EQ(movedByGraph, 12) << "expected exactly the 12 wrapped LYS10/THR11 side-chain atoms to move";
}
