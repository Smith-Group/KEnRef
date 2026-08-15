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
