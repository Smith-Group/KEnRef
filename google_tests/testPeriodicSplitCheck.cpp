/*
 * testPeriodicSplitCheck.cpp
 *
 * Tests for kenref::detail::refuseIfSplitByPeriodicity -- the driver's backstop against refining a
 * molecule that is split across a periodic boundary. See res/PBC-BROKEN.md for why that failure is the
 * one this project cares most about: it does not crash, it silently returns a different energy.
 *
 * All synthetic, deliberately. The point of these cases is to pin the RULE independently of any
 * fixture: a whole molecule must pass no matter how large it is, and a wrapped one must be caught no
 * matter how small the displaced fragment is. The predecessor of this check tested SIZE instead of
 * tearing, which is why it refused a perfectly whole GB3 partway through a 500-step run -- so the
 * "large but whole" case below is a regression test with a specific history.
 */

#include "gtest/gtest.h"

#include <cmath>

#include "core/KEnRefDriver.h"

using kenref::detail::refuseIfSplitByPeriodicity;

namespace {

using Real = double;
using Coords = CoordsMatrixType<Real>;
using Box = Eigen::Matrix<Real, 3, 3>;

//! A cubic box of \p edge nm, in the row-vectors-are-lattice-vectors convention the adapters use.
Box cubicBoxNm(Real edge) {
    Box b = Box::Zero();
    b(0, 0) = b(1, 1) = b(2, 2) = edge;
    return b;
}

/*! \brief A chain of \p n atoms evenly spread along x over \p spanAngstrom, centred on the box.
 *
 * Coordinates are Angstrom and the box is nm, which is the mismatch the real adapters hand over. */
Coords chainAlongX(int n, Real spanAngstrom, Real boxNm) {
    Coords x(n, 3);
    const Real start = boxNm * 10 / 2 - spanAngstrom / 2;
    for (int i = 0; i < n; ++i) {
        x(i, 0) = start + spanAngstrom * i / (n - 1);
        x(i, 1) = boxNm * 10 / 2;
        x(i, 2) = boxNm * 10 / 2;
    }
    return x;
}

} // namespace

//! An ordinary compact set in a roomy box is not a tear.
TEST(PeriodicSplitCheckTest, AcceptsACompactSet) {
    EXPECT_NO_THROW(refuseIfSplitByPeriodicity<Real>(chainAlongX(50, 20.0, 6.0), cubicBoxNm(6.0), "sub"));
}

/*! The regression case: a whole molecule that reaches further from its centroid than half the box.
 *
 * This is exactly what the previous size-based check refused. It compared each atom's distance from the
 * SET'S CENTROID against half the smallest box vector -- and a real protein is not symmetric about its
 * centroid, so a lopsided one reaches further than half its own extent. Measured on the committed GB3
 * sets AFTER the bond-graph repair: the restrained atoms reached 21.657 A against a 21.6561 A half-box,
 * 0.004% over, and a 500-step run aborted at about step 120 on a structure that was never torn.
 *
 * Reproduced here the same way it arises: most of the set packed at one end, which pulls the centroid
 * over, and a tail reaching out from it. The whole thing still fits inside one box length -- it is one
 * piece, in one image, and nothing about it is ambiguous. Extent is not tearing.
 */
TEST(PeriodicSplitCheckTest, AcceptsALopsidedMoleculeReachingPastHalfTheBox) {
    const Real boxNm = 4.33;                      // as in the committed GB3 sets: half-box 21.65 A
    const Real halfBox = boxNm * 10 / 2;
    Coords x(120, 3);
    for (int i = 0; i < 100; ++i) {               // a dense blob, which is where the centroid lands
        x(i, 0) = 5.0 + 0.05 * i;
        x(i, 1) = 20.0;
        x(i, 2) = 20.0;
    }
    for (int i = 100; i < 120; ++i) {             // a tail reaching out to ~x = 33 A
        x(i, 0) = 10.0 + 1.2 * (i - 100);
        x(i, 1) = 20.0;
        x(i, 2) = 20.0;
    }
    const Eigen::RowVector3<Real> centroid = x.colwise().mean();
    Real worstRadius = 0;
    for (Eigen::Index i = 0; i < x.rows(); ++i)
        worstRadius = std::max(worstRadius, (x.row(i) - centroid).norm());
    ASSERT_GT(worstRadius, halfBox)
            << "this case must reproduce the old failure: reach " << worstRadius << " A vs half-box "
            << halfBox << " A";
    ASSERT_LT(x.col(0).maxCoeff() - x.col(0).minCoeff(), boxNm * 10)
            << "...while still fitting inside one box length, i.e. genuinely whole";

    EXPECT_NO_THROW(refuseIfSplitByPeriodicity<Real>(x, cubicBoxNm(boxNm), "sub"));
}

//! A whole set placed OUTSIDE the box (as the make-whole repairs leave it) is still whole.
TEST(PeriodicSplitCheckTest, AcceptsAWholeSetSittingOutsideTheBox) {
    Coords x = chainAlongX(50, 20.0, 6.0);
    x.col(0).array() += Real(6.0) * 10;   // shift a whole box along x; nothing about the set changed
    x.col(2).array() -= Real(6.0) * 10 * 3;
    EXPECT_NO_THROW(refuseIfSplitByPeriodicity<Real>(x, cubicBoxNm(6.0), "sub"));
}

//! The thing the check exists for: one fragment wrapped to the far side of the box.
TEST(PeriodicSplitCheckTest, RefusesAWrappedFragment) {
    const Real boxNm = 6.0;
    Coords x = chainAlongX(50, 20.0, boxNm);
    for (int i = 0; i < 5; ++i)
        x(i, 0) += boxNm * 10;            // five atoms re-imaged to the other end: a tear
    EXPECT_THROW(refuseIfSplitByPeriodicity<Real>(x, cubicBoxNm(boxNm), "sub"), std::runtime_error);
}

//! A single wrapped atom is enough -- the restraint has no cutoff, so one wrong atom is a wrong energy.
TEST(PeriodicSplitCheckTest, RefusesASingleWrappedAtom) {
    const Real boxNm = 6.0;
    Coords x = chainAlongX(50, 20.0, boxNm);
    x(17, 0) -= boxNm * 10;
    EXPECT_THROW(refuseIfSplitByPeriodicity<Real>(x, cubicBoxNm(boxNm), "sub"), std::runtime_error);
}

//! The message must name the axis and say a tear was found, so a user can act on it.
TEST(PeriodicSplitCheckTest, ExplainsWhichBoxVectorIsTorn) {
    const Real boxNm = 6.0;
    Coords x = chainAlongX(50, 20.0, boxNm);
    // Rebuild the chain along z instead, then wrap part of it, so the reported axis is 2 and not 0.
    x.col(2) = x.col(0);
    x.col(0).setConstant(boxNm * 10 / 2);
    for (int i = 0; i < 5; ++i)
        x(i, 2) += boxNm * 10;
    try {
        refuseIfSplitByPeriodicity<Real>(x, cubicBoxNm(boxNm), "guide");
        FAIL() << "expected a refusal";
    } catch (const std::runtime_error &e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("KEnRef refuses to continue"), std::string::npos) << msg;
        EXPECT_NE(msg.find("box vector 2"), std::string::npos) << msg;
        EXPECT_NE(msg.find("'guide'"), std::string::npos) << msg;
    }
}

/*! Triclinic boxes must work too -- the committed GB3 sets use a rhombic dodecahedron.
 *
 * This is why the test is done in fractional coordinates: whatever the cell shape, a periodic image
 * differs by a whole number in exactly one fractional component. A tear along the third lattice vector
 * moves an atom in x and y as well as z, and a naive per-axis test in cartesian space would miss it. */
TEST(PeriodicSplitCheckTest, HandlesATriclinicCell) {
    Box box = Box::Zero();                     // the GB3 dodecahedron, nm
    box(0, 0) = 6.12240;
    box(1, 1) = 6.12240;
    box(2, 0) = 3.06120;
    box(2, 1) = 3.06120;
    box(2, 2) = 4.32919;

    Coords x(60, 3);
    for (int i = 0; i < x.rows(); ++i) {       // a compact blob near the middle of the cell
        x(i, 0) = 30.0 + 0.3 * i;
        x(i, 1) = 30.0 + 0.2 * i;
        x(i, 2) = 20.0 + 0.1 * i;
    }
    EXPECT_NO_THROW(refuseIfSplitByPeriodicity<Real>(x, box, "sub"));

    // Now move six atoms by exactly the third lattice vector: still the same molecule, torn image.
    for (int i = 0; i < 6; ++i)
        for (int k = 0; k < 3; ++k)
            x(i, k) += box(2, k) * 10;
    EXPECT_THROW(refuseIfSplitByPeriodicity<Real>(x, box, "sub"), std::runtime_error);
}

//! No box, no periodicity, no check -- an offline caller can pass a zero box and be left alone.
TEST(PeriodicSplitCheckTest, IgnoresAnUndeclaredBox) {
    EXPECT_NO_THROW(refuseIfSplitByPeriodicity<Real>(chainAlongX(50, 200.0, 6.0), Box::Zero(), "sub"));
}
