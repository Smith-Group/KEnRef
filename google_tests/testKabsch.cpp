#include <Eigen/Geometry>
#include <core/kabsch.h>
#include "gtest/gtest.h"
#include "testHelper.h"

// Just some dummy tests
TEST(KabschTestSuite, TestTranslate) {
    CoordsMatrixType<float> before(8, 3), after(8, 3);
    before << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    const Eigen::Vector3<float> &translationVector = Eigen::Vector3<float>{0, 10000, 0};
    A_applied.translate(translationVector);
    after = before;
    std::cout << "after before translation\n" << after.matrix() << std::endl;
    after = Kabsch_Umeyama<float>::applyTransform(A_applied, after);
    std::cout << "after after translation\n" << after.matrix() << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(after, before + translationVector.transpose().colwise().replicate(8));

    auto A_calculated = Kabsch_Umeyama<float>::find3DAffineTransform(before, after, false, true, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_calculated\n" << A_calculated.matrix() << std::endl;
    std::cout << "A_calculated.matrix().inverse()\n" << A_calculated.matrix().inverse() << std::endl;

    std::cout << "difference\n" << (A_calculated.matrix() - A_applied.matrix()) << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(A_calculated.matrix(), A_applied.matrix(), 1e-12);
}

TEST(KabschTestSuite, TestCenterAtoms) {
    CoordsMatrixType<float> after(8, 3);
    after << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;

    after = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(after);
    TestHelper<float>::EXPECT_MATRIX_NEAR(after.colwise().mean(), Eigen::RowVector3<float>{0, 0, 0});
}

TEST(KabschTestSuite, TestRotateTriangle) {
    CoordsMatrixType<float> before(3, 3), expectedAfter(3, 3);
    before <<
           0, 0, 0,
            1, 0, 0,
            0, -1, 0;
    expectedAfter <<
          0, 0, 0,
            0, 1, 0,
            1, 0, 0;
    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    float angle = M_PI / 2.;  // 90 degrees in radians
    Eigen::Vector3f axis(0, 0, 1);
    Eigen::AngleAxisf rotation(angle, axis);
    A_applied.rotate(rotation);

    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "Q before\n" << before.matrix() << std::endl;
    CoordsMatrixType<float> after = Kabsch_Umeyama<float>::applyTransform(A_applied, before);
    std::cout << "Q after\n" << after.matrix() << std::endl;

    auto A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(before, after, false, true, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_applied.matrix().inverse()\n" << A_applied.matrix().inverse() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;
    std::cout << "difference\n" << (A_predicted.matrix() - A_applied.matrix()) << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(expectedAfter, after);
    TestHelper<float>::EXPECT_MATRIX_NEAR(A_applied.matrix(), A_predicted.matrix());
}



TEST(KabschTestSuite, TestRotateCube) {
    CoordsMatrixType<float> after(8, 3), before(8, 3);
    before << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    before = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(before);
    after = before;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();

    float angle = M_PI / 2.;  // 90 degrees in radians
    Eigen::Vector3f axis(0, 0, 1);
    Eigen::AngleAxisf rotation(angle, axis);

    A_applied.rotate(rotation);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "Q before\n" << after.matrix() << std::endl;
    after = Kabsch_Umeyama<float>::applyTransform(A_applied, after);
    std::cout << "Q after\n" << after.matrix() << std::endl;

    auto A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(before, after, false, true, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
//    std::cout << "A_applied.matrix().inverse()\n" << A_applied.matrix().inverse() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;
    std::cout << "difference\n" << (A_predicted.matrix() - A_applied.matrix()) << std::endl;
    std::cout << "(A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() = "
              << (A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() << std::endl;
    std::cout << std::endl << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(A_predicted.matrix(), A_applied.matrix(), 1e-6);

    std::cout << "----------------- We did the +ve. Let's do the -ve. ------------------" << std::endl;

    after = before;
    A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    angle = -M_PI / 2.;  // 90 degrees in radians
//    Eigen::Vector3f axis(0, 0, 1);
    rotation = Eigen::AngleAxisf(angle, axis);

    A_applied.rotate(rotation);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "Q before\n" << after.matrix() << std::endl;
    after = Kabsch_Umeyama<float>::applyTransform(A_applied, after);
    std::cout << "Q after\n" << after.matrix() << std::endl;

    A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(before, after, false, true, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
//    std::cout << "A_applied.matrix().inverse()\n" << A_applied.matrix().inverse() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;
    std::cout << "difference\n" << (A_predicted.matrix() - A_applied.matrix()) << std::endl;
    std::cout << "(A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() = "
              << (A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() << std::endl;
    std::cout << std::endl << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(A_predicted.matrix(), A_applied.matrix(), 1e-6);
}

TEST(KabschTestSuite, TestScale) {
    CoordsMatrixType<float> after(8, 3), before(8, 3);
    before << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    before = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(before);
    after = before;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    A_applied.scale(55);
    std::cout << "Q before\n" << after.matrix() << std::endl;
    after = Kabsch_Umeyama<float>::applyTransform(A_applied, after);
    std::cout << "Q after\n" << after.matrix() << std::endl;

    auto A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(before, after, false, true, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;
    std::cout << "difference\n" << (A_predicted.matrix() - A_applied.matrix()) << std::endl;
    std::cout << "(A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() = "
              << (A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(A_predicted.matrix(), A_applied.matrix(), 1e-6);
}


TEST(KabschTestSuite, TestComplexTransform){
    CoordsMatrixType<float> after(8, 3), before(8, 3);
    before << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    before = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(before);
    after = before;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    float angle = M_PI / 2.;  // 90 degrees in radians
    Eigen::Vector3f axis(0, 0, 1);
    Eigen::AngleAxisf rotation(angle, axis);

    A_applied.rotate(rotation);
    A_applied.scale(55);
    A_applied.translate(Eigen::Vector3f {400, 200, 7000});

    std::cout << "Q before\n" << after.matrix() << std::endl;
    after = Kabsch_Umeyama<float>::applyTransform(A_applied, after);
    std::cout << "Q after\n" << after.matrix() << std::endl;

    //Find what got you from before to after
    auto A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(before, after, false, true, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;
    std::cout << "difference\n" << (A_predicted.matrix() - A_applied.matrix()) << std::endl;
    std::cout << "(A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() = "
              << (A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(A_predicted.matrix(), A_applied.matrix(), 5e-6);

    //Repeat after centering the reference
    after = before;
    after = Kabsch_Umeyama<float>::applyTransform(A_applied, after);
    before = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(before);
    A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(before, after, true, true, true);
    //TODO complete
}


TEST(KabschTestSuite, TestEigenAffineTransformRoundTrip) {
    CoordsMatrixType<float> after(8, 3), before(8, 3);
    before << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    before = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(before);
    after = before;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    float angle = M_PI / 2.;  // 90 degrees in radians
    Eigen::Vector3f axis(0, 0, 1);
    Eigen::AngleAxisf rotation(angle, axis);

    A_applied.rotate(rotation);
    A_applied.scale(55);
    A_applied.translate(Eigen::Vector3f {400, 200, 7000});

    std::cout << "Q before\n" << after.matrix() << std::endl;
    after = Kabsch_Umeyama<float>::applyTransform(A_applied, after);
    std::cout << "Q after going\n" << after.matrix() << std::endl;
    //TODO TestHelper<float>::EXPECT_MATRIX_NEAR(after, <some precalculated matrix>, 5e-6);

    after = Kabsch_Umeyama<float>::applyInverseOfTransform(A_applied, after);
    std::cout << "Q after round trip\n" << after.matrix() << std::endl;
    std::cout << std::endl;

    std::cout << "difference between after and before\n" << (after - before) << std::endl;
    std::cout << "(after - before).cwiseAbs().maxCoeff() = " << (after - before).cwiseAbs().maxCoeff() << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(after, before, 5e-6);
}

TEST(KabschTestSuite, TestKabschRoundTrip) {
    CoordsMatrixType<float> midWay(8, 3), after2(8, 3), before(8, 3);
    before << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    before = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(before);
    midWay = before;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    float angle = M_PI / 2.;  // 90 degrees in radians
    Eigen::Vector3f axis(0, 0, 1);
    Eigen::AngleAxisf rotation(angle, axis);

    A_applied.rotate(rotation);
    A_applied.scale(55);
    A_applied.translate(Eigen::Vector3f {400, 200, 7000});

    std::cout << "Q before\n" << midWay.matrix() << std::endl;
    midWay = Kabsch_Umeyama<float>::applyTransform(A_applied, midWay);
    std::cout << "Q after going\n" << midWay.matrix() << std::endl;

    auto A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(before, midWay, false, true, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;
    //TODO is next line correct? Any missing comparisons?
    TestHelper<float>::EXPECT_MATRIX_NEAR(A_applied.matrix(), A_predicted.matrix(), 5e-6);

    after2 = before;
    std::cout << "starting the round trip. after2 before transformation:\n" << after2 << std::endl;
    after2 = Kabsch_Umeyama<float>::applyTransform(A_predicted, after2);
    std::cout << "after2 after going\n" << after2.matrix() << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(after2.matrix(), midWay.matrix(), 5e-6);
    after2 = Kabsch_Umeyama<float>::applyInverseOfTransform(A_predicted, after2);
    std::cout << "after2 after round trip\n" << after2.matrix() << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(after2.matrix(), before.matrix(), 5e-6);
    std::cout << std::endl;

    std::cout << "difference between after2 and before\n" << (after2 - before) << std::endl;
    std::cout << "(after - before).cwiseAbs().maxCoeff() = " << (after2 - before).cwiseAbs().maxCoeff() << std::endl;
    TestHelper<float>::EXPECT_MATRIX_NEAR(after2, before, 5e-6);
}

TEST(KabschTestSuite, TestMaxError) {
    // Create datasets with known transform
    CoordsMatrixType<double> in(100, 3), out(100, 3);
    Eigen::Quaternion<double> Q(1, 3, 5, 2);
    Q.normalize();
    Eigen::Matrix3<double> R = Q.toRotationMatrix();
    double scale = 2.0;
    for (int col = 0; col < in.cols(); col++) {
        for (int row = 0; row < in.rows(); row++) {
            in(row, col) = log(2 * col + 10.0) / sqrt(1.0 * row + 4.0) + sqrt(row * 1.0) / (col + 1.0);
        }
    }
    Eigen::Vector3<double> S;
    S << -5, 6, -27;
    for (int row = 0; row < in.rows(); row++)
        out.row(row) = scale * R * in.row(row).transpose() + S;
    // See if we got the transform we expected
    std::cout << "S\n" << S << std::endl;
    std::cout << "scale\n" << scale << std::endl;
    std::cout << "R\n" << R << std::endl;
    Eigen::Transform<double, 3, Eigen::Affine> A_predicted_d = Kabsch_Umeyama<double>::find3DAffineTransform(in, out,
                                                                                                             false,
                                                                                                             true,
                                                                                                             true);
    std::cout << "A_predicted_d\n" << A_predicted_d.matrix() << std::endl;
    TestHelper<double>::EXPECT_MATRIX_NEAR(scale * R.cast<double>(), A_predicted_d.linear(), 1e-13);
    TestHelper<double>::EXPECT_MATRIX_NEAR(S.cast<double>(), A_predicted_d.translation(), 1e-13);

    Eigen::Transform<float, 3, Eigen::Affine> A_f = Kabsch_Umeyama<float>::find3DAffineTransform(
            in.cast<float>(), out.cast<float>(),
            false, true, true);
    TestHelper<float>::EXPECT_MATRIX_NEAR(static_cast<float>(scale) * R.cast<float>(), A_f.linear(), 1e-4);
    TestHelper<float>::EXPECT_MATRIX_NEAR(S.cast<float>(), A_f.translation(), 1e-4);
}