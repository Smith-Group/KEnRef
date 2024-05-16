#include <Eigen/Geometry>
#include <core/kabsch.h>
#include "gtest/gtest.h"
#include "testHelper.h"

// Just some dummy tests
TEST(KabschTestSuite, TestTranslate) {
    CoordsMatrixType<float> query(8, 3), ref(8, 3);
    ref << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    query = ref;

    query.rowwise() += Eigen::RowVector3<float>{0, 10000, 0};

    auto A = Kabsch_Umeyama<float>::find3DAffineTransform(ref, query, false, true);
    std::cout << "Q\n" << query.matrix() << std::endl;
    std::cout << "R\n" << ref.matrix() << std::endl;
    std::cout << "A\n" << A.matrix() << std::endl;
    std::cout << "A.translation().matrix()\n" << A.translation().matrix() << std::endl;
    std::cout << "A.rotation().matrix()\n" << A.rotation().matrix() << std::endl;

    /////////////////////////////////////
    std::cout << "========================\n";
    query = ref;
    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    A_applied.translate(Eigen::Vector3<float>{0, 10000, 0});
    std::cout << "Query before translation\n" << query.matrix() << std::endl;
    query = Kabsch_Umeyama<float>::applyTransform(A_applied, query);
    std::cout << "Query after translation\n" << query.matrix() << std::endl;

    auto A_calculated = Kabsch_Umeyama<float>::find3DAffineTransform(ref, query, false, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_calculated\n" << A_calculated.matrix() << std::endl;
    std::cout << "A_calculated.matrix().inverse()\n" << A_calculated.matrix().inverse() << std::endl;

    std::cout << "difference\n" << (A_calculated.matrix() - A_applied.matrix()) << std::endl;
    EXPECT_LE(((A_calculated.matrix() - A_applied.matrix())).cwiseAbs().maxCoeff(), 1e-12);
    EXPECT_LE(10, 11);
}

TEST(KabschTestSuite, TestRotate) {
    CoordsMatrixType<float> query(8, 3), ref(8, 3);
    ref << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    ref = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(ref);
    query = ref;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();

    float angle = M_PI / 2.;  // 180 degrees in radians
    Eigen::Vector3f axis(0, 0, 1);
    Eigen::AngleAxisf rotation(angle, axis);

    A_applied.rotate(rotation);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "Q before\n" << query.matrix() << std::endl;
    query = Kabsch_Umeyama<float>::applyTransform(A_applied, query);
    std::cout << "Q after\n" << query.matrix() << std::endl;

    auto A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(ref, query, false, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
//    std::cout << "A_applied.matrix().inverse()\n" << A_applied.matrix().inverse() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;
    std::cout << "difference\n" << (A_predicted.matrix() - A_applied.matrix()) << std::endl;
    std::cout << "(A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() = "
              << (A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() << std::endl;
    std::cout << std::endl << std::endl << std::endl;
    EXPECT_LE((A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff(), 1e-6);
}

TEST(KabschTestSuite, TestScale) {
    CoordsMatrixType<float> query(8, 3), ref(8, 3);
    ref << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    ref = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(ref);
    query = ref;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    A_applied.scale(55);
    std::cout << "Q before\n" << query.matrix() << std::endl;
    query = Kabsch_Umeyama<float>::applyTransform(A_applied, query);
    std::cout << "Q after\n" << query.matrix() << std::endl;

    auto A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(ref, query, false, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;
    std::cout << "difference\n" << (A_predicted.matrix() - A_applied.matrix()) << std::endl;
    std::cout << "(A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() = "
              << (A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() << std::endl;
    std::cout << std::endl;
//    GTEST_NONFATAL_FAILURE_("values mismatch");
    EXPECT_LE((A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff(), 1e-6);
}


TEST(KabschTestSuite, TestComplexTransform){
    CoordsMatrixType<float> query(8, 3), ref(8, 3);
    ref << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    ref = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(ref);
    query = ref;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    float angle = M_PI / 2.;  // 180 degrees in radians
    Eigen::Vector3f axis(0, 0, 1);
    Eigen::AngleAxisf rotation(angle, axis);

    A_applied.rotate(rotation);
    A_applied.scale(55);
    A_applied.translate(Eigen::Vector3f {400, 200, 7000});

    std::cout << "Q before\n" << query.matrix() << std::endl;
    query = Kabsch_Umeyama<float>::applyTransform(A_applied, query);
    std::cout << "Q after\n" << query.matrix() << std::endl;

    auto A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(ref, query, false, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;
    std::cout << "difference\n" << (A_predicted.matrix() - A_applied.matrix()) << std::endl;
    std::cout << "(A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() = "
              << (A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff() << std::endl;
    std::cout << std::endl;
    EXPECT_LE((A_predicted.matrix() - A_applied.matrix()).cwiseAbs().maxCoeff(), 5e-6);
}


TEST(KabschTestSuite, TestEigenAffineTransformRoundTrip) {
    CoordsMatrixType<float> query(8, 3), ref(8, 3);
    ref << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    ref = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(ref);
    query = ref;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    float angle = M_PI / 2.;  // 180 degrees in radians
    Eigen::Vector3f axis(0, 0, 1);
    Eigen::AngleAxisf rotation(angle, axis);

    A_applied.rotate(rotation);
    A_applied.scale(55);
    A_applied.translate(Eigen::Vector3f {400, 200, 7000});

    std::cout << "Q before\n" << query.matrix() << std::endl;
    query = Kabsch_Umeyama<float>::applyTransform(A_applied, query);
    std::cout << "Q after going\n" << query.matrix() << std::endl;

    query = Kabsch_Umeyama<float>::applyInverseOfTransform(A_applied, query);
    std::cout << "Q after round trip\n" << query.matrix() << std::endl;
    std::cout << std::endl;

    std::cout << "difference between query and ref\n" << (query - ref) << std::endl;
    std::cout << "(query - ref).cwiseAbs().maxCoeff() = " << (query - ref).cwiseAbs().maxCoeff() << std::endl;
    EXPECT_LE((query - ref).cwiseAbs().maxCoeff(), 5e-6);
}

TEST(KabschTestSuite, TestKabschRoundTrip) {
    CoordsMatrixType<float> query1(8, 3), query2(8, 3), ref(8, 3);
    ref << 0, 0, 0,
            10, 0, 0,
            10, 10, 0,
            0, 10, 0,
            0, 0, 10,
            10, 0, 10,
            10, 10, 10,
            0, 10, 10;
    ref = Kabsch_Umeyama<float>::translateCenterOfMassToOrigin(ref);
    query1 = ref;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    float angle = M_PI / 2.;  // 180 degrees in radians
    Eigen::Vector3f axis(0, 0, 1);
    Eigen::AngleAxisf rotation(angle, axis);

    A_applied.rotate(rotation);
    A_applied.scale(55);
    A_applied.translate(Eigen::Vector3f {400, 200, 7000});

    std::cout << "Q before\n" << query1.matrix() << std::endl;
    query1 = Kabsch_Umeyama<float>::applyTransform(A_applied, query1);
    std::cout << "Q after going\n" << query1.matrix() << std::endl;

    auto A_predicted = Kabsch_Umeyama<float>::find3DAffineTransform(ref, query1, false, true);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;

    query2 = ref;
    std::cout << "starting the round trip. query2 before transformation:\n" << query2 << std::endl;
    query2 = Kabsch_Umeyama<float>::applyTransform(A_predicted, query2);
    std::cout << "Query2 after going\n" << query2.matrix() << std::endl;
    query2 = Kabsch_Umeyama<float>::applyInverseOfTransform(A_predicted, query2);
    std::cout << "Query2 after round trip\n" << query2.matrix() << std::endl;
    std::cout << std::endl;

    std::cout << "difference between query2 and ref\n" << (query2 - ref) << std::endl;
    std::cout << "(query - ref).cwiseAbs().maxCoeff() = " << (query2 - ref).cwiseAbs().maxCoeff() << std::endl;
    EXPECT_LE((query2 - ref).cwiseAbs().maxCoeff(), 5e-6);
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
    Eigen::Transform<double, 3, Eigen::Affine> A_predicted_d = Kabsch_Umeyama<double>::find3DAffineTransform(in, out, false, true);
    std::cout << "A_predicted_d\n" << A_predicted_d.matrix() << std::endl;
    EXPECT_LE((scale * R.cast<double>() - A_predicted_d.linear()).cwiseAbs().maxCoeff(), 1e-13);
    EXPECT_LE((S.cast<double>() - A_predicted_d.translation()).cwiseAbs().maxCoeff(), 1e-13);

    Eigen::Transform<float, 3, Eigen::Affine> A_f = Kabsch_Umeyama<float>::find3DAffineTransform(in.cast<float>(),
                                                                                                 out.cast<float>(),
                                                                                                 false, true);
    EXPECT_LE((static_cast<float>(scale) * R.cast<float>() - A_f.linear()).cwiseAbs().maxCoeff(), 1e-4);
    EXPECT_LE((S.cast<float>() - A_f.translation()).cwiseAbs().maxCoeff(), 1e-4);
}