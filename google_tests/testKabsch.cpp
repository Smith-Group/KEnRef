#include <Eigen/Geometry>
#include <core/kabsch.h>
#include "gtest/gtest.h"

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

    auto A = Kabsch<float>::Find3DAffineTransform(query, ref);
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
    query = Kabsch<float>::applyTransform(A_applied, query);
    std::cout << "Query after translation\n" << query.matrix() << std::endl;

    auto A_calculated = Kabsch<float>::Find3DAffineTransform(query, ref);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_calculated\n" << A_calculated.matrix() << std::endl;
    std::cout << "A_calculated.matrix().inverse()\n" << A_calculated.matrix().inverse() << std::endl;

    std::cout << "difference\n" << (A_calculated.matrix() - A_applied.matrix().inverse()) << std::endl;
    EXPECT_LE(((A_calculated.matrix() - A_applied.matrix().inverse())).cwiseAbs().maxCoeff(), 1e-12);
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
    query = ref;

    auto A_applied = Eigen::Transform<float, 3, Eigen::Affine>::Identity();
    std::cout << "A_applied in the making\n" << A_applied.matrix() << std::endl;

    float angle = M_PI / 2.;  // 180 degrees in radians
    Eigen::Vector3f axis(0, 0, 1);
    Eigen::AngleAxisf rotation(angle, axis);

    A_applied.rotate(rotation);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "Q before\n" << query.matrix() << std::endl;
    query = Kabsch<float>::applyTransform(A_applied, query);
    std::cout << "Q after\n" << query.matrix() << std::endl;

    auto A_predicted = Kabsch<float>::Find3DAffineTransform(query, ref);
    std::cout << "A_applied\n" << A_applied.matrix() << std::endl;
    std::cout << "A_applied.matrix().inverse()\n" << A_applied.matrix().inverse() << std::endl;
    std::cout << "A_predicted\n" << A_predicted.matrix() << std::endl;
    std::cout << "difference\n" << (A_predicted.matrix() - A_applied.matrix().inverse()) << std::endl;
    std::cout << "(A_predicted.matrix() - A_applied.matrix().inverse()).cwiseAbs().maxCoeff() = "
              << (A_predicted.matrix() - A_applied.matrix().inverse()).cwiseAbs().maxCoeff() << std::endl;
    EXPECT_LE((A_predicted.matrix() - A_applied.matrix().inverse()).cwiseAbs().maxCoeff(), 1e-6);
}

TEST(KabschTestSuite, TestScale) {
    GTEST_SKIP_("NOT YET IMPLEMENTED");
    //TODO complete
}

TEST(KabschTestSuite, TestRoundTrip) {
    //TODO complete
}

TEST(KabschTestSuite, TestMaxError) {
    GTEST_SKIP_("Skipping single test until we fix the scale issue");
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
    Eigen::Transform<double, 3, Eigen::Affine> A_d = Kabsch<double>::Find3DAffineTransform(in, out);
    std::cout << "a_d\n" << A_d.matrix() << std::endl;
    EXPECT_LE((scale * R.cast<double>() - A_d.linear()).cwiseAbs().maxCoeff(), 1e-13);
    EXPECT_LE((S.cast<double>() - A_d.translation()).cwiseAbs().maxCoeff(), 1e-13);

    Eigen::Transform<float, 3, Eigen::Affine> A_f = Kabsch<float>::Find3DAffineTransform(in.cast<float>(),
                                                                                         out.cast<float>());
    EXPECT_LE((static_cast<float>(scale) * R.cast<float>() - A_f.linear()).cwiseAbs().maxCoeff(), 1e-4);
    EXPECT_LE((S.cast<float>() - A_f.translation()).cwiseAbs().maxCoeff(), 1e-4);
}