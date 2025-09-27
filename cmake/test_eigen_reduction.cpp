#include <iostream>
#include <Eigen/Dense>
#include <omp.h>

int main() {
    const Eigen::Matrix3d a = Eigen::Matrix3d::Random();
    Eigen::Matrix3d result = Eigen::Matrix3d::Zero();

#pragma omp parallel for reduction(+:result)
    for (int i = 0; i < 4; ++i) {
        result += a;
    }

    const Eigen::Matrix3d expected = 4 * a;
    std::cout << "expected\n" << expected << "\nresult\n" << result << std::endl;
    // std::cout << "difference = " << (result - expected).norm() << std::endl;
    return (result - expected).norm() < 1e-12 ? 0 : 1;
}