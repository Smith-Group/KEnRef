// This code is released in public domain

#ifndef KABSCH_H
#define KABSCH_H

#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <iostream>
#include "core/KEnRef.h"

#define VERBOSE false

// Given two sets of 3D points, find the rotation + translation + scale
// which best transforms the second set to the first (i.e. if applied to the second will result in the first).
//The first and second are also called out and in or after and before.
// Source: http://en.wikipedia.org/wiki/Kabsch_algorithm
// original code obtained from
// https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
// The input 3D points are stored as rows (every point in a row).
//This is the latest modification https://zpl.fi/aligning-point-patterns-with-kabsch-umeyama-algorithm/
template<typename precision>
class Kabsch_Umeyama {
    //TODO check whether there is an available implementation in Eigen.
    //TODO you can check here https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
    //TODO and https://eigen.tuxfamily.org/dox/group__TutorialGeometry.html
public:
    static Eigen::Transform<precision, 3, Eigen::Affine>
    find3DAffineTransform(const CoordsMatrixType<precision> &after, const CoordsMatrixType<precision> &before,
                          bool afterAlreadyCentered = false, bool doTranslate=true, bool doScale = false) {

        // Default output
        Eigen::Transform<precision, 3, Eigen::Affine> A;
        A.linear() = Eigen::Matrix3<precision>::Identity(3, 3);
        A.translation() = Eigen::Vector3<precision>::Zero();

        if (after.rows() != before.rows())
            throw std::runtime_error("find3DAffineTransform(): input data mis-match");

        CoordsMatrixType<precision> after_temp = after.template cast<precision>();
        CoordsMatrixType<precision> before_temp = before.template cast<precision>();

        // Find the centroids then shift to the origin
        Eigen::Vector3<precision> after_centroid = Eigen::Vector3<precision>::Zero();
        Eigen::Vector3<precision> before_centroid = Eigen::Vector3<precision>::Zero();
        if (! afterAlreadyCentered){
            //p_temp = translateCenterOfMassToOrigin(p_temp); //can't use it yet, coz it does not return the COM (yet)
            after_centroid = after_temp.colwise().mean();
            after_temp = (after_temp.rowwise() - after_centroid.transpose()).eval();
        }
        before_centroid = before_temp.colwise().mean();
        before_temp = (before_temp.rowwise() - before_centroid.transpose()).eval();

        // SVD
        Eigen::MatrixX<precision> Cov = after_temp.transpose() * before_temp;
        Eigen::JacobiSVD<Eigen::MatrixX<precision>> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

#if VERBOSE
        std::cout << "Matrix U " << std::endl <<svd.matrixU() << std::endl << std::endl;
        std::cout << "Matrix V " << std::endl << svd.matrixV() << std::endl << std::endl;
#endif
        // Find the rotation

        // this implementaion should be also fine
//        precision d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
//        d = (d > 0) ? 1.0 : -1.0; // improper implementation of sign()

        //This is the Wikipedia implememtation
        precision d = svd.matrixV().determinant() * svd.matrixU().determinant();
        //TODO assert that d == 1 or -1

        Eigen::Matrix3<precision> I = Eigen::Matrix3<precision>::Identity(3, 3);
        I(2, 2) = d;
        Eigen::Matrix3<precision> R = svd.matrixU() * I * svd.matrixV().transpose();

        // Find the scale
#if VERBOSE
        std::cout << "(svd.singularValues().asDiagonal() * I).trace()" << std::endl << (svd.singularValues().asDiagonal() * I).trace() << std::endl;
        std::cout << "(after_temp.transpose() * after_temp).trace()" << std::endl << (after_temp.transpose() * after_temp).trace() << std::endl;
        std::cout << "calculateVariance(after_temp)" << std::endl << calculateVariance(after_temp) << std::endl;
        std::cout << "(svd.singularValues().asDiagonal() * I).trace()" << std::endl << (svd.singularValues().asDiagonal() * I).trace() << std::endl;
#endif
//        precision scale = (svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixU().transpose()).trace() / (p_temp * p_temp.transpose()).trace();
//        precision scale = calculateVariance(p_temp) / (svd.singularValues().asDiagonal() * I).trace();
//        if(scale == 0) scale = 1;
        precision scale = doScale ? (after_temp * after_temp.transpose()).trace() / (svd.singularValues().asDiagonal() * I).trace() : 1.;

        // The final transform
        A.linear() = scale *  R;
        if (doTranslate)
            A.translation() = after_centroid - scale * R * before_centroid;
        return A;
    }

    // Translate the center of mass of a 3D point cloud to the origin. This could save some processing time if you
    // consistently fit some points (atoms) to the same reference structure everytime.
    static inline CoordsMatrixType<precision> translateCenterOfMassToOrigin(const CoordsMatrixType<precision>& points) {
        // Compute the center of mass
        Eigen::RowVector3<precision> center_of_mass = points.colwise().mean();
        // Translate the points so that the center of mass is at the origin
        CoordsMatrixType<precision> translated_points = points.rowwise() - center_of_mass;
        return translated_points;
    }

    static inline precision calculateVariance(const CoordsMatrixType<precision>& A) {
        // Compute the centroid
        Eigen::RowVector3<precision> centroid = A.colwise().mean();
        // Compute the distances from the centroid
        Eigen::VectorX<precision> distances = (A.rowwise() - centroid).rowwise().norm();
        // Compute the mean of the distances
        precision mean_distance = distances.mean();
        // Compute the variance of the distances
        precision variance = (distances.array() - mean_distance).square().mean();
        return variance;
    }

    static CoordsMatrixType<precision>
    applyTransform(Eigen::Transform<precision, 3, Eigen::Affine> affine, CoordsMatrixType<precision> coords) {
        return (coords.rowwise().homogeneous() * affine.matrix().transpose()).leftCols(3).eval();
    }

    static CoordsMatrixType<precision>
    applyInverseOfTransform(Eigen::Transform<precision, 3, Eigen::Affine> affine, CoordsMatrixType<precision> coords) {
        return (coords.rowwise().homogeneous() * affine.inverse().matrix().transpose()).leftCols(3).eval();
    }


};


#undef VERBOSE


template
class Kabsch_Umeyama<float>;

template
class Kabsch_Umeyama<double>;

#endif
