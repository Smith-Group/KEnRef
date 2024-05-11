// This code is released in public domain

#ifndef KABSCH_H
#define KABSCH_H

#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <iostream>
#include "core/KEnRef.h"

#define VERBOSE false

// Given two sets of 3D points, find the rotation + translation + scale
// which best transforms the first set to the second (i.e. if applied to the first will result in the second).
//The first and second are also called in and out or before and after.
// Source: http://en.wikipedia.org/wiki/Kabsch_algorithm
// original code obtained from
// https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
// The input 3D points are stored as rows (every point in a row).
template<typename precision>
class Kabsch {
    //TODO check whether there is an available implementation in Eigen.
    //TODO you can check here https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
    //TODO and https://eigen.tuxfamily.org/dox/group__TutorialGeometry.html
public:
    static Eigen::Transform<precision, 3, Eigen::Affine>
    Find3DAffineTransform(const CoordsMatrixType<precision> &p, const CoordsMatrixType<precision> &q,
                          bool secondAlreadyCentered = false) {

        // Default output
		Eigen::Transform<precision,3,Eigen::Affine> A;
		A.linear() = Eigen::Matrix3<precision>::Identity(3, 3);
		A.translation() = Eigen::Vector3<precision>::Zero();

		if (p.rows() != q.rows())
            throw std::runtime_error("Find3DAffineTransform(): input data mis-match");

		Eigen::MatrixX3<precision> p_temp = p.template cast <precision> ();
		Eigen::MatrixX3<precision> q_temp = q.template cast <precision> ();

        // Find the centroids then shift to the origin
        Eigen::Vector3<precision> p_ctr = Eigen::Vector3<precision>::Zero();
        Eigen::Vector3<precision> q_ctr = Eigen::Vector3<precision>::Zero();
        p_ctr = p_temp.colwise().mean();
        p_temp = (p_temp.rowwise() - p_ctr.transpose()).eval();
        if (! secondAlreadyCentered){
            //q_temp = translateCenterOfMassToOrigin(q_temp); //can't use it, coz it does not return the COM
            q_ctr = q_temp.colwise().mean();
            q_temp = (q_temp.rowwise() - q_ctr.transpose()).eval();
        }

		// SVD
		//  Eigen::MatrixXf Cov = p_temp * q_temp.transpose();
		Eigen::MatrixX<precision> Cov = p_temp.transpose() * q_temp; //FIXME I am not sure about this line
		Eigen::JacobiSVD<Eigen::MatrixX<precision>> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

#if VERBOSE
			std::cout << "Matrix U " << std::endl <<svd.matrixU() << std::endl << std::endl;
			std::cout << "Matrix V " << std::endl << svd.matrixV() << std::endl << std::endl;
#endif

        // Find the rotation
        precision d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
        if (d > 0)
            d = 1.0;
        else
            d = -1.0;
        Eigen::Matrix3<precision> I = Eigen::Matrix3<precision>::Identity(3, 3);
        I(2, 2) = d;
        Eigen::Matrix3<precision> R = svd.matrixV() * I * svd.matrixU().transpose();

        // The final transform
        A.linear() = R;
        A.translation() = q_ctr - R * p_ctr;
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
class Kabsch<float>;

template
class Kabsch<double>;

#endif
