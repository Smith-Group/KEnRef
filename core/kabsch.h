// This code is released in public domain

#ifndef KABSCH_H
#define KABSCH_H

#include <Eigen/Geometry>
#include <iostream> //FIXME remove this import line later

// Given two sets of 3D points, find the rotation + translation + scale
// which best maps the first set to the second.
// Source: http://en.wikipedia.org/wiki/Kabsch_algorithm

// The input 3D points are stored as rows (every point in a row).
Eigen::Affine3f Find3DAffineTransform(const Eigen::MatrixX3f p, const Eigen::MatrixX3f q) {

  // Default output
  Eigen::Affine3f A;
  A.linear() = Eigen::Matrix3f::Identity(3, 3);
  A.translation() = Eigen::Vector3f::Zero();

  if (p.rows() != q.rows())
    throw "Find3DAffineTransform(): input data mis-match";

  // First find the scale, by finding the ratio of sums of some distances,
  // then bring the datasets to the same scale.
  double dist_p = 0, dist_q = 0;
  for (int row = 0; row < p.rows()-1; row++) {
    dist_p  += (p.row(row+1) - p.row(row)).norm();
    dist_q += (q.row(row+1) - q.row(row)).norm();
  }
  if (dist_p <= 0 || dist_q <= 0)
    return A;

  Eigen::MatrixX3f p_temp = p;
  Eigen::MatrixX3f q_temp = q;

  double scale = dist_q/dist_p;
  q_temp = q / scale;

  // Find the centroids then shift to the origin
  Eigen::Vector3f p_ctr = Eigen::Vector3f::Zero();
  Eigen::Vector3f q_ctr = Eigen::Vector3f::Zero();
  for (int row = 0; row < p_temp.rows(); row++) {
    p_ctr += p_temp.row(row);
    q_ctr += q_temp.row(row);
  }
  p_ctr /= p_temp.rows();
  q_ctr /= q_temp.rows();
  for (int row = 0; row < p_temp.rows(); row++) {
	  p_temp.row(row) -= p_ctr;
	  q_temp.row(row) -= q_ctr;
  }

  // SVD
//  Eigen::MatrixXf Cov = p_temp * q_temp.transpose();
  Eigen::MatrixXf Cov = p_temp.transpose() * q_temp; //FIXME I am not sure about this line
  Eigen::JacobiSVD<Eigen::MatrixXf> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

  std::cout << "Matrix U " << std::endl <<svd.matrixU() << std::endl << std::endl;
  std::cout << "Matrix V " << std::endl << svd.matrixV() << std::endl << std::endl;

  // Find the rotation
  float d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
  if (d > 0)
    d = 1.0;
  else
    d = -1.0;
  Eigen::Matrix3f I = Eigen::Matrix3f::Identity(3, 3);
  I(2, 2) = d;
  Eigen::Matrix3f R = svd.matrixV() * I * svd.matrixU().transpose();

  // The final transform
  A.linear() = scale * R;
  A.translation() = scale*(q_ctr - R*p_ctr);

  return A;
}

// A function to test Find3DAffineTransform()
void TestFind3DAffineTransform(){

  // Create datasets with known transform
  Eigen::MatrixX3f in(100, 3), out(100, 3);
  Eigen::Quaternion<float> Q(1, 3, 5, 2);
  Q.normalize();
  Eigen::Matrix3f R = Q.toRotationMatrix();
  double scale = 2.0;
  for (int col = 0; col < in.cols(); col++) {
	  for (int row = 0; row < in.rows(); row++) {
		  in(row, col) = log(2*col + 10.0)/sqrt(1.0*row+ 4.0) + sqrt(row*1.0)/(col + 1.0);
	  }
  }
  Eigen::Vector3f S;
  S << -5, 6, -27;
  for (int row = 0; row < in.rows(); row++)
    out.row(row) = scale*R*in.row(row).transpose() + S;

  Eigen::Affine3f A = Find3DAffineTransform(in, out);

  // See if we got the transform we expected
  if ( (scale * R - A.linear()).cwiseAbs().maxCoeff() > 1e-5 ||	//1e-13)  // (for double)
       (S - A.translation()).cwiseAbs().maxCoeff() > 1e-5)	//1e-13)  // (for double)
    throw "Could not determine the affine transform accurately enough";
}

#endif
