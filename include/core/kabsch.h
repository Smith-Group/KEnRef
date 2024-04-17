// This code is released in public domain

#ifndef KABSCH_H
#define KABSCH_H

#include <Eigen/Geometry>
#include <iostream> //FIXME remove this import line later

#define VERBOSE false

// Given two sets of 3D points, find the rotation + translation + scale
// which best maps the first set to the second.
// Source: http://en.wikipedia.org/wiki/Kabsch_algorithm
// original code obtained from
// https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp

// The input 3D points are stored as rows (every point in a row).
template <typename precision>
class Kabsch{
    //TODO check whether there is an available implementation in Eigen.
    //TODO you can check here https://eigen.tuxfamily.org/dox/group__TutorialLinearAlgebra.html
    //TODO and https://eigen.tuxfamily.org/dox/group__TutorialGeometry.html
public:
    static Eigen::Transform<precision,3,Eigen::Affine> Find3DAffineTransform(const CoordsMatrixType<precision>& p, const CoordsMatrixType<precision>& q) {

        // Default output
		Eigen::Transform<precision,3,Eigen::Affine> A;
		A.linear() = Eigen::Matrix3<precision>::Identity(3, 3);
		A.translation() = Eigen::Vector3<precision>::Zero();

		if (p.rows() != q.rows())
            throw std::runtime_error("Find3DAffineTransform(): input data mis-match");

		// First find the scale, by finding the ratio of sums of some distances,
		// then bring the datasets to the same scale.
		precision dist_p = 0, dist_q = 0;
		for (int row = 0; row < p.rows()-1; row++) {
			dist_p  += (p.row(row+1) - p.row(row)).norm();
			dist_q += (q.row(row+1) - q.row(row)).norm();
		}
		if (dist_p <= 0 || dist_q <= 0)
			return A;

		Eigen::MatrixX3<precision> p_temp = p.template cast <precision> ();
		Eigen::MatrixX3<precision> q_temp = q.template cast <precision> ();

		precision scale = dist_q/dist_p;
		q_temp /= scale;

		// Find the centroids then shift to the origin
		Eigen::Vector3<precision> p_ctr = Eigen::Vector3<precision>::Zero();
		Eigen::Vector3<precision> q_ctr = Eigen::Vector3<precision>::Zero();
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
		Eigen::MatrixX<precision> Cov = p_temp.transpose() * q_temp; //FIXME I am not sure about this line
		Eigen::JacobiSVD<Eigen::MatrixX<precision>> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

#if VERBOSE
			std::cout << "Matrix U " << std::endl <<svd.matrixU() << std::endl << std::endl;
			std::cout << "Matrix V " << std::endl << svd.matrixV() << std::endl << std::endl;
#endif

		// Find the rotation
		float d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
		if (d > 0)
			d = 1.0;
		else
			d = -1.0;
		Eigen::Matrix3<precision> I = Eigen::Matrix3<precision>::Identity(3, 3);
		I(2, 2) = d;
		Eigen::Matrix3<precision> R = svd.matrixV() * I * svd.matrixU().transpose();

		// The final transform
		A.linear() = scale * R;
		A.translation() = scale*(q_ctr - R*p_ctr);
		return A;
	}
};

// A function to test Find3DAffineTransform()
void TestFind3DAffineTransform(){

  // Create datasets with known transform
  Eigen::MatrixX3<double> in(100, 3), out(100, 3);
  Eigen::Quaternion<double> Q(1, 3, 5, 2);
  Q.normalize();
  Eigen::Matrix3<double> R = Q.toRotationMatrix();
  double scale = 2.0;
  for (int col = 0; col < in.cols(); col++) {
	  for (int row = 0; row < in.rows(); row++) {
		  in(row, col) = log(2*col + 10.0)/sqrt(1.0*row+ 4.0) + sqrt(row*1.0)/(col + 1.0);
	  }
  }
  Eigen::Vector3<double> S;
  S << -5, 6, -27;
  for (int row = 0; row < in.rows(); row++)
    out.row(row) = scale*R*in.row(row).transpose() + S;

  Eigen::Transform<double,3,Eigen::Affine> A = Kabsch<double>::Find3DAffineTransform(in, out);

  // See if we got the transform we expected
    if ((scale * R.cast<double>() - A.linear()).cwiseAbs().maxCoeff() > 1e-5 ||    //1e-13)  // (for double)
        (S.cast<double>() - A.translation()).cwiseAbs().maxCoeff() > 1e-5) {    //1e-13)  // (for double)
        throw std::runtime_error("Could not determine the affine transform accurately enough");
    }
}

#undef VERBOSE


template class Kabsch<float>;
template class Kabsch<double>;

#endif
