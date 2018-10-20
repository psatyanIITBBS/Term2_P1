#include "kalman_filter.h"
#include <iostream>
#include "Eigen/Dense"
#include <cmath>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
	cout << " Present x in KF = " << x_ << endl;
	x_ = F_ * x_;
	cout << " Predicted x in KF = " << x_ << endl;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
	cout << " P in KF = " << P_ << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	cout << " Measured data z in Update = " << z << endl;
	VectorXd y_ = z - H_ * x_;
	cout << " Error in predicted y_ in KF = " << y_ << endl;
	MatrixXd Ht_ = H_.transpose();
	cout << " Ht_ in KF = " << Ht_ << endl;
	MatrixXd S_ = H_ * P_ * Ht_ + R_;
	MatrixXd Si_ = S_.inverse();
	MatrixXd K_ =  P_ * Ht_ * Si_;
	cout << " K_ in KF = " << K_ << endl;

	MatrixXd I = MatrixXd::Identity(4, 4);
	cout << " I in KF = " << I << endl;
	//new state
	x_ = x_ + (K_ * y_);
	cout << " Corrected state x in Update = " << x_ << endl;
	P_ = (I - K_ * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
	cout << " z in EKF = " << z << endl;
	double a = sqrt(x_[0]*x_[0] + x_[1]*x_[1]);
	cout << " a in EKF = " << a << endl;
	double b = atan2(x_[1],x_[0]);
	b = z[1]- b;
	while ( b > M_PI || b < -M_PI ) {
	    if ( b > M_PI ) {

	      b -= M_PI;
	    } else {
	      b += M_PI;
	    }
	  }

	cout << " b in EKF = " << b << endl;
	double c = (x_[0]*x_[2] + x_[1]*x_[3])/a;
	cout << " c in EKF = " << c << endl;

	VectorXd y_ = z - H_ * x_;
	y_[0] = z[0]- a;
	y_[1] =  b;
	y_[2] = z[2]- c;
	cout << " y_ in EKF = " << y_ << endl;
	MatrixXd Ht_ = H_.transpose();
	//cout << " Ht_ in EKF = " << Ht_ << endl;
	MatrixXd S_ = H_ * P_ * Ht_ + R_;
	MatrixXd Si_ = S_.inverse();
	MatrixXd K_ =  P_ * Ht_ * Si_;
	//cout << " K_ in EKF = " << K_ << endl;

	MatrixXd I = MatrixXd::Identity(4, 4);
	//cout << " I in EKF = " << I << endl;
	//new state
	x_ = x_ + (K_ * y_);
	P_ = (I - K_ * H_) * P_;
}
