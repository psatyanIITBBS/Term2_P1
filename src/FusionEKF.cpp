#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
		  	  0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
		  	  0, 0.0009, 0,
			  0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  Q_ = MatrixXd(4,4);
  F_ = MatrixXd(4,4);
  P_ = MatrixXd(4,4);

  Q_ << 1,0,1,0,
		0,1,0,1,
		1,0,1,0,
		0,1,0,1;

  F_ << 1,0,1,0,
		0,1,0,1,
		0,0,1,0,
		0,0,0,1;

  P_ << 1,0,0,0,
		0,1,0,0,
		0,0,10,0,
		0,0,0,10;

  H_laser_ << 1,0,0,0,
		  	  0,1,0,0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  double noise_ax = 9.0;
  double noise_ay = 9.0;
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF initializing: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    	cout<<"First data - RADAR data - ";
    	float rho_ = measurement_pack.raw_measurements_[0];
    	float phi_ = measurement_pack.raw_measurements_[1];
    	float rhoDot_ = measurement_pack.raw_measurements_[2];

    	ekf_.x_[0] = rho_*cos(phi_);
    	if (ekf_.x_[0]<0.0001){
    		ekf_.x_[0] = 0.0001;
    	}
    	ekf_.x_[1] = rho_*sin(phi_);
    	if (ekf_.x_[1]<0.0001){
    	    		ekf_.x_[1] = 0.0001;
    	    	}
    	ekf_.x_[2] = rhoDot_*cos(phi_);
    	if (ekf_.x_[2]<0.00001){
    	    		ekf_.x_[2] = 0.00001;
    	    	}
    	ekf_.x_[3] = rhoDot_*sin(phi_);
    	if (ekf_.x_[3]<0.00001){
    	    		ekf_.x_[3] = 0.00001;
    	    	}

    	cout << "Initial State vector: "<<ekf_.x_ << endl;
    	ekf_.F_ = F_;
    	ekf_.Q_ = Q_;
    	Hj_ = tools.CalculateJacobian(ekf_.x_);
    	ekf_.H_ = Hj_;
    	ekf_.R_ = R_radar_;
    	ekf_.P_ = P_;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
    	cout<<"First data - LASER data - ";
    	ekf_.x_[0] = measurement_pack.raw_measurements_[0];
    	ekf_.x_[1] = measurement_pack.raw_measurements_[1];

    	//cout << "Initial State vector: "<<ekf_.x_ << endl;
    	ekf_.F_ = F_;
    	ekf_.Q_ = Q_;
    	ekf_.H_ = H_laser_;
    	ekf_.R_ = R_laser_;
    	ekf_.P_ = P_;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt;
  dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  double dtsq = dt*dt;
  double dtcube = dtsq*dt/2.0;
  double dtfour = dtsq*dtsq/4.0;

  cout << "dt = " << dt << endl;
  F_ << 1,0,dt,0,
		0,1,0,dt,
		0,0,1,0,
		0,0,0,1;
  Q_ << dtfour*noise_ax,0,dtcube*noise_ax,0,
		0,dtfour*noise_ay,0,dtcube*noise_ay,
		dtcube*noise_ax,0,dtsq*noise_ax,0,
		0,dtcube*noise_ay,0,dtsq*noise_ay;

  ekf_.F_ = F_;
  ekf_.Q_ = Q_;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
	//return;
    // Radar updates
	z = VectorXd(3);
	z[0] = measurement_pack.raw_measurements_[0];
	z[1] = measurement_pack.raw_measurements_[1];
	z[2] = measurement_pack.raw_measurements_[2];

  	Hj_ = tools.CalculateJacobian(ekf_.x_);
  	ekf_.H_ = Hj_;
  	ekf_.R_ = R_radar_;
	ekf_.UpdateEKF(z);
  }

  else
  {
    // Laser updates
	z = VectorXd(2);
	z[0] = measurement_pack.raw_measurements_[0];
	z[1] = measurement_pack.raw_measurements_[1];
  	ekf_.H_ = H_laser_;
  	ekf_.R_ = R_laser_;
	ekf_.Update(z);
  }

  previous_timestamp_ = measurement_pack.timestamp_;
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
