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
 
  float noise_ax = 9.0;
  float noise_ay = 9.0;

  Tools tools;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


	if (measurement_pack.sensor_type_ == measurement_pack.LASER) {
		cout << "Laser " << endl;
	} else {
		cout << "Radar " << endl;
	}

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
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		float rho = measurement_pack.raw_measurements_(0);
		float phi = measurement_pack.raw_measurements_(1);
		float rho_dot = measurement_pack.raw_measurements_[2];

		float x = rho * cos(phi);
		float y = rho * sin(phi);
		float vx = rho_dot * cos(phi);
		float vy = rho_dot * sin(phi);

		ekf_.x_ << x, y, vx, vy;

	}
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0, 0.0;

    }

	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
			   0, 1, 0, 0,
		       0, 0, 1000, 0,
		       0, 0, 0, 1000;

	H_laser_ << 1, 0, 0, 0,
				0, 1, 0, 0;


    // done initializing, no need to predict or update
	previous_timestamp_ = measurement_pack.timestamp_;

    is_initialized_ = true;

	cout << "Init x_ = " << ekf_.x_ << endl;
	cout << "Init P_ = " << ekf_.P_ << endl;

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
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  cout << "dt " << dt << endl;
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
	         0, 1, 0, dt,
	         0, 0, 1, 0,
	         0, 0, 0, 1;

 
  //2. Set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << std::pow(dt, 4) / 4 * noise_ax, 0, std::pow(dt, 3) / 2 * noise_ax, 0,
	  0, std::pow(dt, 4) / 4 * noise_ay, 0, std::pow(dt, 3) / 2 * noise_ay,
	  std::pow(dt, 3) / 2 * noise_ax, 0, std::pow(dt, 2) * noise_ax, 0,
	  0, std::pow(dt, 3) / 2 * noise_ay, 0, std::pow(dt, 2) * noise_ay;

  ekf_.Predict();

  cout << "After Predict F " << ekf_.F_ << endl;
  cout << "After Predict Q " << ekf_.Q_ << endl;

  cout << "After Predict x_ = " << ekf_.x_ << endl;
  cout << "After Predict P_ = " << ekf_.P_ << endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  ekf_.R_ = this->R_radar_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
	  ekf_.R_ = this->R_laser_;
	  ekf_.H_ = this->H_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "After Update x_ = " << ekf_.x_ << endl;
  cout << "After Update P_ = " << ekf_.P_ << endl;

}
