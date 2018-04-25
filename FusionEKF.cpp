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
  Hj_.setZero();
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
  H_laser_ << 1,0,0,0,
              0,1,0,0;
  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 8, 6, 5, 2;
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rhodot = measurement_pack.raw_measurements_[2];
      ekf_.x_ << rho * cos(phi), rho * sin(phi), 0, 0;
      ekf_.R_ = MatrixXd(3, 3);
      ekf_.R_ = R_radar_; // R_ is measurement noise is constant for sensor
      ekf_.H_ = MatrixXd(3, 4);
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);  // H_ for Radar going to be changed sinse each time it need to be calculated
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      double px = measurement_pack.raw_measurements_[0];
      double py = measurement_pack.raw_measurements_[1];
      ekf_.x_ << px, py, 0, 0;
      ekf_.R_ = MatrixXd(2, 2);
      ekf_.R_= R_laser_;
      ekf_.H_ = MatrixXd(2, 4);
      ekf_.H_ = H_laser_;   //H_ for LAser is going to be the same
    }

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    cout << " Initialization done " << endl;
    return;
  }
  cout << "calculating dt"<<endl;
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.;
  cout << " dt coming as " << dt << endl;
  float noise_ax = 9;
  float noise_ay = 9;
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1,0,dt,0,
             0,1,0,dt,
             0,0,1,0,
             0,0,0,1;

  float dt_2 = dt * dt;
	float dt_3 = dt_2 * dt;
	float dt_4 = dt_3 * dt;


	//set the process covariance matrix Q
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


  ekf_.Predict();
  cout << " Prediction done " << endl;
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
    ekf_.R_ = MatrixXd(3, 3);
    ekf_.R_ = R_radar_;
    ekf_.H_ = MatrixXd(3, 4);
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.R_ = MatrixXd(2, 2);
    ekf_.R_= R_laser_;
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  previous_timestamp_ = measurement_pack.timestamp_;
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;



}
