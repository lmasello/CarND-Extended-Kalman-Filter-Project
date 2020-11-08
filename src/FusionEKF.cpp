#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing the kalman filter
  ekf_.x_ = VectorXd(4);
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.Q_ = MatrixXd(4, 4);

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  P_ = MatrixXd(4, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement projection matrix - lidar
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // process covariance matrix
  P_ <<  1, 0, 0, 0,
         0, 1, 0, 0,
         0, 0, 1000, 0,
         0, 0, 0, 1000;

  // Initialize the acceleration noise components according to values suggested by Udacity
  noise_ax_ = 9;
  noise_ay_ = 9;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Udacity note: Although radar gives velocity data in the form of the range rate rho_dot,
      // a radar measurement does not contain enough information to determine the state variable velocities
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      ekf_.x_ << rho * cos(phi), rho * sin(phi), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1],
                  0, 0;
    }

    ekf_.P_ << P_;
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  double dt2 = dt * dt;
  double dt3 = dt2 * dt;
  double dt4 = dt3 * dt;
  previous_timestamp_ = measurement_pack.timestamp_;
  // Update the transition matrix with dt
  ekf_.F_ <<  1, 0, dt, 0,
              0, 1, 0, dt,
              0, 0, 1, 0,
              0, 0, 0, 1;
  // Set the process covariance matrix
  ekf_.Q_ <<  dt4 / 4 * noise_ax_, 0, dt3 / 2 * noise_ax_, 0,
              0, dt4 / 4 * noise_ay_, 0, dt3 / 2 * noise_ay_,
              dt3 / 2 * noise_ax_, 0, dt2 * noise_ax_, 0,
              0, dt3 / 2 * noise_ay_, 0, dt2 * noise_ay_;
  ekf_.Predict();

  /**
   * Update
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.R_ = MatrixXd(3, 3);
    ekf_.R_ << R_radar_;
    ekf_.H_ = MatrixXd(3, 4);
    ekf_.H_ << tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateWithRadar(measurement_pack.raw_measurements_);
  } else {
    ekf_.R_ = MatrixXd(2, 2);
    ekf_.R_ << R_laser_;
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << H_laser_;
    ekf_.UpdateWithLidar(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
