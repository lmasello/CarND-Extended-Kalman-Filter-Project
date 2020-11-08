#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &y) {
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  // New estimate
  x_ = x_ + (K * y);
  P_ -= K * H_ * P_;
}

void KalmanFilter::UpdateWithLidar(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  Update(y);
}

void KalmanFilter::UpdateWithRadar(const VectorXd &z) {
   double px = x_(0);
   double py = x_(1);
   double vx = x_(2);
   double vy = x_(3);
   double rho = sqrt(px * px + py * py);
   double phi = atan2(py, px);  // values between -pi and pi.
   double rho_dot = (px * vx + py * vy) / rho;

   VectorXd h = VectorXd(3);
   h << rho, phi, rho_dot;
   VectorXd y = z - h;
   // Adjust the phi angle in y so that it is between -pi and pi
   y(1) = tools.WrapMinMax(y(1), -M_PI, M_PI);
   Update(y);
}
