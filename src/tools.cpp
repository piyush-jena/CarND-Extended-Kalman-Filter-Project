#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;
using std::endl;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */

   VectorXd rmse(4);
   rmse << 0,0,0,0;

   int N1 = estimations.size();
   int N2 = ground_truth.size();

   if (N1 == 0 || N1 != N2){
      cout << "Dimensional error in estimation vector or ground truth vector or both" << endl;
      return rmse;
   }

   for (int i = 0 ; i < N1 ; i++)
   {
      VectorXd diff = estimations[i] - ground_truth[i];
      rmse += (diff.array() * diff.array()).matrix();
   }

   rmse = rmse / N1;
   rmse = rmse.array().sqrt();

   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
   MatrixXd H(3,4);

   float px = x_state(0);
   float py = x_state(1);
   float vx = x_state(2);
   float vy = x_state(3);

   //Set the Jacobian to zero
   H << 0,0,0,0,
      0,0,0,0,
      0,0,0,0;

   float rho = pow((pow(px,2) + pow(py,2)), 0.5);
   if( rho < 0.0001){
    cout << "Value of rho too small - possible div by 0. Reassigning rho = 0.0005";
    rho = 0.0001;
   }

   float inv_rho = pow(rho, -1);
   H(0,0) = px * inv_rho;
   H(1,0) = -py * pow(inv_rho,2);
   H(2,0) = py * (vx*py - vy*px) * pow(inv_rho, 3);
   H(0,1) = py * inv_rho;
   H(1,1) = px * pow(inv_rho, 2);
   H(2,1) = px * (vy*px - vx*py) * pow(inv_rho,3);
   H(2,2) = H(0,0);
   H(2,3) = H(0,1);

   return H;
}
