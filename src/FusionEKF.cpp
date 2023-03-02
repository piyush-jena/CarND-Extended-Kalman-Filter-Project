#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_      = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
                0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0,      0,
                0,    0.0009, 0,
                0,    0,      0.09;

    //Initalize the F, P & Q matrices. 
    ekf_.F_ = MatrixXd(4,4);
    ekf_.P_ = MatrixXd(4,4);
    ekf_.Q_ = MatrixXd(4,4);
    
    //Assign the values for the laser sensor matrix, H_laser
    H_laser_ << 1,0,0,0,
                    0,1,0,0;
    
    //Set the values for acceleration noise
    noise_ax = 9;
    noise_ay = 9; 
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
        // Initializing the state vector and setting initial state
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 0.5, 0.5;
    
    
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            // RADAR gives measurements in Polar coordinate system

            float rho = measurement_pack.raw_measurements_[0];
            float theta = measurement_pack.raw_measurements_[1];
            float rho_dot = measurement_pack.raw_measurements_[2]; 

            ekf_.x_(0) = rho*cos(theta);
            ekf_.x_(1) = rho*sin(theta);
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            float px = measurement_pack.raw_measurements_[0];
            float py = measurement_pack.raw_measurements_[1];
        
            ekf_.x_(0) = px;
            ekf_.x_(1) = py;
        }

        previous_timestamp_ = measurement_pack.timestamp_;
        //assign initial values to the state transition matrix, F
        ekf_.F_ << 1,0,1,0,
                0,1,0,1,
                0,0,1,0,
                0,0,0,1;    
        //assign initial values to the covariance matrix, P. Adjust the variance values to reflect uncertainty in initial state
        ekf_.P_  << 1,0,  0,0,
                    0,1,  0,0,
                    0,0,500,0,
                    0,0,  0,500;

        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    float dt2 = dt*dt;
    float dt3 = dt2*dt;
    float dt4 = dt3*dt;

    // Update the state transition matrix, F
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;

    // Set the process covariance matrix, Q
    ekf_.Q_ <<  dt4/4*noise_ax, 0, dt3/2*noise_ax, 0,
	            0, dt4/4*noise_ay, 0, dt3/2*noise_ay,
	            dt3/2*noise_ax, 0, dt2*noise_ax, 0,
	            0, dt3/2*noise_ay, 0, dt2*noise_ay;

    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        //Calculate the Jacobian matrix about the current predicted state and set the EKF state transition matrix, H
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = Hj_;

        //Initialize the EKF object measurement covariance matrix, R, to the right size and assign the correct values
        ekf_.R_ = MatrixXd(3,3);
        ekf_.R_ = R_radar_; 
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        ekf_.H_ = H_laser_;

        //Initialize the EKF object measurement covariance matrix, R, to the right size and assign the correct values
        ekf_.R_ = MatrixXd(2,2);
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
