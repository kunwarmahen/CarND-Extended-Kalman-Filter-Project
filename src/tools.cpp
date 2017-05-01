#include <iostream>
#include "tools.h"
#include <cmath>
#include <stdlib.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if (estimations.size() != ground_truth.size()
		|| estimations.size() == 0) {
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); ++i) {

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse / estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	MatrixXd Hj(3, 4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px + py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);
	float c4 = py*(vx*py - vy*px);
	float c5 = px*(px*vy - py*vx);

	//check division by zero
	if (fabs(c1) < 0.001) {
		//cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		c1 = 0.001;
	}

	if (fabs(c3) < 0.001) {
		//cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		c3 = 0.001;
	}


	//compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		 (c4/ c3), (c5 / c3), (px / c2), (py / c2);

	return Hj;
}

VectorXd Tools::CalculateHofX(const VectorXd& x_state) {

	VectorXd h(3);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero

	if (px == 0) {
		//cout << "px is close to zero" << endl;
		px = 0.1;
	} 

	if (py == 0) {
		//cout << "py is close to zero" << endl;
		py = 0.1;
	}


//compute the Jacobian matrix
	float px2 = px * px;
	float py2 = py * py;
	float rho = sqrt(px2 + py2);
	float mult = (px * vx) + (py * vy);
	float phi = checkPIValue(atan2(py, px));
	float roh_dot = mult / rho;

	if (px2 <= 0.001) {
		//cout << "px2 is close to zero" <<endl;
		px2 = 0.001;
	}

	if (py2 <= 0.001) {
		//cout << "py2 is close to zero" << endl;
		py2 = 0.001;
	}

	if (rho <= 0.0001) {
		//cout << "sqrt of px2 and py2 is close to zero" << endl;
		rho = 0.0001;

	}

	h << rho, phi, roh_dot;

	return h;
}

void Tools::printDimension(MatrixXd x) {
	cout << x.rows() << "x" << x.cols() << endl;
}

void Tools::printMatrix(MatrixXd x) {
	cout << x << endl;
}

double Tools::checkPIValue(double x) {
	while (x < -M_PI) {
		//cout << "checkPIValue" << endl;
		x += 2 * M_PI;
	}
	while (x > M_PI) {
		//cout << "checkPIValue" << endl;
		x -= 2 * M_PI;
	}

	if (x > M_PI || x < -M_PI) {
		cout << "Error" << endl;
	}

	return x;
}