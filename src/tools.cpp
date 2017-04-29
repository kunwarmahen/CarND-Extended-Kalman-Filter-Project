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

	//check division by zero
	if (fabs(c1) < 0.0001) {
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px / c2), (py / c2), 0, 0,
		-(py / c1), (px / c1), 0, 0,
		py*(vx*py - vy*px) / c3, px*(px*vy - py*vx) / c3, px / c2, py / c2;

	return Hj;
}

VectorXd Tools::CalculateHofX(const VectorXd& x_state) {

	MatrixXd h(3, 1);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero

	if (px == 0 || py == 0) {
		cout << "px or py cannot be zero";
		return h;
	}

	//compute the Jacobian matrix
	float px2 = px * px;
	float py2 = py * py;

	h << sqrt(px2 + py2),
		constrainAngle(atan2(py,px)),
		(px * vx) + (py * vy) / sqrt(px2 + py2);

	return h;
}

void Tools::printDimension(MatrixXd x) {
	cout << x.rows() << "x" << x.cols() << endl;
}

void Tools::printMatrix(MatrixXd x) {
	cout << x << endl;
}

double Tools::constrainAngle(double x) {
/*	if (x > M_PI)
		x = fmod(x - M_PI, 2 * M_PI) - M_PI;
	if (x < -M_PI)
		x = fmod(x + M_PI, 2 * M_PI) + M_PI;
		*/
	while (x < -M_PI) {
		x += 2 * M_PI;
	}
	while (x > M_PI) {
		x -= 2 * M_PI;
	}

	if (x > M_PI || x < -M_PI)
		cout << "Error" << endl;

	return x;
}