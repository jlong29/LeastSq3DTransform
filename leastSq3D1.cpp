/*
Recreation of Python Example of:

Least-Squares Fitting of Two 3-D Point Sets

by Arun, Huang, and Blostein 1987

from:

https://gist.github.com/scturtle/c3037529098338eccc403a9842870273

Evaluated here in cpp using the Eigen Library

by: John D. Long II, PhD

*/

// STL
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>

//Project Utilities
#include "utilities.h"

//Other Projects
#include "eigen_mvn.h"
#include "lieAlgebras.h"

using namespace Lsq3D;
using namespace lieAlgebras;

// UTILS
float degrees2radians(float deg)
{
  return deg*(PI/180.0f);
}

//////////
// MAIN //
//////////
int main()
{
	//Load Lsq3D configuration
	Lsq3Dconfig config("leastSqs3Dparams.yaml");

	//Copy in configuration data for ease of use
	int N = config.N;
	std::vector<float> sigB(config.sigB.begin(), config.sigB.end());
	std::vector<float> corB(config.corB.begin(), config.corB.end());
	std::vector<float> sigO(config.sigO.begin(), config.sigO.end());

	//Transform Euler Angles from degrees to radians
	std::transform(config.EA.begin(), config.EA.end(), config.EA.begin(), degrees2radians);
	Eigen::Vector3f EA;
	EA << config.EA[0], config.EA[1], config.EA[2];
	std::cout << "Euler Angles (radians): " << EA.transpose() << std::endl;

	// Define the covariance matrix and the mean
	Eigen::VectorXf mu(3);
	Eigen::MatrixXf covB(3, 3);
	Eigen::MatrixXf covO(3, 3);

	//All noise is zero mean
	mu << 0.0f, 0.0f, 0.0f;

	//Construct Baseline Covariance matrix from variance and correlation coefficients
	covB << sigB[0], corB[0]*sqrt(sigB[0])*sqrt(sigB[1]), corB[1]*sqrt(sigB[0])*sqrt(sigB[2]),
			corB[0]*sqrt(sigB[0])*sqrt(sigB[1]), sigB[1], corB[2]*sqrt(sigB[1])*sqrt(sigB[2]),
			corB[1]*sqrt(sigB[0])*sqrt(sigB[2]), corB[2]*sqrt(sigB[1])*sqrt(sigB[2]), sigB[2];

	std::cout << "\ncovB:\n" << covB << std::endl;

	//Construct Outlier Covariance matrix
	covO << sigO[0], 0.0f, 0.0f,
			0.0f, sigO[1], 0.0f,
			0.0f, 0.0f, sigO[2];

	//Create Baseline Noise Sampler
	Mvn mvnB(3, N);
	mvnB.setParameters(mu, covB);

	//Create Outlier Noise Sampler
	Mvn mvnO(3, N);
	mvnO.setParameters(mu, covO);

	//Create Python Sample Data
	std::vector< Eigen::Vector3f > data;
	data.reserve(N);
	for (int ii = 0; ii < N; ii++)
	{
		Eigen::Vector3f vec;
		vec(0) = 1.0f + ii*(3.0f-1.0f)/(float)N;
		vec(1) = sin(vec(0));
		vec(2) = cos(vec(0));
		data.push_back(vec);
	}

	// //Create Rotation matrix from Euler Angles
	Eigen::Matrix3f R;
	composeRotationFromEuler(EA, R);

	std::cout << "\nR:\n" << R << std::endl;

	// for (int ii = 0; ii < config.N; ii++)
	// {
	// 	std::cout << data[ii].transpose() << std::endl;
	// }
	return 0;
}