/*
Recreation of Python Example of:

Least-Squares Fitting of Two 3-D Point Sets

by Arun, Huang, and Blostein 1987

from:

https://gist.github.com/scturtle/c3037529098338eccc403a9842870273

Evaluated here in cpp using the Eigen Library

by: John D. Long II, PhD

*/

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "eigen_mvn.h"
#include "utilities.h"

using namespace Lsq3D;

//////////
// MAIN //
//////////
int main()
{	
	//Load Lsq3D configuration
	Lsq3Dconfig config("leastSqs3Dparams.yaml");

	// Define the covariance matrix and the mean
	Eigen::VectorXf mu(3);
	Eigen::MatrixXf covB(3, 3);
	Eigen::MatrixXf covO(3, 3);

	//All noise is zero mean
	mu << 0.0f, 0.0f, 0.0f;

	//Construct Baseline Covariance matrix from variance and correlation coefficients
	covB << config.sigB[0], config.corB[0]*sqrt(config.sigB[0])*sqrt(config.sigB[1]), config.corB[1]*sqrt(config.sigB[0])*sqrt(config.sigB[2]),
	   		config.corB[0]*sqrt(config.sigB[0])*sqrt(config.sigB[1]), config.sigB[1], config.corB[2]*sqrt(config.sigB[1])*sqrt(config.sigB[2]),
			config.corB[1]*sqrt(config.sigB[0])*sqrt(config.sigB[2]), config.corB[2]*sqrt(config.sigB[1])*sqrt(config.sigB[2]), config.sigB[2];

	//Construct Outlier Covariance matrix
	covO << config.sigO[0], 0.0f, 0.0f,
			0.0f, config.sigO[1], 0.0f,
			0.0f, 0.0f, config.sigO[2];

	//Create Baseline Noise Sampler
	Mvn mvnB(3, config.N);
	mvnB.setParameters(mu, covB);

	//Create Outlier Noise Sampler
	Mvn mvnO(3, config.N);
	mvnO.setParameters(mu, covO);


	//Create Python Sample Data
	std::vector< Eigen::Vector3f > data;
	data.reserve(config.N);
	for (int ii = 0; ii < config.N; ii++)
	{
		Eigen::Vector3f vec;
		vec(0) = 1.0f + ii*(3.0f-1.0f)/(float)config.N;
		vec(1) = sin(vec(0));
		vec(2) = cos(vec(0));
		data.push_back(vec);
	}


	for (int ii = 0; ii < config.N; ii++)
	{
		std::cout << data[ii].transpose() << std::endl;
	}
	return 0;
}