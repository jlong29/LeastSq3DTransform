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
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

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
float radians2degrees(float rad)
{
	return rad*(180.0f/PI);
}

// Forward Declarations
//Mean over each dimension of 3D vector
Eigen::Vector3f meanVectorEigenV3f( std::vector< Eigen::Vector3f > const& data);

template<typename T>
std::vector<float> linspace(T start_in, T end_in, int num_in);

template<typename T>
void print_vector(const std::string& vecName, std::vector<T> vec);

// Function to return the next random number
int getNum(std::vector<int>& v);

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

	// 3D Transformation
	Eigen::Vector3f EA, T;
	//Transform Euler Angles from degrees to radians
	std::transform(config.EA.begin(), config.EA.end(), config.EA.begin(), degrees2radians);
	EA << config.EA[0], config.EA[1], config.EA[2];
	T << config.T[0], config.T[1], config.T[2];

	//Create Rotation matrix from Euler Angles
	Eigen::Matrix3f R;
	composeRotationFromEuler(EA, R);

	//Create Python Sample Data (p1) and Transformed Data (p2)
	std::vector<float> seed = linspace(1,3,N);

	std::vector< Eigen::Vector3f > p1;
	p1.reserve(N);
	for (int ii = 0; ii < N; ii++)
	{
		Eigen::Vector3f vec;
		vec(0) = seed[ii];
		vec(1) = sin(vec(0));
		vec(2) = cos(vec(0));
		p1.push_back(vec);
	}

	//Create Transformed Sample Data
	std::vector< Eigen::Vector3f > p2;
	p2.reserve(N);

	for (int ii = 0; ii < N; ii++)
	{
		Eigen::Vector3f vec = R*p1[ii] + T;
		p2.push_back(vec);
	}

	//Check whether to add baseline and outlier noise to the transformed data (p2)
	if (config.addNoise)
	{
		//All noise, baseline and outliers is zero mean
		Eigen::VectorXf mu(3);
		mu << 0.0f, 0.0f, 0.0f;

		//Create Baseline Noise Sampler
		Mvn mvnB(3, N);

		// Define the covariance matrix
		Eigen::MatrixXf covB(3, 3);


		//Construct Baseline Covariance matrix from variance and correlation coefficients
		covB << sigB[0], corB[0]*sqrt(sigB[0])*sqrt(sigB[1]), corB[1]*sqrt(sigB[0])*sqrt(sigB[2]),
				corB[0]*sqrt(sigB[0])*sqrt(sigB[1]), sigB[1], corB[2]*sqrt(sigB[1])*sqrt(sigB[2]),
				corB[1]*sqrt(sigB[0])*sqrt(sigB[2]), corB[2]*sqrt(sigB[1])*sqrt(sigB[2]), sigB[2];

		std::cout << "covB:\n" << covB << std::endl << std::endl;

		//Set Parameters
		mvnB.setParameters(mu, covB);
		mvnB.sample();

		//Add Baseline Noise to p1
		for (int ii = 0; ii < N; ii++)
		{
			Eigen::Vector3f vec;
			vec << mvnB.Z[3*ii + 0], mvnB.Z[3*ii + 1], mvnB.Z[3*ii + 2];
			p2[ii] += vec;
		}

		if (config.addOutliers)
		{
			//Set the number of outliers
			int Noutliers = (int)((float)N*config.ProbO);
			std::cout << "Outliers: " << Noutliers << " out of " << N << std::endl;

			//Create Outlier Noise Sampler
			Mvn mvnO(3, Noutliers);
			//Construct Outlier Covariance matrix
			Eigen::MatrixXf covO(3, 3);
			covO << sigO[0], 0.0f, 0.0f,
				0.0f, sigO[1], 0.0f,
				0.0f, 0.0f, sigO[2];
			
			std::cout << "covO:\n" << covO << std::endl << std::endl;

			//Set parameters
			mvnO.setParameters(mu, covO);

			/* initialize random seed: */
			srand (time(NULL));

			//Set outliers
			mvnO.sample();
			std::vector<int> indices;
			for (int ii = 0; ii < N; ii++)
				indices.push_back(ii);

			for (int ii = 0; ii < Noutliers; ii++)
			{
				Eigen::Vector3f vec;
				vec << mvnO.Z[3*ii + 0], mvnO.Z[3*ii + 1], mvnO.Z[3*ii + 2];

				//Generate non-repeating random index between 0 and N
				int idx = getNum(indices);
				p2[idx] += vec;
			}
		}
	}

	if (config.inspectData)
	{
		print_vector<float>("seed", seed);
		std::cout << "Euler Angles (radians): " << EA.transpose() << std::endl;
		std::cout << "Translation: " << T.transpose() << std::endl;

		std::cout << "Inspect p1 and Transformed p2:" << std::endl;
		for (int ii = 0; ii < N; ii++)
		{
			std::cout << "\tp1_" << ii << ": " << p1[ii].transpose() << std::endl;
			std::cout << "\tp2_" << ii << ": " << p2[ii].transpose() << std::endl<<std::endl;
		}
	}

	//Compute the mean vector of data and p2
	Eigen::Vector3f p1c = meanVectorEigenV3f(p1);
	Eigen::Vector3f p2c = meanVectorEigenV3f(p2);

	std::cout << "p1 mean: " << p1c.transpose() << std::endl;
	std::cout << "p2 mean: " << p2c.transpose() << std::endl;

	//Subtract the means from each respective set of data
	std::vector< Eigen::Vector3f > q1;
	std::vector< Eigen::Vector3f > q2;
	q1.reserve(N); q2.reserve(N);

	for (int ii = 0; ii < N; ii++)
	{
		q1.push_back(p1[ii] - p1c);
		q2.push_back(p2[ii] - p2c);
	}

	if (config.inspectData)
	{
		std::cout << "Mean subtracted q1 and Transformed q2:" << std::endl;
		for (int ii = 0; ii < N; ii++)
		{
			std::cout << "\tq1_" << ii << ": " << q1[ii].transpose() << std::endl;
			std::cout << "\tq2_" << ii << ": " << q2[ii].transpose() << std::endl<<std::endl;
		}
	}
	//Compute H matrix as input to SVD
	Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
	for (int ii = 0; ii < N; ii++)
	{
		H += q1[ii]*q2[ii].transpose();
	}

	//SVD
	Eigen::Matrix3f U, Ut, V;   //Orthogonal matrices
	Eigen::JacobiSVD<Eigen::Matrix3f> svd(H, Eigen::ComputeFullV | Eigen::ComputeFullU);
	U  = svd.matrixU();
	V  = svd.matrixV();

	//Getting the combo transposes right took some guessing
	Ut = U.transpose();
	Eigen::Matrix3f R2 = V*Ut;

	Eigen::Vector3f EA2;
	getEulerAngles(R2, EA2);

	std::cout << "Inspect H and Compare Initial to Estimated Transformation:\n" << std::endl;
	std::cout << "H:\n" << H << std::endl;
	std::cout << "R:\n" << R << std::endl;
	std::cout << "\nR2:\n" << R2 << std::endl;
	std::cout << "det(R2) = " << R2.determinant() << std::endl << std::endl;

	std::cout << "EA:  " << radians2degrees(EA(0)) << ", " << radians2degrees(EA(1)) << ", " << radians2degrees(EA(2)) << std::endl;
	std::cout << "EA2: " << radians2degrees(EA2(0)) << ", " << radians2degrees(EA2(1)) << ", " << radians2degrees(EA2(2)) << std::endl << std::endl;

	std::cout << "T  = " << T.transpose() << std::endl;
	Eigen::Vector3f T2 = p2c - R2*p1c;
	std::cout << "T2 = " << T2.transpose() << std::endl << std::endl;

	//Create Estimate of Transformed Sample Data
	if (config.inspectData)
	{
		std::vector< Eigen::Vector3f > p3;
		p3.reserve(N);

		for (int ii = 0; ii < N; ii++)
		{
			Eigen::Vector3f vec = R2*p1[ii] + T2;
			p3.push_back(vec);
		}
		std::cout << "Inspect p1, Transformed p2, and Estimated p3:" << std::endl;
		for (int ii = 0; ii < N; ii++)
		{
			std::cout << "\tp1_" << ii << ": " << p1[ii].transpose() << std::endl;
			std::cout << "\tp2_" << ii << ": " << p2[ii].transpose() << std::endl;
			std::cout << "\tp3_" << ii << ": " << p3[ii].transpose() << std::endl<<std::endl;
		}
	}

	return 0;
}

//Mean over each dimension of 3D vector
Eigen::Vector3f meanVectorEigenV3f( std::vector< Eigen::Vector3f > const& data)
{
	Eigen::Vector3f mean;
	mean << 0.0f, 0.0f, 0.0f;
	for (size_t ii = 0; ii < data.size(); ii++)
	{
		mean(0) += data[ii](0);
		mean(1) += data[ii](1);
		mean(2) += data[ii](2);
	}
	mean(0) /= data.size();
	mean(1) /= data.size();
	mean(2) /= data.size();

	return mean;
}

//Like Numpy's linspace
template<typename T>
std::vector<float> linspace(T start_in, T end_in, int num_in)
{

  std::vector<float> linspaced;

  float start = static_cast<float>(start_in);
  float end = static_cast<float>(end_in);
  float num = static_cast<float>(num_in);

  if (num == 0) { return linspaced; }
  if (num == 1) 
	{
	  linspaced.push_back(start);
	  return linspaced;
	}

  float delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
	{
	  linspaced.push_back(start + delta * i);
	}
  linspaced.push_back(end); // I want to ensure that start and end
							// are exactly the same as the input
  return linspaced;
}

template<typename T>
void print_vector(const std::string& vecName, std::vector<T> vec)
{
  std::cout << vecName << " of size " << vec.size() << std::endl;
  for (T d : vec)
	std::cout << d << " ";
  std::cout << std::endl;
}

// Function to return the next random number
int getNum(std::vector<int>& v)
{
 
    // Size of the vector
    int n = v.size();
 
    // Generate a random number
    srand(time(NULL));
 
    // Make sure the number is within
    // the index range
    int index = rand() % n;
 
    // Get random number from the vector
    int num = v[index];
 
    // Remove the number from the vector
    std::swap(v[index], v[n - 1]);
    v.pop_back();
 
    // Return the removed number
    return num;
}