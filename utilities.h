#ifndef UTILITIES_H__
#define UTILITIES_H__

#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <chrono>

//CHRONO
typedef std::chrono::steady_clock::time_point timepoint;
typedef std::chrono::duration<double, std::milli> durMilli;

#include "yaml-cpp/yaml.h"

namespace Lsq3D
{

//Class for reading in Least Squares 3D Configuration
class Lsq3Dconfig
{
public:

	//Data Directory for NFK
	std::string dataDir;

	//The number of datapoints
	int N;

	//Baseline Noise
	std::vector<float> sigB;
	//Baseline Coorelation: cor12, cor13, cor23
	std::vector<float> corB;

	//Outlier Probability
	float ProbO;
	//Outlier Noise
	std::vector<float> sigO;

	// 3D Transformation
	//Euler Angles for Rotation Matrix (degrees)
	std::vector<float> EA;
	//Translation
	std::vector<float> T;

	//Constructors
	Lsq3Dconfig():initialized(false){};

	Lsq3Dconfig(std::string const& filename):
		initialized(false)
	{
		loadConfiguration(filename);
	};

	//Deconstructor
	~Lsq3Dconfig(){};

	//Checks if initialized and if not, loads configuration
	int loadConfiguration(std::string const& filename);
private:
	bool initialized;

	//Object-specific yaml file parser
	bool parseYaml(std::string const& filename);
	//Checks a string for Enviroment Variable references, ${VarName} and expands
	int expandEnvVar(std::string& text);

	//Helpers
	void printKeyValue(std::string key, uint value)
	{
		std::cout << key << ": " << value << std::endl;
	}
	void printKeyValue(std::string key, int value)
	{
		std::cout << key << ": " << value << std::endl;
	}
	void printKeyValue(std::string key, float value)
	{
		std::cout << key << ": " << value << std::endl;
	}
	void printKeyValue(std::string key, std::vector<float> value)
	{
		std::cout << key << ": ";
		for (auto val : value){ std::cout << val << ", "; }; std::cout << std::endl;
	}
	void printKeyValue(std::string key, std::string value)
	{
		std::cout << key << ": " << value << std::endl;
	}
	
	bool getItKeyValueStrExpandEnvVar(YAML::const_iterator it, const char* key, std::string& returnStr)
	{
		if (it->first.as<std::string>() == key)
		{
			returnStr = it->second.as<std::string>();
			expandEnvVar(returnStr);
			printKeyValue(it->first.as<std::string>(), returnStr);
			return true;
		}
		return false;
	}
};

} //namespace Lsq3D

#endif  //UTILITIES_H_