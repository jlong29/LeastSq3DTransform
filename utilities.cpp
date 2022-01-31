#include "utilities.h"

namespace Lsq3D
{

int Lsq3Dconfig::loadConfiguration(std::string const& filename)
{
	if (initialized)
	{
		fprintf(stderr, "Lsq3Dconfig::loadConfiguration: configuration already initialzed\n");
		return -1;
	} else
	{
		initialized = parseYaml(filename);
		if (!initialized)
		{
			fprintf(stderr, "Lsq3Dconfig::loadConfiguration: Error parsing yaml file\n");
		return -1;
		}
	}
	return 0;
}

bool Lsq3Dconfig::parseYaml(std::string const& filename)
{

	fprintf(stdout, "\nLoading Configuration File...\n");
	YAML::Node config = YAML::LoadFile(filename);

	if (config["Lsq3Dconfig"])
	{
		std::cout << "LEAST SQUARES 3D PARAMS:" << std::endl;
		for(YAML::const_iterator it=config["Lsq3Dconfig"].begin();it != config["Lsq3Dconfig"].end();++it)
		{
			//Data Directory
			if (getItKeyValueStrExpandEnvVar(it, "dataDir", dataDir)){
				continue;
			}
			//Get N
			if (it->first.as<std::string>() == "N")
			{
				N = it->second.as<int>();
				printKeyValue(it->first.as<std::string>(), N);
				continue;
			}
			//Get sigB
			if (it->first.as<std::string>() == "sigB")
			{
				for (YAML::const_iterator vit=it->second.begin(); vit!=it->second.end();++vit)
				{
					sigB.push_back(vit->as<float>());
				}
				printKeyValue(it->first.as<std::string>(), sigB);
			}
			//Get corB
			if (it->first.as<std::string>() == "corB")
			{
				for (YAML::const_iterator vit=it->second.begin(); vit!=it->second.end();++vit)
				{
					corB.push_back(vit->as<float>());
				}
				printKeyValue(it->first.as<std::string>(), corB);
			}
			//Get ProbO
			if (it->first.as<std::string>() == "ProbO")
			{
				ProbO = it->second.as<float>();
				printKeyValue(it->first.as<std::string>(), ProbO);
				continue;
			}
			//Get sigO
			if (it->first.as<std::string>() == "sigO")
			{
				for (YAML::const_iterator vit=it->second.begin(); vit!=it->second.end();++vit)
				{
					sigO.push_back(vit->as<float>());
				}
				printKeyValue(it->first.as<std::string>(), sigO);
			}
			//Get Euler Angles
			if (it->first.as<std::string>() == "EA")
			{
				for (YAML::const_iterator vit=it->second.begin(); vit!=it->second.end();++vit)
				{
					EA.push_back(vit->as<float>());
				}
				printKeyValue(it->first.as<std::string>(), EA);
			}
		}
	}
	std::cout << std::endl;

	return true;
}

//Checks a string for Enviroment Variable references, ${VarName} and expands
int Lsq3Dconfig::expandEnvVar(std::string& text)
{
	//Search tags
	std::string startTag = "${";
	std::string endTag   = "}";
	
	//Look for startTag
	size_t startPos = text.find(startTag);
	
	//Process Environment Variable Reference from left-to-right
	while (startPos!=std::string::npos) {
		
		//Look for endTag
		size_t endPos = text.find(endTag);
		if (endPos == std::string::npos){
			fprintf(stderr, "Lsq3D:expandEnvVar: No endTag \"}\" match...aborting\n");
			return -1;
		}
		//Catch Order Error
		if (endPos < startPos)
		{
			fprintf(stderr, "Lsq3D:expandEnvVar: SynTax Error! endTag %s before startTag %s ...aborting\n", endTag.c_str(), startTag.c_str());
			return -1;
		}

		//VARIABLE_NAME starts +2 after startPos and goes until endPos-1
		std::string varName = text.substr(startPos+2, (endPos-1) - (startPos+1));
		const char* var = getenv(varName.c_str());
		if (var==NULL){
			fprintf(stderr, "Lsq3D:expandEnvVar: Could not find Env Var: %s ...aborting\n", varName.c_str());
			return -1;
		}

		//Construct tmp string to contain update
		std::string tmp = text.substr(0,startPos) + var + text.substr(endPos+1);

		//Update text with tmp
		text  = tmp;
		
		//Look for next Environment Variable
		startPos = text.find("${");
	}
	return 0;
}

} //namespace kalmanFilter