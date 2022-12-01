#include <iostream>
#include <fstream>
#include <cstdlib>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <string>
#include <ctime>
#include <gsl/gsl_rng.h>

#include <DI_ComputeLengthDynamic.h>

int main(int argc, char* argv[])
{
	//Check number of parameters
	if (argc != 11 + 1)
	{
		std::cout << "Usage: " << argv[0] << " #Run, start_time, Max_time, ParameterSetIndex (kres=0.0001*ParameterSetIndex), ka, kd, kcat, dtwmax, load, savetime, ssTime" << std::endl;
		//std::cout << "Usage: " << argv[0] << " #Run, start_time, Max_time, firstParameterSetIndex1, finalParameterSetIndex1, firstParameterSetIndex2, finalParameterSetIndex2, firstParameterSetIndex3, finalParameterSetIndex3, IncrementIndex1, IncrementIndex2, IncrementIndex3, dtw1, dtw2, tm0, tmf, dtm, ts0, tsf, dts, load, savetime" << std::endl;
		return 1;
	}

	//Getting the date/time of execution
	//time_t now = time(0);
	//char* dt = ctime(&now);
	//std::cout << "The local date and time is: " << dt << std::endl;

	//Input parameters
	int run = std::stoi(argv[1]);
	double t = std::stod(argv[2]);
	double T = std::stod(argv[3]);
	int paramSet = std::stoi(argv[4]);
	double k_a = std::stod(argv[5]);
	double k_d = std::stod(argv[6]);
	double k_cat = std::stod(argv[7]);
	double k_res = double(paramSet)*0.0001;
	int dtwmax = std::stoi(argv[8]);
	int load = std::stoi(argv[9]);
	double savetime = std::stod(argv[10]);
	double ssTime = std::stod(argv[11]);

	double ratio = (k_a*k_res)/(k_d*k_cat);
	std::cout << "k_a = " << k_a << " k_d = " << k_d << " k_cat = " << k_cat << " k_res = " << k_res << std::endl;
	std::cout << "Control Parameter = " << ratio << std::endl;

	//Saving the parameters log in a file.
	std::string dirName("/home/joel/Documents/DI/Data/ParameterSet_" + std::to_string(paramSet) + "/Run_" + std::to_string(run));
	std::filesystem::create_directories(dirName); //create folder if doesn't exist.

	std::string parameters(dirName + "/Parameters.dat");
	std::ofstream parametersOut;
	parametersOut.open(parameters.c_str(), std::ofstream::app);
	parametersOut << "#k_a " << "\t" << "k_d" << "\t" << "k_cat" << "\t" << "k_res" << std::endl;
	parametersOut <<  k_a << "\t" << k_d << "\t" << k_cat << "\t" << k_res << std::endl;
	parametersOut.close();

	std::cout << "Run# " << run << std::endl;

	if(DI_ComputeLengthDynamic(run, k_a, k_d, k_cat, k_res, t, T, dtwmax, dirName, load, savetime, ssTime) == 1)
	{
		std::cout << "Error: exit programm." << std::endl;
		return 1;
	}

	return 0;
}
