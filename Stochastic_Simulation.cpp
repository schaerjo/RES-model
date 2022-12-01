#include <iostream>
#include <fstream>
#include <cstdlib>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <string>
#include <ctime>
#include <gsl/gsl_rng.h>

#include <ComputeLengthDynamic7.h>

int main(int argc, char* argv[])
{
	//Check number of parameters
	if (argc != 15 + 1)
	{
		std::cout << "Usage: " << argv[0] << " #Run, start_time, Max_time, ParameterSetIndex1, ParameterSetIndex2, ParameterSetIndex3, dtwmax, tm0, dtm, ts0, dts, load, savetime, ssTime (in units of Max_time), ParameterSet" << std::endl;
		//std::cout << "Usage: " << argv[0] << " #Run, start_time, Max_time, firstParameterSetIndex1, finalParameterSetIndex1, firstParameterSetIndex2, finalParameterSetIndex2, firstParameterSetIndex3, finalParameterSetIndex3, IncrementIndex1, IncrementIndex2, IncrementIndex3, dtw1, dtw2, tm0, tmf, dtm, ts0, tsf, dts, load, savetime" << std::endl;
		return 1;
	}

	//Getting the date/time of execution
	time_t now = time(0);
	char* dt = ctime(&now);
	//std::cout << "The local date and time is: " << dt << std::endl;

	//Input parameters
	//int start_run = std::stoi(argv[1]);
	//int end_run = start_run + std::stoi(argv[2]);
	int run = std::stoi(argv[1]);
	double t = std::stod(argv[2]);
	double T = std::stod(argv[3]);
	int paramSet1 = std::stoi(argv[4]);
	//int firstParamSet1 = std::stoi(argv[4]);
	//int finalParamSet1 = std::stoi(argv[5]);
	int paramSet2 = std::stoi(argv[5]);
	//int firstParamSet2 = std::stoi(argv[6]);
	//int finalParamSet2 = std::stoi(argv[7]);
	int paramSet3 = std::stoi(argv[6]);
	//int firstParamSet3 = std::stoi(argv[8]);
	//int finalParamSet3 = std::stoi(argv[9]);
	//int Increment1 = std::stoi(argv[10]);
	//int Increment2 = std::stoi(argv[11]);
	//int Increment3 = std::stoi(argv[12]);
	int dtwmax = std::stoi(argv[7]);
	int tm0 = std::stoi(argv[8]);
	int tmf = (int) T;
	//int tmf = std::stoi(argv[10]);
	int dtm = std::stoi(argv[9]);
	int ts0 = std::stoi(argv[10]);
	int tsf = (int) T;
	//int tsf = std::stoi(argv[19]);
	int dts = std::stoi(argv[11]);
	int load = std::stoi(argv[12]);
	double savetime = std::stod(argv[13]);
	double ssTime = std::stod(argv[14]);
	int set = std::stoi(argv[15]);
	double motorTime = 0;
	//int i=0;

	//Parameter Set
	double k_a = 0; //attachment rate
	double k_d = 0; //detachment rate
	double k_cat = 0; //catastroph rate
	double p_res = 0; //rescue probability if tip state is GTPX
	double k_hyd = 0; //hydrolization rate for +end attached tubulin
	double k_hydx = 0;//hydrolization rate for repaired tubulin
	double k_rep = 0;//rate of tubulin repairement
	double lambda = 0;//Motors concentration is given by c=lambda*c0 where c0=1nM
	double k_ma0 = 0;//Motors attachment rate for 1 nM concentration
	double k_ma = 0;//Motors attachment rate
	double k_md = 0;//Motors detachment rate
	double k_mhop = 0;//Motors hopping rate
	long double p_x = 0;//Probability of repair upon hopping
	int Lmax = 0;

	/*for(int paramSet1 = firstParamSet1; paramSet1 <= finalParamSet1; paramSet1 += Increment1)
	{
	for(int paramSet2 = firstParamSet2; paramSet2 <= finalParamSet2; paramSet2 += Increment2)
	{
	for(int paramSet3 = firstParamSet3; paramSet3 <= finalParamSet3; paramSet3 += Increment3)
	{*/
		//std::cout << "paramSet1 = " << paramSet1 << std::endl;
		//std::cout << "paramSet2 = " << paramSet2 << std::endl;
		//std::cout << "paramSet3 = " << paramSet3 << std::endl;
		//Determine the set of parameter
		k_a = 2; //per s
		k_d = 27; //per s
		k_cat = 0.005; //per s
		//p_res = 0.15;
		p_res = 0.15;
		//k_hyd = 0.05; //per s
		k_hyd = 0.125;
		k_hydx = 0.000125; //per s
		//k_hydx = 0.0000001*paramSet3; //per s
		//k_rep = 0.0;//per s
		k_rep = 0.00001;//per s
		lambda = 0.1*double(paramSet1);
		k_ma0 = 0.000033;//Motors attachment rate for 1 nM concentration per s
		k_ma = lambda*k_ma0;//per s
		k_md = 0.09;//per s
		k_md = 0.001*double(paramSet3);
		k_mhop = 16;//per s
		p_x = 0.00001*paramSet2;
		//p_x = 0.0018;
		Lmax=-1;
		//Lmax = paramSet2;

		//Saving the parameters log in a file.
		std::string dirName("/home/joel/Documents/LDMIX/Data/set" + std::to_string(set) + "/ParameterSet_" + std::to_string(paramSet1) + "_" + std::to_string(paramSet2) + "_" + std::to_string(paramSet3) + "/Run_" + std::to_string(run));
		std::filesystem::create_directories(dirName); //create folder if doesn't exist.
		std::string parameters(dirName + "/Parameters.dat");
		std::ofstream parametersOut;
		parametersOut.open(parameters.c_str(), std::ofstream::app);
		parametersOut << "Date: " << dt << std::endl << " #run: " << run << std::endl << " ka= " << k_a << std::endl << " kd= " << k_d << std::endl << " kcat= " << k_cat << std::endl << " pres= " << p_res << std::endl << " khyd= " << k_hyd << std::endl << " khydx= " << k_hydx << std::endl << " k_rep= " << k_rep << std::endl << " lambda= " << lambda << std::endl << " kma= " << k_ma << std::endl << " kmd= " << k_md << std::endl << " kmhop= " << k_mhop << std::endl << " p_x= " << std::endl << p_x << " dtwmax= " << dtwmax << std::endl << " tm0= " << tm0 << std::endl << " tmf= " << tmf << std::endl << " dtm= " << dtm << std::endl << " ts0= " << ts0 << std::endl << " tsf= " << tsf << std::endl << " dts= " << dts << std::endl << " T= " << T << std::endl << " ssT= " << ssTime << std::endl << "motorTime= " << motorTime << std::endl << "Lmax= " << Lmax << std::endl;

		parametersOut.close();

		//Counting elapsed time of running
		//auto start = std::chrono::steady_clock::now();

		std::cout << "Run# " << run << std::endl;

		//std::string dirName2(dirName + "/Run_" + std::to_string(run));
		//std::filesystem::create_directories(dirName2);

		if(ComputeLengthDynamic7(run, k_a, k_d, k_cat, p_res, k_hyd, k_hydx, k_rep, k_ma, k_md, k_mhop, p_x, t, T, dtwmax, dirName, tm0, tmf, dtm, ts0, tsf, dts, load, savetime, ssTime, motorTime, Lmax) == 1)
		{
			std::cout << "Error: exit programm." << std::endl;
			return 1;
		}

		//Real elapsed time during running
		//auto end = std::chrono::steady_clock::now();
		//double elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds> (end - start).count();
		//std::cout << "Elapsed time: " << elapsed_seconds << " seconds" << std::endl;

		//std::cout << " #run: " << run << " ka= " << k_a << " kd= " << k_d << " kcat= " << k_cat << " pres= " << p_res << " khyd= " << k_hyd << " khydx= " << k_hydx << " k_rep= " << k_rep << " lambda= " << lambda << " kma= " << k_ma << " kmd= " << k_md << " kmhop= " << k_mhop << " p_x= " << p_x << " dtw= " << dtw << " dtw2= " << dtw2 << std::endl;
	//}
	//}
	//}

	return 0;
}
