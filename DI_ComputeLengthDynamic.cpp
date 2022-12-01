#include "DI_ComputeLengthDynamic.h"

/*
 * C++ function that compute the length dynamic of a microtubule linear model by using an adapted version of the Gillespie algorithm.
 *
 * The MT switches between a phase in which it grows by adding subunits at rate k_a at the plus end and another phase in which it shrinks by detaching subunits at rate k_d at the plus end.
 * The catastroph (growing to shrinking) events happen stochastically at rate k_cat.
 * The rescue (shrinking to growing) events happen stochastically at rate k_res.
 * The length at definite times given by dtw, the growing/shrinking times/lengths and the number of Catastroph/rescue/total shrinkages are saved in files.
 *
 * The simulation starts at time t and ends at time T
 * */

 typedef std::vector<int> ivec;
 typedef std::vector<ivec> imat;
 typedef std::vector<double> dvec;
 typedef std::vector<dvec> dmat;


int DI_ComputeLengthDynamic(int run, double k_a, double k_d, double k_cat, double k_res, long double t, double T, double dtwmax, std::string const dirName, int load, double savetime, double ssTime)
{
	std::cout << std::setprecision(12);
	std::cout << std::fixed;

  //Random Number Generator using GSL library.
  const gsl_rng_type * randT;
  gsl_rng * rng;

  gsl_rng_env_setup();

  randT = gsl_rng_default;
  rng = gsl_rng_alloc (randT);

	//printf ("generator type: %s\n", gsl_rng_name (rng));
  //printf ("seed = %lu\n", gsl_rng_default_seed);

	//Counting elapsed time of running
	auto start = std::chrono::steady_clock::now();

	//Seed for the RNG taken from the std random device
	//pcg_extras::seed_seq_from<std::random_device> seed_source;

	//Make RNG
	//pcg32 rng(seed_source);

	//Make a copy of RNG state to use later
	//pcg32 rng_checkpoint = rng;

  //Create a real uniform distribution
	//std::uniform_real_distribution<long double> uniform_dist(0.0, 1.0);

	/*std::random_device rd; //Will be used to obtain a seed for the random number engine
	std::mt19937 rng(rd()); //Standard mersenne_twister_engine seeded with rd()*/

	//Number of events
	int nSubunitAttachment = 0;
	int nSubunitDetachment = 0;
	int nCatastroph = 0;
	int nRescue = 0;
	int nTotalShrinkage = 0;
	int nTotalEvent = 0;
	double timeAtNucleation = 0;

	//Some variables
	double unif_rn = 0.0; //RN drawn from a uniform distribution
	long double dt = 0.0; // time step
	double k_tot = 0.0; //total events rate
	double tw = (double)t; //writing time
	double timeAtChange = t;//time at change of global phase.
	int lengthAtChange = 1;//MT length at change of global phase
	double growingTime = 0.0;
	double shrinkingTime = 0.0;
	int growingLength = 0;
	int shrinkingLength = 0;
	double avgdt = 0.0;
	int loop = 0;
	int maxL = 1;
	int i=0;
	int Nss = 0; //Number of save state. If 0, all files are overwritten
	double dtw = 1;
	int subunitDetachmentReturn = 0; //integer returned by the motorsDetachment method: -1:Error // 0:No Motor at the tip // 1: Motor at the tip

  //Vectors for collected data that will be written to file.
  dvec lifetimeList;
  imat MTLvst, polyLengthList;
  dmat polyTimeList;

  //Define filenames where to write results.
	std::string const numbers(dirName + "/numbers.dat");
	std::string const events(dirName + "/events.dat"); //Store the number of events
	std::string const exeState(dirName + "/exestate.dat");// time/length at change/nucleation + Nss
	std::string const length(dirName + "/lengthVStime.dat");
	std::string const polyLength(dirName + "/polyLength.dat");
	std::string const Lifetime(dirName + "/lifetime.dat");

  //output stream for writing resuls in files.
	std::ofstream numbersOut;
  std::ofstream eventsOut;
	std::ofstream exeStateOut;
	std::ofstream lengthOut;
	std::ofstream polyLengthOut;
	std::ofstream LifetimeOut;

  //Variables and vectors used for MT overloaded constructor
	int L, globalState;
	std::vector<int> Sub(0);
	std::vector<int> Mot(0);

  //Reading results files in case of loading a savestate
	if(load == 1 || load == 2)
	{
		//Load a savestate
		std::ifstream infile;
		infile.open(numbers.c_str());
		if(infile.fail())
		{
			std::cout << "Error: can't open the file: " << numbers << std::endl;
			return 1;
		}
		std::string line;
		while(std::getline(infile, line))
		{
			std::stringstream ss(line);
			ss >> L >> globalState;
		}
		infile.close();

		infile.open(events.c_str());
		if(infile.fail())
		{
			std::cout << "Error: can't open the file: " << events << std::endl;
			return 1;
		}
		while(std::getline(infile, line))
		{
			std::stringstream ss(line);
			ss >> nCatastroph >> nRescue >> nTotalShrinkage >> nSubunitAttachment >> nSubunitDetachment;
		}
		infile.close();

		infile.open(exeState.c_str());
		if(infile.fail())
		{
			std::cout << "Error: can't open the file: " << exeState << std::endl;
			return 1;
		}
		while(std::getline(infile, line))
		{
			std::stringstream ss(line);
			ss >> Nss >> timeAtChange >> lengthAtChange >> timeAtNucleation >> maxL >> t >> tw;
		}
		infile.close();
	}
	else
	{
		//No savestate loaded
		L = 1;
		globalState = 0;
	}

	for(i=0; i<L; i++)
	{
		Sub.push_back(0);
	}

	//Overloaded Constructor
	Filament MT(L, globalState, L, 0, 0, Sub, Mot, 0);
	Sub.clear();
	Mot.clear();

	//Writing initial state.
	MTLvst.push_back({(int)tw, MT.getLength()});

	//Main Loop
	while (t < T)
	{
		auto end = std::chrono::steady_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds> (end - start).count();
		//std::cout << "Elapsed time in seconds: " << elapsed_seconds << std::endl;
		if(elapsed_seconds >= savetime && savetime != -1)
		{
	    //Opening the results files
	  	if(load == 1)
	  	{	numbersOut.open(numbers, std::ofstream::app);
	  		eventsOut.open(events, std::ofstream::app);
	  		exeStateOut.open(exeState, std::ofstream::app);

	  		lengthOut.open(length, std::ofstream::app);
	  		polyLengthOut.open(polyLength, std::ofstream::app);
	  		LifetimeOut.open(Lifetime, std::ofstream::app);

	  	}
	  	else if (load ==2)
	  	{
				numbersOut.open(numbers);
				numbersOut << "#L" << "\t" << "Global State" << std::endl;
	  		eventsOut.open(events);
	  		eventsOut << "#nCat" << "\t" << "nRes" << "\t" << "nTotShrink" << "\t" << "nSubunitAttachment" << "\t" << "nSubunitDetachment" << std::endl;
	  		exeStateOut.open(exeState, std::ofstream::app);

	  		lengthOut.open(length, std::ofstream::app);
	  		polyLengthOut.open(polyLength);
	  		polyLengthOut << "#growing length" << "\t" << "MT length at the catastrophe" << "\t" << "shrinking length" << "\t" << "#growing time" << "\t" << "shrinking time" << std::endl;
	  		LifetimeOut.open(Lifetime);
	  		LifetimeOut << "#MT lifetime" << std::endl;
	  	}
	  	else
	  	{
				numbersOut.open(numbers);
	  		eventsOut.open(events);
	  		exeStateOut.open(exeState);

	  		eventsOut << "#nCat" << "\t" << "nRes" << "\t" << "nTotShrink" << "\t" << "nSubunitAttachment" << "\t" << "nSubunitDetachment" << std::endl;
	  		exeStateOut << "#Nss" << "\t" << "timeAtChange" << "\t" << "lengthAtChange" << "\t" << "timeAtNucleation" << "\t" << "maxL" << "\t" << "t" << "\t" << "tw" << std::endl;
				numbersOut << "#L" << "\t" << "Global State" << std::endl;

	  		lengthOut.open(length);
	  		lengthOut << "#time (s)" << "\t" << "Length (ntub)" << std::endl;
	  		polyLengthOut.open(polyLength);
	  		polyLengthOut << "#growing length" << "\t" << "MT length at the catastrophe" << "\t" << "shrinking length" << "\t" << "#growing time" << "\t" << "shrinking time" << std::endl;
	  		LifetimeOut.open(Lifetime);
	  		LifetimeOut << "#MT lifetime" << std::endl;
	  	}
			//Writing the results in files.
			for(i=0; i<(int)MTLvst.size(); i++)
	  	{
	    	lengthOut << MTLvst[i][0] << "\t" << MTLvst[i][1] << std::endl;
	  	}
	  	for(i=0; i<(int)polyLengthList.size(); i++)
	  	{
	    	polyLengthOut << polyLengthList[i][0] << "\t" << polyLengthList[i][1] << "\t" << polyLengthList[i][2] << "\t" << polyTimeList[i][0] << "\t" << polyTimeList[i][1] << std::endl;
	  	}
	  	for(i=0; i<(int)lifetimeList.size(); i++)
	  	{
	    	LifetimeOut << lifetimeList[i] << std::endl;
	  	}

			Nss += 1;
			numbersOut << std::fixed << std::setprecision(12) << MT.getLength() << "\t" << MT.getGlobalState() << std::endl;
			exeStateOut << std::fixed << std::setprecision(12) << Nss << "\t" << timeAtChange << "\t" << lengthAtChange << "\t" << timeAtNucleation << "\t" << maxL << "\t" << t << "\t" << tw << std::endl;
			eventsOut << std::fixed << std::setprecision(12) << nCatastroph << "\t" << nRescue << "\t" << nTotalShrinkage << "\t" << nSubunitAttachment << "\t" << nSubunitDetachment << std::endl;

			std::cout << "#################################" << std::endl;
			std::cout << "Run: " << run << std::endl;
			std::cout << "Savestate at time: " << t << std::endl;
			std::cout << "#################################" << std::endl;

			//Closing results files
			lengthOut.close();
			exeStateOut.close();
			eventsOut.close();
			polyLengthOut.close();
			LifetimeOut.close();
			numbersOut.close();
    	gsl_rng_free (rng);
			return 0;
		}

		//Checks
		if(MT.getLength() != MT.getNtp() + MT.getNtr() + MT.getNd())
		{
			std::cout << "Error: The total number of subunit in any state doesn't agree with the MT length." << std::endl;
			return 1;
		}

		nTotalEvent += 1;

		if(MT.getGlobalState() == 0)//if growing
		{
			k_tot = k_a + (MT.getLength()>0)*k_cat;

			//Random draw of time until next event
			unif_rn = gsl_rng_uniform_pos(rng);
			dt = -log(unif_rn)/k_tot; //when next reaction occurs
			t += dt;
			avgdt += dt;
			loop += 1;

			//Writing data at definite intervals
			if(tw < 10*dtw || dtw == dtwmax)
			{
				if(t >= tw + dtw)
				{
					tw += dtw;
          MTLvst.push_back({(int)tw, MT.getLength()});
				}
			}
			else
			{
				dtw = 10*dtw;
				if(t >= tw + dtw)
				{
					tw += dtw;
          MTLvst.push_back({(int)tw, MT.getLength()});
				}
			}

			unif_rn = gsl_rng_uniform_pos(rng)*k_tot;//Determine which reaction occurs next

			if (unif_rn < k_a)//Subunit attachment
			{
				MT.plusEndAttachment();
				if(MT.getLength() > maxL)
				{
					maxL = MT.getLength();
				}
				if(t > ssTime)
				{
					nSubunitAttachment += 1;
				}
			}
			else
			{
				//std::cout << "Catastroph: MT length is: " << MT.getLength() << std::endl;
				MT.catastroph();
				growingTime = t - timeAtChange;
				growingLength = MT.getLength() - lengthAtChange;
				timeAtChange = t;
				lengthAtChange = MT.getLength();
				if(t > ssTime)
				{
					nCatastroph += 1;
				}
			}
		}
		else //if shrinking
		{
			//std::cout << "Shrinking" << std::endl;
			k_tot = k_d + k_res;

			//Random draw of time until next event
			//unif_rn = uniform_dist(rng);
      unif_rn = gsl_rng_uniform_pos(rng);
			dt = -log(unif_rn)/k_tot; //when next reaction occurs
			t += dt;

			avgdt += dt;
			loop += 1;

			//Writing data at definite intervals
			if(tw < 10*dtw || dtw == dtwmax)
			{
				if(t >= tw +dtw)
				{
					tw += dtw;
          MTLvst.push_back({(int)tw, -MT.getLength()});
				}
			}
			else
			{
				dtw = 10*dtw;
				if(t >= tw +dtw)
				{
					tw += dtw;
          MTLvst.push_back({(int)tw, -MT.getLength()});
				}
			}

			unif_rn = gsl_rng_uniform_pos(rng)*k_tot;//Determine which reaction occurs next

			if (unif_rn < k_d)//Subunit Detachment
			{
				//std::cout << "Detachment: MT length is: " << MT.getLength() << std::endl;
				subunitDetachmentReturn = MT.plusEndDetachment();
				if(subunitDetachmentReturn == -1)
				{
					std::cout << "Error in plusEndDetachment method // Phase: Shrinking -> exit" << std::endl;
					return 1;
				}

				if(t > ssTime)
				{
					nSubunitDetachment += 1;
				}

				if(MT.getLength() == 0)
				{
					shrinkingTime = t - timeAtChange;
					shrinkingLength = -MT.getLength() + lengthAtChange;
					if(t > ssTime)
					{
						nTotalShrinkage += 1;
            polyLengthList.push_back({growingLength, lengthAtChange, shrinkingLength});
            polyTimeList.push_back({growingTime, shrinkingTime});
            lifetimeList.push_back(t - timeAtNucleation);
					}
					timeAtChange = t;
					lengthAtChange = 1;
					timeAtNucleation = t;
					MT.rescue();
					MT.plusEndAttachment();
				}
			}
			else //Rescue event
			{
				shrinkingTime = t - timeAtChange;
				shrinkingLength = -MT.getLength() + lengthAtChange;
				if(t > ssTime)
				{
					nRescue += 1;
					polyLengthList.push_back({growingLength, lengthAtChange, shrinkingLength});
					polyTimeList.push_back({growingTime, shrinkingTime});
				}
				timeAtChange = t;
				lengthAtChange = MT.getLength();
				MT.rescue();
			}
		}
	}

  lifetimeList.push_back(-t + timeAtNucleation);

  //Opening the results files
	if(load == 1)
	{
		numbersOut.open(numbers, std::ofstream::app);
		eventsOut.open(events, std::ofstream::app);
		exeStateOut.open(exeState, std::ofstream::app);

		lengthOut.open(length, std::ofstream::app);
		polyLengthOut.open(polyLength, std::ofstream::app);
		LifetimeOut.open(Lifetime, std::ofstream::app);

	}
	else if (load ==2)
	{
		numbersOut.open(numbers);
		numbersOut << "#L" << "\t" << "Global State" << std::endl;
		eventsOut.open(events);
		eventsOut << "#nCat" << "\t" << "nRes" << "\t" << "nTotShrink" << "\t" << "nSubunitAttachment" << "\t" << "nSubunitDetachment" << std::endl;
		exeStateOut.open(exeState, std::ofstream::app);

		lengthOut.open(length, std::ofstream::app);
		polyLengthOut.open(polyLength);
		polyLengthOut << "#growing length" << "\t" << "MT length at the catastrophe" << "\t" << "shrinking length" << "\t" << "#growing time" << "\t" << "shrinking time" << std::endl;
		LifetimeOut.open(Lifetime);
		LifetimeOut << "#MT lifetime" << std::endl;
	}
	else
	{
		numbersOut.open(numbers);
		eventsOut.open(events);
		exeStateOut.open(exeState);

		eventsOut << "#nCat" << "\t" << "nRes" << "\t" << "nTotShrink" << "\t" << "nSubunitAttachment" << "\t" << "nSubunitDetachment" << std::endl;
		exeStateOut << "#Nss" << "\t" << "timeAtChange" << "\t" << "lengthAtChange" << "\t" << "timeAtNucleation" << "\t" << "maxL" << "\t" << "t" << "\t" << "tw" << std::endl;
		numbersOut << "#L" << "\t" << "Global State" << std::endl;

		lengthOut.open(length);
		lengthOut << "#time (s)" << "\t" << "Length (ntub)" << std::endl;
		polyLengthOut.open(polyLength);
		polyLengthOut << "#growing length" << "\t" << "MT length at the catastrophe" << "\t" << "shrinking length" << "\t" << "#growing time" << "\t" << "shrinking time" << std::endl;
		LifetimeOut.open(Lifetime);
		LifetimeOut << "#MT lifetime" << std::endl;
	}

  //Writing the results in the files.
  for(i=0; i<(int)MTLvst.size(); i++)
  {
    lengthOut << MTLvst[i][0] << "\t" << MTLvst[i][1] << std::endl;
  }
  for(i=0; i<(int)polyLengthList.size(); i++)
  {
    polyLengthOut << polyLengthList[i][0] << "\t" << polyLengthList[i][1] << "\t" << polyLengthList[i][2] << "\t" << polyTimeList[i][0] << "\t" << polyTimeList[i][1] << std::endl;
  }
  for(i=0; i<(int)lifetimeList.size(); i++)
  {
    LifetimeOut << lifetimeList[i] << std::endl;
  }

  //Saving the final state.
	Nss += 1;
	numbersOut << std::fixed << std::setprecision(12) << MT.getLength() << "\t" << MT.getGlobalState() << std::endl;
	exeStateOut << std::fixed << std::setprecision(12) << Nss << "\t" << timeAtChange << "\t" << lengthAtChange << "\t" << timeAtNucleation << "\t" << maxL << "\t" << t << "\t" << tw << std::endl;
	eventsOut << std::fixed << std::setprecision(12) << nCatastroph << "\t" << nRescue << "\t" << nTotalShrinkage << "\t" << nSubunitAttachment << "\t" << nSubunitDetachment << std::endl;

	//Closing results files
	lengthOut.close();
	polyLengthOut.close();
	LifetimeOut.close();
	exeStateOut.close();
	eventsOut.close();
	numbersOut.close();

  gsl_rng_free (rng);

	return 0;
}
