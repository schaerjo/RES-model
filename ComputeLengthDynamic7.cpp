#include "ComputeLengthDynamic7.h"

/*
 * C++ function that compute the length dynamic of a microtubule linear model in which the subunits can be in three different state (GTP, GTPX, GDP) with molecular motors
 * by using an adapted version of the Gillespie algorithm.
 *
 * The MT switches between a phase in which it grows by adding subunits at rate k_a at the plus end and another phase in which it shrinks by detaching subunits at rate k_d at the
 * plus end.
 * The catastroph (growing to shrinking) events happen stochastically at rate k_cat.
 * A rescue can happen with probability p_res if the subunit at the tip is an exchanged one (GTPX)
 *
 * Subunits are added in the T state. They hydrolize (T to D) at rate k_hyd. D-Subunits can be exchanged spontaneously with T-Subunits at rate k_rep. These subunits are defined to be
 * in the state TX. TX-Subunits hydrolize at a rate k_hydx.
 *
 * Molecular motors can attach anywhere on the lattice at rate k_ma. They walk toward the plus end by hopping to the nearest neighbor at rate k_mhop if it's not already occupied
 * by another motor. When they hop, motors can induce a subunit exchange at the site they are hopping from with a probability p_x. Motors detach from the lattice at rate k_md.
 *
 * The length at definite times given by dtw, the growing/shrinking times/lengths and the number of Catastroph/rescue/total shrinkages/ spontaneous exchanges/ motor induced
 * exchanges are saved in files.
 *
 * The average number of motors and the average motor distribution in the steady-state are computed.
 *
 * The simulation starts at time t and ends at time T
 * */

 typedef std::vector<int> ivec;
 typedef std::vector<ivec> imat;
 typedef std::vector<double> dvec;
 typedef std::vector<dvec> dmat;


int ComputeLengthDynamic7(int run, double k_a, double k_d, double k_cat, double p_res, double k_hyd, double k_hydx, double k_rep, double k_ma, double k_md, double k_mhop, double p_x, long double t, double T, double dtwmax, std::string const dirName, int tm0, int tmf, int dtm, int ts0, int tsf, int dts, int load, double savetime, double ssTime, double motorTime, int Lmax)
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
	int nHydro = 0;
	int nHydroX = 0;
	int nSpontaneousRepair = 0;
	int nMotorRepair = 0;
	int nMotorAttachment = 0;
	int nMotorDetachment = 0;
	int nMotorHops = 0;
	int nMotorFallingOffTheTip = 0;
	int nMotorHoppingOffTheTip = 0;
	int nTotalShrinkage = 0;
	int nTotalEvent = 0;
	double timeAtNucleation = 0;

	//Rates of events depending on number of subunits in a given state or number of occupied/unoccupied sites by a motor
	double w1 = 0.0;
	double w2 = 0.0;
	double w3 = 0.0;
	double w4 = 0.0;
	double w5 = 0.0;
	double w6 = 0.0;

	//Some variables
	double unif_rn = 0.0; //RN drawn from a uniform distribution
	int num = 0; //number of a subunit or number of a motor
	int pos = 0; //Lattice position of a subunit or of a motor
	int draw = 0; //increment variable used to find position of a given motor
	int motorsHoppingReturn = 0; //integer return by the motorsHopping method: -1:Error // 0:motor at the tip or blocked // 1:hopping // 2:Hopping+Repair // 3 Hopping off the tip // 4 Hopping off the tip + repair
	int subunitDetachmentReturn = 0; //integer returned by the motorsDetachment method: -1:Error // 0:No Motor at the tip // 1: Motor at the tip
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
	int j=0;
	int i=0;
	int tm = tm0-dtm;
	int ts = ts0-dts;
	int Nss = 0; //Number of save state. If 0, all files are overwritten
	double dtw = 1;

  //Vectors for collected data that will be written to file.
  dvec lifetimeList;
  imat MTLvst, Nstatesvst, polyLengthList;
  dmat polyTimeList, Nmvst;

  //Define filenames where to write results.
	std::string const numbers(dirName + "/numbers.dat");
	std::string const subunits(dirName + "/subunits.dat");
	std::string const motors(dirName + "/motors.dat");

	std::string const events(dirName + "/events.dat"); //Store the number of events
	std::string const exeState(dirName + "/exestate.dat");// time/length at change/nucleation + Nss
	std::string const length(dirName + "/lengthVStime.dat");
	std::string const nMotors(dirName + "/nMotorsVStime.dat");
	std::string const nStates(dirName + "/nStatesVStime.dat");
	std::string const polyLength(dirName + "/polyLength.dat");
	//std::string const polyTime(dirName + "/polyTime.dat");
	std::string const Lifetime(dirName + "/lifetime.dat");

  //output stream for writing resuls in files.
  std::ofstream eventsOut;
	std::ofstream exeStateOut;
	std::ofstream lengthOut;
	std::ofstream nMotorsOut;
	std::ofstream nStatesOut;
	std::ofstream polyLengthOut;
	//std::ofstream polyTimeOut;
	std::ofstream LifetimeOut;

  //Variables and vectors used for MT overloaded constructor
	int L, globalState, Ntp, Ntr, Nd, Nm;
	std::vector<int> Sub(0);
	std::vector<int> Mot(0);

  //Reading results files in case of loading a savestate
	if(load == 1)
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
			ss >> L >> globalState >> Ntp >> Ntr >> Nd >> Nm;
		}
		infile.close();

		infile.open(subunits.c_str());
		if(infile.fail())
		{
			std::cout << "Error: can't open the file: " << subunits << std::endl;
			return 1;
		}
		while(std::getline(infile, line))
		{
			int substate;
			std::stringstream ss(line);
			ss >> substate;
			Sub.push_back(substate);
		}
		infile.close();

		infile.open(motors.c_str());
		if(infile.fail())
		{
			std::cout << "Error: can't open the file: " << motors << std::endl;
			return 1;
		}
		while(std::getline(infile, line))
		{
			int motor;
			std::stringstream ss(line);
			ss >> motor;
			Mot.push_back(motor);
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
			ss >> nCatastroph >> nRescue >> nTotalShrinkage >> nSpontaneousRepair >> nMotorRepair >> nMotorHoppingOffTheTip >> nMotorFallingOffTheTip >> nMotorAttachment >> nMotorDetachment >> nMotorHops;
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
			ss >> Nss >> timeAtChange >> lengthAtChange >> timeAtNucleation >> maxL >> t >> tw >> tm >> ts;
		}
		infile.close();
	}
	else
	{
		//No savestate loaded
		L = 1;
		globalState = 0;
		Ntp = 1;
		Ntr = 0;
		Nd = 0;
		Nm = 0;
		Sub.push_back(0);
	}

	//Overloaded Constructor
	Filament MT(L, globalState, Ntp, Ntr, Nd, Sub, Mot, Nm);
	Sub.clear();
	Mot.clear();

	//dataOut << "#avgNcat \t avgNres \t avgNtotalshrinkage \t avgNspontaneousrepair \t avgNmotorrepair \t avgPolylength \t avgPolytime \t avgDepolylength \t avgDepolytime \t CorrL \t CorrT \t avgMTlifetime" << std::endl;

	//Writing initial state.
	MTLvst.push_back({(int)tw, MT.getLength()});
	Nmvst.push_back({tw, (double)MT.getNm()/(double) MT.getLength()});
	Nstatesvst.push_back({(int)tw, MT.getNtp(), MT.getNd(), MT.getNtr()});

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
	  	{
	  		eventsOut.open(events, std::ofstream::app);
	  		exeStateOut.open(exeState, std::ofstream::app);

	  		lengthOut.open(length, std::ofstream::app);
	  		nMotorsOut.open(nMotors, std::ofstream::app);
	  		nStatesOut.open(nStates, std::ofstream::app);
	  		polyLengthOut.open(polyLength, std::ofstream::app);
	  		LifetimeOut.open(Lifetime, std::ofstream::app);

	  	}
	  	else if (load ==2)
	  	{
	  		eventsOut.open(events);
	  		eventsOut << "#nCat" << "\t" << "nRes" << "\t" << "nTotShrink" << "\t" << "nSpontEx" << "\t" << "nMotEx" << "\t" << "nMotHopOffTip" << "\t" << "nMotFallOffTip" << "\t" << "nMotAttachment" << "\t" << "nMotDetachment" << "\t" << "nMotorHops" << std::endl;
	  		exeStateOut.open(exeState, std::ofstream::app);

	  		lengthOut.open(length, std::ofstream::app);
	  		nMotorsOut.open(nMotors, std::ofstream::app);
	  		nStatesOut.open(nStates, std::ofstream::app);
	  		polyLengthOut.open(polyLength);
	  		polyLengthOut << "#growing length" << "\t" << "MT length at the catastrophe" << "\t" << "shrinking length" << "\t" << "#growing time" << "\t" << "shrinking time" << std::endl;
	  		LifetimeOut.open(Lifetime);
	  		LifetimeOut << "#MT lifetime" << std::endl;
	  	}
	  	else
	  	{
	  		eventsOut.open(events);
	  		exeStateOut.open(exeState);

	  		eventsOut << "#nCat" << "\t" << "nRes" << "\t" << "nTotShrink" << "\t" << "nSpontEx" << "\t" << "nMotEx" << "\t" << "nMotHopOffTip" << "\t" << "nMotFallOffTip" << "\t" << "nMotAttachment" << "\t" << "nMotDetachment" << "\t" << "nMotorHops" << std::endl;
	  		exeStateOut << "#Nss" << "\t" << "timeAtChange" << "\t" << "lengthAtChange" << "\t" << "timeAtNucleation" << "\t" << "maxL" << "\t" << "t" << "\t" << "tw" << "\t" << "tm" << "\t" << "ts" << std::endl;

	  		lengthOut.open(length);
	  		lengthOut << "#time (s)" << "\t" << "Length (ntub)" << std::endl;
	  		nMotorsOut.open(nMotors);
	  		nMotorsOut << "#time (s)" << "\t" << "nMotor per polymerized tubulin" << std::endl;
	  		nStatesOut.open(nStates);
	  		nStatesOut << "#time (s)" << "\t" << "nGTP" << "\t" << "nGDP" << "\t" << "nGTPx"  << std::endl;
	  		polyLengthOut.open(polyLength);
	  		polyLengthOut << "#growing length" << "\t" << "MT length at the catastrophe" << "\t" << "shrinking length" << "\t" << "#growing time" << "\t" << "shrinking time" << std::endl;
	  		//polyTimeOut.open(polyTime);
	  		//polyTimeOut << "#growing time" << "\t" << "Simulation time at the catastrophe" << "\t" << "shrinking time" << std::endl;
	  		LifetimeOut.open(Lifetime);
	  		LifetimeOut << "#MT lifetime" << std::endl;
	  	}
			//Writing the results in files.
			for(i=0; i<(int)MTLvst.size(); i++)
	  	{
	    	lengthOut << MTLvst[i][0] << "\t" << MTLvst[i][1] << std::endl;
	    	nMotorsOut << Nmvst[i][0] << "\t" << Nmvst[i][1] << std::endl;
	    	nStatesOut << Nstatesvst[i][0] << "\t" << Nstatesvst[i][1] << "\t" << Nstatesvst[i][2] << "\t" << Nstatesvst[i][3]  << std::endl;
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
			MT.savestate(numbers, subunits, motors);
			exeStateOut << std::fixed << std::setprecision(12) << Nss << "\t" << timeAtChange << "\t" << lengthAtChange << "\t" << timeAtNucleation << "\t" << maxL << "\t" << t << "\t" << tw << "\t" << tm << "\t" << ts << std::endl;
			eventsOut << std::fixed << std::setprecision(12) << nCatastroph << "\t" << nRescue << "\t" << nTotalShrinkage << "\t" << nSpontaneousRepair << "\t" << nMotorRepair << "\t" << nMotorHoppingOffTheTip << "\t" << nMotorFallingOffTheTip << "\t" << nMotorAttachment << "\t" << nMotorDetachment << "\t" << nMotorHops << std::endl;

			std::cout << "#################################" << std::endl;
			std::cout << "Run: " << run << std::endl;
			std::cout << "Savestate at time: " << t << std::endl;
			std::cout << "#################################" << std::endl;

			//Closing results files
			lengthOut.close();
			nMotorsOut.close();
			nStatesOut.close();
			exeStateOut.close();
			eventsOut.close();
			polyLengthOut.close();
			//polyTimeOut.close();
			LifetimeOut.close();
    	gsl_rng_free (rng);
			return 0;
		}

		draw = 1;
		pos = 0;
		num = 1;

		//Checks
		if(MT.getLength() != MT.getNtp() + MT.getNtr() + MT.getNd())
		{
			std::cout << "Error: The total number of subunit in any state doesn't agree with the MT length." << std::endl;
			return 1;
		}

		nTotalEvent += 1;

		if(MT.getGlobalState() == 0)//if growing
		{
			w1 = MT.getNtp()*k_hyd;
			w2 = MT.getNtr()*k_hydx;
			w3 = MT.getNd()*k_rep;
			if(t < motorTime)
			{
				w4 = 0;
				w5 = 0;
				w6 = 0;
			}
			else
			{
				w4 = (MT.getLength() - MT.getNm())*k_ma;
				w5 = MT.getNm()*k_md;
				w6 = MT.getNm()*k_mhop;
			}
			k_tot = k_a + (MT.getLength()>0)*k_cat + w1 + w2 + w3 + w4 + w5 + w6;

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
					//std::cout << "tw= " << tw << " and dtw= " << dtw << " and 10*dtw = " << 10*dtw << std::endl;
					//std::cout << "(int)tw = " << (int)tw << std::endl;
          MTLvst.push_back({(int)tw, MT.getLength()});
          Nmvst.push_back({tw, (double)MT.getNm()/(double) MT.getLength()});
          Nstatesvst.push_back({(int)tw, MT.getNtp(), MT.getNd(), MT.getNtr()});
				}
			}
			else
			{
				dtw = 10*dtw;
				if(t >= tw + dtw)
				{
					tw += dtw;
					//std::cout << "tw= " << tw << " and dtw= " << dtw << " and 10*dtw = " << 10*dtw << std::endl;
					//std::cout << "(int)tw = " << (int)tw << std::endl;
          MTLvst.push_back({(int)tw, MT.getLength()});
          Nmvst.push_back({tw, (double)MT.getNm()/(double) MT.getLength()});
          Nstatesvst.push_back({(int)tw, MT.getNtp(), MT.getNd(), MT.getNtr()});
				}
			}

			//Counting the motors and computing the motor and the subunits states distributions
			if((t >= double(tm) + double(dtm)) && tm <= tmf)
			{
				tm += dtm;

				//writing motor dist in a file
				std::string mdist(dirName + "/mdist_t_" + std::to_string(tm) + ".dat");
				std::ofstream mdistOut;
				mdistOut.open(mdist.c_str());
				for(j=0; j<MT.getNm(); j++)
				{
					mdistOut << std::fixed << MT.getMotorPosition(j) + 1 << "\t" << MT.getLength() - MT.getMotorPosition(j) - 1 << std::endl;
				}
				mdistOut.close();
			}
			if((t >= double(ts) + double(dts)) && ts <= tsf)
			{
				ts += dts;

				//Writing the subunits states distribution
				std::string statedist(dirName + "/statedist_t_" + std::to_string(ts) + ".dat");
				std::ofstream statedistOut;
				statedistOut.open(statedist.c_str());
				for(j=0; j<MT.getLength(); j++)
				{
					statedistOut << std::fixed << j+1 << "\t" << MT.getLength() - j - 1 << "\t" << MT.getSingleSubunitState(j) << std::endl;
				}
				statedistOut.close();
			}

			unif_rn = gsl_rng_uniform_pos(rng)*k_tot;//Determine which reaction occurs next

			if (unif_rn < k_a)//Subunit attachment
			{
				MT.plusEndAttachment();
				nSubunitAttachment += 1;
				if(MT.getLength() > maxL)
				{
					maxL = MT.getLength();
				}
				if(MT.getLength() == Lmax)
				{
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
			else{
				unif_rn -= k_a;
				if(unif_rn < k_cat)//Catastroph event
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
				else{
					unif_rn -= k_cat;
					if(unif_rn < w1)//Hydrolization of GTP+ tubulin
					{
						//std::uniform_int_distribution<> uniform_int_dist(1,MT.getNtp());
						//num = uniform_int_dist(rng);
            num = gsl_rng_uniform_int(rng, MT.getNtp()); //Sample integer in [0,n-1]
						//std::cout << "Test: num= " << num << " MT length= " << MT.getLength() << std::endl;
						pos = MT.getSubunit(0, num+1); //finding the position of changed subunits
						//std::cout << "Test: pos= " << pos << std::endl;
						if(pos == -1 || pos > MT.getLength())
						{
							std::cout << "Error during hyd in the getSubunit method: No valid position returned // Phase: Growing -> exit" << std::endl;
							return 1;
						}
						if(MT.changeSubunit(pos, 2) != 0)
						{
							std::cout << "Error in the changeSubunit method // Phase: Growing -> exit" << std::endl;
							return 1;
						}
						nHydro += 1;
					}
					else{
						unif_rn -= w1;
						if(unif_rn < w2)//Hydrolization of GTPX tubulin
						{
							//std::uniform_int_distribution<> uniform_int_dist(1,MT.getNtr());
							//num = uniform_int_dist(rng);
              num = gsl_rng_uniform_int(rng, MT.getNtr()); //Sample integer in [0,n-1]
							pos = MT.getSubunit(1, num+1); //finding the position of changed subunits
							if(pos == -1 || pos > MT.getLength())
							{
								std::cout << "Error during hydx in the getSubunit method: No valid position returned // Phase: Growing -> exit" << std::endl;
								return 1;
							}
							if(MT.changeSubunit(pos, 2) != 0)
							{
								std::cout << "Error in the changeSubunit method // Phase: Growing -> exit" << std::endl;
								return 1;
							}
							nHydroX += 1;
						}
						else{
							unif_rn -= w2;
							if(unif_rn < w3)//Exchange of GDP tubulin
							{
								//std::uniform_int_distribution<> uniform_int_dist(1,MT.getNd());
								//num = uniform_int_dist(rng);
                num = gsl_rng_uniform_int(rng, MT.getNd()); //Sample integer in [0,n-1]
								pos = MT.getSubunit(2, num+1); //finding the position of changed subunits
								if(MT.getOccupancySite(pos) == 0)
								{
									if(pos == -1 || pos > MT.getLength())
									{
										std::cout << "Error during exchange in the getSubunit method: No valid position returned // Phase: Growing -> exit" << std::endl;
										return 1;
									}
									if(MT.changeSubunit(pos, 1) != 0)
									{
										std::cout << "Error in the changeSubunit method // Phase: Growing -> exit" << std::endl;
										return 1;
									}
									if(t > ssTime)
									{
										nSpontaneousRepair += 1;
									}
								}
							}
							else{
								unif_rn -= w3;
								if(unif_rn <  w4)//Motor attachment
								{
									//std::uniform_int_distribution<> uniform_int_dist(1,MT.getLength()-MT.getNm());
									//num = uniform_int_dist(rng);//find randomly an empty site
                  num = gsl_rng_uniform_int(rng, MT.getLength()-MT.getNm()); //Sample integer in [0,n-1]
									while(pos < MT.getLength())
									{
										if(MT.getOccupancySite(pos) == 0)
										{
											if(draw == num+1)
											{
												break;
											}
											else
											{
												draw += 1;
											}
										}
										pos += 1;
									}
									if(draw != num+1)
									{
										std::cout << "Error in finding an empty site for motor attachment // Phase: Growing -> exit" << std::endl;
										return 1;
									}
									if(MT.motorsAttachment(pos) != 0)
									{
										std::cout << "Error in the motorAttachment method // Phase: Growing -> exit" << std::endl;
										return 1;
									}
									if(t > ssTime)
									{
										nMotorAttachment += 1;
									}
								}
								else{
									unif_rn -= w4;
									if(unif_rn < w5)//Motor detachment
									{
										//std::uniform_int_distribution<> uniform_int_dist(0,MT.getNm()-1);
										//num = uniform_int_dist(rng);
                    num = gsl_rng_uniform_int(rng, MT.getNm()); //Sample integer in [0,n-1]
										if(MT.motorsDetachment(num) != 0)
										{
											std::cout << "Error in the motorDetachment method // Phase: Growing -> exit" << std::endl;
											return 1;
										}
										if(t > ssTime)
										{
											nMotorDetachment += 1;
										}
									}
									else//Motor hopping
									{
										//cout << "motor hopping // Growing phase" << endl;
										//std::uniform_int_distribution<> uniform_int_dist(0,MT.getNm()-1);
										//num = uniform_int_dist(rng);
                    num = gsl_rng_uniform_int(rng, MT.getNm()); //Sample integer in [0,n-1]
										motorsHoppingReturn = MT.motorsHopping(num, gsl_rng_uniform_pos(rng), p_x);
										if(motorsHoppingReturn == -1)
										{
											std::cout << "Error in the motorHopping method // Phase: Growing -> exit" << std::endl;
											return 1;
										}
										else if(motorsHoppingReturn == 1)
										{
											if(t > ssTime)
											{
												nMotorHops += 1;
											}
										}
										else if(motorsHoppingReturn == 2)
										{
											if(t > ssTime)
											{
												nMotorHops += 1;
												nMotorRepair += 1;
											}
										}
										else if(motorsHoppingReturn == 3)
										{
											if(t > ssTime)
											{
												nMotorHops += 1;
												nMotorHoppingOffTheTip += 1;
											}
										}
										else if(motorsHoppingReturn == 4)
										{
											if(t > ssTime)
											{
												nMotorHops += 1;
												nMotorHoppingOffTheTip += 1;
												nMotorRepair += 1;
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		else //if shrinking
		{
			//std::cout << "Shrinking" << std::endl;
			w1 = MT.getNtp()*k_hyd;
			w2 = MT.getNtr()*k_hydx;
			w3 = MT.getNd()*k_rep;
			if(t < motorTime)
			{
				w4 = 0;
				w5 = 0;
				w6 = 0;
			}
			else
			{
				w4 = (MT.getLength() - MT.getNm())*k_ma;
				w5 = MT.getNm()*k_md;
				w6 = MT.getNm()*k_mhop;
			}

			k_tot = k_d + w1 + w2 + w3 + w4 + w5 + w6;

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
          Nmvst.push_back({tw, (double)MT.getNm()/(double) MT.getLength()});
          Nstatesvst.push_back({(int)tw, MT.getNtp(), MT.getNd(), MT.getNtr()});
				}
			}
			else
			{
				dtw = 10*dtw;
				if(t >= tw +dtw)
				{
					tw += dtw;
          MTLvst.push_back({(int)tw, -MT.getLength()});
          Nmvst.push_back({tw, (double)MT.getNm()/(double) MT.getLength()});
          Nstatesvst.push_back({(int)tw, MT.getNtp(), MT.getNd(), MT.getNtr()});
				}
			}

			//Counting the motors and computing the motor distribution
			if((t >= double(tm) + double(dtm)) && tm <= tmf)
			{
				tm += dtm;
				//writing motor dist in a file
				std::string mdist(dirName + "/mdist_t_" + std::to_string(tm) + ".dat");
				std::ofstream mdistOut;
				mdistOut.open(mdist.c_str());
				for(j=0; j<MT.getNm(); j++)
				{
					mdistOut << std::fixed << MT.getMotorPosition(j) + 1 << "\t" << MT.getLength() - MT.getMotorPosition(j) - 1 << std::endl;
				}
				mdistOut.close();
			}
			if((t >= double(ts) + double(dts)) && ts <= tsf)
			{
				ts += dts;

				//Writing the subunits states distribution
				std::string statedist(dirName + "/statedist_t_" + std::to_string(ts) + ".dat");
				std::ofstream statedistOut;
				statedistOut.open(statedist.c_str());
				for(j=0; j<MT.getLength(); j++)
				{
					statedistOut << std::fixed << j+1 << "\t" << MT.getLength() - j - 1 << "\t" << MT.getSingleSubunitState(j) << std::endl;
				}
				statedistOut.close();
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
				else if(subunitDetachmentReturn == 1)
				{
					if(t > ssTime)
					{
						nMotorFallingOffTheTip += 1;
					}
				}
				nSubunitDetachment += 1;
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
						/*polyLengthOut << growingLength << "\t" << lengthAtChange << "\t" << shrinkingLength << std::endl;
						polyTimeOut << growingTime << "\t" << timeAtChange << "\t" << shrinkingTime << std::endl;
						LifetimeOut << t - timeAtNucleation << std::endl;*/
					}
					timeAtChange = t;
					lengthAtChange = 1;
					timeAtNucleation = t;
					MT.rescue();
					//std::cout << "Complete depolymerization at time " << t << std::endl;
					MT.plusEndAttachment();
				}
				else //Possible rescue ?
				{
					if(MT.getSingleSubunitState(MT.getLength()-1) == 1)
					{
						if(MT.getSingleSubunitState(MT.getLength()-1) == -1)
						{
							std::cout << "Error: invalid state returned by getSingleSubunitState // Phase: Shrinking -> exit" << std::endl;
							return 1;
						}
						unif_rn = gsl_rng_uniform_pos(rng);
						if(unif_rn < p_res)
						{
							MT.rescue();
							shrinkingTime = t - timeAtChange;
							shrinkingLength = -MT.getLength() + lengthAtChange;
							if(t > ssTime)
							{
								nRescue += 1;
                polyLengthList.push_back({growingLength, lengthAtChange, shrinkingLength});
                polyTimeList.push_back({growingTime, shrinkingTime});
								/*polyLengthOut << growingLength << "\t" << lengthAtChange << "\t" << shrinkingLength << std::endl;
								polyTimeOut << growingTime << "\t" << timeAtChange << "\t" << shrinkingTime << std::endl;
								LifetimeOut << t - timeAtNucleation << std::endl;*/
							}
							timeAtChange = t;
							lengthAtChange = MT.getLength();
						}
					}
				}
			}
			else{
				unif_rn -= k_d;
				if(unif_rn <  w1)//Hydrolization of GTP+ tubulin
				{
					//std::uniform_int_distribution<> uniform_int_dist(1,MT.getNtp());
					//num = uniform_int_dist(rng);
          num = gsl_rng_uniform_int(rng, MT.getNtp()); //Sample integer in [0,n-1]
					pos = MT.getSubunit(0, num+1);//finding the position of changed subunits
					if(pos == -1 || pos > MT.getLength())
					{
						std::cout << "Error during hyd in the getSubunit method: No valid position returned // Phase: Shrinking -> exit" << std::endl;
						return 1;
					}
					if(MT.changeSubunit(pos, 2) != 0)
					{
						std::cout << "Error in the changeSubunit method // Phase: Shrinking -> exit" << std::endl;
						return 1;
					}
					nHydro += 1;
				}
				else{
					unif_rn -= w1;
					if(unif_rn < w2)//Hydrolization of GTPX tubulin
					{
						//std::uniform_int_distribution<> uniform_int_dist(1,MT.getNtr());
						//num = uniform_int_dist(rng);
            num = gsl_rng_uniform_int(rng, MT.getNtr()); //Sample integer in [0,n-1]
						pos = MT.getSubunit(1, num+1);//finding the position of changed subunits
						if(pos == -1 || pos > MT.getLength())
						{
							std::cout << "Error during hydx in the getSubunit method: No valid position returned // Phase: Shrinking -> exit" << std::endl;
							return 1;
						}
						if(MT.changeSubunit(pos, 2) != 0)
						{
							std::cout << "Error in the changeSubunit method // Phase: Shrinking -> exit" << std::endl;
							return 1;
						}
						nHydroX += 1;
					}
					else{
						unif_rn -= w2;
						if(unif_rn < w3)//Repair of GDP tubulin
						{
							//std::uniform_int_distribution<> uniform_int_dist(1,MT.getNd());
							//num = uniform_int_dist(rng);
              num = gsl_rng_uniform_int(rng, MT.getNd()); //Sample integer in [0,n-1]
							pos = MT.getSubunit(2, num+1);//finding the position of changed subunits
							if(MT.getOccupancySite(pos) == 0)
							{
								if(pos == -1 || pos > MT.getLength())
								{
									std::cout << "Error during exchange in the getSubunit method: No valid position returned // Phase: Shrinking -> exit" << std::endl;
									return 1;
								}
								if(MT.changeSubunit(pos, 1) != 0)
								{
									std::cout << "Error in the changeSubunit method // Phase: Shrinking -> exit" << std::endl;
									return 1;
								}
								if(t > ssTime)
								{
									nSpontaneousRepair += 1;
								}
							}
						}
						else{
							unif_rn -= w3;
							if(unif_rn < w4)//Motor Attachment
							{
								//std::uniform_int_distribution<> uniform_int_dist(1,MT.getLength()-MT.getNm());
								//num = uniform_int_dist(rng);//find randomly an empty site
                num = gsl_rng_uniform_int(rng, MT.getLength()-MT.getNm()); //Sample integer in [0,n-1]
								while(pos < MT.getLength())
								{
									if(MT.getOccupancySite(pos) == 0)
									{
										if(draw == num+1)
										{
											break;
										}
										else
										{
											draw += 1;
										}
									}
									pos += 1;
								}
								if(draw != num+1)
								{
									std::cout << "Error in finding an empty site for motor attachment // Phase:Shrinking -> exit" << std::endl;
									return 1;
								}
								if(MT.motorsAttachment(pos) != 0)
								{
									std::cout << "Error in the motorAttachment method // Phase: Shrinking -> exit" << std::endl;
									return 1;
								}
								if(t > ssTime)
								{
									nMotorAttachment += 1;
								}
							}
							else{
								unif_rn -= w4;
								if(unif_rn < w5)//motor detachment
								{
									//std::uniform_int_distribution<> uniform_int_dist(0,MT.getNm()-1);
									//num = uniform_int_dist(rng);
                  num = gsl_rng_uniform_int(rng, MT.getNm()); //Sample integer in [0,n-1]
									if(MT.motorsDetachment(num) != 0)
									{
										std::cout << "Error in the motorDetachment method // Phase: Shrinking -> exit" << std::endl;
										return 1;
									}
									if(t > ssTime)
									{
										nMotorDetachment += 1;
									}
								}
								else//motor hopping
								{
									//std::uniform_int_distribution<> uniform_int_dist(0,MT.getNm()-1);
									//num = uniform_int_dist(rng);
                  num = gsl_rng_uniform_int(rng, MT.getNm()); //Sample integer in [0,n-1]
									motorsHoppingReturn = MT.motorsHopping(num, gsl_rng_uniform_pos(rng), p_x);
									if(motorsHoppingReturn == -1)
									{
										std::cout << "Error in the motorHopping method -> exit" << std::endl;
										return 1;
									}
									else if(motorsHoppingReturn == 1)
									{
										if(t > ssTime)
										{
											nMotorHops += 1;
										}
									}
									else if(motorsHoppingReturn == 2)
									{
										if(t > ssTime)
										{
											nMotorHops += 1;
											nMotorRepair += 1;
										}
									}
									else if(motorsHoppingReturn == 3)
									{
										if(t > ssTime)
										{
											nMotorHops += 1;
											nMotorHoppingOffTheTip += 1;
										}
									}
									else if(motorsHoppingReturn == 4)
									{
										if(t > ssTime)
										{
											nMotorHops += 1;
											nMotorHoppingOffTheTip += 1;
											nMotorRepair += 1;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

  lifetimeList.push_back(-t + timeAtNucleation);
	//LifetimeOut << -t + timeAtNucleation << std::endl;

  //Opening the results files
	if(load == 1)
	{
		eventsOut.open(events, std::ofstream::app);
		exeStateOut.open(exeState, std::ofstream::app);

		lengthOut.open(length, std::ofstream::app);
		nMotorsOut.open(nMotors, std::ofstream::app);
		nStatesOut.open(nStates, std::ofstream::app);
		polyLengthOut.open(polyLength, std::ofstream::app);
		//polyTimeOut.open(polyTime, std::ofstream::app);
		//LifetimeOut.open(Lifetime, std::ofstream::app);

	}
	else if (load ==2)
	{
		eventsOut.open(events);
		eventsOut << "#nCat" << "\t" << "nRes" << "\t" << "nTotShrink" << "\t" << "nSpontEx" << "\t" << "nMotEx" << "\t" << "nMotHopOffTip" << "\t" << "nMotFallOffTip" << "\t" << "nMotAttachment" << "\t" << "nMotDetachment" << "\t" << "nMotorHops" << std::endl;
		exeStateOut.open(exeState, std::ofstream::app);

		lengthOut.open(length, std::ofstream::app);
		nMotorsOut.open(nMotors, std::ofstream::app);
		nStatesOut.open(nStates, std::ofstream::app);
		polyLengthOut.open(polyLength);
		polyLengthOut << "#growing length" << "\t" << "MT length at the catastrophe" << "\t" << "shrinking length" << "\t" << "#growing time" << "\t" << "shrinking time" << std::endl;
		//LifetimeOut.open(Lifetime);
		//LifetimeOut << "#MT lifetime" << std::endl;
	}
	else
	{
		eventsOut.open(events);
		exeStateOut.open(exeState);

		eventsOut << "#nCat" << "\t" << "nRes" << "\t" << "nTotShrink" << "\t" << "nSpontEx" << "\t" << "nMotEx" << "\t" << "nMotHopOffTip" << "\t" << "nMotFallOffTip" << "\t" << "nMotAttachment" << "\t" << "nMotDetachment" << "\t" << "nMotorHops" << std::endl;
		exeStateOut << "#Nss" << "\t" << "timeAtChange" << "\t" << "lengthAtChange" << "\t" << "timeAtNucleation" << "\t" << "maxL" << "\t" << "t" << "\t" << "tw" << "\t" << "tm" << "\t" << "ts" << std::endl;

		lengthOut.open(length);
		lengthOut << "#time (s)" << "\t" << "Length (ntub)" << std::endl;
		nMotorsOut.open(nMotors);
		nMotorsOut << "#time (s)" << "\t" << "nMotor per polymerized tubulin" << std::endl;
		nStatesOut.open(nStates);
		nStatesOut << "#time (s)" << "\t" << "nGTP" << "\t" << "nGDP" << "\t" << "nGTPx"  << std::endl;
		polyLengthOut.open(polyLength);
		polyLengthOut << "#growing length" << "\t" << "MT length at the catastrophe" << "\t" << "shrinking length" << "\t" << "#growing time" << "\t" << "shrinking time" << std::endl;
		//LifetimeOut.open(Lifetime);
		//LifetimeOut << "#MT lifetime" << std::endl;
	}

  //Writing the results in the files.
  for(i=0; i<(int)MTLvst.size(); i++)
  {
    lengthOut << MTLvst[i][0] << "\t" << MTLvst[i][1] << std::endl;
    nMotorsOut << Nmvst[i][0] << "\t" << Nmvst[i][1] << std::endl;
    nStatesOut << Nstatesvst[i][0] << "\t" << Nstatesvst[i][1] << "\t" << Nstatesvst[i][2] << "\t" << Nstatesvst[i][3]  << std::endl;
  }
  for(i=0; i<(int)polyLengthList.size(); i++)
  {
    polyLengthOut << polyLengthList[i][0] << "\t" << polyLengthList[i][1] << "\t" << polyLengthList[i][2] << "\t" << polyTimeList[i][0] << "\t" << polyTimeList[i][1] << std::endl;
  }
  /*for(i=0; i<(int)lifetimeList.size(); i++)
  {
    LifetimeOut << lifetimeList[i] << std::endl;
  }*/

  //Saving the final state.
	Nss += 1;
	MT.savestate(numbers, subunits, motors);
	exeStateOut << std::fixed << std::setprecision(12) << Nss << "\t" << timeAtChange << "\t" << lengthAtChange << "\t" << timeAtNucleation << "\t" << maxL << "\t" << t << "\t" << tw << "\t" << tm << "\t" << ts << std::endl;
	eventsOut << std::fixed << std::setprecision(12) << nCatastroph << "\t" << nRescue << "\t" << nTotalShrinkage << "\t" << nSpontaneousRepair << "\t" << nMotorRepair << "\t" << nMotorHoppingOffTheTip << "\t" << nMotorFallingOffTheTip << "\t" << nMotorAttachment << "\t" << nMotorDetachment << "\t" << nMotorHops << std::endl;

	//Closing results files
	lengthOut.close();
	nMotorsOut.close();
	nStatesOut.close();
	polyLengthOut.close();
	//polyTimeOut.close();
	//LifetimeOut.close();
	exeStateOut.close();
	eventsOut.close();

  gsl_rng_free (rng);

	return 0;
}
