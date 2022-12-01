#ifndef DEF_FILAMENT
#define DEF_FILAMENT

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <iomanip>

class Filament
{
	public:

		Filament(); //Default Constructor
		Filament(int length, bool globalState, int Ntp, int Ntr, int Nd, std::vector<int> Subunits, std::vector<int> Motors_Occupancy, int Nm); //Overloaded operator
		~Filament(); //Destructor
		void savestate(std::string numbers, std::string subunits, std::string motors);
		void plusEndAttachment();
		int plusEndDetachment();
		void catastroph();
		void rescue();
		int changeSubunit(int pos, int state);//change state of given subunit
		int motorsAttachment(int pos);
		int motorsDetachment(int num);
		int motorsHopping(int num, double rng, double repairProbability);
		int getNm() const;
		int getOccupancySite(int pos) const;//Get position of nth occupied/unoccupied site
		void printOccupiedIndices() const;//Print the positions of all motors on lattice in the terminal
		int getSubunit(int state, int num) const;//Get position of nth subunit in given state
		int getSingleSubunitState(int pos) const;//get state of a given subunits
		std::vector<int> getSubunitsStates() const;// get string of subunits states
		bool getGlobalState() const;
		int getLength() const;
		int getNtp() const;
		int getNtr() const;
		int getNd() const;
		int getMotorPosition(int j) const;
	
	private:

		int a_length;
		bool a_globalState; //0 -> growing state 1 -> shrinking state
		int a_Ntp; //Number of subunits in GTP state from plus end attachment
		int a_Ntr; //Number of subunits in GTP state from repair
		int a_Nd; //Number of subunits in GDP state
		std::vector<int> a_Subunits;//vector of int representing subunits state in filament, 0 is GTP +end, 1 is GTP repair and 2 is gdp state
		std::vector<int> a_Motors_Occupancy;//Vector of size a_Nm containig the lattice indices of sites occupied by a motor
		int a_Nm;//Number of attached motors
};

#endif
