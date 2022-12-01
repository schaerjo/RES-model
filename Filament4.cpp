#include "Filament4.h"

//Constructor
Filament::Filament() : a_length(1), a_globalState(0), a_Ntp(1), a_Ntr(0), a_Nd(0), a_Subunits(1,0), a_Motors_Occupancy(0), a_Nm(0)
{
	//std::cout << "Filament constructor: a_Subunits.size() = " << a_Subunits.size() << " and a_Subunits[0] = " << a_Subunits[0] << std::endl;
}

//Overloaded Constructor
Filament::Filament(int length, bool globalState, int Ntp, int Ntr, int Nd, std::vector<int> Subunits, std::vector<int> Motors_Occupancy, int Nm) : a_length(length), a_globalState(globalState), a_Ntp(Ntp), a_Ntr(Ntr), a_Nd(Nd), a_Subunits(Subunits), a_Motors_Occupancy(Motors_Occupancy), a_Nm(Nm)
{
	
}

//Destructor
Filament::~Filament()
{

}

//Function that save the state of the object to be reloaded in a later execution.
void Filament::savestate(std::string numbers, std::string subunits, std::string motors)
{
	std::ofstream numbersOut; //Write out: length, globalstate, Ntp, Ntr, Nd, Nm 
	std::ofstream subunitsOut;
	std::ofstream motorsOut;
	std::cout << std::setprecision(12);

	//open the state files and writing the data
	numbersOut.open(numbers.c_str());
	numbersOut << std::fixed  << a_length << "\t" << a_globalState << "\t" << a_Ntp << "\t" << a_Ntr << "\t" << a_Nd << "\t" << a_Nm << std::endl;
	numbersOut.close();

	subunitsOut.open(subunits.c_str());
	for(int j=0; j<a_length; j++)
	{
		subunitsOut << std::fixed << a_Subunits[j] << std::endl;
	}
	subunitsOut.close();

	motorsOut.open(motors.c_str());
	for(int i=0; i<a_Nm; i++)
	{
		motorsOut << std::fixed << a_Motors_Occupancy[i] << std::endl;
	}
	motorsOut.close();
}

void Filament::plusEndAttachment()
{
	//std::cout << "Inside plusEndAttachment(): before: MT length = " << a_length << " a_Subunits.size() = " << a_Subunits.size() << std::endl;
	a_length +=1;
	a_Ntp += 1;
	a_Subunits.push_back(0);
	//std::cout << "Inside plusEndAttachment(): after: MT length = " << a_length << " a_Subunits.size() = " << a_Subunits.size() << std::endl;
}

int Filament::plusEndDetachment()
{
	int pop = 0;
	int j = 0;
	if(a_length == 0)
	{
		std::cout << "Error: MT length is already 0: wrong utilization of plusEndDetachment() -> exit" << std::endl;
		return -1;
	}
	else
	{
		if(a_Subunits.back() == 0){
			a_Ntp -= 1;
		}
		else if(a_Subunits.back() == 1){
			a_Ntr -= 1;
		}
		else{
			a_Nd -= 1;
		}
		
		if(getOccupancySite(a_length-1) == 1)
		{
			while(pop == 0 && j < a_Nm)
			{
				if(a_Motors_Occupancy[j] == a_length-1)
				{
					a_Motors_Occupancy.erase(a_Motors_Occupancy.begin()+j);
					pop = 1;
				}
				j += 1;
			}
			if(pop == 0)
			{
				std::cout << "Error in plusEndDetachment: Tip Motor update has failed." << std::endl;
				return -1;
			}
			a_Nm -=1;
			if(a_Nm != int(a_Motors_Occupancy.size()))
			{
				std::cout << "Error in plusEndDetachment: number of motors and size of motor occupancy vector don't agree." << std::endl;
				return -1;
			}
		}
		a_length -=1;
		a_Subunits.pop_back();
		return pop;
	}
}

void Filament::catastroph()
{
	a_globalState = 1;
}

void Filament::rescue()
{
	a_globalState = 0;
}

//Change the state of a given subunit
int Filament::changeSubunit(int pos, int state)
{
	//std::cout << "Inside changeSubunit(): before: a_Subunits[" << pos << "] = " << a_Subunits[pos] << " and new state = " << state << std::endl;
	//std::cout << "Inside changeSubunit(): before: a_Ntp = " << a_Ntp << " a_Ntr = " << a_Ntr << " a_Nd = " << a_Nd << std::endl;
	if(pos < 0 || pos >= a_length || (state != 0 && state != 1 && state != 2))
	{
		std::cout << "Error: wrong utilization of changeSubunit(pos, state) -> exit" << std::endl;
		return 1;
	}
	else
	{
		if(a_Subunits[pos] == 0){
			a_Ntp -= 1;
		}
		else if(a_Subunits[pos] == 1){
			a_Ntr -= 1;
		}
		else{
			a_Nd -= 1;
		}

		if(state == 0){
			a_Ntp += 1;
		}
		else if(state == 1){
			a_Ntr += 1;
		}
		else{
			a_Nd += 1;
		}

		a_Subunits[pos] = state;
		//std::cout << "Inside changeSubunit(): after: a_Subunits[" << pos << "] = " << a_Subunits[pos] << " and new state = " << state << std::endl;
		//std::cout << "Inside changeSubunit(): after: a_Ntp = " << a_Ntp << " a_Ntr = " << a_Ntr << " a_Nd = " << a_Nd << std::endl;
		return 0;
	}
}

int Filament::motorsAttachment(int pos)
{
	if(getOccupancySite(pos) == 0)
	{
		a_Motors_Occupancy.push_back (pos);
		a_Nm += 1;
		
		//Check
		if(a_Nm != int(a_Motors_Occupancy.size()))
		{
			std::cout << "Error in motorAttachment: number of motors and size of motor occupancy vector don't agree." << std::endl;
			return 1;
		}
		return 0;
	}
	else
	{
		std::cout << "Error wrong utilization of motorsAttachment." << std::endl;
		return 1;
	}
}

int Filament::motorsDetachment(int num)
{
	if(num < 0 || num >= int(a_Motors_Occupancy.size()))
	{
		std::cout << "Error in motorDetachment: invalid motor number" << std::endl;
		return 1;
	}
	
	a_Motors_Occupancy.erase (a_Motors_Occupancy.begin()+num);
	a_Nm -= 1;
	
	//Check
	if(a_Nm != int(a_Motors_Occupancy.size()))
	{
		std::cout << "Error in motorDetachment: number of motors and size of motor occupancy vector don't agree." << std::endl;
		return 1;
	}
	return 0;
}

int Filament::motorsHopping(int num, double rng, double repairProbability)
{
	if(num < 0 || num >= int(a_Motors_Occupancy.size()))
	{
		std::cout << "Error in motorsHopping: invalid motor number" << std::endl;
		return -1;
	}
	else
	{
		//Test if motors is at the tip
		if(a_Motors_Occupancy[num] + 1 == a_length)
		{
			//Possible repair event
			if(a_Subunits[a_Motors_Occupancy[num]] == 2 && rng <= repairProbability)
			{
				if(changeSubunit(a_Motors_Occupancy[num],1) != 0)
				{
					std::cout << "Error in motorsHopping: changeSubunits returned an error." << std::endl;
					std::cout << "MotorOccupancy = " << a_Motors_Occupancy[num] << " MT length= " << a_length << std::endl;
					return -1;
				}
				//Hopping of the MT.
				a_Motors_Occupancy.erase(a_Motors_Occupancy.begin()+num);
				a_Nm -= 1;
				return 4;
			}
			//Hopping of the MT.
			a_Motors_Occupancy.erase(a_Motors_Occupancy.begin()+num);
			a_Nm -= 1;
			return 3;
		}
		//Test if next site is empty
		else if(getOccupancySite(a_Motors_Occupancy[num] + 1) == 0)
		{
			//Possible repair event
			if(a_Subunits[a_Motors_Occupancy[num]] == 2 && rng <= repairProbability)
			{
				if(changeSubunit(a_Motors_Occupancy[num],1) != 0)
				{
					std::cout << "Error in motorsHopping: changeSubunits returned an error." << std::endl;
					return -1;
				}
				a_Motors_Occupancy[num] += 1;
				return 2;
			}
			//Hopping
			a_Motors_Occupancy[num] += 1;
			return 1;
		}
		else
		{
			return 0;
		}
	}
}

int Filament::getNm() const
{
	//Check
	if(a_Nm != int(a_Motors_Occupancy.size()))
	{
		std::cout << "Error in getNm: number of motors and size of motor occupancy vector don't agree." << std::endl;
	}
	return int(a_Motors_Occupancy.size());
}

//Find if the nth lattice site is occupied by a motor: return 1 if occupied or 0 if unoccupied
int Filament::getOccupancySite(int pos) const
{
	int Occ = 0;
	for(int j=0; j<int(a_Motors_Occupancy.size()); j++)
	{
		if(a_Motors_Occupancy[j] == pos)
		{
			Occ = 1;
		}
	}
	return Occ;
}

int Filament::getMotorPosition(int j) const
{
	return a_Motors_Occupancy[j];
}

void Filament::printOccupiedIndices() const
{
	for(int j=0; j < int(a_Motors_Occupancy.size()); j++)
	{
		std::cout << "a_Motors_Occupancy[" << j << "] = " << a_Motors_Occupancy[j] << std::endl;
	}
}

//Find the position of the nth subunit in a given state from the begining of the string, i.e. the minus end
int Filament::getSubunit(int state, int num) const
{
	int count = 1;
	for(int i=0; i < a_length; i++){
		if(state == a_Subunits[i]){
			//std::cout << "Inside getSubunit(): a_Subunits[" << i << "] = " << a_Subunits[i] << std::endl;
			//std::cout << "Inside getSubunit(): num= " << num << " count = " << count << std::endl;
			if(count == num){
				return i;
			}
			else{
				count += 1;
			}
		}
	}
	//std::cout <<"Inside getSubunit(): count= " << count << std::endl;
	std::cout << "Error in getSubunit(state, num) method -> exit" << std::endl;
	return -1;
}

int Filament::getSingleSubunitState(int pos) const
{
	if(pos < 0 || pos >= a_length)
	{
		std::cout << "Error: invalid position given to getSingleSubunitState(pos) -> exit" << std::endl;
		return -1;
	}
	else
	{
		return a_Subunits[pos];
	}
}

std::vector<int> Filament::getSubunitsStates() const
{
	return a_Subunits;	
}

bool Filament::getGlobalState() const
{
	return a_globalState; //0 = growing; 1 = shrinking
}

int Filament::getLength() const
{
	return a_length;
}

int Filament::getNtp() const
{
	return a_Ntp;
}

int Filament::getNtr() const
{
	return a_Ntr;
}

int Filament::getNd() const
{
	return a_Nd;
}
