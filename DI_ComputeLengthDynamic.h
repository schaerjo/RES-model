#ifndef DEF_COMPUTE_LENGTH_DYNAMIC
#define DEF_COMPUTE_LENGTH_DYNAMIC

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <chrono>
#include <iomanip>
#include <string>
#include <random>
#include <cmath>
//#include "pcg_random.hpp"
#include <gsl/gsl_rng.h>
#include <vector>
#include "Filament4.h"

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
 * The simulation starts at time t and ends at time T
 * */
 int DI_ComputeLengthDynamic(int run, double k_a, double k_d, double k_cat, double k_res, long double t, double T, double dtwmax, std::string const dirName, int load, double savetime, double ssTime);

#endif
