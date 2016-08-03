//==============================================//
// Comment notes:                               //
//   // or /*     - Keith's comments.           //
//   //@ or /*@   - Syrian's comments.          //
//   /*@@         - Syrian's comment questions. //
//==============================================//

#include <iostream>
#include <fstream>
#include <cmath>

#include <MT/mtrand.h>


int isteps; /*@ Sandvik writes equilibriation steps (consider Sandvik's in == 0 for our code)*/

//@ equilibriation steps

int msteps; //@ Sandvik writes msteps

//@ measurement steps 

int nd; //@ Sandvik writes nd

// number of runs

double beta; //@ Sandvik uses beta to write T = beta then properly defines beta = 1 / beta

// inverse temperature

double ht; //@ Sandvik writes ht

// transverse field strength

void mcstep();

void checkl(int step);

void init();

//@ initializes the system: initconf, lattice, pvectors, isingpart, zerodat

//@Part of init()-----------------------------------------------------

void initconf();

/*@ generates initial random spin state and sets initial values for l (Hamiltonian string length) and nh (number of non identity Hamiltonians) */

void lattice();

/*@@ creates lattice for coordinates for the lattice sites (1D in this case though) and vectors for distances that take into account the periodic boundary conditions (gives shorter of two possible distances, although this is probably simpler because we're just 1d) */

void pvectors();

/*@@ calculate the acceptance probabilities for the n-changing update */

void isingpart();

/*@@ not entirely sure what this does */

void zerodat();

//@ sets measurement data variables to 0.0

//@End of init---------------------------------------------------------

void openlog();

/*@@ looks like writes information to a (12) log.txt (note a 10?), note the closelog function so probably writing the equilibrium step thing to log.txt */

void closelog();

//@ close the (12) log.txt

void errcheck();

/*@@ looks like checks if certain spins weren't possible in the first place? uses str(1,i) & str(2,i) arrays */

void writeconf();

/*@@ writes length of lattice in x,y,z directions (we'll just use x), beta, L (hamiltonian string length?), spn(i) spin array values, str(1,i) (bond array indicates site left for our case?), str(2,i) (bond array indicates site right for our case?) uses (10) (20)*/ 

//@====================================================================

int main() {

unsigned long int time_seed = time(NULL);
MTRand_open rand_num(time_seed);

//@ random number generator (0,1). Note: omitted RAN(), INITRANDOM, and WRITERAND(rndf)

init();

//@ consider Sandvik's (in == 0)

for (int i=0; i<isteps; i++) {

	mcstep();

  checkl(i);

	if (((i+1)%(isteps/10))==0) {
	
		openlog();

		cout << "Done equilibriation step: " << i << "\n\n";	
		
		closelog();

		errcheck();

	}

}
  
writeconf();

for (int j=0; j<nd; j++) {
	







  return 0;


}
