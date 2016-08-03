//@================================================================//
//@ Comment notes:                                                 //
//@   // or /*     - Keith's comments.                             //
//@   //@ or /*@   - Syrian's comments.                            //
//@   //& or /*&   - Syrian's comments + questions.                //
//@   //## or /*## - Syrian's comments + flagged for removal/mod.  //
//@================================================================//

#include <iostream>
#include <fstream>
#include <cmath>

#include <MT/mtrand.h>

//@----------------------------------------------------------------

//@ is.h information (variables)

//@----------------------------------------------------------------

int nx = 8;

//@ size in x-direction (since 1d, it's only direction)

int ny = 1;

//## size in y-direction (remove later if easily translate into just 1d considerations)

int nz = 1;

//## size in z-direction (remove later if easily translate into just 1d considerations)

int n3 = nx*ny*nz;

//## total number of spins (remove later if easily translate into just 1d considerations)

int ll = 10000;

//@ max value for expansion truncation L (hamiltonian string cut off max)

int lls = 500;

/*& maximum length of subsequence, what does this mean */

int l;

//@ expansion truncation, or hamiltonian string cut off max

double ht;

// transverse field strength

double beta;

// inverse temperature

double eshift;

/*& sum of all constants added to the hamiltonian (see hamiltonian equations perhaps?) */

int nh;

//@ current expansion order of n, or number of non-identity hamiltonians in hamiltonian string

int mlls;

//@ maximum substring length that has been reached

double ar1, ar2, ar3, ar4;

//@ acceptance rates

//@----------------------------------------------------------------

//@----------------------------------------------------------------

//@ other global (variables)

//@----------------------------------------------------------------

int isteps;

//@ equilibriation steps

int msteps;

//@ measurement steps 

int nd;

// number of runs

//@----------------------------------------------------------------

//@----------------------------------------------------------------

//@ initializing arrays from is.h

//@----------------------------------------------------------------

int spn[n3] = {0};

//@ spin array

int stra[ll] = {0};
int strb[ll] = {0};

/*& (NOTE reason for names, see str1 below) operator string(s)? uses ll which is the max expansion truncation */

int xyz1[n3] = {0};
int xyz2[n3] = {0};
int xyz3[n3] = {0};

/*@ corrdinates for given spin i, where x-coordinate is xyz1, y-coordinate is xyz2, z-coordinate is xyz3. note: Sandvik (I believe) makes a typo when he typed that y-coordinate is equal to the same array as x-coordinate. */

int xyzi[n3] = {0};

/*## this one gets hairy. Sandvik has it originally as xyzi(nx,ny,nz) which is a 3 dimensional array in fortran77 representing spin number given by coordinates x,y,z = xyzi(x,y,z). so although I am not exactly following him at this point, it's likely best to switch to a 1D chosen from here on out, hence why n3 still works at the moment. Anyway, flagged for modification. */

int disx[nx*nx] = {0};
int disy[ny*ny] = {0};
int disz[nz*nz] = {0};

/*## this one also could get hairy. although this seems like the best route, just watch for the math. represents the periodic distances in x,y, and z directions between two points (remove later if easily translate into just 1d considerations). possibly still needs modifications as well since original is disx(nx,nx) and etc. */

int pfrst[n3+1] = {0};

/*& first position to search in probability list given INT(R*n3), what is R (random number generator?) and is this similar to frstspinop array? */

int plast[n3+1] = {0};

/*& last position to search in probability list given INT(R*n3), what is R (random number generator?) and is this similar to lastspinop array? */

bool lis[3][3][nx/2 + 1][ny/2 + 1][nz/2 + 1] = {0};

/*& //## flag for allowed spin configuration given distance between sites. what even is this. ok so original is basically the equivalent of a boolean 5 dimensional array it seems: LOGICAL lis(-1:1,-1:1,0:nx/2,0:ny/2,0:nz/2) */

int ris[nx/2 + 1][ny/2 + 1][nz/2 + 1] = {0};

/*& //## Sandvik has this as essentially floating numbers in an array, is it necessary (find out in rest of the code), can I even attempt that in c++? represents ising interaction with constants added given distance */

int pint[n3 + 1] = {0};

/*& list of integrated probabilities lists, Sandvik uses floating point */

int ap1[ll + 1] = {0};

/*& probability for increasing n by 1, given current n. expansion order n (number of non identity hamiltonians in string). Sandvik uses floating point */

int dp1[ll + 1] = {0};

/*& probability for decreasing n by 1, given current n. expansion order n (number of non identity hamiltonians in string). */

int lsub[n3] = {0};

/*& length of subsequences?*/

int pos1[lls*n3] = {0};

/*& positions?, watch the math */

int str1[lls*n3] = {0};

//& substring operators?, watch the math (NOTE above name stra and strb arrays)

int con1[(lls+1)*n3] = {0};

//& constraints? watch the math

//@----------------------------------------------------------------

//@----------------------------------------------------------------

//@ initializing functions directly in main

//@----------------------------------------------------------------

void mcstep();

//@ monte carlo step

void checkl(int step);

/*& check the hamiltonian string length perhaps? has pvectors, openlog, and closelog */

void init();

//@ initializes the system: initconf, lattice, pvectors, isingpart, zerodat

void openlog();

/*& looks like writes information to a (12) log.txt (note a 10?), note the closelog function so probably writing the equilibrium step thing to log.txt */

void closelog();

//@ close the (12) log.txt

void errcheck();

/*& looks like checks if certain spins weren't possible in the first place? uses str(1,i) & str(2,i) arrays */

void writeconf();

/*& writes length of lattice in x,y,z directions (we'll just use x), beta, L (hamiltonian string length?), spn(i) spin array values, str(1,i) (bond array indicates site left for our case?), str(2,i) (bond array indicates site right for our case?) uses (10) (20)*/ 

void measure();

/*& measuring observables looks like. note: uses a common which seems like fortran77's version of a sort of pass by reference across different fortran .f files*/

void results();

/*@ writes calculated observables into (10) mag.dat for su and xu, (10) enr.dat for e, c, avni, avnt, and avnu */

void writeacc(int msteps);

/*& writes diagonal acceptance 1, diagonal acceptance 2, off-diagonal substitutions, spin flips / site, max lsub (max hamiltonian string length?) */

//@----------------------------------------------------------------

//@----------------------------------------------------------------

//@ initializing functions not directly in main

//@----------------------------------------------------------------

void initconf();

/*@ generates initial random spin state and sets initial values for l (Hamiltonian string length) and nh (number of non identity Hamiltonians) */

void lattice();

/*& creates lattice for coordinates for the lattice sites (1D in this case though) and vectors for distances that take into account the periodic boundary conditions (gives shorter of two possible distances, although this is probably simpler because we're just 1d) */

void pvectors();

/*& calculate the acceptance probabilities for the n-changing update */

void isingpart();

/*& not entirely sure what this does */

void zerodat();

//@ sets measurement data variables to 0.0

//@----------------------------------------------------------------

//@----------------------------------------------------------------

//@ main function

//@----------------------------------------------------------------

int main() {

	unsigned long int time_seed = time(NULL);
	MTRand_open rand_num(time_seed);

	//@ random number generator (0,1). Note: omitted RAN(), INITRANDOM, and WRITERAND(rndf)

	init();

	//@ consider Sandvik's (in == 0)

	for (int i=0; i<isteps; i++) {

		mcstep();

		checkl(i);

		if (((i+1)%(isteps/10))==0) { //@ to record every 10 steps
	
			openlog();

			cout << "Done equilibriation step: " << i << "\n\n";	
		
			closelog();

			errcheck();

		}

	}
		
	writeconf();

	for (int j=0; j<nd; j++) { /*@ so in Sandvik's code the 10 label in the loop is basically another way to have "enddo" in the form of "10 continue" (probably has utility for goto 10 or something)*/
	
		double ar1 = 0.0;
		double ar2 = 0.0;
		double ar3 = 0.0;
		double ar4 = 0.0;
		for (int k=0; k<msteps; k++) {

			mcstep();

			if (((k+1)%2)==0) {

				measure();

			}

			if (((k+1)%(msteps/10))==0) { //@ to record every 10 steps

				openlog();

				cout << "Done measurement step: " << j << ", " << k << "\n\n";

				closelog();

				errcheck();

			}

		}

		results();

		writeacc(msteps);

		writeconf();

	}


	return 0;


}

//@----------------------------------------------------------------

//@----------------------------------------------------------------

//@ Helper functions

//@----------------------------------------------------------------

void init() {

	//@ omitted readconf function
	initconf();
	lattice();
	pvectors();
	isingpart();
	zerodat();

}

//@----------------------------------------------------------------

void initconf() {

	
