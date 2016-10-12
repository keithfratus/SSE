//@==================================================================//
//@ Comment notes:                                                   //
//@   // or /*        - Keith's comments.                            //
//@   crtl+f: @       - Syrian's comments.                           //
//@   crtl+f: $       - Syrian's comments + questions.               //
//@   crtl+f: ~       - Syrian's comments + flagged for removal/mod. //
//@   crtl+f: RELOOK  - Relook at before compiling                   //
//@==================================================================//

//@==================================================================//
//@ Modification notes:
//@		1. Started changing n3's to nx's, xyz1 and xyzi to x1, and 
//@			omitting ny, nz, xyz2, xyz3, disy, disz
//@		2. Removed random number generator parts and using MTRand.
//@		3. openlog(), closelog() functions.
//@   4. readconf() not used because doesn't seem like we'd want to?
//@   5. removed quite a bit
//@   6. reorganized, temp loop, reset to 0 start of each temp loop,
//@      put rng to reset start of each temp loop, <M^2>, <M^4>, <E>,
//@			 p=10, isteps, msteps, and EnVsisteps (calcenrcool) 
//@		7. LF rng line fix: random number from 0 to nx 
//@			instead of 0 to (nx-1)
//@==================================================================//

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <string>

#include "MT/mtrand.h"

using namespace std;



//@----------------------------------------------------------------
//@ Declare (and initialize some) variables
//@----------------------------------------------------------------



const int nx = 16;

//@ size in x-direction (since 1d, it's only direction)

const double ht = 1.5;

// transverse field strength

const double start_temp = 0.50;

const int num_temp_step = 5;

const int ll = 496;

//@ max value for expansion truncation L (hamiltonian string cut off max)



const double power = 1.5;

//@ power of ising interactions decay 

const long int msteps = 64*nx*12500;

//@ measurement steps

//@const double avn1ctemp = 1.0;

//@ for first "extra" idea

const int lls = 500;

//$ maximum length of subsequence

const double temp_step = 0.01;

const double final_temp = start_temp + temp_step*(num_temp_step - 1);

const long int isteps = 2000;

//@ equilibriation steps

int l;

//@ expansion truncation, or hamiltonian string cut off max

int nh;

//@ current expansion order of n, or number of non-identity hamiltonians in hamiltonian string

int mlls;

//@ maximum substring length that has been reached

int nmsr;

// number of measurements counted (for getting averages)

double temp;

double beta;

// inverse temperature

long double ar1;

long double ar2;

long double ar3;

long double ar4;

//@ acceptance rates

long double eshft;

/*$ sum of all constants added to the hamiltonian (see hamiltonian equations perhaps?) */

long double su;

//$

long double xu;

//$

long double avnu;

//$

long double avni;

//$

long double avnt;

//$

long double avn1; 

//@ average energy perhaps

long double avn2;

//$


long double av_magpow2;

long double av_magpow4;

//@ <M^2> and <M^4>

long double q;

//@ average # of [i,i] operators + # of [i,0] operators

long double q_2;

//@ average of the square of # of [i,i] operators + # of [i,0] operators

long double mx;

//@ average transverse field

long double dmx_dh;

//@ average transverse field susceptibility



//@----------------------------------------------------------------
//@ Declare/initialize arrays
//@----------------------------------------------------------------



int spn[nx+1] = {0};

//@ spin array

int stra[ll+1] = {0};

int strb[ll+1] = {0};

/*~ (NOTE reason for names, see str1 below) operator string(s)? uses ll which is the max expansion truncation. go down p list decides which operator it is. i and j can be much farther apart, but there are still two things to specify hence the 1,2 part of the multidimensions first array gives i second array gives j put together to get Hij. */

int x1[nx+1] = {0}; //@ changed xyz1 to x1

/*~ corrdinates for given spin i*/

int disx[nx+1][nx+1] = {0};

/*~ represents the periodic distances in x directions between two points*/

int pfrst[nx+1] = {0};

/*$ first position to search in probability list given INT(R*n3), what is R (random number generator?) and is this similar to frstspinop array? */

int plast[nx+1] = {0};

/*$ last position to search in probability list given INT(R*n3), what is R (random number generator?) and is this similar to lastspinop array? */

bool lis[3][3][nx/2 + 1] = {0};

/*$ //~ For allowed spin configuration given distance between sites original was LOGICAL lis(-1:1,-1:1,0:nx/2,0:ny/2,0:nz/2)*/

double ris[nx/2 + 1] = {0};

/*$ //~ ising interaction with constants added given distance */

double pint[nx + 1] = {0};

/*$ //~ list of integrated probabilities lists */

double ap1[ll + 1] = {0};

/*$ probability for increasing n by 1, given current n. expansion order n (number of non identity hamiltonians in string). */

double dp1[ll + 1] = {0};

/*$ probability for decreasing n by 1, given current n. expansion order n (number of non identity hamiltonians in string). */

int lsub[nx+1] = {0};

/*$ //~ length of subsequences?*/

int pos1[lls+1][nx+1] = {0};

/*$ //~ positions?*/

int str1[lls+1][nx+1] = {0};

//$ //~ substring operators? (NOTE above name stra and strb arrays)

bool con1[lls+1][nx+1] = {0};

//$ //~ constraints? watch the math

int pos[lls+1] = {0};

//$

int opr[lls+1] = {0};

//$

int tmp1[ll+1] = {0};

//$

int tmp2[ll+1] = {0};

//$

bool lstr[ll+1] = {0};

//$

int spn1[nx+1] = {0};

//$



//@----------------------------------------------------------------
//@ Declare functions within main
//@----------------------------------------------------------------



void init(MTRand_open& rand_num);

//@ initializes the system: initconf, lattice, pvectors, isingpart, zerodat

void mcstep(MTRand_open& rand_num);

//@ monte carlo step

void checkl(long int step, ofstream& myfile1, MTRand_open& rand_num);

/*$ check the hamiltonian string length perhaps? has pvectors, openlog, and closelog */

void errcheck();

//@ void writeconf(ofstream& myfile2);

/*$ writes length of lattice in x, beta, L (hamiltonian string length?), spn(i) spin array values, str(1,i) (bond array indicates site left for our case?), str(2,i) (bond array indicates site right for our case?)*/ 

void measure();

/*$ measuring observables looks like. note: uses a common which seems like fortran77's version of a sort of pass by reference across different fortran .f files*/

void results(ofstream& myfile3, ofstream& myfile4, ofstream& myfile6, ofstream& myfile7, ofstream& myfile8, ofstream& myfile9, ofstream& myfile10, ofstream& myfile11, ofstream& myfile12, ofstream& myfile13);

/*@ writes calculated observables into (10) mag.dat for su and xu, (10) enr.dat for e, c, avni, avnt, and avnu */

void writeacc(ofstream& myfile5);

/*$ writes diagonal acceptance 1, diagonal acceptance 2, off-diagonal substitutions, spin flips / site, max lsub (max hamiltonian string length?) */

//void calcenrcool (int step, ofstream& myfile8);

//@



//@----------------------------------------------------------------
//@ Declare functions within functions
//@----------------------------------------------------------------



void initconf(MTRand_open& rand_num);

/*@ generates initial random spin state and sets initial values for l (Hamiltonian string length) and nh (number of non identity Hamiltonians) */

void lattice();

/*$ creates lattice for coordinates for the lattice sites (1D in this case though) and vectors for distances that take into account the periodic boundary conditions (gives shorter of two possible distances, although this is probably simpler because we're just 1d) */

void pvectors();

/*$ calculate the acceptance probabilities for the n-changing update */

void isingpart();

/*$ not entirely sure what this does */

void zerodat();

//@ sets measurement data variables to 0.0

void diaupdate(MTRand_open& rand_num);

//@ [-1,i] = flips i; [0,0] = 1; [i,i] = h; [i,j] = ising operator

void partition();

//$ ?

void offupdate(int s, MTRand_open& rand_num);

//$


void writeopen(ofstream& myfile1, ofstream& myfile3, ofstream& myfile4, ofstream& myfile5, ofstream& myfile6, ofstream& myfile7, ofstream& myfile8, ofstream& myfile9, ofstream& myfile10, ofstream& myfile11, ofstream& myfile12, ofstream& myfile13);

//NOTE: REMOVED MYFILE2

void zero_between_temps();

/*@ Important: setting variables that need to be 0 (or 0.0 if a double) to that between temperatures!!
		setting arrays that need to be 0 (or 0.0 if a double) to that between temperatures!!*/



//@----------------------------------------------------------------
//@ Main function
//@----------------------------------------------------------------



int main() {

	ofstream writelog;

	//@ ofstream writeconfig;

	ofstream writemag;

	ofstream writeenr;

	ofstream writeac;

	ofstream writec_temp;

	ofstream writemag_2vstemp;

	ofstream writeevstemp;

	ofstream writebindercumulantvstemp;

	ofstream writeq;

	ofstream writeq_2;

	ofstream writemx;

	ofstream writedmx_dh;

	writeopen(writelog, writemag, writeenr, writeac, writec_temp, writemag_2vstemp, writeevstemp, writebindercumulantvstemp, writeq, writeq_2, writemx, writedmx_dh);

	//@ Initialize the write-data-to-file.
	//@ Name and begin the write-data-to-file.
	
	for (int y=0; y<num_temp_step; y++) {

  	temp = start_temp + y*temp_step;
		
		//@ start temp loop

    beta = 1.0/temp;

    cout << "Now processing T = " << temp << ", beta = " << beta << "\n\n";

    writelog << "Now processing T = " << temp << ", beta = " << beta << "\n\n";

   	writec_temp << temp << " ";

		writemag_2vstemp << temp << " ";

		writebindercumulantvstemp << temp << " ";

		writeevstemp << temp << " ";

   	writeq << temp << " ";

		writeq_2 << temp << " ";

		writemx << temp << " ";

		writedmx_dh << temp << " ";

		unsigned long int time_seed = time(NULL);

		MTRand_open ran(time_seed);

		//@ random number generator (0,1). Note: omitted RAN(), INITRANDOM, and WRITERAND(rndf)

		zero_between_temps();

		/*@ Important: setting variables that need to be 0 (or 0.0 if a double) to that between temperatures!!
		setting arrays that need to be 0 (or 0.0 if a double) to that between temperatures!!*/

		init(ran);

		for (long int i=1; i<(isteps+1); i++) {

			mcstep(ran);

			checkl(i, writelog, ran);

			//@ if (temp==avn1ctemp){

				//@ calcenrcool(i, writeavn1c_cool);

			//@ }
			
			if ((i%(isteps/10))==0) { //@ to record every 10 steps

				writelog << "Done equilibriation step: " << i << "\n\n";	

				errcheck();

			}

		}

		//@ writeconf(writeconfig);

		ar1 = 0.0;
		
		ar2 = 0.0;

		ar3 = 0.0;

		ar4 = 0.0;

		for (long int k=1; k<(msteps+1); k++) {

			mcstep(ran);

			if ((k%2)==0) {

				measure();

			}

			if ((k%(msteps/10))==0) { //@ to record every 10 steps

				writelog << "Done measurement step " << k << "\n\n";

				errcheck();

			}

		}

		results(writemag, writeenr, writec_temp, writemag_2vstemp, writeevstemp, writebindercumulantvstemp, writeq, writeq_2, writemx, writedmx_dh);

		writeacc(writeac);

		//@ writeconf(writeconfig);
	
		cout << "Completed temp step of temp = " << temp << ", beta = " << beta << "\n\n";

		writelog << "Completed temp step of temp = " << temp << ", beta = " << beta << "\n\n";

	}

	writelog.close();

	//@ writeconfig.close();

	writemag.close();

	writeenr.close();

	writeac.close();

	writec_temp.close();

	writemag_2vstemp.close();

	writeevstemp.close();

	writebindercumulantvstemp.close();

	writeq.close();

	writeq_2.close();

	writemx.close();

	writedmx_dh.close();

	cout << "Done.\n\n";

	return 0;


}



//@----------------------------------------------------------------
//@ Helper functions
//@----------------------------------------------------------------



void init(MTRand_open& rand_num) {

	initconf(rand_num);

	lattice();

	pvectors();

	isingpart();

	zerodat();


}



//@----------------------------------------------------------------



void initconf(MTRand_open& rand_num) {

	/*@ Generates the initial random spin state, and sets initial values
	for l (expansion truncation) and nh (current expansion order n). */

	int c;
		
	l = 20;

	nh = 0;

	for (int i=1; i<(ll+1); i++) {

		stra[i] = 0;

		strb[i] = 0;

	}
	
	for (int j=1; j<(nx+1); j++) {

		spn[j] = -1;

	}

	//@ setting all spins to -1 first, could use the fill_n thing, but I'll leave it as is at least for now

	for (int k=1; k<(nx/2 + 1); k++) { 

		c = min((int(rand_num()*nx) + 1), nx);
		
		while (true) {

			if (spn[c] == -1) {

				break;

			}

			c = min ((int(rand_num()*nx) + 1), nx);

		}

		// the result of this is to randomly select half of the spins and flip them to +1, seems like an overly complicated way to initialize the spins, but whatever

		spn[c] = 1;

	}

	
}



//@----------------------------------------------------------------



void lattice() {

		/*@Creates maps with coordinates for the lattice sites "(xyz, xyzi)", 
	and vectors disx,"disy,disz" for distances that take into account the 
	periodic boundary conditions (gives the shorter of the two possible
	distances).*/

	// many of these data structures may no longer be necessary/useful in 1D code, but possibly useful to hold on to for higher dimensional generalization	
	
	int i = 0;

	for (int ix=1; ix<(nx+1); ix++) {

		i = i + 1;		

		x1[i] = ix;

	}
	
	// object below encodes distances on a periodic 1D lattice

	for (int j=1; j<(nx+1); j++) {

		for (int k=1; k<(nx+1); k++) {

			disx[k][j] = abs(k - j);

			if ((disx[k][j])>(nx/2)) {

				disx[k][j] = nx - disx[k][j]; /*@ think 1D chain with PBC so it's like a ring which means max disx is nx/2 */

			}

		}

	}


}



//@----------------------------------------------------------------



void pvectors() {

	//@ Calculates the acceptance probabilities for the n-changing update.

	for (int i=0; i<l; i++ ) {

		ap1[i] = ht * beta * double(nx) / (double(l - i));

	} 

	for (int j=1; j<(l+1); j++) {

		dp1[j] = double(l - j + 1) / (ht * beta * double(nx));

	}


}



//@----------------------------------------------------------------



void isingpart() {

	int ix; //@note ix used for inside below loop as well
	
	double r1, r2, p1, p2;

	for (ix=0; ix<(nx/2 + 1); ix++) { /*$ what do these exactly accomplish*/

		if (ix==0) {

			ris[ix] = ht;

			lis[2][2][ix] = true;

			lis[2][0][ix] = false;

			lis[0][2][ix] = false;

			lis[0][0][ix] = true;

		} else {

			r1 = double(ix);

			r2 = double(nx - ix);

			ris[ix] = (-1.0)/pow(r1, power) - 1.0/pow(r2, power); //$~ TWEAK THIS FOR PROPER "DISTANCES OR THIS THING LOOK AT THOSE r^2's"

			if (ris[ix]<=0) {

				lis[2][2][ix] = true;

				lis[2][0][ix] = false;

				lis[0][2][ix] = false;

				lis[0][0][ix] = true;

			} else {

				lis[2][2][ix] = false;

				lis[2][0][ix] = true;

				lis[0][2][ix] = true;

				lis[0][0][ix] = false;

			}

			ris[ix] = abs(ris[ix]);

		}

	}

	pint[0] = 0.0;

	pint[1] = ris[0];

	eshft = 0.0;

	for (int i=2; i<(nx+1); i++) {

		ix = disx[1][ x1[i] ]; /*@ to be clear, this is finding the distance between the first site (site x1[0]) and every other site (site x1[i]), note: 1 to nx*/

		pint[i] = pint[i-1] + ris[ix];

		eshft = eshft + (ris[ix]/2);

	}
	
	for (int j=1; j<(nx+1); j++) {

		pint[j] = pint[j]/pint[nx];

	}

	for (int k=0; k<(nx+1); k++) {

		p1 = double(k)/double(nx);

		p2 = double(k+1)/double(nx);

		if (k==0) {

			pfrst[k] = 1;

		} else {

			for (int c=pfrst[k-1]; c<(nx+1); c++) {

				if (pint[c]>=p1) {

					pfrst[k] = c;

					break; //@ breaks and goes to below "if (k==nx)" statement

				}

			}

		}

		if (k==nx) {

			plast[k] = nx;

		} else {

			for (int d=pfrst[k]; d<(nx+1); d++) {

				if (pint[d]>=p2) {

					plast[k] = d;

					break; /*@ breaks and continues "for (int k=0; k<(nx+1); k++)" loop */
				}

			}

		}

	}


}



//@----------------------------------------------------------------



void mcstep(MTRand_open& rand_num) {

	diaupdate(rand_num);

	partition();

	for (int i=1; i<(nx+1); i++) {

		offupdate(i, rand_num);

	}


}



//@----------------------------------------------------------------



void diaupdate(MTRand_open& rand_num) {

	/*@ [-1,i] = flips i
			[0,0]  = 1
			[i,i]  = h
			[i,j]  = Ising operator */

	int s1, s2, p0, p1, p2, ix, nac1, nac2, ntr2;

	double p;

	bool failed_all_ifs = false;

	nac1 = 0;

	nac2 = 0;

	ntr2 = 0;

	for (int i=1; i<(l+1); i++) {

		s1 = stra[i];

		s2 = strb[i];

		if (s1==0) { //@ H00

			p = ap1[nh];

			if (p>=1) {

				stra[i] = min(int(rand_num()*nx)+1, nx);

				strb[i] = stra[i];

				nh = nh + 1;

				nac1 = nac1 + 1;

			} else if (rand_num()<p) {

				stra[i] = min(int(rand_num()*nx)+1, nx);

				strb[i] = stra[i];

				nh = nh + 1;

				nac1 = nac1 + 1;

			}
			
		} else if (s1==s2) { //@ Hii

			p = dp1[nh];

			if (p>=1) {

				stra[i] = 0;

				strb[i] = 0;

				nh = nh - 1;

				nac1 = nac1 + 1;

			} else if (rand_num()<p) {

				stra[i] = 0;

				strb[i] = 0;

				nh = nh - 1;

				nac1 = nac1 + 1;

			}
					
		} else if (s1==-1) { //@ H(-1)j

			spn[s2] = -spn[s2];

		}
					
		s1 = stra[i];
		/*~ Can likely omit this considering he already defined s1 to this and didn't change s1 inbetween this. */						

		if (s1>0) { //@ doesn't this work on Hii as well if i>0

			while (true) {

				if (!failed_all_ifs) {

					ntr2 = ntr2 + 1;

					p = rand_num();

					p0 = int(p*(nx+1));

					p1 = pfrst[p0];

					p2 = plast[p0];

				} else {

				failed_all_ifs = false;

				//@ for last branch with an if, else statement

				}

				if (p1==p2) {

					s2 = ((s1 + p1 - 2) % nx) + 1;

					ix = disx[ x1[s1] ][ x1[s2] ];

					if ((lis[spn[s1]+1][spn[s2]+1][ix])) {

						strb[i] = s2;

						nac2 = nac2 + 1;

						break;

					} else {

						continue;

					}

				} else if (p2==(p1+1)) {

					if (pint[p1]>=p) {

						s2 = ((s1 + p1 - 2) % nx) + 1;

					} else {

						s2 = ((s1 + p2 - 2) % nx) + 1;

					}

					ix = disx[ x1[s1] ][ x1[s2] ];		

					if (lis[ spn[s1]+1 ][ spn[s2]+1 ][ix]) {

						strb[i] = s2;

						nac2 = nac2 + 1;

						break;

					} else {

						continue;

					}

				} else {

					p0 = (p1 + p2)/2;

					if ((pint[p0-1]<p)&&(pint[p0]>=p)) {

						s2 = ((s1 + p0 - 2) % nx) + 1;

						ix = disx[ x1[s1] ][ x1[s2] ];

						if (lis[ spn[s1]+1 ][ spn[s2]+1 ][ix]) {

							strb[i] = s2;

							nac2 = nac2 + 1;

							break;

						} else {

							continue;

						}

					} else {

						if ((pint[p0])>=p) {

							p2 = p0;

						} else {

							p1 = p0;

						}

						failed_all_ifs = true;

						continue;

						//@ see "near" beginning of while(true) statement

					}

				}	

			}

		}

	}

	ar1 = ar1 + (long double)(nac1)/((long double)(l));

	if (ntr2!=0) {

		ar2 = ar2 + (long double)(nac2)/((long double)(ntr2));

	}


}



//@----------------------------------------------------------------



void partition() {

	int s1, s2;
	
	for (int i=1; i<(nx+1); i++) {

		lsub[i] =0;

		con1[0][i] = false;

	}

	for (int j=1; j<(l+1); j++) {

		s1 = stra[j];

		if (s1!=0) {

			s2 = strb[j];

			if (s1==s2) {

				lsub[s1] = lsub[s1] + 1;

				pos1[ lsub[s1] ][s1] = j;

				str1[ lsub[s1] ][s1] = 1;

				con1[ lsub[s1] ][s1] = false;

			} else if (s1==(-1)) {

				lsub[s2] = lsub[s2] + 1;

				pos1[ lsub[s2] ][s2] = j;

				str1[ lsub[s2] ][s2] = 2;

				con1[ lsub[s2] ][s2] = false;

			} else {

				con1[ lsub[s1] ][s1] = true;

				con1[ lsub[s2] ][s2] = true;

			}

		}

	}
	
	for (int k=1; k<(nx+1); k++) {

		if (con1[0][k]) {

			con1[ lsub[k] ][k] = true;

		}
			
		if (lsub[k]>mlls) {

			mlls = lsub[k];

		}

	}


}



//@----------------------------------------------------------------



void offupdate(int s, MTRand_open& rand_num) {
	
	int ii, p1, p2, lens, flp, top, nac, nupd;

	nupd = 0;

	lens = lsub[s];

	for (int i=1; i<(lens+1); i++) {

		if ((!con1[i][s])) {

			nupd = nupd + 1;

			pos[nupd] = i;

		}

		opr[i] = str1[i][s];

	}

	if (lens<2) {

		if ((nupd == lens)&&((!con1[0][s]))) { //$ ?

			if (rand_num()<0.75) {

				spn[s] = -spn[s];

				ar4 = ar4 + 1.0;

				return;

			} else {

				return;

			}

		} else {

			return;

		}

	}

	nac = 0;

	flp = 1;

	for (int j=1; j<(2*nupd); j++) {

		p1 = min(int(rand_num()*nupd)+1, nupd);

		p1 = pos[p1];

		p2 = (p1 % lens) + 1;

		if (opr[p1]==opr[p2]) {

			opr[p1] = 3 - opr[p1];

			opr[p2] = 3 - opr[p2];

			nac = nac + 1;

		} else {

			top = opr[p1];

			opr[p1] = opr[p2];

			opr[p2] = top; //@ why even have top if this is all its used for?

		}
		
		if (p2<p1) {

		 flp = -flp;

		}

	}

	spn[s] = spn[s] * flp;

	for (int k=1; k<(lens+1); k++) {

		ii = pos1[k][s];

		if (opr[k]==1) {

			stra[ii] = strb[ii];

		} else {

			stra[ii] = -1;

		}

	}

	ar3 = ar3 + (long double)(nac);

	if ((nupd == lens)&&((!con1[0][s]))) {

		if (rand_num()<0.75) {

			spn[s] = -spn[s];

			ar4 = ar4 + 1.0;

		}

	}


}



//@----------------------------------------------------------------



void checkl(long int step, ofstream& myfile1, MTRand_open& rand_num) {

	int p, dl, l1;

	//@if ((step%500)==0){

	//@	cout << "step = " << step << " has l = " << l <<  " and has nh = " << nh << "\n\n";

	//@}

	dl = l/10 + 2;

	if (nh<(l - dl/2)) {

		return;

	}

	l1 = l + dl;

	for (int i=1; i<(l1+1); i++) {

		lstr[i] = true;

	}

	for (int j=1; j<(dl+1); j++) {

		while (true) {

			p = min(int(rand_num()*l1) + 1, l1);

			if (lstr[p]) {

				lstr[p] = false;

				break;

			} else {

				continue;

			}

		}

	}

	int k = 0;

	for (int c=1; c<(l1+1); c++) {

		if (lstr[c]) {

			k = k + 1;

			tmp1[c] = stra[k];

			tmp2[c] = strb[k];

		} else {

			tmp1[c] = 0;

			tmp2[c] = 0;

		}

	}

	l = l1;

	for (int d=1; d<(l+1); d++) {

		stra[d] = tmp1[d];

		strb[d] = tmp2[d];

	}

	pvectors();

	myfile1 << step << " increased l to " << l << "\n\n";

	//@ cout << "\n\n step " << step << " l is " << l << " ";  


}



//@----------------------------------------------------------------



void errcheck() {

	int s1, s2, ix;

	for (int i=1; i<(nx+1); i++) {

		spn1[i] = spn[i];

	}
	
	for (int j=1; j<(l+1); j++) {

		s1 = stra[j];

		s2 = strb[j];

		if ((s1<-1)||(s1>nx)) {

 			cout << "Illegal s1 operator on an errcheck() j = " << j << " where s1 = " << s1 << "\n\n";

			return;

		}
		
		if ((s2<0)||(s2>nx)) {

			cout << "Illegal s2 operator on an errcheck() j = " << j << " where s2 = " << s2 << "\n\n";

			return;

		}

		if (s1>0) {

			ix = disx[ x1[s1] ][ x1[s2] ];			

			if (!lis[ spn1[s1]+1 ][ spn1[s2]+1 ][ix]) {

				cout << "Illegal Ising bond on an errcheck() j = " << j << " where Ising bond = " << lis[ spn[s1]+1 ][ spn[s2]+1 ][ix] << ", ix = " << ix << "\n\n";

				return;

			}

		} else if (s1==-1) {

			spn1[s2] = -spn1[s2];

		}

	}

	for (int k=1; k<(nx+1); k++) {

		if (spn[k]!=spn1[k]) {

			cout << "Wrong state propagion, spin index = " << k << "\n\n"; /*@ Unless he meant the french word, propagion is probably a typo for propagation*/

			return;

		}

	}


}



//@----------------------------------------------------------------



/*
void writeconf(ofstream& myfile2) {
	
	myfile2 << "Lx   = " << nx << "\n\n";

	myfile2 << "Ht   = " << ht << "\n\n";

	myfile2 << "beta = " << beta << "\n\n";

	myfile2 << "L    = " << l << "\n\n";

	myfile2 << "----------------" << "\n\n\n";

	for (int i=1; i<(nx+1); i++) {

		myfile2 << i << " " << spn[i] << "\n";

	}

	myfile2 << "\n" << "----------------" << "\n\n\n";

	for (int j=1; j<(l+1); j++) {

		myfile2 << stra[j] << " " << strb[j] << "\n";

	}

	myfile2 << "\n" << "----------------" << "\n\n";

}
*/



//@----------------------------------------------------------------



void measure() {

	int s1, s2, mu, nu, ni, nt, nh1, ssum, last;
	long double su1, xu1, magpow2, magpow4;

	nu = 0;

	ni = 0;

	nt = 0;

	mu = 0;

	ssum = 0;

	last = 0;

	for (int i=1; i<(nx+1); i++) {

		mu = mu + spn[i];

	}

	magpow2 = pow(((long double)(mu)/((long double)(nx))), 2.0);

	//@ changed to long double

	magpow4 = pow(((long double)(mu)/((long double)(nx))), 4.0);

	//@ changed to long double

	//@ <M^2> and <M^4>

	su1 = (long double)(pow((long double)mu, 2.0));

	nh1 = 0;

	for (int j=1; j<(l+1); j++) {

		s1 = stra[j];

		if (s1!=0) {

			nh1 = nh1 + 1; //@ Hij where i=/=0, non identity

			s2 = strb[j];

			if (s1==-1) {

				nt = nt + 1; //@ H(-1)j, "spin flip" with flip and h strength?

				ssum = ssum + (nh1 - last)*mu;

				last = nh1;

				spn[s2] = -spn[s2];

				mu = mu + 2*spn[s2];

			} else if (s1==s2) {

				nu = nu + 1; //@ Hii, "spin flip" just h strength?

			} else {

				ni = ni + 1;

			}

			su1 = su1 + (long double)(mu*mu);

		}

	}

	ssum = ssum + (nh1 - last)*mu;

	su1 = su1/((long double)((nh1 + 1)*nx));

	if (nh1!=0) {

		xu1 = (long double)(pow((long double)ssum,2.0))/((long double)(nh1 * nx) * (long double)(nh1 + 1));

	} else {

		xu1 = 0.0;

	}

	xu1 = (long double)(beta)*(xu1 + su1/((long double)(nh1+1)));

	su = su + su1;

	xu = xu + xu1;

	avnu = avnu + (long double)(nu);

	avni = avni + (long double)(ni);

	avnt = avnt + (long double)(nt);

	avn1 = avn1 + (long double)(nh1);

	avn2 = avn2 + (long double)(pow( nh1, 2.0 ));

	nmsr = nmsr + 1;


	av_magpow2 = av_magpow2 + magpow2;

	av_magpow4 = av_magpow4 + magpow4;

	//@ <M^2> and <M^4> before dividing by nmsr

	
	q_2 = q_2 + (long double)pow((long double)(nt + nu) , 2.0);

	// <q^2> before dividing by nmsr

}



//@----------------------------------------------------------------



void results(ofstream& myfile3, ofstream& myfile4, ofstream& myfile6, ofstream& myfile7, ofstream& myfile8, ofstream& myfile9, ofstream& myfile10, ofstream& myfile11, ofstream& myfile12, ofstream& myfile13) {
	//~ ************ CHANGING


	long double c, e;

	xu = xu/((long double)(nmsr));

	su = su/((long double)(nmsr));


	q = (avnu + avnt)/((long double)(nmsr));

	//@ average # of [i,i] operators + # of [i,0] operators

	q_2 = q_2/((long double)(nmsr));

	//@ average of the square of # of [i,i] operators + # of [i,0] operators

	avnu = avnu/((long double)(nmsr));

	avni = avni/((long double)(nmsr));

	avnt = avnt/((long double)(nmsr));

	avn1 = avn1/((long double)(nmsr));

	avn2 = avn2/((long double)(nmsr));

	c = (avn2 - pow(avn1,2.0) - avn1)/((long double)(nx));

	e = eshft - (avn1 - avnu)/((long double)(beta) * (long double)(nx));

	avnu = avnu/(ht * (long double)(beta) * (long double)(nx));

	avni = eshft/((long double)(nx)) - avni/((long double)(beta) * (long double)(nx));

	avnt = -avnt/((long double)(beta) * (long double)(nx));

	av_magpow2 = av_magpow2/((long double)(nmsr));

	av_magpow4 = av_magpow4/((long double)(nmsr));

	//@ <M^2> and <M^4>

	long double Binder_cumulant = ((long double) 1) - ( ((long double)av_magpow4) / (((long double) 3) * ((long double)(pow(av_magpow2, 2.0)) )) );

	//@*** CHANGED TO 2.0 IN POW

	//@ Binder_cumulant

	mx = ( (long double)(temp) / ((long double)(ht))) * q - (long double)(nx);

	//@ transverse field: <Mx> = (T/h)*<q> - N

	dmx_dh = ( (long double)(temp) / ((long double)(pow((long double)(ht), 2.0)))) * ( q_2 - pow(q, 2.0) - q);

	//@ transverse field susceptibility: d<Mx>/dh = (T/h^2)*( <q^2> - <q>^2 - <q> )

	myfile3 << su << " " << xu << "\n";

	myfile4 << e << " " << c << " " << avni << " " << avnt << " " << avnu << "\n";

	if (temp == final_temp){

		myfile6 << c;

		myfile7 << av_magpow2;

		myfile8 << e;

		myfile9 << Binder_cumulant;

		myfile10 << q;

		myfile11 << q_2;

		myfile12 << mx;

		myfile13 << dmx_dh;

	} else {

		myfile6 << c << "\n";

		myfile7 << av_magpow2 << "\n";

		myfile8 << e << "\n";

		myfile9 << Binder_cumulant << "\n";

		myfile10 << q << "\n";

		myfile11 << q_2 << "\n";

		myfile12 << mx << "\n";

		myfile13 << dmx_dh << "\n";

	}

	zerodat(); 

}



//@----------------------------------------------------------------



void zerodat() {
	
	su = 0.0;

	xu = 0.0;

	avnu = 0.0;

	avni = 0.0;

	avnt = 0.0;

	avn1 = 0.0;

	avn2 = 0.0;

	nmsr = 0;


	av_magpow2 = 0.0;

	av_magpow4 = 0.0;

	q = 0.0;

	mx = 0.0;

	dmx_dh = 0.0;

	q_2 = 0.0;

}



//@----------------------------------------------------------------



void writeacc(ofstream& myfile5) {

	ar1 = ar1/((long double)(msteps));

	ar2 = ar2/((long double)(msteps));

	ar3 = ar3/((long double)(msteps));

	ar4 = ar4/((long double)(msteps));

	myfile5 << "Diagonal acceptance 1     : " << ar1 << "\n"; 

	myfile5 << "Diagonal acceptance 2     : " << ar2 << "\n"; 

	myfile5 << "Off-Diagonal substitions  : " << ar3 << "\n"; 

	myfile5 << "Spin flips / site         : " << ar4 << "\n"; 

	myfile5 << "Max lsub                  : " << mlls << "\n\n";
 
}



//@----------------------------------------------------------------



//void calcenrcool(int step, ofstream& myfile8) {

	//int nh_c = 0;

	//for (int i=1; i<(l+1); i++) {

		//if (stra[i]!=0) {

			//nh_c = nh_c + 1;

		//}

	//}

	//myfile8 << step << " " << nh_c << "\n";


//}



//@----------------------------------------------------------------



void writeopen(ofstream& myfile1, ofstream& myfile3, ofstream& myfile4, ofstream& myfile5, ofstream& myfile6, ofstream& myfile7, ofstream& myfile8, ofstream& myfile9, ofstream& myfile10, ofstream& myfile11, ofstream& myfile12, ofstream& myfile13) {

	string nxstring;

	string htstring;

	string starttempstring;

	string finaltempstring;

	if (nx==1024) {

		nxstring = "1024";

	} else if (nx==512) {

		nxstring = "512";

	} else if (nx==256) {

		nxstring = "256";

	} else if (nx==128) {

		nxstring = "128";

	} else if (nx==64) {

		nxstring = "64";

	} else if (nx==32) {

		nxstring = "32";

	} else if (nx==16) {

		nxstring = "16";

	}

	if (ht==4.0) {

		htstring = "4.0";

	} else if (ht==3.5) {

		htstring = "3.5";

	} else if (ht==3.0) {

		htstring = "3.0";

	} else if (ht==2.5) {

		htstring = "2.5";

	} else if (ht==2.0) {

		htstring = "2.0";

	} else if (ht==1.5) {

		htstring = "1.5";

	} else if (ht==1.0) {

		htstring = "1.0";

	} else if (ht==0.5) {

		htstring = "0.5";

	}


	if (start_temp==0.50) {

		starttempstring = "0.50";

		if (final_temp==0.54) {
		
			finaltempstring = "0.54";

		} else if (final_temp==2.16) {

			finaltempstring = "2.16";

		} else if (final_temp==5.49) {

			finaltempstring = "5.49";

		}

  } else if (start_temp==0.55) {

    starttempstring = "0.55";

    finaltempstring = "0.59";

  } else if (start_temp==0.60) {

    starttempstring = "0.60";

    finaltempstring = "0.64";

  } else if (start_temp==0.65) {

    starttempstring = "0.65";

    finaltempstring = "0.69";

  } else if (start_temp==0.70) {

    starttempstring = "0.70";

    finaltempstring = "0.74";

  } else if (start_temp==0.75) {

    starttempstring = "0.75";

    finaltempstring = "0.79";

  } else if (start_temp==0.80) {

    starttempstring = "0.80";

    finaltempstring = "0.84";

  } else if (start_temp==0.85) {

    starttempstring = "0.85";

    finaltempstring = "0.89";

  } else if (start_temp==0.90) {

    starttempstring = "0.90";

    finaltempstring = "0.94";

  } else if (start_temp==0.95) {

    starttempstring = "0.95";

    finaltempstring = "0.99";

  } else if (start_temp==1.00) {

    starttempstring = "1.00";

    finaltempstring = "1.04";

  } else if (start_temp==1.05) {

    starttempstring = "1.05";

    finaltempstring = "1.09";

  } else if (start_temp==1.10) {

    starttempstring = "1.10";

    finaltempstring = "1.14";

  } else if (start_temp==1.15) {

    starttempstring = "1.15";

    finaltempstring = "1.19";

  } else if (start_temp==1.20) {

    starttempstring = "1.20";

    finaltempstring = "1.24";

  } else if (start_temp==1.25) {

    starttempstring = "1.25";

    finaltempstring = "1.29";

  } else if (start_temp==1.30) {

    starttempstring = "1.30";

    finaltempstring = "1.34";

  } else if (start_temp==1.35) {

    starttempstring = "1.35";

    finaltempstring = "1.39";

  } else if (start_temp==1.40) {

    starttempstring = "1.40";

    finaltempstring = "1.44";

  } else if (start_temp==1.45) {

    starttempstring = "1.45";

    finaltempstring = "1.49";

  } else if (start_temp==1.50) {

    starttempstring = "1.50";

    finaltempstring = "1.54";

  } else if (start_temp==1.55) {

    starttempstring = "1.55";

    finaltempstring = "1.59";

  } else if (start_temp==1.60) {

    starttempstring = "1.60";

    finaltempstring = "1.64";

  } else if (start_temp==1.65) {

    starttempstring = "1.65";

    finaltempstring = "1.69";

  } else if (start_temp==1.70) {

    starttempstring = "1.70";

    finaltempstring = "1.74";

  } else if (start_temp==1.75) {

    starttempstring = "1.75";

    finaltempstring = "1.79";

  } else if (start_temp==1.80) {

    starttempstring = "1.80";

    finaltempstring = "1.84";

  } else if (start_temp==1.85) {

    starttempstring = "1.85";

    finaltempstring = "1.89";

  } else if (start_temp==1.90) {

    starttempstring = "1.90";

    finaltempstring = "1.94";

  } else if (start_temp==1.95) {

    starttempstring = "1.95";

    finaltempstring = "1.99";

  } else if (start_temp==2.00) {

    starttempstring = "2.00";

    finaltempstring = "2.04";

  } else if (start_temp==2.05) {

    starttempstring = "2.05";

    finaltempstring = "2.09";

  } else if (start_temp==2.10) {

    starttempstring = "2.10";

    finaltempstring = "2.14";

  } else if (start_temp==2.15) {

    starttempstring = "2.15";

    finaltempstring = "2.24";

  } else if (start_temp==2.25) {

    starttempstring = "2.25";

    finaltempstring = "2.34";

  } else if (start_temp==2.35) {

    starttempstring = "2.35";

    finaltempstring = "2.44";

  } else if (start_temp==2.45) {

    starttempstring = "2.45";

    finaltempstring = "2.54";

  } else if (start_temp==2.55) {

    starttempstring = "2.55";

    finaltempstring = "2.64";

  } else if (start_temp==2.65) {

    starttempstring = "2.65";

    finaltempstring = "2.74";

  } else if (start_temp==2.75) {

    starttempstring = "2.75";

    finaltempstring = "2.84";

  } else if (start_temp==2.85) {

    starttempstring = "2.85";

    finaltempstring = "2.94";

  } else if (start_temp==2.95) {

    starttempstring = "2.95";

    finaltempstring = "3.04";

  } else if (start_temp==3.05) {

    starttempstring = "3.05";

    finaltempstring = "3.14";

  } else if (start_temp==3.15) {

    starttempstring = "3.15";

    finaltempstring = "3.24";

  } else if (start_temp==3.25) {

    starttempstring = "3.25";

    finaltempstring = "3.34";

  } else if (start_temp==3.35) {

    starttempstring = "3.35";

    finaltempstring = "3.44";

  } else if (start_temp==3.45) {

    starttempstring = "3.45";

    finaltempstring = "3.54";

  } else if (start_temp==3.55) {

    starttempstring = "3.55";

    finaltempstring = "3.64";

  } else if (start_temp==3.65) {

    starttempstring = "3.65";

    finaltempstring = "3.74";

  } else if (start_temp==3.75) {

    starttempstring = "3.75";

    finaltempstring = "3.84";

  } else if (start_temp==3.85) {

    starttempstring = "3.85";

    finaltempstring = "3.99";

  } else if (start_temp==4.00) {

    starttempstring = "4.00";

    finaltempstring = "4.14";

  } else if (start_temp==4.15) {

    starttempstring = "4.15";

    finaltempstring = "4.29";

  } else if (start_temp==4.30) {

    starttempstring = "4.30";

    finaltempstring = "4.44";

  } else if (start_temp==4.45) {

    starttempstring = "4.45";

    finaltempstring = "4.59";

  } else if (start_temp==4.60) {

    starttempstring = "4.60";

    finaltempstring = "4.74";

  } else if (start_temp==4.75) {

    starttempstring = "4.75";

    finaltempstring = "4.89";

  } else if (start_temp==4.90) {

    starttempstring = "4.90";

    finaltempstring = "5.04";

  } else if (start_temp==5.05) {

    starttempstring = "5.05";

    finaltempstring = "5.19";

  } else if (start_temp==5.20) {

    starttempstring = "5.20";

    finaltempstring = "5.34";

  } else if (start_temp==5.35) {

    starttempstring = "5.35";

    finaltempstring = "5.49";

	} else if (start_temp==2.17) {

		starttempstring = "2.17";

		finaltempstring = "3.82";

	} else if (start_temp==3.83) {

		starttempstring = "3.83";

		finaltempstring = "5.49";

	}

	string log_txt = "log_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	//@ string config_txt = "config_nx=" + to_string(nx) + "_task5_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + ".txt";

	string mag_txt = "mag_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	string enr_txt = "enr_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	string acc_txt = "acc_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	string c_txt = "sp_heat_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	string mag_pow2_txt = "mag_pow2_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	string avg_en_txt = "avg_en_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	string binder_txt = "binder_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	string q_txt = "q_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	string q_2_txt = "q_2_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	string mx_txt = "mx_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	string dmx_dh_txt = "dmx_dh_nx=" + nxstring + "_ht=" + htstring + "_temps" + starttempstring + "to" + finaltempstring + "_TFIMSSE.txt";

	myfile1.open (log_txt.c_str());

	//@ myfile2.open(config_txt);

	myfile3.open (mag_txt.c_str());

	myfile4.open (enr_txt.c_str());

	myfile5.open (acc_txt.c_str());

	myfile6.open (c_txt.c_str());

	myfile7.open (mag_pow2_txt.c_str());

	myfile8.open (avg_en_txt.c_str());

	myfile9.open (binder_txt.c_str());

	myfile10.open (q_txt.c_str());

	myfile11.open (q_2_txt.c_str());

	myfile12.open (mx_txt.c_str());

	myfile13.open (dmx_dh_txt.c_str());

}



//@----------------------------------------------------------------



void zero_between_temps () {

	mlls = 0;

	ar1 = 0.0;

	ar2 = 0.0;

	ar3 = 0.0;

	ar4 = 0.0;

	for (int c=0; c<(ll+1); c++) {

		ap1[c] = 0.0;

		dp1[c] = 0.0;

		tmp1[c] = 0;

		tmp2[c] = 0;

		lstr[c] = 0;

	}

	for (int d=0; d<(nx+1); d++) {

		pfrst[d] = 0;

		plast[d] = 0;

		lsub[d] = 0;

		spn1[d] = 0;

		for (int f=0; f<(lls+1); f++) {

			pos1[f][d] = 0;

			str1[f][d] = 0;

			con1[f][d] = 0;

			pos[f] = 0;

			opr[f] = 0;

		}

	}


}



//@----------------------------------------------------------------
