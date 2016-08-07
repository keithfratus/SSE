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
//@		2. Removed random number generator parts and using MTRand
//@		3. openlog(), closelog() functions.
//@   4. readconf() not used because doesn't seem like we'd want to?
//@==================================================================//

/*@ So obviously I'm using quite a bit of global variables. This is mostly because it's
immediately easier, and seemed to be what you started to do (since this isn't some code
for a big company, etc.). However, I can certainly switch to lessen the amount of 
global variables should you suggest so. */  

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "MT/mtrand.h"

using namespace std;

//@----------------------------------------------------------------

//@ variables (is.h)

//@----------------------------------------------------------------

const int nx = 8;

//@ size in x-direction (since 1d, it's only direction)

//~ int ny = 1;

//~ size in y-direction (remove later if easily translate into just 1d considerations)

//~ int nz = 1;

//~ size in z-direction (remove later if easily translate into just 1d considerations)

//~ int n3 = nx*ny*nz;

//~ total number of spins (remove later if easily translate into just 1d considerations)

const int ll = 10000;

//@ max value for expansion truncation L (hamiltonian string cut off max)

const int lls = 500;

/*$ maximum length of subsequence, what does this mean */

int l;

//@ expansion truncation, or hamiltonian string cut off max, note: hard-coded to be 20 in initconf()

double ht = 0.5;

// transverse field strength

double beta = 1.0;

// inverse temperature

double eshft;

/*$ sum of all constants added to the hamiltonian (see hamiltonian equations perhaps?) */

int nh;

//@ current expansion order of n, or number of non-identity hamiltonians in hamiltonian string

int mlls;

//@ maximum substring length that has been reached

double ar1, ar2, ar3, ar4;

//@ acceptance rates

//@----------------------------------------------------------------

int isteps = 100000;

//@ equilibriation steps

int msteps = 10000;

//@ measurement steps 

int nd = 1;

// number of runs

int nmsr;

//$

double su, xu, avnu, avni, avnt, avn1, avn2;

//$

//@----------------------------------------------------------------

//@----------------------------------------------------------------

//@ initializing arrays (is.h)

//@----------------------------------------------------------------

int spn[nx+1] = {0};

//@ //~ spin array, flag n3

int stra[ll+1] = {0};
int strb[ll+1] = {0};

/*$ (NOTE reason for names, see str1 below) operator string(s)? uses ll which is the max expansion truncation. go down p list decides which operator it is. i and j can be much farther apart, but there are still two things to specify hence the 1,2 part of the multidimensions first array gives i second array gives j put together to get Hij*/

int x1[nx+1] = {0}; //~ changed xyz1 to x1
//~ int xyz2[n3] = {0};
//~ int xyz3[n3] = {0};

/*@ //~ corrdinates for given spin i, where x-coordinate is xyz1, y-coordinate is xyz2, z-coordinate is xyz3. note: Sandvik (I believe) makes a typo when he typed that y-coordinate is equal to the same array as x-coordinate. */

//~ int xyzi[n3] = {0};

/*~ this one gets hairy. Sandvik has it originally as xyzi(nx,ny,nz) which is a 3 dimensional array in fortran77 representing spin number given by coordinates x,y,z = xyzi(x,y,z). so although I am not exactly following him at this point, it's likely best to switch to a 1D chosen from here on out, hence why n3 still works at the moment. Anyway, flagged for modification. */

int disx[nx+1][nx+1] = {0};
//~ int disy[ny*ny] = {0};
//~ int disz[nz*nz] = {0};

/*~ this one also could get hairy. although this seems like the best route, just watch for the math. represents the periodic distances in x,y, and z directions between two points (remove later if easily translate into just 1d considerations). possibly still needs modifications as well since original is disx(nx,nx) and etc. */

int pfrst[nx+1] = {0};

/*$ //~ first position to search in probability list given INT(R*n3), what is R (random number generator?) and is this similar to frstspinop array? */

int plast[nx+1] = {0};

/*$ //~ last position to search in probability list given INT(R*n3), what is R (random number generator?) and is this similar to lastspinop array? */

bool lis[3][3][nx/2 + 1] = {0};

/*$ //~ For allowed spin configuration given distance between sites. ok so original is basically the equivalent of a boolean 5 dimensional array it seems: LOGICAL lis(-1:1,-1:1,0:nx/2,0:ny/2,0:nz/2). I just split into 4 arrays. */

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

/*$ //~ positions?, watch the math */

int str1[lls+1][nx+1] = {0};

//$ //~ substring operators?, watch the math (NOTE above name stra and strb arrays)

bool con1[lls+1][nx+1] = {0};

//$ //~ constraints? watch the math

//@----------------------------------------------------------------

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

//@----------------------------------------------------------------

//@ initializing functions directly in main

//@----------------------------------------------------------------

void mcstep(MTRand_open& rand_num);

//@ monte carlo step

void checkl(int step, ofstream& myfile1, MTRand_open& rand_num);

/*$ check the hamiltonian string length perhaps? has pvectors, openlog, and closelog */

void init(MTRand_open& rand_num);

//@ initializes the system: initconf, lattice, pvectors, isingpart, zerodat

//@ void openlog();

/*$ looks like writes information to a (12) log.txt (note a 10?), note the closelog function so probably writing the equilibrium step thing to log.txt */

//@ void closelog();

//@ close the (12) log.txt

void errcheck();

/*$ looks like checks if certain spins weren't possible in the first place? uses str(1,i) and str(2,i) arrays */

void writeconf(ofstream& myfile2);

/*$ writes length of lattice in x,y,z directions (we'll just use x), beta, L (hamiltonian string length?), spn(i) spin array values, str(1,i) (bond array indicates site left for our case?), str(2,i) (bond array indicates site right for our case?) uses (10) (20)*/ 

void measure();

/*$ measuring observables looks like. note: uses a common which seems like fortran77's version of a sort of pass by reference across different fortran .f files*/

void results(ofstream& myfile3, ofstream& myfile4);

/*@ writes calculated observables into (10) mag.dat for su and xu, (10) enr.dat for e, c, avni, avnt, and avnu */

void writeacc(int msteps, ofstream& myfile5);

/*$ writes diagonal acceptance 1, diagonal acceptance 2, off-diagonal substitutions, spin flips / site, max lsub (max hamiltonian string length?) */

//@----------------------------------------------------------------

//@----------------------------------------------------------------

//@ initializing functions not directly in main

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

//@----------------------------------------------------------------

//@----------------------------------------------------------------

//@ main function

//@----------------------------------------------------------------

int main() {

	ofstream myfile1; //@ Initialize the write-data-to-file.
	myfile1.open ("log.txt"); //@ Name and begin the write-data-to-file.
	ofstream myfile2;
	myfile2.open ("writeconf.txt");
	ofstream myfile3;
	myfile3.open ("mag.txt");
	ofstream myfile4;
	myfile4.open ("enr.txt");
	ofstream myfile5;
	myfile5.open ("acc.txt");

	unsigned long int time_seed = time(NULL);
	MTRand_open rand_num(time_seed);

	//@ random number generator (0,1). Note: omitted RAN(), INITRANDOM, and WRITERAND(rndf)

	init(rand_num);

	//@ consider Sandvik's (in == 0)

	cout << "I have made a configuration and have started first mcstep loop" << "\n\n";

	for (int i=1; i<(isteps+1); i++) {

		mcstep(rand_num);

		checkl(i, myfile1, rand_num);

		if ((i%(isteps/10))==0) { //@ to record every 10 steps
	
			//@ openlog();

			myfile1 << "Done equilibriation step: " << i << "\n\n";	
		
			//@ closelog();

			errcheck();

			//cout << i << " " << stra[1] << " " << strb[1] << endl;

		}

	}

	cout << "I have finished the first mcstep loops, error checked, and etc." << "\n\n";

	writeconf(myfile2);

	cout << "I have written down the configuration for that loop" << "\n\n";

	for (int j=1; j<(nd+1); j++) {

	//cout << "Loop test " << j << endl;

		ar1 = 0.0;
		ar2 = 0.0;
		ar3 = 0.0;
		ar4 = 0.0;

		for (int k=1; k<(msteps+1); k++) {

			mcstep(rand_num);

			if ((k%2)==0) {

				measure();

			}

	//cout << "Main test b" << endl;

			if ((k%(msteps/10))==0) { //@ to record every 10 steps

				//@ openlog();

				myfile1 << "Done measurement step: " << j << ", " << k << "\n\n";

				//@ closelog();

				errcheck();

			}

		}

		results(myfile3, myfile4);

		//cout << "Main res test" << endl;

		writeacc(msteps, myfile5);

		writeconf(myfile2);

	}
	cout << "I have finished the second loop, wrote results, and etc." << "\n\n";
	myfile1.close(); //@ Finish the write-data-to-file.
	myfile2.close();
	myfile3.close();
	myfile4.close();
	myfile5.close();

	return 0;


}

//@----------------------------------------------------------------

//@----------------------------------------------------------------

//@ Helper functions

//@----------------------------------------------------------------

void init(MTRand_open& rand_num) {

	//@ omitted readconf function
	initconf(rand_num);
	lattice();
	pvectors();
	isingpart();
	zerodat();

}

//@----------------------------------------------------------------

void initconf(MTRand_open& rand_num) {
	
	int c;
	
	l = 20;
	nh = 0;
	for (int i=1; i<(ll+1); i++) {
		stra[i] = 0; //@ gives i of Hij, these are relatively uneeded since already initialized all to 0
		strb[i] = 0; //@ gives j of Hij, these are relatively uneeded since already initialized all to 0
	}

	for (int j=1; j<(nx+1); j++) {
		spn[j] = -1;
	}

	//@ setting all spins to -1 first

	for (int k=1; k<(nx/2 + 1); k++) { /*$ note confirm how this is supposed to play out if n3/2=integer or n3/2=floating point=rounded integer in fortran77 version at least */
		c = min((int(rand_num()*nx) + 1), nx);
		
		while (true) {
			if (spn[c] == -1) {
				break;
			}
			c = min ((int(rand_num()*nx) + 1), nx);
		}

/*$ NOTE: so this one may be tricky or I may be over thinking it. at first I just thought it was a simple goto continue basically, but he doesn't place to "goto 10"'s 10 at the end of the loop, he puts it at the first min(,) statement. which I believe means the code wants to repeat finding another random site that IS spn[c] = -1 and flip that to spn[c] 1 so that EXACTLY (if n3 is an even number, roughly half if n3 is odd) half of the spins are +1 (spin up) and half the spins are (-1) after randomly choosing spins and flipping or leaving the same +1 (spin up), until this is true. So I'm going to simply implement a forever loop that accomplishes that for now (apparently the internet hates goto statements) unless I think of something better. */

		spn[c] = 1;
	}
	
}

//@----------------------------------------------------------------

void lattice() { //~ lots of omits due to original using 3 dimensions
	int i = 0;
	for (int ix=1; ix<(nx+1); ix++) {
		i += 1;		
		x1[i] = ix;
	}
	
	for (int j=1; j<(nx+1); j++) { //@ (think rows)
		for (int k=1; k<(nx+1); k++) { //@ (think cols)
			disx[k][j] = abs(k - j);
			if ((disx[k][j])>(nx/2)) { 
				disx[k][j] = nx - disx[k][j]; /*@ think 1D chain with PBC so it's like a ring which means max disx is nx/2 */
			}
		}
	}

}

//@----------------------------------------------------------------

void pvectors() {

	for (int i=0; i<l; i++ ) { //@ see above comment
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

	for (ix=0; ix<(nx/2 + 1); ix++) { /*$ //~ note because see lis..[max], test different values perhaps, it's because he introduces a 0 index when fortran arrays are usually 1 to number inclusive both ends, whereas c++ is 0 to number inclusive for 0 and exclusive for number */
		if (ix==0) {
			ris[ix] = ht;
			lis[2][2][ix] = true;
			lis[2][0][ix] = false;
			lis[0][2][ix] = false;
			lis[0][0][ix] = true;
		} else {
			r1 = double(ix);
			r2 = double(nx - ix);
			ris[ix] = (-1.0)/(r1*r1) - 1.0/(r2*r2);
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
	pint[1] = ris[0]; //$ //~ perhaps should switch with max ris value instead, but unsure
	eshft = 0.0;
	for (int i=2; i<(nx+1); i++) {
		ix = disx[1][x1[i]]; /*@ to be clear, this is finding the distance between the first site (site x1[0]) and every other site (site x1[i]), note: 1 to nx*/
		pint[i] = pint[i-1] + ris[ix];
		eshft += (ris[ix]/2);
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
					break; //@ go to next if statement the if (k==nx) one
				}
			}
		}

		if (k==nx) {
			plast[k] = nx;
		} else {
			for (int d=pfrst[k]; d<(nx+1); d++) {
				if (pint[d]>=p2) {
					plast[k] = d;
					break; /*@ go to next statement which is beginning of the for (int k=0; k<(nx+1); k++) loop */
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

	int s1, s2, p0, p1, p2, ix, nac1, nac2, ntr2;
	double p;
		//cout << "Begin diagonal update" << "\n\n";
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
				nh += 1;
				nac1 += 1;
			} else if ((rand_num())<p) {
				stra[i] = min(int(rand_num()*nx)+1, nx);
				strb[i] = stra[i];
				nh += 1;
				nac1 += 1;
			}
			
		} else if (s1==s2) { //@ Hii
			p = dp1[nh];
			if (p>=1) {
				stra[i] = 0;
				strb[i] = 0;
				nh -= 1;
				nac1 += 1;
			} else if (rand_num()<p) {
				stra[i] = 0;
				strb[i] = 0;
				nh -= 1;
				nac1 += 1;
			}
					
		} else if (s1==-1) { //@ H(-1)j
			spn[s2] = -spn[s2];
		}
					
		s1 = stra[i]; //$ why did he redefine s1
						
		if (s1>0) { //@ doesn't this work on Hii as well if i>0
					
		//cout << "enter if statement" << endl;

			while (true) {
				ntr2 += 1;
				p = rand_num();
				p0 = int(p*nx);
				p1 = pfrst[p0];
				p2 = plast[p0];

				fromlastelse:
				//@ goto from far below

				if (p1==p2) {
					s2 = ((s1 + p1 - 2) % nx) + 1;
					ix = disx[x1[s1]][x1[s2]];
					if ((lis[spn[s1]+1][spn[s2]+1][ix])) {
						strb[i] = s2;
						nac2 += 1;
						//cout << "b1" << endl;
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
					ix = disx[x1[s1]][x1[s2]];		
					if ((lis[spn[s1]+1][spn[s2]+1][ix])) {
						strb[i] = s2;
						nac2 += 1;
						//cout << "b2" << endl;
						break;
					} else {
						continue;
					}

				} else {
					p0 = (p1 + p2)/2;
					if (((pint[p0-1])<p)&&((pint[p0])>=p)) {
						s2 = ((s1 + p0 - 2) % nx) + 1;
						ix = disx[x1[s1]][x1[s2]];
						if ((lis[spn[s1]+1][spn[s2]+1][ix])) {
							strb[i] = s2;
							nac2 += 1;
							//cout << "b3" << endl;
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
						goto fromlastelse;
						//@ see "near" beginning of while(true) statement
					}
				}	
			}
		}
	}
			//cout << "End diagonal update" << "\n\n";
	ar1 += double(nac1)/double(l);
	if (ntr2!=0) {
		ar2 += double(nac2)/double(ntr2);
	}

}

//@----------------------------------------------------------------

void partition() {
		//cout << "I have this part of the mcstep and am on the partition function" << "\n\n";
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
				lsub[s1] += 1;
				pos1[lsub[s1]][s1] = j;
				str1[lsub[s1]][s1] = 1;
				con1[lsub[s1]][s1] = false;
			} else if (s1==(-1)) {
				lsub[s2] += 1;
				pos1[lsub[s2]][s2] = j;
				str1[lsub[s2]][s2] = 2;
				con1[lsub[s2]][s2] = false;
			} else {
				con1[lsub[s1]][s1] = true;
				con1[lsub[s2]][s2] = true;
			}
		}
	}
	
	for (int k=1; k<(nx+1); k++) {
		if (con1[0][k]) {
			con1[lsub[k]][k] = true;
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
			nupd += 1;
			pos[nupd] = i;
		}
		opr[i] = str1[i][s];
	}

	if (lens<2) {
		if ((nupd == lens)&&((!con1[0][s]))) { //$ ?
			if (rand_num()<0.75) {
				spn[s] = -spn[s];
				ar4 += 1.0;
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
			nac += 1;
		} else {
			top = opr[p1];
			opr[p1] = opr[p2];
			opr[p2] = top; //@ why even have top if this is all its used for?
		}
		
		if (p2<p1) {
		 flp = -flp;
		}
	}

	spn[s] *= flp;

	for (int k=1; k<(lens+1); k++) {
		ii = pos1[k][s];
		if (opr[k]==1) {
			stra[ii] = strb[ii];
		} else {
			stra[ii] = -1;
		}
	}

	ar3 += double(nac);

	if ((nupd == lens)&&((!con1[0][s]))) {
		if (rand_num()<0.75) {
			spn[s] = -spn[s];
			ar4 += 1.0;
		}
	}

}

//@----------------------------------------------------------------

void checkl(int step, ofstream& myfile1, MTRand_open& rand_num) {

	//myfile1 << "MF1 test \n\n";

	int p, dl, l1;

	dl = l/10 + 2;
	if (nh<(l - dl/2)) {
		return;
	}

//cout << "test b" << endl;

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
			k += 1;
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

	//cout << "A screen test" << endl;
	
	myfile1 << step << " increased l to " << l << "\n\n";
	//myfile1 << "MF2 test \n\n";

	//cout << "MF2 screen test" << endl;


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
			//$ writes "illegal s1 operator" to somewhere seemingly not log.txt
 			cout << "Illegal s1 operator on an errcheck() j = " << j << " where s1 = " << s1 << "\n\n";
			return;
		}
		
		if ((s2<0)||(s2>nx)) {
			cout << "Illegal s2 operator on an errcheck() j = " << j << " where s2 = " << s2 << "\n\n";
			return;
		}

		if (s1>0) {
			ix = disx[x1[s1]][x1[s2]];			
			if ((!lis[spn1[s1]+1][spn1[s2]+1][ix])) {
				cout << "Illegal Ising bond on an errcheck() j = " << j << " where Ising bond = " << lis[spn[s1]+1][spn[s2]+1][ix] << ", ix = " << ix << "\n\n";
				return;
			}
		} else if (s1==-1) {
			spn1[s2] = -spn1[s2];
		}
	}

	for (int k=1; k<(nx+1); k++) {
		if (spn[k]!=spn1[k]) {
			cout << "Wrong state propagion, spin index = " << k << "\n\n";
			return;
		}
	}

}

//@----------------------------------------------------------------

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

//@----------------------------------------------------------------

void measure() {

	int s1, s2, mu, nu, ni, nt, nh1, ssum, last;
	double su1, xu1;

	nu = 0;
	ni = 0;
	nt = 0;
	mu = 0;
	ssum = 0;
	last = 0;
	for (int i=1; i<(nx+1); i++) {
		mu += spn[i];
	}

	su1 = double(pow(mu,2));
	nh1 = 0;
	for (int j=1; j<(l+1); j++) {
		s1 = stra[j];
		if (s1!=0) {
			nh1 += 1;
			s2 = strb[j];
			if (s1==-1) {
				nt += 1;
				ssum += (nh1 - last)*mu;
				last = nh1;
				spn[s2] = -spn[s2];
				mu += 2*spn[s2];
			} else if (s1==s2) {
				nu += 1;
			} else {
				ni += 1;
			}

			su1 += double(mu*mu);
		}
	}

	ssum += (nh1 - last)*mu;
	su1 = su1/double((nh1 + 1)*nx);
	if (nh1!=0) {
		xu1 = double(pow(ssum,2))/(double(nh1*nx)*double(nh1 + 1));
	} else {
		xu1 = 0.0;
	}

	xu1 = double(beta)*(xu1 + su1/double(nh1+1));
	su += su1;
	xu += xu1;
	avnu += double(nu);
	avni += double(ni);
	avnt += double(nt);
	avn1 += double(nh1);
	avn2 += double(pow(nh1,2));
	nmsr += 1;


}

//@----------------------------------------------------------------

void results(ofstream& myfile3, ofstream& myfile4) {

	//cout << "res test A" << endl;

	//cout << "su test " << su << endl;

	double c, e;

	xu = xu/double(nmsr);
	su = su/double(nmsr);
	avnu = avnu/double(nmsr);
	avni = avni/double(nmsr);
	avnt = avnt/double(nmsr);
	avn1 = avn1/double(nmsr);
	avn2 = avn2/double(nmsr);
	c = (avn2 - pow(avn1,2) - avn1)/double(nx);
	e = eshft - (avn1 - avnu)/(double(beta)*double(nx));
	avnu = avnu/(ht*double(beta)*double(nx));
	avni = eshft/double(nx) - avni/(double(beta)*double(nx));
	avnt = -avnt/(double(beta)*double(nx));

	myfile3 << su << " " << xu << "\n";
	myfile4 << e << " " << c << " " << avni << " " << avnt << " " << avnu << "\n";

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

}

//@----------------------------------------------------------------

void writeacc(int msteps, ofstream& myfile5) {

	ar1 = ar1/double(msteps);
	ar2 = ar2/double(msteps);
	ar3 = ar3/double(msteps);
	ar4 = ar4/double(msteps);
	myfile5 << "Diagonal acceptance 1     : " << ar1 << "\n"; 
	myfile5 << "Diagonal acceptance 2     : " << ar2 << "\n"; 
	myfile5 << "Off-Diagonal substitions  : " << ar3 << "\n"; 
	myfile5 << "Spin flips / site         : " << ar4 << "\n"; 
	myfile5 << "Max lsub                  : " << mlls << "\n\n";
 
}

//@----------------------------------------------------------------
