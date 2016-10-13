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



//@----------------------------------------------------------------
//@ Main function
//@----------------------------------------------------------------



int main() {

	ifstream writelog;

	//@ ofstream writeconfig;

	ifstream writemag;

	ifstream writeenr;

	ifstream writeac;

	ifstream writec_temp;

	ifstream writemag_2vstemp;

	ifstream writeevstemp;

	ifstream writebindercumulantvstemp;

	ifstream writeq;

	ifstream writeq_2;

	ifstream writemx;

	ifstream writedmx_dh;

//-----

	ifstream writelog;

  stringstream writelog_ss;

  writelog_ss << "log_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writelog_s = writelog_s.str();

  char * writelog_str = new char [writelog_s.length()+1];

  strcpy(writelog_str, writelog_s.c_str());

  writelog.open(writelog_str);

//-----

	ifstream writemag;

  stringstream writemag_ss;

  writemag_ss << "mag_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writemag_s = writemag_s.str();

  char * writemag_str = new char [writemag_s.length()+1];

  strcpy(writemag_str, writemag_s.c_str());

  writemag.open(writemag_str);

//-----

	ifstream writeenr;

  stringstream writeenr_ss;

  writeenr_ss << "enr_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writeenr_s = writeenr_s.str();

  char * writeenr_str = new char [writeenr_s.length()+1];

  strcpy(writeenr_str, writeenr_s.c_str());

  writeenr.open(writeenr_str);

//-----

	ifstream writeac;

  stringstream writeac_ss;

  writeac_ss << "ac_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writeac_s = writeac_s.str();

  char * writeac_str = new char [writeac_s.length()+1];

  strcpy(writeac_str, writeac_s.c_str());

  writeac.open(writeac_str);

//-----

	ifstream writec_temp;

  stringstream writec_temp_ss;

  writec_temp_ss << "sp_heat_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writec_temp_s = writec_temp_s.str();

  char * writec_temp_str = new char [writec_temp_s.length()+1];

  strcpy(writec_temp_str, writec_temp_s.c_str());

  writec_temp.open(writec_temp_str);

//-----

	ifstream writemag_2vstemp;

  stringstream writemag_2vstemp_ss;

  writemag_2vstemp_ss << "mag_2_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writemag_2vstemp_s = writemag_2vstemp_s.str();

  char * writemag_2vstemp_str = new char [writemag_2vstemp_s.length()+1];

  strcpy(writemag_2vstemp_str, writemag_2vstemp_s.c_str());

  writemag_2vstemp.open(writemag_2vstemp_str);

//-----

	ifstream writeevstemp;

  stringstream writeevstemp_ss;

  writeevstemp_ss << "e_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writeevstemp_s = writeevstemp_s.str();

  char * writeevstemp_str = new char [writeevstemp_s.length()+1];

  strcpy(writeevstemp_str, writeevstemp_s.c_str());

  writeevstemp.open(writeevstemp_str);

//-----

	ifstream writebindercumulantvstemp;

  stringstream writebindercumulantvstemp_ss;

  writebindercumulantvstemp_ss << "binder_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writebindercumulantvstemp_s = writebindercumulantvstemp_s.str();

  char * writebindercumulantvstemp_str = new char [writebindercumulantvstemp_s.length()+1];

  strcpy(writebindercumulantvstemp_str, writebindercumulantvstemp_s.c_str());

  writebindercumulantvstemp.open(writebindercumulantvstemp_str);

//-----

	ifstream writeq;

  stringstream writeq_ss;

  writeq_ss << "q_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writeq_s = writeq_s.str();

  char * writeq_str = new char [writeq_s.length()+1];

  strcpy(writeq_str, writeq_s.c_str());

  writeq.open(writeq_str);

//-----

	ifstream writeq_2;

  stringstream writeq_2_ss;

  writeq_2_ss << "q_2_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writeq_2_s = writeq_2_s.str();

  char * writeq_2_str = new char [writeq_2_s.length()+1];

  strcpy(writeq_2_str, writeq_2_s.c_str());

  writeq_2.open(writeq_2_str);

//-----

	ifstream writemx;

  stringstream writemx_ss;

  writemx_ss << "mx_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writemx_s = writemx_s.str();

  char * writemx_str = new char [writemx_s.length()+1];

  strcpy(writemx_str, writemx_s.c_str());

  writemx.open(writemx_str);

//-----

	ifstream writedmx_dh;

  stringstream writedmx_dh_ss;

  writedmx_dh_ss << "dmx_dh_nx" << nx << "ht" + ht << "temp" << start_temp << "to" << final_temp << "_TFIMSSE.txt";

  string writedmx_dh_s = writedmx_dh_s.str();

  char * writedmx_dh_str = new char [writedmx_dh_s.length()+1];

  strcpy(writedmx_dh_str, writedmx_dh_s.c_str());

  writedmx_dh.open(writedmx_dh_str);

//-----

	//@ Initialize the write-data-to-file.
	//@ Name and begin the write-data-to-file.
	
	for (int y=0; y<num_temp_step; y++) {

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

		MTRand_open rand_num(time_seed);

		//@ random number generator (0,1). Note: omitted RAN(), INITRANDOM, and WRITERAND(rndf)

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



//---------------------------

//initconf

//---------------------------



		/*@ Important: setting variables that need to be 0 (or 0.0 if a double) to that between temperatures!!
		setting arrays that need to be 0 (or 0.0 if a double) to that between temperatures!!*/

		/*@ Generates the initial random spin state, and sets initial values
		for l (expansion truncation) and nh (current expansion order n). */

		int cc;
		
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

			cc = min((int(rand_num()*nx) + 1), nx);
		
			while (true) {

				if (spn[c] == -1) {

					break;

				}

				cc = min ((int(rand_num()*nx) + 1), nx);

			}

			// the result of this is to randomly select half of the spins and flip them to +1, seems like an overly complicated way to initialize the spins, but whatever

			spn[c] = 1;

		}



//---------------------------

//lattice

//---------------------------



			/*@Creates maps with coordinates for the lattice sites "(xyz, xyzi)", 
		and vectors disx,"disy,disz" for distances that take into account the 
		periodic boundary conditions (gives the shorter of the two possible
		distances).*/

		// many of these data structures may no longer be necessary/useful in 1D code, but possibly useful to hold on to for higher dimensional generalization	
	
		int z = 0;

		for (int ix=1; ix<(nx+1); ix++) {

			z = z + 1;		

			x1[z] = ix;

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



//---------------------------

//pvectors

//---------------------------



		for (int i=0; i<l; i++ ) {

			ap1[i] = ht * beta * double(nx) / (double(l - i));

		} 

		for (int j=1; j<(l+1); j++) {

			dp1[j] = double(l - j + 1) / (ht * beta * double(nx));

		}



//---------------------------

//isingpart

//---------------------------



		int ix1; //@note ix used for inside below loop as well
	
		double r1, r2, p1, p2;

		for (ix1=0; ix1<(nx/2 + 1); ix1++) { /*$ what do these exactly accomplish*/

			if (ix1==0) {

				ris[ix1] = ht;

				lis[2][2][ix1] = true;

				lis[2][0][ix1] = false;

				lis[0][2][ix1] = false;

				lis[0][0][ix1] = true;

			} else {

				r1 = double(ix1);

				r2 = double(nx - ix1);

				ris[ix1] = (-1.0)/pow(r1, power) - 1.0/pow(r2, power); //$~ TWEAK THIS FOR PROPER "DISTANCES OR THIS THING LOOK AT THOSE r^2's"

				if (ris[ix1]<=0) {

					lis[2][2][ix1] = true;

					lis[2][0][ix1] = false;

					lis[0][2][ix1] = false;

					lis[0][0][ix1] = true;

				} else {

					lis[2][2][ix1] = false;

					lis[2][0][ix1] = true;

					lis[0][2][ix1] = true;

					lis[0][0][ix1] = false;

				}

				ris[ix1] = abs(ris[ix1]);

			}

		}

		pint[0] = 0.0;

		pint[1] = ris[0];

		eshft = 0.0;

		for (int i=2; i<(nx+1); i++) {

			ix1 = disx[1][ x1[i] ]; /*@ to be clear, this is finding the distance between the first site (site x1[0]) and every other site (site x1[i]), note: 1 to nx*/

			pint[i] = pint[i-1] + ris[ix1];

			eshft = eshft + (ris[ix1]/2);

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



//---------------------------

//zerodat

//---------------------------



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

//---------------------------

		for (long int u=1; u<(isteps+1); u++) {

			//---------------------------

			//diaupdate

			//---------------------------



				/*@ [-1,i] = flips i
						[0,0]  = 1
						[i,i]  = h
						[i,j]  = Ising operator */

				int s1, s2, p0, p3, p4, ix2, nac1, nac2, ntr2;

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

								p3 = pfrst[p0];

								p4 = plast[p0];

							} else {

							failed_all_ifs = false;

							//@ for last branch with an if, else statement

							}

							if (p3==p4) {

								s2 = ((s1 + p3 - 2) % nx) + 1;

								ix2 = disx[ x1[s1] ][ x1[s2] ];

								if ((lis[spn[s1]+1][spn[s2]+1][ix2])) {

									strb[i] = s2;

									nac2 = nac2 + 1;

									break;

								} else {

									continue;

								}

							} else if (p4==(p3+1)) {

								if (pint[p3]>=p) {

									s2 = ((s1 + p3 - 2) % nx) + 1;

								} else {

									s2 = ((s1 + p4 - 2) % nx) + 1;

								}

								ix2 = disx[ x1[s1] ][ x1[s2] ];		

								if (lis[ spn[s1]+1 ][ spn[s2]+1 ][ix2]) {

									strb[i] = s2;

									nac2 = nac2 + 1;

									break;

								} else {

									continue;

								}

							} else {

								p0 = (p3 + p4)/2;

								if ((pint[p0-1]<p)&&(pint[p0]>=p)) {

									s2 = ((s1 + p0 - 2) % nx) + 1;

									ix2 = disx[ x1[s1] ][ x1[s2] ];

									if (lis[ spn[s1]+1 ][ spn[s2]+1 ][ix2]) {

										strb[i] = s2;

										nac2 = nac2 + 1;

										break;

									} else {

										continue;

									}

								} else {

									if ((pint[p0])>=p) {

										p4 = p0;

									} else {

										p3 = p0;

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



			//---------------------------

			//partition

			//---------------------------



				int s3, s4;
	
				for (int i=1; i<(nx+1); i++) {

					lsub[i] =0;

					con1[0][i] = false;

				}

				for (int j=1; j<(l+1); j++) {

					s3 = stra[j];

					if (s3!=0) {

						s4 = strb[j];

						if (s3==s4) {

							lsub[s3] = lsub[s3] + 1;

							pos3[ lsub[s3] ][s3] = j;

							str1[ lsub[s3] ][s3] = 1;

							con1[ lsub[s3] ][s3] = false;

						} else if (s3==(-1)) {

							lsub[s4] = lsub[s4] + 1;

							pos3[ lsub[s4] ][s4] = j;

							str1[ lsub[s4] ][s4] = 2;

							con1[ lsub[s4] ][s4] = false;

						} else {

							con1[ lsub[s3] ][s3] = true;

							con1[ lsub[s4] ][s4] = true;

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



//---------------------------

//offupdate

//---------------------------



				for (int g=1; g<(nx+1); g++) {

						int ii, p5, p6, lens, flp, top, nac, nupd;

	nupd = 0;

	lens = lsub[g];

	for (int i=1; i<(lens+1); i++) {

		if ((!con1[i][g])) {

			nupd = nupd + 1;

			pos[nupd] = i;

		}

		opr[i] = str1[i][g];

	}

	if (lens<2) {

		if ((nupd == lens)&&((!con1[0][g]))) { //$ ?

			if (rand_num()<0.75) {

				spn[g] = -spn[g];

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

		p5 = min(int(rand_num()*nupd)+1, nupd);

		p5 = pos[p5];

		p6 = (p5 % lens) + 1;

		if (opr[p5]==opr[p6]) {

			opr[p5] = 3 - opr[p5];

			opr[p6] = 3 - opr[p6];

			nac = nac + 1;

		} else {

			top = opr[p5];

			opr[p5] = opr[p6];

			opr[p6] = top; //@ why even have top if this is all its used for?

		}
		
		if (p6<p5) {

		 flp = -flp;

		}

	}

	spn[g] = spn[g] * flp;

	for (int k=1; k<(lens+1); k++) {

		ii = pos1[k][g];

		if (opr[k]==1) {

			stra[ii] = strb[ii];

		} else {

			stra[ii] = -1;

		}

	}

	ar3 = ar3 + (long double)(nac);

	if ((nupd == lens)&&((!con1[0][g]))) {

		if (rand_num()<0.75) {

			spn[g] = -spn[g];

			ar4 = ar4 + 1.0;

		}

	}

		}



//---------------------------

//checkl

//---------------------------

	int o, dl, l1;

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

			o = min(int(rand_num()*l1) + 1, l1);

			if (lstr[o]) {

				lstr[o] = false;

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



//---------------------------

//pvectors

//---------------------------



		for (int i=0; i<l; i++ ) {

			ap1[i] = ht * beta * double(nx) / (double(l - i));

		} 

		for (int j=1; j<(l+1); j++) {

			dp1[j] = double(l - j + 1) / (ht * beta * double(nx));

		}



//---------------------------

//back to checkl

//---------------------------



	myfile1 << u << " increased l to " << l << "\n\n";

	//@ cout << "\n\n step " << step << " l is " << l << " ";  



//---------------------------



			//@ if (temp==avn1ctemp){

				//@ calcenrcool(i, writeavn1c_cool);

			//@ }
			
			if ((u%(isteps/10))==0) { //@ to record every 10 steps

				writelog << "Done equilibriation step: " << u << "\n\n";	



//---------------------------

//errcheck

//---------------------------



	int s5, s6, ixa;

	for (int i=1; i<(nx+1); i++) {

		spn1[i] = spn[i];

	}
	
	for (int j=1; j<(l+1); j++) {

		s5 = stra[j];

		s6 = strb[j];

		if ((s5<-1)||(s5>nx)) {

 			cout << "Illegal s5 operator on an errcheck() j = " << j << " where s5 = " << s5 << "\n\n";

			return;

		}
		
		if ((s6<0)||(s6>nx)) {

			cout << "Illegal s6 operator on an errcheck() j = " << j << " where s6 = " << s6 << "\n\n";

			return;

		}

		if (s5>0) {

			ixa = disx[ x1[s5] ][ x1[s6] ];			

			if (!lis[ spn1[s5]+1 ][ spn1[s6]+1 ][ixa]) {

				cout << "Illegal Ising bond on an errcheck() j = " << j << " where Ising bond = " << lis[ spn[s5]+1 ][ spn[s2]+1 ][ix] << ", ix = " << ix << "\n\n";

				return;

			}

		} else if (s5==-1) {

			spn1[s6] = -spn1[s6];

		}

	}

	for (int k=1; k<(nx+1); k++) {

		if (spn[k]!=spn1[k]) {

			cout << "Wrong state propagion, spin index = " << k << "\n\n"; /*@ Unless he meant the french word, propagion is probably a typo for propagation*/

			return;

		}

	}



//---------------------------



			}

		}

		//@ writeconf(writeconfig);

		ar1 = 0.0;
		
		ar2 = 0.0;

		ar3 = 0.0;

		ar4 = 0.0;

		for (long int u=1; u<(msteps+1); u++) {

			//---------------------------

			//diaupdate

			//---------------------------



				/*@ [-1,i] = flips i
						[0,0]  = 1
						[i,i]  = h
						[i,j]  = Ising operator */

				int s1, s2, p0, p3, p4, ix2, nac1, nac2, ntr2;

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

								p3 = pfrst[p0];

								p4 = plast[p0];

							} else {

							failed_all_ifs = false;

							//@ for last branch with an if, else statement

							}

							if (p3==p4) {

								s2 = ((s1 + p3 - 2) % nx) + 1;

								ix2 = disx[ x1[s1] ][ x1[s2] ];

								if ((lis[spn[s1]+1][spn[s2]+1][ix2])) {

									strb[i] = s2;

									nac2 = nac2 + 1;

									break;

								} else {

									continue;

								}

							} else if (p4==(p3+1)) {

								if (pint[p3]>=p) {

									s2 = ((s1 + p3 - 2) % nx) + 1;

								} else {

									s2 = ((s1 + p4 - 2) % nx) + 1;

								}

								ix2 = disx[ x1[s1] ][ x1[s2] ];		

								if (lis[ spn[s1]+1 ][ spn[s2]+1 ][ix2]) {

									strb[i] = s2;

									nac2 = nac2 + 1;

									break;

								} else {

									continue;

								}

							} else {

								p0 = (p3 + p4)/2;

								if ((pint[p0-1]<p)&&(pint[p0]>=p)) {

									s2 = ((s1 + p0 - 2) % nx) + 1;

									ix2 = disx[ x1[s1] ][ x1[s2] ];

									if (lis[ spn[s1]+1 ][ spn[s2]+1 ][ix2]) {

										strb[i] = s2;

										nac2 = nac2 + 1;

										break;

									} else {

										continue;

									}

								} else {

									if ((pint[p0])>=p) {

										p4 = p0;

									} else {

										p3 = p0;

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



			//---------------------------

			//partition

			//---------------------------



				int s3, s4;
	
				for (int i=1; i<(nx+1); i++) {

					lsub[i] =0;

					con1[0][i] = false;

				}

				for (int j=1; j<(l+1); j++) {

					s3 = stra[j];

					if (s3!=0) {

						s4 = strb[j];

						if (s3==s4) {

							lsub[s3] = lsub[s3] + 1;

							pos3[ lsub[s3] ][s3] = j;

							str1[ lsub[s3] ][s3] = 1;

							con1[ lsub[s3] ][s3] = false;

						} else if (s3==(-1)) {

							lsub[s4] = lsub[s4] + 1;

							pos3[ lsub[s4] ][s4] = j;

							str1[ lsub[s4] ][s4] = 2;

							con1[ lsub[s4] ][s4] = false;

						} else {

							con1[ lsub[s3] ][s3] = true;

							con1[ lsub[s4] ][s4] = true;

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



//---------------------------

//offupdate

//---------------------------



				for (int g=1; g<(nx+1); g++) {

						int ii, p5, p6, lens, flp, top, nac, nupd;

	nupd = 0;

	lens = lsub[g];

	for (int i=1; i<(lens+1); i++) {

		if ((!con1[i][g])) {

			nupd = nupd + 1;

			pos[nupd] = i;

		}

		opr[i] = str1[i][g];

	}

	if (lens<2) {

		if ((nupd == lens)&&((!con1[0][g]))) { //$ ?

			if (rand_num()<0.75) {

				spn[g] = -spn[g];

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

		p5 = min(int(rand_num()*nupd)+1, nupd);

		p5 = pos[p5];

		p6 = (p5 % lens) + 1;

		if (opr[p5]==opr[p6]) {

			opr[p5] = 3 - opr[p5];

			opr[p6] = 3 - opr[p6];

			nac = nac + 1;

		} else {

			top = opr[p5];

			opr[p5] = opr[p6];

			opr[p6] = top; //@ why even have top if this is all its used for?

		}
		
		if (p6<p5) {

		 flp = -flp;

		}

	}

	spn[g] = spn[g] * flp;

	for (int k=1; k<(lens+1); k++) {

		ii = pos1[k][g];

		if (opr[k]==1) {

			stra[ii] = strb[ii];

		} else {

			stra[ii] = -1;

		}

	}

	ar3 = ar3 + (long double)(nac);

	if ((nupd == lens)&&((!con1[0][g]))) {

		if (rand_num()<0.75) {

			spn[g] = -spn[g];

			ar4 = ar4 + 1.0;

		}

	}

		}



//---------------------------



			if ((u%2)==0) {



//---------------------------

//measure fxn

//---------------------------



	int s7, s8, mu, nu, ni, nt, nh1, ssum, last;
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

		s7 = stra[j];

		if (s7!=0) {

			nh1 = nh1 + 1; //@ Hij where i=/=0, non identity

			s8 = strb[j];

			if (s7==-1) {

				nt = nt + 1; //@ H(-1)j, "spin flip" with flip and h strength?

				ssum = ssum + (nh1 - last)*mu;

				last = nh1;

				spn[s8] = -spn[s8];

				mu = mu + 2*spn[s8];

			} else if (s7==s8) {

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

//--------------------------

			}



			if ((u%(msteps/10))==0) { //@ to record every 10 steps

				writelog << "Done measurement step " << u << "\n\n";

//---------------------------

//errcheck

//---------------------------



	int s5, s6, ixa;

	for (int i=1; i<(nx+1); i++) {

		spn1[i] = spn[i];

	}
	
	for (int j=1; j<(l+1); j++) {

		s5 = stra[j];

		s6 = strb[j];

		if ((s5<-1)||(s5>nx)) {

 			cout << "Illegal s5 operator on an errcheck() j = " << j << " where s5 = " << s5 << "\n\n";

			return;

		}
		
		if ((s6<0)||(s6>nx)) {

			cout << "Illegal s6 operator on an errcheck() j = " << j << " where s6 = " << s6 << "\n\n";

			return;

		}

		if (s5>0) {

			ixa = disx[ x1[s5] ][ x1[s6] ];			

			if (!lis[ spn1[s5]+1 ][ spn1[s6]+1 ][ixa]) {

				cout << "Illegal Ising bond on an errcheck() j = " << j << " where Ising bond = " << lis[ spn[s5]+1 ][ spn[s2]+1 ][ix] << ", ix = " << ix << "\n\n";

				return;

			}

		} else if (s5==-1) {

			spn1[s6] = -spn1[s6];

		}

	}

	for (int k=1; k<(nx+1); k++) {

		if (spn[k]!=spn1[k]) {

			cout << "Wrong state propagion, spin index = " << k << "\n\n"; /*@ Unless he meant the french word, propagion is probably a typo for propagation*/

			return;

		}

	}



//---------------------------

			}

		}



//---------------------------

//results

//---------------------------



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

//---------------------------

//zerodat

//---------------------------




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




//--------------------------- 

//back to results

//--------------------------- 


//--------------------------- 

//--------------------------- 

//writeacc

//---------------------------



	ar1 = ar1/((long double)(msteps));

	ar2 = ar2/((long double)(msteps));

	ar3 = ar3/((long double)(msteps));

	ar4 = ar4/((long double)(msteps));

	myfile5 << "Diagonal acceptance 1     : " << ar1 << "\n"; 

	myfile5 << "Diagonal acceptance 2     : " << ar2 << "\n"; 

	myfile5 << "Off-Diagonal substitions  : " << ar3 << "\n"; 

	myfile5 << "Spin flips / site         : " << ar4 << "\n"; 

	myfile5 << "Max lsub                  : " << mlls << "\n\n";



//---------------------------

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
