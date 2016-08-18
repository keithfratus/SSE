#include <string.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <iomanip>
#include <cmath>

#include "mtrand.h"

using namespace std;

/*sizes NEW LAPTOP

char: 1
char16_t: 2
char32_t: 4
wchar_t: 4
signed char: 1
signed short int: 2
signed int: 4
signed long int: 8
signed long long int: 8
unsigned char: 1
unsigned short int: 2
unsigned int: 4
unsigned long int: 8
unsigned long long int: 8
float: 4
double: 8
long double: 16
bool: 1

*/

//@ Manipulate these 12 variables below:

bool OBCtrue_PBCfalse = false; //@ Toggle true for Open Boundary Conditions or false for Periodic Boundary Conditions

const double start_temp = 0.5;

const double temp_step = 0.01;

const int num_temp_step = 500; //@ Works as [0,num_temp_step) so first step is 0 start_temp, last is num_temp_step - 1

const double final_temp = start_temp + temp_step*(num_temp_step - 1); //@ For when writing final to get rid of comma

const long int num_steps = 2000000;

// const int num_columns = 15;

// const int num_rows = 15;

const long long int num_sites = 16;

const int skip_step = 2;

const long long int cut_off_step = 100000;

int spin_config[num_sites] = {0};

int hor_bond[num_sites] = {0};

// int ver_bond[num_sites] = {0};

const double catchingcool = 1.0;

int calculateEnergy();

int calculateEnergyDifference(int site);

void monteCarloStep(double beta, int& configEn, long int& curMag, MTRand_open& rand_site);

void calcenrcool(int step, ofstream& myfile2, int configEn);



//@----------------------------------------------------------------


//@ Main function


//@----------------------------------------------------------------



int main() {

    ofstream myfile1; //@ Writing file

    ofstream myfile2; //@ Writing file

    ofstream myfile3; //@ Writing file

    if (OBCtrue_PBCfalse){

        myfile1.open ("1DTempVsAvMag2OBC.txt");

        myfile2.open ("1DTempis1CDVsEnerg.txt");

        myfile3.open ("1DTempVsAvEnerg.txt");

    } else {

        myfile1.open ("1DTempVsAvMag2PBC.txt");

        myfile2.open ("1DTempis1CDVsAvEnerg.txt");

        myfile3.open ("1DTempVsAvEnerg.txt");

    }

    // loop over different beta values

    for (int a=0; a<num_temp_step; a++) {

        double temp = start_temp + a*temp_step;

        double beta = 1.0/temp;

        cout << "Now processing T = " << temp << ", beta = " << beta << "\n\n";

        myfile1 << temp << " ";

        // myfile2 << temp << " ";

        myfile3 << temp << " ";

        long int mag_sum = 0;

        long int energy_sum = 0;

        //unsigned long int energy_sum_pow2 = 0;

        unsigned long int mag_sum_pow2 = 0;

        //unsigned long int mag_sum_pow4 = 0;

        unsigned long int num_good = 0;

        unsigned long time_seed = time(NULL);

		    // now we initialize a random spin configuration

        MTRand_open site_gen(time_seed);

        for (int b=0; b<num_sites; b++) {

		      if (site_gen() < 0.5) {

						spin_config[b] = 0;

		      } else {

		      	spin_config[b] = 1;

		      }

        }

				// compute energy of initial spin configuration

        int configEn = calculateEnergy();

				// compute mag of initial spin configuration

        long int curMag = 0;

        for (int c=0; c<num_sites; c++) {

            if (spin_config[c]==0) {

                curMag--;

            } else {

                curMag++;

            }

        }

				// now that initial spin configuration is set, we must perform MC analysis

				// start with cool down cycle

				cout << "Performing initial cool down cycle...\n\n";

				for (int d=0; d<cut_off_step; d++) {

			    monteCarloStep(beta, configEn, curMag, site_gen);

					if (temp==catchingcool) {

						calcenrcool(d, myfile2, configEn);
					
					}

				}

				// now we begin actual averaging sequence

				cout << "Performing averaging procedure...\n\n";

				for (int e=0; e<num_steps; e++) {

						monteCarloStep(beta, configEn, curMag, site_gen);

					 // sample the configuration every N steps

						if ((e%skip_step)==0) {

					 // we want to add this mag to the sum

								mag_sum += curMag;

								energy_sum += configEn;

								//energy_sum_pow2 += pow(configEn, 2);

								mag_sum_pow2 += pow(curMag, 2);

								//mag_sum_pow4 += pow(curMag, 4);

								num_good++;

						}

				}

				// now save/display the relevant information at the end of the MC analysis

				//@For specific heat

				long double avMag = ( (long double) mag_sum ) / ( (long double) num_good );

				long double avEn = ( (long double) energy_sum) / ( (long double) num_good );

				//long double avEn_pow2 = ( (long double) energy_sum_pow2) / ( (long double) num_good );

				//@Below shows Cv = (<E^2> - <E>^2)/(NT^2)

				//long double specific_heat = ( ( (long double) avEn_pow2) - ((long double) pow(avEn, 2)) ) / ( ((long double)num_sites) * ((long double) temp) * ((long double) temp) );

				//@For Binder cumulant

				//long double avMag_pow4 = ( (long double) mag_sum_pow4) / ( (long double) num_good );

				long double avMag_pow2 = ( (long double) mag_sum_pow2) / ( (long double) num_good );

				//@Below shows U(T,C) = 1 - <M^4>/(3<M^2>^2)

				//long double Binder_cumulant = ((long double) 1) - ( ((long double)avMag_pow4) / (((long double) 3) * ((long double)(pow(avMag_pow2, 2)) )) );

				cout << "The final magnetization is " << curMag << "\n";

				cout << "The average magnetization is " << avMag << "\n";

				if (temp == final_temp) {

					myfile1 << avMag_pow2;

					//myfile2 << Binder_cumulant;

					myfile3 << abs(avEn);

				} else {

					myfile1 << avMag_pow2 << "\n";

					//myfile2 << Binder_cumulant << "\n";

					myfile3 << abs(avEn) << "\n";

				}

			} // end of temp loop

			myfile1.close();
			myfile2.close();
			myfile3.close();
			return 0;


}



//@----------------------------------------------------------------


//@      																		 helper functions below


//@----------------------------------------------------------------



int calculateEnergy() {

	int energy = 0;
	
	int right_site;

	for (int i = 0; i < num_sites; ++i) {

		right_site = i + 1;

		if (right_site==num_sites) {

			right_site = 0;

		}

		if (spin_config[i]==spin_config[right_site]) {

			hor_bond[i] = -1;

		} else {

			hor_bond[i] = 1;
		}
	
		energy += hor_bond[i];
	
	}

	return energy;


}



//@----------------------------------------------------------------



int calculateEnergyDifference(int site) {

	int en_diff = 0;

	int left_site = site - 1;
	
	if (left_site==(-1)) {

		left_site = num_sites - 1;

	}

	en_diff = (-2)*(hor_bond[site] + hor_bond[left_site]);
	
	return en_diff;


}



//@----------------------------------------------------------------



void monteCarloStep(double beta, int& configEn, long int& curMag, MTRand_open& rand_site){

    long double new_site_seed = rand_site();

    new_site_seed = num_sites*new_site_seed; //@ We can do this because it's multiplying some number between 0-1

    int new_site = (int) new_site_seed;

    //@ Let's compute the energy of the configuration if we were to flip the above "randomly" chosen spin

    int en_diff = calculateEnergyDifference(new_site);

    //@ now decide whether to accept this proposed state

    bool accept = false;

    if (en_diff < 0) {

        accept = true;

    } else {

        double therm_rand = rand_site();

        double exp_en = exp( ( (-1.0)*beta*en_diff));

        //@ Instead of changing J strength, we could manipulate the above line (multiply inside exp), otherwise assume J=1

    	if (exp_en > therm_rand) {

        	accept = true;

    	}

    }

    if (accept) {

    	//@ actually flip the spin configuration
			
		  spin_config[new_site] = ( (spin_config[new_site]) + 1 )%2;

		  // update energy

		  configEn += en_diff;

			int left_site = new_site - 1;
	
			if (left_site==(-1)) {

				left_site = num_sites - 1;

			}

			hor_bond[new_site] = hor_bond[new_site]*(-1);

			hor_bond[left_site] = hor_bond[left_site]*(-1);

		  // we also must update magnetization

		  if ( spin_config[new_site] == 0 ) {

		  // spin was 1, and is now 0, spin has decreased by two

		      curMag = curMag-2;

		  } else {

		      curMag = curMag+2;

		  }

    }

		//@ if wasn't accepted, energy is not updated, and neither is mag

    // Do NOT start summing the mag yet!!!

		// Do NOT start summing the mag yet!!!

}



//@----------------------------------------------------------------



void calcenrcool(int step, ofstream& myfile2, int configEn) {

	myfile2 << step << " " << configEn << "\n";

}



//@----------------------------------------------------------------
