#include <iostream>
#include <fstream>
#include <cmath>

#include <MT/mtrand.h>


int isteps;

int msteps;

int nd;

// number of runs

double beta;

// inverse temperature

double ht;

// transverse field strength

void mcstep();

void checkl(int step);



int main() {

 
  for (int i=0; i<isteps; i++) {

    mcstep();

    checkl(i);



  }
  











  return 0;


}
