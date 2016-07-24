#include <iostream>
#include <fstream>
#include "math.h"

#include <MT/mtrand.h>


const int lx = 5;

const int ly = 5;

const int nn = 25;

// number of sites, equals lx times ly

const int nb = 50;

// number of bonds, equal to twice the number of sites

int nh;

// number of non-identity operators

int mm = 20;

// cut-off for Taylor expansion, set to initial value above

const int nbins;

const int msteps;

const int isteps = 10000;

// number of cool-down steps

const double beta;

double aprob;

double dprob;




// -----------------------------------------------------------------------------


int spin[nn] = {0};

int bsites[2][nb];

int * opstring;

int frstspinop[nn];

int lastspinop[nn];

int * vertexlist;

int cutofflengths[isteps] = {0};


void makelattice();

void initconfig(MTRand_open& rand_gen);

void diagonalupdate(MTRand_open& rand_gen);

void loopupdate(MTRand_open& rand_gen);

void adjustcutoff(int step);

void measureobservables();

void writeresults(int msteps, int bins);


int main() {

  
  unsigned long time_seed = time(NULL);

  MTRand_open site_gen(time_seed);



  aprob = 0.5*beta*nb;
  
  dprob = 1.0/aprob;


  makelattice();

  initconfig(site_gen);


  for (int i=0; i < isteps; i++) {

    diagonalupdate(site_gen);

    loopupdate(site_gen);

    adjustcutoff(i);

  }


  cout << "Finished equilibration, M = " << mm << endl;


  for (int j=0; j<nbins; j++) {

    for (int i=0; i<msteps; i++) {

      diagonalupdate(site_gen);

      loopupdate(site_gen);

      measureobservables();

    }

    writeresults(msteps,j);

  }

  return 0;

}



// ------------------------------------------------------



void makelattice() {

  int s;

  int x2;

  int y2;

  for (int y1=0; y1<ly; y1++) {
    
    for (int x1=0; x1<lx; x1++) {

      s = 1+x1+y1*lx;

      x2 = (x1+1)%lx;

      y2=y1;

      bsites[0][(s-1)] = s;

      bsites[1][(s-1)] = 1+x2+y2*lx;

      x2=x1;

      y2 = (y1+1)%ly;

      bsites[0][(s+nn-1)] = s;

      bsites[1][(s+nn-1)] = 1+x2+y2*lx;

    }

  }


}



// ------------------------------------------------------


void initconfig(MTRand_open& rand_gen) {



  for (int i=0; i<nn; i++) {

    if (rand_gen() < 0.5) {

      spin[i] = -1;

    } else {

      spin[i] = 1;

    }

  }


  opstring = new int[mm];

  for (int i=0; i<mm; i++) {

    opstring[i] = 0;

  }
  
  nh = 0;
  
  vertexlist = new int[(4*mm)];

}





// -----------------------------------------------------------

void diagonalupdate(MTRand_open& rand_gen) {

  int op;

  int b;

  for (int i=0; i<mm; i++) {

    op = opstring[i];

    if ( op == 0) {

      double new_bond_seed = rand_gen();

      new_bond_seed = nb*new_bond_seed;

      int new_bond = (int) new_bond_seed;

      b = new_bond+1;


      if ( spin[(bsites[0][(b-1)]-1)] != spin[(bsites[1][(b-1)]-1)] ) {

	double diff_term = ( (double) (mm-nh) );

	if ( ( aprob >= diff_term ) || ( aprob >= diff_term*rand_gen() )  ) {

	  opstring[i] = 2*b;

	  nh = nh+1;

	}

      }


    } else if ( op%2 == 0 ) {

      double p = dprob*(mm-nh+1);

      if ( ( p >= 1 ) || ( p >= rand_gen() ) ) {

	opstring[i] = 0;

	nh = nh-1;

      }

    } else {

      b = op/2;

      spin[(bsites[0][(b-1)]-1)] = -spin[(bsites[0][(b-1)]-1)];

      spin[(bsites[1][(b-1)]-1)] = -spin[(bsites[1][(b-1)]-1)];

    }

  }

}


// -----------------------------------------------------------


void loopupdate(MTRand_open& rand_gen) {

  for (int i=0; i<nn; i++) {

    frstspinop[i] = -1;

    lastspinop[i] = -1;

  }
  
  
  for (int v0=0; v0<(4*mm); v0 += 4) {

    int op = opstring[(v0/4)];

    if (op != 0) {

      int b = op/2;

      int s1 = bsites[0][(b-1)];

      int s2 = bsites[1][(b-1)];

      int v1 = lastspinop[(s1-1)];

      int v2 = lastspinop[(s2-1)];

      if ( v1 != -1 ) {

	vertexlist[v1] = v0;

	vertexlist[v0] = v1;

      } else {

	frstspinop[(s1-1)] = v0;

      }

      if ( v2 != -1 ) {

	vertexlist[v2] = v0+1;

	vertexlist[(v0+1)] = v2;

      } else {

	frstspinop[(s2-1)] = v0+1;

      }

      lastspinop[(s1-1)] = v0+2;

      lastspinop[(s2-1)] = v0+3;

    } else {

      for (int q=v0; q<(v0+4); q++) {

	vertexlist[q] = 0;

      }

    }



  }

  // finished loop over v0


  for (int s1=1; s1<(nn+1); s1++) {

    int v1 = frstspinop[(s1-1)];

    if ( v1 != -1 ) {

      int v2 = lastspinop[(s1-1)];

      vertexlist[v2] = v1;

      vertexlist[v1] = v2;

    }

  }

  // finished loop over s1

  for (int v0 = 0; v0 < (4*mm); v0 += 2) {

    if (vertexlist[v0] < 1 ) {

      continue;

    }

    int v1 = v0;

    if ( rand_gen() < 0.5 ) {

      





    }










  }












}







