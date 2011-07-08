#pragma once
#include <vector>
#include <complex>
#include <algorithm>
#include <fstream>
#include <istream>
#include <iterator>
#include <pthread.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <omp.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <unistd.h>

#include "defines.h"
#include "elements.h"
#include <grace_np.h>

using namespace std;
class zpoint 
//zk, nk
{
        public:
               double zk;
               complex<double> nk,aki;
               zpoint(double, complex<double>);
};


class multilayer
{
    //mkdensity parameter
    //double gene_ra,gene_rb,gene_rc,gene_ba,gene_bb,gene_bc;
public:
        int convoluted;
    vector<zpoint> zprofile,zp_raw;
    vector<compound> cpmd0;
    vector<vector<double> > zpi_save,rf_save;
    void tozprofile(vector<complex<double> >,vector<vector<double> >);//read formulae and zprofiles, convert to nk;
    void tozprofile(vector<complex<double> >,vector<vector<double> >,vector<compound>);//read formulae and zprofiles, convert to nk;
    double lambda,k0;
    double sigma0;
    int by_rF,born,outputdata;
    void convolution(double);
    multilayer(double );
    double rf(double);
    double born_rf(double);
    void rfvector(double,double,double);
    double fresnel(double);
    void plot(void);
};

