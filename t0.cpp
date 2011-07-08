#include<iostream>
#include "multilayer.h"
#include "elements.h"

using namespace std;
int main()
{
        double lambda=1.239;
string s0("Au H2O 0.998");
vector<complex<double> > nbulki;
istringstream iss (s0);
vector <string > vs;
    std::copy (istream_iterator < string >(iss),
               istream_iterator < string >(), back_inserter (vs));
if(vs.size()<1) exit(0);
vector<compound> vcpmd0;
for(vector<string>::iterator pvs=vs.begin(); pvs != vs.end();pvs++){
        compound::compound cpmd0(*pvs,lambda,1.0);
        if(!cpmd0.els.size()) continue;
        if(!(cpmd0.els.size()==1 && (!cpmd0.n0.size() || ! cpmd0.n0.front()))) {
                double rho0=1.0;
                if(isdigit( (pvs+1)->at(0))) {
                        pvs++;
                        rho0=strtod(& ( pvs->at(0)),NULL);
                }
              //  cout<<"rho= "<<rho0<<endl;
                cpmd0.rho_scale(rho0);
        }
        nbulki.push_back(cpmd0.nk);
        vcpmd0.push_back(cpmd0);
}
vector<vector<double> > rz;
double z0=-5.,dz=0.5,d0=26.,d1=45;

{
        vector<double> va;
        va.push_back(0.);va.push_back(0.);va.push_back(0.);
        rz.push_back(va);
        va.clear();
        va.push_back(5.);va.push_back(0.);va.push_back(0.25);
        rz.push_back(va);
        va.clear();
        va.push_back(30.);va.push_back(1.);va.push_back(0);
        rz.push_back(va);
}
/*
for(unsigned int i=0;i<400;i++){
        vector<double> va;
        va.push_back(z0);
        double rho0=z0/d0 - 1.0;
        if(fabs(rho0)<= 1.0) {
                rho0= 0.5*M_PI/(2.*sqrt(3.))*(1.0 - rho0*rho0);  //sphere
        }else rho0=0.;
        va.push_back(rho0); //Au
        if(z0<d1-d0) va.push_back(0.);//water
        else 
        if (z0>d1) va.push_back(1.);
        else va.push_back(1.-(z0-d1)/d0*(z0-d1)/d0);
        rz.push_back(va);
        z0 += dz;
}
*/
multilayer::multilayer ml0(lambda);
ml0.tozprofile(nbulki,rz,vcpmd0);
ml0.convolution(1.0);//gaussian roughness
ml0.rfvector(0.01,1.,0.005);
ml0.plot();
return 0;
ofstream out1("rf0.txt");
for(double q=0.01;q<2.;q+=0.0005){
        double r0=ml0.rf(q)/ml0.fresnel(q);
        //cout<<q<<' '<<r0<<endl;
        out1<<q<<' '<<r0<<endl;
}
out1.close();
ml0.convolution(6.0);//gaussian roughness
out1.open("rf1.txt");
for(double q=0.01;q<2.;q+=0.0005){
        double r0=ml0.rf(q)/ml0.fresnel(q);
        //cout<<q<<' '<<r0<<endl;
        out1<<q<<' '<<r0<<endl;
}
out1.close();
out1.open("rho0.txt");
for(unsigned int i=0;i<ml0.zprofile.size();i++){
        //cout<<q<<' '<<r0<<endl;
        out1<<ml0.zprofile.at(i).zk<<' '<<ml0.zprofile.at(i).nk.real()<<endl;
}
out1.close();
out1.open("rho1.txt");
for(unsigned int i=0;i<ml0.zp_raw.size();i++){
        //cout<<q<<' '<<r0<<endl;
        out1<<ml0.zp_raw.at(i).zk<<' '<<ml0.zp_raw.at(i).nk.real()<<endl;
}
out1.close();


return 0;
}
