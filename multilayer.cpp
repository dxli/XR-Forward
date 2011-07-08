#include<iostream>
#include "multilayer.h"
using namespace std;

zpoint::zpoint(double zk0,complex<double> nk0)
{
    zk=zk0;
    nk=nk0;
}
void multilayer::tozprofile(vector<complex<double> > nbulki,vector<vector<double> > zpi)
//sum up zprofiles into nk
{
        zpi_save=zpi;
    for (vector<vector<double> >::iterator pf=zpi.begin();pf != zpi.end();pf++)
    {
        complex<double> nk0(0.,0.);
        for (unsigned int i=1;i<pf->size();i++)
            nk0 += nbulki.at(i-1)*pf->at(i);
        zprofile.push_back(zpoint(pf->at(0),nk0));
    }
    for (vector<zpoint>::iterator pf=zprofile.begin()+1;pf+1 != zprofile.end();pf++) {
        pf->aki=complex<double>(0.,2.)*k0*((pf+1)->zk - pf->zk);
        //      cout<<pf->zk<<' '<<pf->nk.real()<<endl;
    }
}

void multilayer::tozprofile(vector<complex<double> > nbulki,vector<vector<double> > zpi,vector<compound> c0)
{
        cpmd0=c0;
        tozprofile(nbulki,zpi);
}

multilayer::multilayer( double l)
{
    lambda=l;
    k0 = 2.*M_PI/lambda;
    convoluted=0;
}

double multilayer::fresnel(double q)
//fresnel reflectivity
{
    double t=q/(4.*M_PI)*lambda;
//    double cosk0=t*t;
//   cosk0 =  cosk0*(1./2.-cosk0*(1./24.-1./720.*cosk0)); // 1 - cos(t)
    //double cosk0=1-a;
    //c=sqrt(a+a);
    //double sink0=sqrt(cosk0*(2.0-cosk0));
    t *=t;
    complex<double> sink02 =zprofile.front().nk+t;
    complex<double> sink12 =sink02 - zprofile.back().nk;
    complex<double> sink0=sqrt(sink02);
    complex<double>  sink1=sqrt(sink12);
    //cosk0=1-cosk0;
    //sink1=(sink1 - zprofile.back().nk)*(1.+zprofile.back().nk); // (1 -a0) /(1-nk)=(1-a0)(1+nk+nk^2)=1 -(a0 +a0nk-nk-nk^2)
    //complex<double> cosk1=(1. - sink1)/cosk0;
    //a0 = sqrt(b0+b0);
    //sink1= sqrt(sink1*(2.0-sink1))/sink0;
    /*
    complex<double> ans=zprofile.front().nk*sink0+zprofile.back().nk*sink1;
    ans=( ( zprofile.front().nk + zprofile.back().nk)*t - zprofile.back().nk*zprofile.back().nk)*(zprofile.back().nk-zprofile.front().nk)/(ans*ans);
    return(ans.real()*ans.real()+ans.imag()*ans.imag());
    */
    sink0=(sink0-sink1)/(sink0+sink1);
    return(sink0.real()*sink0.real()+sink0.imag()*sink0.imag());
}

double multilayer::born_rf(double q)
// reflectivity for a given angle, within Born Approximation
{
complex<double> ans(0.,0.),iq(0.,q);
for(unsigned int i=1;i<zprofile.size();i++){
        ans += (zprofile.at(i).nk - zprofile.at(i-1).nk)*exp( iq*zprofile.at(i).zk);
}
ans /= zprofile.back().nk;
return(ans.real()*ans.real()+ans.imag()*ans.imag());
}

double multilayer::rf(double q)
// reflectivity for a given angle
{

    // sink.at(0)=pref0->sink0;
    //cosk.at(0)=pref0->cosk0;
    double t=q/(4.*M_PI)*lambda;
    complex<double> a0=zprofile.front().nk + t*t,xs;
    complex<double> r0(0,0),sink1,sink0;
    //xs=(a0 - nk[i])*(1.+nk[i]); // cos\theta/(1-nk) =1 - ( a0 +a0 nk -nk - nk^2)
    //sink1= sqrt(xs*(2.-xs)); //sin\theta = sqrt(1-cos^2\theta)
    //sink1 -= sink1*nk[i]; // sin\theta *(1- nk)

    complex<double> aki(1,0);
    //ak.at(--i)=complex<double>(1,0);
    int i=zprofile.size()-1;
    sink1 = sqrt(a0 - zprofile.back().nk);
    while (--i>=0) {
        //xs=(a0 - nk[i])*(1.+nk[i]); // cos\theta/(1-nk) =1 - ( a0 +a0 nk -nk - nk^2)
        //sink0= sqrt(xs*(2.-xs)); //sin\theta = sqrt(1-cos^2\theta)
        //sink0 -= sink0*nk[i]; // sin\theta *(1- nk)
        sink0 = sqrt(a0 - zprofile.at(i).nk);
        xs=(sink0-sink1)/(sink0+sink1);
        sink1=sink0;
        r0 *= aki;
        r0 = (r0+xs)/(1.+xs*r0);
        aki=exp(zprofile.at(i).aki*sink1);
    }
    return( r0.real()*r0.real()+r0.imag()*r0.imag());
}

bool operator <(const zpoint & a, const zpoint & b) {
        return a.zk < b.zk;
}
void multilayer::convolution(double sigma)
// convolute the given zprofile with gaussian       
{
        convoluted=1;
        sigma0=sigma;
        double zrange=3.*sigma;
        double isigmas2=sqrt(0.5)/sigma;
        std::sort(zprofile.begin(),zprofile.end());

        double dz=0.01*sigma;
        double dz2=(zprofile.at(1).zk-zprofile.front().zk);
        for(unsigned int i=1;i<zprofile.size();i++) {
double a=zprofile.at(i).zk -zprofile.at(i-1).zk;
                if (a >0. && a<dz2) dz2=a;
        }
        dz2 *= 0.025;
        if(dz > dz2) dz=dz2;
        vector<complex<double> > nk0,dnk0;
        zprofile.push_back(zpoint(zprofile.back().zk+zrange,zprofile.back().nk));

         //   dnk0.push_back( complex<double>(0.,0.));
    for(unsigned int i=1;i<zprofile.size();i++) dnk0.push_back( 0.5*(zprofile.at(i).nk - zprofile.at(i-1).nk));
    dnk0.push_back(complex<double>(0.,0.));
vector<int> n0i,n1i;//slab numbers before and after convolution
    for(unsigned int i=1;i<zprofile.size();i++) n0i.push_back((int)( (zrange + zprofile.at(i).zk-zprofile.front().zk)/dz)-1);
    //int jmax=(int)( (zprofile.back().zk+zrange+zrange)/dz+0.5);
    n1i.resize((int)((zprofile.back().zk-zprofile.front().zk+zrange+zrange)/dz)+1);
    for(unsigned int i=0,j=0;j<n1i.size();j++) {
            if(i<n0i.size()) 
                    if(j>=n0i.at(i)) i++;
            n1i.at(j)=i;
            //cout<<j<<' '<<i<<endl;
    }
    int jmax=n1i.size();
    int dside = (int) ( zrange/dz)+1;
vector<zpoint> zp2;
int ic=0;
    for (int i=0;i<jmax;i++){
       //     cout<<i<<endl;
        ic=n1i.at(i); // center gene
        int il=i-dside;
        il=il<0?0:n1i.at(il);
        int ir=i+dside;
        if (ir>=n1i.size()) ir=n1i.size() -1;
        ir= n1i.at(ir);
        double z = (i+0.5)*dz+zprofile.front().zk-zrange ;
        complex<double> d0 = zprofile.at(il).nk;
        //std::cout<<i<<"("<<jmax<<"): "<<il<<' '<<ic<<' '<<ir<<endl;
        //std::cerr<<i<<"("<<jmax<<"): "<<il<<' '<<ic<<' '<<ir<<endl;
        for (unsigned int j=il;j<=ir;j++) {
            if(j+1<zprofile.size()) d0 += (1.+gsl_sf_erf ((z - zprofile.at(j+1).zk)*isigmas2))*dnk0.at(j);
        }
        zp2.push_back(zpoint(z,d0));
    }
    //cout<<"1\n";
    zp_raw=zprofile;
    zprofile=zp2;
   for (vector<zpoint>::iterator pf=zprofile.begin()+1;pf+1 != zprofile.end();pf++) {
        pf->aki=complex<double>(0.,2.)*k0*((pf+1)->zk - pf->zk);
        //      cout<<pf->zk<<' '<<pf->nk.real()<<endl;
    }
    //cout<<"2\n";

}


#ifndef EXIT_SUCCESS
#  define EXIT_SUCCESS 0
#endif

#ifndef EXIT_FAILURE
#  define EXIT_FAILURE -1
#endif

void my_error_function(const char *msg)
{
    fprintf(stderr, "library message: \"%s\"\n", msg);
}

void multilayer::rfvector(double q0,double q1,double dq)
        //
{
        rf_save.clear();
        q0=fabs(q0);
        q1=fabs(q1);
        if(q0>q1) {
                double a=q0;
                q0=q1;
                q1=a;
        }
        if(q0<0.001) q0=0.001;
        if(q0==q1) {
                q0=0.01;
                q1=1.;
        }
        if( !( (q1-q0)/dq < 2000.|| (q1-q0)/dq<2.)) dq = (q1-q0)/2000.;

        for(double q=q0;q<=q1;q+=dq){
                vector<double> rf0;
                rf0.push_back(q);
                if(born) {
                if(by_rF) rf0.push_back(born_rf(q));
                else rf0.push_back(fresnel(q)*born_rf(q));

                }else{
                if(by_rF) rf0.push_back(rf(q)/fresnel(q));
                else rf0.push_back(rf(q));
                }
                rf_save.push_back(rf0);
        }
        //cout<<rf_save.size()<<endl;
}
void multilayer::plot(void)
        //plotout
{
    //if (GraceOpen(3072)){
    if (GraceOpenVA("gracebat", 16384, "-nosafe", "-noask",NULL)==-1){
        cerr<<"Can't run Grace. \n";
        exit(EXIT_FAILURE);
    }
    GracePrintf("page size 612, 792");
    GracePrintf("ARRANGE(2, 1, 0.15, 0.15,0.40)");
    //plotting rho
    GracePrintf("with g1");
    GracePrintf("legend on");
    GracePrintf("legend 0.75,0.55");
    int si=0;
    double ymax=zpi_save.at(0).at(1);
    double ymin=zpi_save.at(0).at(1);
    //std::cout<<"<br>"<<zpi_save.at(0).size()<<"<br>";
for(unsigned int i=1;i<zpi_save.at(0).size();i++){
    GracePrintf("s%d on",si);
    GracePrintf("s%d type xy",si);
    GracePrintf("s%d symbol 0",si);
    GracePrintf("s%d line linestyle 1",si);
    {
            ostringstream oss;
    //std::cout<<"<br>cpmd0.size()="<<cpmd0.size()<<"<br>\n";
    //std::cout<<"<br>cpmd0="<<cpmd0.at(si).els.size()<<"<br>";
                for(unsigned int j=0;j<cpmd0.at(si).els.size();j++){
                oss<<cpmd0.at(si).els.at(j).symbol;
                if(j<cpmd0.at(si).n0.size()) 
                        if(cpmd0.at(si).n0.at(j) != 1) oss <<std::string("\\s")<<cpmd0.at(si).n0.at(j)<<"\\N";
        }
     //           std::cout<<oss.str().c_str();
    GracePrintf("s%d legend \"%s\"",si,oss.str().c_str());
    }
    //std::cout<<"<br>1<br>";
        for(unsigned int j=0;j<zpi_save.size();j++){
        if(ymin>zpi_save.at(j).at(i)) ymin=zpi_save.at(j).at(i);
        if(ymax<zpi_save.at(j).at(i)) ymax=zpi_save.at(j).at(i);
                GracePrintf("s%d point %g,%g",si,zpi_save.at(j).at(0),zpi_save.at(j).at(i));
                if(j+1<zpi_save.size()) GracePrintf("s%d point %g,%g",si,zpi_save.at(j+1).at(0),zpi_save.at(j).at(i));
        }
                GracePrintf("s%d point %g,%g",si,2*zpi_save.back().at(0)-zpi_save.at(zpi_save.size()-2).at(0),zpi_save.back().at(i));
                GracePrintf("s%d point %g,%g",si,3*zpi_save.back().at(0)-2*zpi_save.at(zpi_save.size()-2).at(0),zpi_save.back().at(i));
                si++;
}
//convoluted
    //cout<<"3\n";
{
    GracePrintf("s%d on",si);
    GracePrintf("s%d type xy",si);
    GracePrintf("s%d symbol 0",si);
    GracePrintf("s%d line linestyle 1",si);
if(convoluted==1)
    GracePrintf("s%d legend \"\\f{Symbol}r\\f{}\\sel\\N,\\f{Symbol}s\\f{}=%g\\#{c5}\"",si,sigma0);
else
    GracePrintf("s%d legend \"\\f{Symbol}r\\f{}\\sel\"",si,sigma0);
    unsigned int j;
        for(j=0;j<zprofile.size();j++){
        if(ymin>zprofile.at(j).nk.real()/zprofile.back().nk.real()) ymin=zprofile.at(j).nk.real()/zprofile.back().nk.real();
        if(ymax<zprofile.at(j).nk.real()/zprofile.back().nk.real()) ymax=zprofile.at(j).nk.real()/zprofile.back().nk.real();
                GracePrintf("s%d point %g,%g",si,zprofile.at(j).zk,zprofile.at(j).nk.real()/zprofile.back().nk.real());
               if(j+1<zprofile.size()) GracePrintf("s%d point %g,%g",si,zprofile.at(j+1).zk,zprofile.at(j).nk.real()/zprofile.back().nk.real());
        }
                GracePrintf("s%d point %g,%g",si,2*zprofile.back().zk-zprofile.at(j-2).zk,1.);
                GracePrintf("s%d point %g,%g",si,3*zprofile.back().zk-2.*zprofile.at(j-2).zk,1.);
                si++;
}
   GracePrintf("xaxis  ticklabel char size 1.5");
    GracePrintf("yaxis  ticklabel char size 1.5");

    GracePrintf("autoscale");
    GracePrintf("world ymax %g",ymax*1.25);
    GracePrintf("autoticks");
        GracePrintf("xaxis label \"\\+\\+z (\\f{Times-Roman}%c\\f{})\"",(unsigned char) 197);
    GracePrintf("yaxis label \"\\+\\+\\f{Symbol} r\\f{}(z)/\\f{Symbol}r\\f{}(\\f{Symbol}%c\\f{})\"",(unsigned char) 0xa5);

   GracePrintf("redraw");
    GraceFlush();
    GracePrintf("with g0");
    GracePrintf("YAXES SCALE LOGARITHMIC");
   GracePrintf("xaxis  ticklabel char size 1.5");
    GracePrintf("yaxis  ticklabel char size 1.5");
    GracePrintf("yaxis  tick major 10");
    GracePrintf("yaxis  tick minor 9");
    GracePrintf("yaxis  ticklabel format power");
    GracePrintf("yaxis  ticklabel prec 0");

    
    si=0;
    GracePrintf("s%d on",si);
    GracePrintf("s%d type xy",si);
    GracePrintf("s%d symbol 0",si);
    GracePrintf("s%d line linestyle 1",si);
    ymax=rf_save.at(0).at(1);
    ymin=rf_save.at(0).at(1);
for(unsigned int j=0;j<rf_save.size();j++) {
        if(ymin>rf_save.at(j).at(1)) ymin=rf_save.at(j).at(1);
        if(ymax<rf_save.at(j).at(1)) ymax=rf_save.at(j).at(1);
                GracePrintf("s%d point %g,%g",si,rf_save.at(j).at(0),rf_save.at(j).at(1));
}
if( outputdata){
        for(unsigned int j=0;j<rf_save.size();j++) {
                std::cout<<rf_save.at(j).at(0)<<' '<<rf_save.at(j).at(1)<<endl;
}
   GraceFlush();
   if(GraceIsOpen()) GraceClose();
return;
}
    GracePrintf("autoscale");
    //cout<<ymin<<" - "<<ymax<<endl;
    GracePrintf("world ymin %g",ymin*0.75);
    GracePrintf("world ymax %g",ymax*1.25);
    GracePrintf("autoticks");
    if(born) GracePrintf("subtitle \"Born Approximation\"");
    else GracePrintf("subtitle \"Parrett\"");
  GracePrintf("xaxis label \"\\+\\+q (%c\\S-1\\N\\+\\+)\"",(unsigned char) 197);
    if(by_rF) GracePrintf("yaxis label \"\\+\\+R(q\\sz\\N\\+\\+)/R\\sF\\N\"");
    else GracePrintf("yaxis label \"\\+\\+R(q\\sz\\N\\+\\+)\"");
   GracePrintf("redraw");
   pid_t pid0=getpid();
   ostringstream sfn0;
   sfn0<<pid0;
   std::string fn0=std::string("/tmp/")+sfn0.str();
   std::string fnin=fn0+std::string(".eps");
   std::string fnout=fn0+std::string(".png");
   GracePrintf("HARDCOPY DEVICE \"EPS\"");
   GracePrintf("print to \"%s\"",fnin.c_str());
   GracePrintf("print");
   GracePrintf("exit");
   GraceFlush();
   if(GraceIsOpen()) GraceClose();
    waitpid(-1,NULL,0);
    pid0=fork();
    if(pid0 == 0) execl("/usr/bin/convert","convert", fnin.c_str(),fnout.c_str(),NULL);
    waitpid(pid0,NULL,0);
    struct stat results;
    stat(fnout.c_str(),&results);
    //std::cout<<"size1="<<results.st_size;
    char buffer[results.st_size+2];
    ifstream in1(fnout.c_str(),ios::in | ios::binary);
    in1.read(buffer,results.st_size);
    in1.close();
    unlink(fnin.c_str());
    unlink(fnout.c_str());
    std::cout.write(buffer,results.st_size);
}
