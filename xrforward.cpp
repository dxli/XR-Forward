// calculate x-ray absorption length, critical angles, etc.
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>
#include<cstdlib>
#include<unistd.h>
#include <sys/resource.h>
#include "multilayer.h"
#include "elements.h"

struct cgiinput {
    std::string name,field;
};

int main()
{

    /* */
    /* strcasecmp() is not supported in Windows-- use strcmpi() instead */
    //std::cout<<"Content-type: image/png\n\n";
    //std::cout<<"<html><head></head><body>\n";
    std::string s0;
    /*
    if ( !(content_length = atoi(getenv("CONTENT_LENGTH"))) ) {
        printf("Content-Type: text/plain\n\n") ;
        printf("getcgivars(): No Content-Length was sent with the POST request.\n") ;
        exit(1) ;
    }
    if (content_length > 406) content_length=406;
    */
    getline(std::cin,s0);
    std::vector<struct cgiinput> cgi0;
    struct cgiinput cgiinput0;
    //process stdin for name&field values
    for (unsigned int i=0,j=0;j<s0.size();j++) {

        if (s0.at(j) == '=' || s0.at(j) == '&' || j==s0.size() -1) {
            switch (s0.at(j)) {
            case '=':
                cgiinput0.name=s0.substr(i,j-i);
                i=j+1;
                break;
            default:
                cgiinput0.field=s0.substr(i,j-i);
                i=j+1;
                cgi0.push_back(cgiinput0);
            }
        }
    }
    double en=-1,lambda=-1,sigma=-1,qmin=0.01,qmax=1.,dq=0.0025;
    std::string f0,r0;
    vector<std::string> buffer0;
    int by_rF=0,born=0,outputdata=0;

    //read energy wavelength
    for (unsigned int i=0;i<cgi0.size();i++) {
    //std::cout<<"<br>"<<cgi0.at(i).name<<" "<<cgi0.at(i).field<<endl;
        if (cgi0.at(i).name == std::string("xenergy") ||cgi0.at(i).name == std::string("xlambda")||cgi0.at(i).name == std::string("sigma")||cgi0.at(i).name == std::string("qmin") ||cgi0.at(i).name == std::string("qmax") ||cgi0.at(i).name == std::string("dq") || cgi0.at(i).name == std::string("outputdata") ) {
            std::istringstream iss(cgi0.at(i).field);
            std::vector<double> va;
            std::copy(istream_iterator<double>(iss),istream_iterator<double>(), back_inserter (va));
    //std::cout<<va.size()<<" "<<cgi0.at(i).name<<" en="<<en<<" lambda="<<lambda<<endl;
            if (va.size()>0) {
                if (cgi0.at(i).name == std::string("xenergy")) {
                    en=va.at(0);
                } 
                if (cgi0.at(i).name == std::string("xlambda")) {
                    lambda=va.at(0);
                } 
                if (cgi0.at(i).name == std::string("sigma")) {
                    sigma=va.at(0);
                } 
                if (cgi0.at(i).name == std::string("qmin")) {
                    qmin=va.at(0);
                } 
                if (cgi0.at(i).name == std::string("qmax")) {
                    qmax=va.at(0);
                } 
                if (cgi0.at(i).name == std::string("dq")) {
                    dq=va.at(0);
                } 
            } else {
                if (cgi0.at(i).name == std::string("xenergy")) {
                    en=-1;
                } else {
                    lambda=-1;
                }
            }
        }
        if (cgi0.at(i).name == std::string("outputdata") ) outputdata=1; //enable by_rF
        if (cgi0.at(i).name == std::string("fresnel") ) by_rF=1; //enable by_rF
        if (cgi0.at(i).name == std::string("born") ) born=1; //use Born Approximation
        if (cgi0.at(i).name == std::string("formulae") ) {//processing formulae
            for (unsigned int j=0;j<cgi0.at(i).field.size();j++) {
                if (cgi0.at(i).field.at(j) !='+') f0.push_back(cgi0.at(i).field.at(j));
                else f0.push_back(' ');
            }
        }
                if (cgi0.at(i).name == std::string("rho") ) {//processing density profile
            for (unsigned int j=0;j<cgi0.at(i).field.size();j++) {
                if (cgi0.at(i).field.at(j) =='+') {
                        r0.push_back(' ');
                        continue;
                }
                if (cgi0.at(i).field.at(j) =='%'){
                        if (cgi0.at(i).field.substr(j,3)==std::string("%0D") ||cgi0.at(i).field.substr(j,3)==std::string("%0d")  ) {
                                buffer0.push_back(r0);
                                r0.clear();
                        }
                        if (cgi0.at(i).field.substr(j,3)==std::string("%2D") ||cgi0.at(i).field.substr(j,3)==std::string("%2d")  ) r0.push_back('-');
                        if (cgi0.at(i).field.substr(j,3)==std::string("%2B") ||cgi0.at(i).field.substr(j,3)==std::string("%2b")  ) r0.push_back('+');
                        j+=2;
                        continue;
                }
                        r0.push_back(cgi0.at(i).field.at(j));
            }
            if(r0.size()) buffer0.push_back(r0);
  // std::cout<<"<br>"<<r0;
        }
   //std::cout<<"<br>"<<r0<<" en="<<en<<" lambda="<<lambda<<endl;

    }
   //std::cout<<"<br>en="<<en<<" lambda="<<lambda;
   //std::cout<<"<br>"<<r0;
       char *pc= getenv("CONTENT_TYPE");
    if ( pc == NULL) {
        std::cout<<"Content-Type: text/plain\n\n";
        std::cout<<"getcgivars(): Unsupported Content-Type."<<s0<<std::endl;
        exit(1) ;
    }
    s0=string(pc);
    if (  s0 != std::string("application/x-www-form-urlencoded")) {
        std::cout<<"Content-Type: text/plain\n\n";
        std::cout<<"getcgivars(): Unsupported Content-Type."<<s0<<std::endl;
        exit(1) ;
    } else {
            if(outputdata){
                    //output text for xr
        std::cout<<"Content-Type: text/plain\n\n";
            } else {
                    //plot image
    std::cout<<"Content-type: image/png\n\n";
            }
    }

    if (  en>=0.02 &&en<=123.9) {
        lambda = E_to_l(en);
    }else en= E_to_l(lambda);
    if (! ( lambda>=0.02 &&lambda<=12.39)) {
        std::cout<<"lambda = "<<lambda<<" &Aring; ? Invalid input</body></html>";
        return 0;
    }

    istringstream iss(f0);
    vector <string > vs;
    std::copy (istream_iterator < string >(iss), istream_iterator < string >(), back_inserter (vs));
    if (vs.size()<1) {
        return 0;
    }
    vector<complex<double> > nbulki;

    //std::cout<<vs.size()<<" ";
    //std::cout<<"<table border=\"2\"><tbody>\n";
   // std::cout<<"<tr><td>X-ray energy = "<<en<<" KeV</td><td> &lambda; = "<<lambda<<" &Aring; </td></tr>\n";
   vector<compound> vcpmd0;
    for (vector<string>::iterator pvs=vs.begin(); pvs != vs.end();pvs++) {
        //std::cout<<*pvs<<endl;
        compound::compound cpmd0(*pvs,lambda,1.0);
        if (!cpmd0.els.size()) continue;
        if (!( pvs + 1 ==vs.end() )) {
            double rho0=1e-24; // in g angstrom^-3
            if (isdigit( (pvs+1)->at(0))) {
                pvs++;
                rho0=strtod(& ( pvs->at(0)),NULL)*1e-24;
                cpmd0.rho_scale(rho0);
            }else rho0=cpmd0.rho;
            //std::cout<<"rho= "<<rho0<<endl;
        }
        nbulki.push_back(cpmd0.nk);
        vcpmd0.push_back(cpmd0);
    }
    vector<vector<double> > rhoprofile;
    for(unsigned int i=0;i<buffer0.size();i++){
        istringstream issl(buffer0.at(i));
        vector <double > va;
    std::copy (istream_iterator < double >(issl),
               istream_iterator < double >(), back_inserter (va));
            //std::cout<<"<br>rho= "<<va.size()<<endl;
    if(va.size()<nbulki.size()+1) continue;
    rhoprofile.push_back(va);
    }
    //std::cout<<"size="<<rhoprofile.size();

    
multilayer::multilayer ml0(lambda);
ml0.tozprofile(nbulki,rhoprofile,vcpmd0);
if(sigma>=0.05) ml0.convolution(sigma);//gaussian roughness
ml0.by_rF=by_rF;
ml0.outputdata=outputdata;
ml0.born=born;
ml0.rfvector(qmin,qmax,dq);
ml0.plot();
    return 0;
}

