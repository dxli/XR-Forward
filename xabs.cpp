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
    std::cout<<"Content-type: text/html\n\n";
    std::string s0;
    char *pc= getenv("CONTENT_TYPE");
    if (  pc == NULL) {
        std::cout<<"Content-Type: text/plain\n\n";
        std::cout<<"getcgivars(): Unsupported Content-Type."<<s0<<std::endl;
        exit(1) ;
    }
    s0=string(pc);
    if (  s0 != std::string("application/x-www-form-urlencoded")) {
        std::cout<<"Content-Type: text/plain\n\n";
        std::cout<<"getcgivars(): Unsupported Content-Type."<<s0<<std::endl;
        exit(1) ;
    }
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
    double en=-1,lambda=-1;
    std::string f0;
    //read energy wavelength
    for (unsigned int i=0;i<cgi0.size();i++) {
  //  std::cout<<cgi0.at(i).name<<" "<<cgi0.at(i).field<<endl;
        if (cgi0.at(i).name == std::string("xenergy") ||cgi0.at(i).name == std::string("xlambda") ) {
            std::istringstream iss(cgi0.at(i).field);
            std::vector<double> va;
            std::copy(istream_iterator<double>(iss),istream_iterator<double>(), back_inserter (va));
    //std::cout<<va.size()<<" "<<cgi0.at(i).name<<" en="<<en<<" lambda="<<lambda<<endl;
            if (va.size()>0) {
                if (cgi0.at(i).name == std::string("xenergy")) {
                    en=va.at(0);
                } else {
                    lambda=va.at(0);
                }
            } else {
                if (cgi0.at(i).name == std::string("xenergy")) {
                    en=-1;
                } else {
                    lambda=-1;
                }
            }
        }
        if (cgi0.at(i).name == std::string("formulae") ) {//processing formulae
            for (unsigned int j=0;j<cgi0.at(i).field.size();j++) {
                if (cgi0.at(i).field.at(j) !='+') f0.push_back(cgi0.at(i).field.at(j));
                else f0.push_back(' ');
            }
        }
    }
    std::cout<<"<html><head></head><body><center>\n";
   // std::cout<<"en="<<en<<" lambda="<<lambda;
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
        std::cout<<"no formula found</body></html>";
        return 0;
    }
    //std::cout<<vs.size()<<" ";
    std::cout<<"<table border=\"2\"><tbody>\n";
    std::cout<<"<tr><td>X-ray energy = "<<en<<" KeV</td><td> &lambda; = "<<lambda<<" &Aring; </td></tr>\n";
    for (vector<string>::iterator pvs=vs.begin(); pvs != vs.end();pvs++) {
        //std::cout<<*pvs<<endl;
        compound cpmd0(*pvs,lambda,1.0);
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
        std::cout<<"<tr><td>";
        for(unsigned int j=0;j<cpmd0.els.size();j++){
                std::cout<<cpmd0.els.at(j).symbol;
                if(j<cpmd0.n0.size()) 
                        if(cpmd0.n0.at(j) != 1) std::cout<<"<sub>"<<cpmd0.n0.at(j)<<"</sub>";
        }
        std::cout<<"</td><tr>\n";
        std::cout<<"<tr><td> &rho; = "<<cpmd0.rho*1e24<<" g cm<sup>-3</sup></td><td> &rho; <sub>el</sub> ="<<cpmd0.rho_el<<" &Aring; <sup>-3</sup></td></tr>\n";
        double delta=re_classical*lambda*lambda*cpmd0.rho_el/(2.*M_PI);
        std::cout<<"<tr><td>critical angle = "<<sqrt(re_classical*cpmd0.rho_el/M_PI)*lambda*180./M_PI<<" &deg; </td><td>critical Q ="<<4.*sqrt(M_PI*re_classical*cpmd0.rho_el)<<" &Aring; <sup>-1</sup></td></tr>\n";
        std::cout<<"<tr><td>&beta; = "<<cpmd0.beta<<"</td><td> Absorption Length ="<<1e-7/cpmd0.mu<<" mm</td></tr>";
    }

    std::cout<<"</tbody></table></body></html>";
    return 0;
}

