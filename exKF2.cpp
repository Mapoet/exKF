
#include <iostream>
#include <random>
#include <numeric>
#include"exKF.hpp"
int main(int argc,char**argv){
    typedef exKF::exKF<double>    exKF;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> randn;
    double lambda=argc>=2?atof(argv[1]):0.0001,beta=50;
    exKF::Array  x0={0.0,0.0,0.0,0.0,0.0,0.0},p={0.0},theory={0.0,0.0},measure={0.0,0.0};
    exKF::Matrix P(6,6),Q(6,6),R(2,2),xc;
    P(0,0)=1.0;
    P(1,1)=1.0;
    P(2,2)=1.0;
    P(3,3)=1.0;
    P(4,4)=1.0;
    P(5,5)=1.0;
    Q=lambda*P;
    R(0,0)=beta*beta;
    R(1,1)=beta*beta;
    exKF kf(x0,P,p,[](const exKF::DiffArray&argsin,const exKF::Array&parameter){
        exKF::DiffArray argsout(6);
        argsout[0]=argsin[0]+argsin[2]*parameter[1]+argsin[4]*parameter[1]*parameter[1]/2;
        argsout[1]=argsin[1]+argsin[3]*parameter[1]+argsin[5]*parameter[1]*parameter[1]/2;
        argsout[2]=argsin[2]+argsin[4]*parameter[1];
        argsout[3]=argsin[3]+argsin[5]*parameter[1];
        argsout[4]=argsin[4];
        argsout[5]=argsin[5];
        return argsout;
    });
    for (auto i=0;i<10000;i++){
        p[0]=0.1*i;
        kf.predict(p,Q);
        theory[0]=0.1*i*sin(p[0]/300*3.14159265)*(1+randn(gen)*0.02)+0.1*i;
        theory[1]=0.2*i*sin(p[0]/500*3.14159265)*(1+randn(gen)*0.02)+0.2*i;
        measure[0]=theory[0]+randn(gen)*beta;
        measure[1]=theory[1]+randn(gen)*beta;
        kf.update_linear(measure,R,[](const exKF::DiffArray&argsin,const exKF::Array&parameter){
        exKF::DiffArray argsout(2);
        argsout[0]=argsin[0];
        argsout[1]=argsin[1];
        return argsout;
    });
        xc=kf.getState();
        std::cout<<theory[0]<<"\t"<<theory[1]<<"\t"<<measure[0]<<"\t"<<measure[1]<<"\t"<<xc(0,0)<<"\t"<<xc(1,0)<<"\n";
    }
    return 0;
}