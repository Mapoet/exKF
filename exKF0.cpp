
#include <iostream>
#include <random>
#include <numeric>
#include"exKF.hpp"
int main(int argc,char**argv){
    typedef exKF::exKF<double>    exKF;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> randn;
    double lambda=argc>=2?atof(argv[1]):0.1,beta=50;
    exKF::Array  x0={0.0,0.0},p={0.0},theory={0.0,0.0},measure={0.0,0.0};
    exKF::Matrix P(2,2),Q(2,2),R(2,2),xc;
    P(0,0)=1.0;
    P(1,1)=1.0;
    Q=lambda*P;
    R(0,0)=beta*beta;
    R(1,1)=beta*beta;
    exKF kf(x0,P,p,[](const exKF::DiffArray&argsin,const exKF::Array&parameter){
        exKF::DiffArray argsout(2);
        argsout[0]=argsin[0];
        argsout[1]=argsin[1];
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