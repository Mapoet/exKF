
#include <iostream>
#include <random>
#include <numeric>
#include"exKF.hpp"
int main(){
    typedef std::valarray<double> Array;
    typedef Eigen::MatrixXd       Matrix;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> randn;
    double lambda=0.1,beta=5;
    Array  x0={0.0,0.0,0.1,0.2},p={0.0},theory={0.0,0.0},measure={0.0,0.0};
    Matrix P(4,4),Q(4,4),R(2,2),xc;
    P(0,0)=1.0;
    P(1,1)=1.0;
    P(2,2)=1.0;
    P(3,3)=1.0;
    Q=lambda*P;
    exKF::exKF<double> kf(x0,P,p,[](const std::valarray<exKF::Diff<double>>&argsin,const std::valarray<double>&parameter){
        std::valarray<exKF::Diff<double>> argsout(4);
        argsout[0]=argsin[0]+argsin[2]*parameter[1];
        argsout[1]=argsin[1]+argsin[3]*parameter[1];
        argsout[2]=argsin[2];
        argsout[3]=argsin[3];
        return argsout;
    });
    for (auto i=0;i<10000;i++){
        p[0]=0.1*i;
        kf.predict(p,Q);
        theory[0]=0.1*i*sin(p[0]/300*3.14159265)+0.1*i;
        theory[1]=0.2*i*sin(p[0]/500*3.14159265)+0.2*i;
        measure[0]=theory[0]+randn(gen)*beta;
        measure[1]=theory[1]+randn(gen)*beta;
        R(0,0)=beta*beta;
        R(1,1)=beta*beta;
        kf.update(measure,R,[](const std::valarray<exKF::Diff<double>>&argsin,const std::valarray<double>&parameter){
        std::valarray<exKF::Diff<double>> argsout(2);
        argsout[0]=argsin[0];
        argsout[1]=argsin[1];
        return argsout;
    });
        xc=kf.getState();
        std::cout<<theory[0]<<"\t"<<theory[1]<<"\t"<<measure[0]<<"\t"<<measure[1]<<"\t"<<xc(0,0)<<"\t"<<xc(1,0)<<"\n";
    }
    return 0;
}