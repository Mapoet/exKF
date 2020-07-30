
#include <iostream>
#include <random>
#include"exKF.hpp"
int main(){
    typedef std::valarray<double> Array;
    typedef Eigen::MatrixXd       Matrix;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> randn;
    double lambda=0.1,beta=0.5;
    Array  x0={0.0,0.1},p={0.0},theory={0.0},measure={0.0};
    Matrix P(2,2),Q(2,2),R(1,1),xc;
    P(0,0)=1.0;
    P(0,1)=0.0;
    P(1,0)=0.0;
    P(1,1)=1.0;
    Q=lambda*P;
    exKF::exKF<double> kf(x0,P,p,[](const std::valarray<exKF::Diff<double>>&argsin,const std::valarray<double>&parameter){
        std::valarray<exKF::Diff<double>> argsout(2);
        argsout[0]=argsin[0]+argsin[1]*parameter[1];
        argsout[1]=argsin[1];
        return argsout;
    });
    for (auto i=0;i<10000;i++){
        p[0]=0.1*i;
        kf.predict(p,Q);
        theory[0]=0.1*i;
        measure[0]=theory[0]+randn(gen)*beta;
        R(0,0)=beta*beta;
        kf.update(measure,R,[](const std::valarray<exKF::Diff<double>>&argsin,const std::valarray<double>&parameter){
        std::valarray<exKF::Diff<double>> argsout(1);
        argsout[0]=argsin[0];
        return argsout;
    });
        xc=kf.getState();
        std::cout<<theory[0]<<"\t"<<measure[0]<<"\t"<<xc(0,0)<<"\n";
    }
    return 0;
}