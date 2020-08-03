#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <random>
#include <numeric>
#include <vector>
#include <map>
#include <eigen3/Eigen/Eigen>
#include "exKF.hpp"
const double CLIGHT=3e8;
typedef exKF::exKF<double> EXKF;
typedef std::vector<std::size_t>   IArray; 
EXKF::Matrix getR(const std::size_t &n, const double &nu0)
{
    EXKF::Matrix R(n, n);
    for (std::size_t i = 0; i < n; i++)
        for (std::size_t j = 0; j < n; j++)
            R(i, j) = (i == j ? 2.0 : 1.0) * nu0 * nu0;
    return R;
}
int main(int argc, char **argv)
{
    EXKF::Matrix TS; 
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> randn;
    std::uniform_int_distribution<> randu(4,8);
    std::uniform_real_distribution<> randx(0,1500),randy(0,2000);
    double lambda=argc>=2?atof(argv[1]):0.01,beta=5e-9;
    EXKF::Array  x0={0.0,0.0,0.0,0.0,0.0,0.0},p={0.0},theory={0.0,0.0,0.0},measure;
    IArray obsidx,iobs;
    EXKF::Matrix P(6,6),Q,xc;
    P(0,0)=1.0;
    P(1,1)=1.0;
    P(2,2)=1.0/CLIGHT;
    P(3,3)=1.0/10;
    P(4,4)=1.0/10;
    P(5,5)=1.0/10/CLIGHT;
    Q=lambda*P;
    EXKF kf(x0,P,p,[](const EXKF::DiffArray&argsin,const EXKF::Array&parameter){
        EXKF::DiffArray argsout(6);
        argsout[0]=argsin[0]+argsin[3]*parameter[1];
        argsout[1]=argsin[1]+argsin[4]*parameter[1];
        argsout[2]=argsin[2]+argsin[5]*parameter[1];
        argsout[3]=argsin[3];
        argsout[4]=argsin[4];
        argsout[5]=argsin[5];
        return argsout;
    });
    TS.resize(20,2);
    obsidx.resize(TS.rows());
    for(std::size_t i=0;i<TS.rows();i++)
    {
        TS(i,0)=randx(gen);
        TS(i,1)=randy(gen);
        obsidx[i]=i;
        std::cout<<i<<":("<<TS(i,0)<<","<<TS(i,1)<<").\n";
    }
    for (auto i=0;i<10000;i++){
        p[0]=0.1*i;
        kf.predict(p,Q);
        theory[0]=0.1*i*sin(p[0]/300*3.14159265)+0.1*i;
        theory[1]=0.2*i*sin(p[0]/500*3.14159265)+0.2*i;
        theory[2]=0.3e-8+0.04e-12*i;
        std::random_shuffle(obsidx.begin(),obsidx.end());
        auto k=randu(gen);
        iobs.resize(k);
        measure.resize(k);
        for(std::size_t j=0;j<k;j++)
        { 
            iobs[j]=obsidx[j];
            auto t=0.0,dx=theory[0]-TS(iobs[j],0),dy=theory[1]-TS(iobs[j],1);
            t=sqrt(dx*dx+dy*dy);
            t=t/CLIGHT-theory[2]+randn(gen)*beta;
            measure[j]=t;
        }
        //argsin:Xms,Yms,dTms
        kf.update_nolinear(measure,getR(k,beta),[iobs,measure,TS](const EXKF::DiffArray &argsin, const EXKF::Array &parameter) {
        size_t nv=measure.size();
        EXKF::DiffArray argsout(nv);
        for(auto j=0; j<nv;j++){
            auto Lis=exKF::sqrt((argsin[0]-TS(iobs[j],0))*(argsin[0]-TS(iobs[j],0))+
            (argsin[1]-TS(iobs[j],1))*(argsin[1]-TS(iobs[j],1)));
            argsout[j]=Lis/CLIGHT-argsin[2];
        }
        return argsout;
    });
        xc=kf.getState();
        for(auto j=0;j<k;j++)
        std::cout<<(j==0?"measure:(":",(")<<iobs[j]<<","<<measure[j]<<")";
        std::cout<<".\nstatus:\t"<<theory[0]<<"\t"<<theory[1]<<"\t"<<theory[2]<<"\t"<<xc(0,0)<<"\t"<<xc(1,0)<<"\t"<<xc(2,0)<<"\n";
    }
    return 0;    
}