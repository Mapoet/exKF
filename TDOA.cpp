//
//  TDOA.cpp
//  TDOA
//
//  Created by Mapoet Niphy on 2020/7/29.
//  Copyright © 2020年 Mapoet Niphy. All rights reserved.
//
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <eigen3/Eigen/Eigen>
#include "exKF.hpp"
const double CLIGHT=3e8;
typedef exKF::exKF<double> EXKF;
typedef std::map<std::size_t,double> TDOA; 
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
    TDOA obs;
    //TS:nx2
    EXKF::Matrix TS,xs; 
    //argsin:Xms,Yms,dTms
    auto meansure = [obs,TS](const EXKF::DiffArray &argsin, const EXKF::Array &parameter) {
        size_t nv=obs.size(),i=0;
        EXKF::DiffArray argsout(nv);
        for(auto it=obs.begin(); it!=obs.end();it++,i++){
            auto Lis=exKF::sqrt((argsin[0]-TS(it->first,0))*(argsin[0]-TS(it->first,0))+
            (argsin[1]-TS(it->first,1))*(argsin[1]-TS(it->first,1)));
            argsout[i]=Lis/CLIGHT-argsin[2];
        }
        return argsout;
    };
}