//
//  exkf_tnna.h
//  EXKF
//
//  Created by Mapoet Niphy on 2020/7/29.
//  Copyright © 2020年 Mapoet Niphy. All rights reserved.
//
#define EIGEN_USE_BLAS
#include <memory>
#include <chrono>
//#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
//#include "Diff_tnna.h"
#include "Diff.hpp"
namespace exKF{
    
    template<typename Cell>
    class exKF{
        public:
        typedef std::valarray<Diff<Cell>> (*Apply)(const std::valarray<Diff<Cell>>& argsin,const std::valarray<Cell>&parameter);
        private:
        Apply _theory;
        Eigen::Matrix<Cell,-1,-1> _P;
        std::valarray<Diff<Cell>> _X;
        std::valarray<Cell> _pdp;
        public:
        template<typename Theory>
        exKF(const std::valarray<Cell> &state,const Eigen::Matrix<Cell,-1,-1> &P,const std::valarray<Cell> &p,const Theory&theory){
            auto n=state.size();
            assert(n==P.rows()&&state.size()==P.cols());
            _X.resize(n);
            for(std::size_t i=0;i<n;i++)
                _X[i]._val=state[i];
            _pdp=p;
            _P=P;
            _pdp.resize(p.size()*2);
            _theory=theory;
        }
        void predict(std::valarray<Cell>&p,const Eigen::Matrix<Cell,-1,-1> &Q){
            Eigen::Matrix<Cell,-1,-1> F;
            auto np=p.size();
            auto n=_X.size();
            assert(np*2==_pdp.size()&&n==Q.rows()&&n==Q.cols());
            for(std::size_t i=0;i<np;i++)
            _pdp[np+i]=p[i]-_pdp[i];
            // update _X
            for(std::size_t i=0;i<n;i++)
            {
                _X[i]._dval.resize(n,Cell(0.0));
                _X[i]._dval[i]=Cell(1.0);
            }
            _X=_theory(_X,_pdp);
            assert(_X.size()==n);
            // get F
            F.resize(n,n);
            for(std::size_t i=0;i<n*n;i++)F(i/n,i%n)=_X[i/n]._dval[i%n];
            // update _P
            _P=F*_P*F.transpose()+Q;
        }
        template<typename Measure>
        void update(const std::valarray<Cell> &value,const Eigen::Matrix<Cell,-1,-1> &R,const Measure & measure){
            Eigen::Matrix<Cell,-1,-1> H,K,dX;
            Eigen::Matrix<Cell,-1,-1> dY;
            auto n=_X.size();
            auto np=_pdp.size()/2;
            auto nv=value.size();
            assert(nv==R.rows()&&nv==R.cols());
            for(std::size_t i=0;i<n;i++)
            {
                _X[i]._dval.resize(n,Cell(0.0));
                _X[i]._dval[i]=Cell(1.0);
            }
            auto Y=measure(_X,_pdp);
            assert(Y.size()==nv);
            dY.resize(nv,1);
            // get \delta Y
            for(std::size_t i=0;i<nv;i++)dY(i)=value[i]-Y[i]._val;
            H.resize(nv,n);
            // get H
            for(std::size_t i=0;i<nv*n;i++)H(i/n,i%n)=_X[i/n]._dval[i%n];
            // get S[K]
            K=H*_P*H.transpose()+R;
            // get K (gain matrix)
            K=K.inverse();
            K=_P*H.transpose()*K;
            // update _X
            dX=K*dY;
            for(std::size_t i=0;i<n;i++)
            {
                _X[i]._val=_X[i]._val+dX(i);
                _X[i]._dval.resize(n,Cell(0.0));
                _X[i]._dval[i]=Cell(1.0);
            }
            // update _P
            _P=_P-K*H*_P;
            for(std::size_t i=0;i<np;i++)
            _pdp[i]=_pdp[i]+_pdp[i+np];
        }
        Eigen::Matrix<Cell,-1,-1> getState()const{
            Eigen::Matrix<Cell,-1,-1> state;
            auto n=_X.size();
            state.resize(n,n+1);
            for(std::size_t i=0;i<n;i++)
            for(std::size_t j=0;j<n+1;j++)
            state(i,j)=(j==0?_X[i]._val:_P(i,j-1));
            return state;
        }
    };



}