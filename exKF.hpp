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
namespace exKF
{

    template <typename Cell>
    class exKF
    {
    public:
        typedef std::valarray<Cell> Array;
        typedef std::valarray<Diff<Cell>> DiffArray;
        typedef Eigen::Matrix<Cell, -1, -1> Matrix;
        typedef DiffArray (*Apply)(const DiffArray &argsin, const Array &parameter);

    private:
        Apply _theory;
        Matrix _P;
        DiffArray _X;
        Array _pdp;
        Cell _eps2;

    public:
        template <typename Theory>
        exKF(const Array &state, const Matrix &P, const Array &p, const Theory &theory, const Cell &eps = Cell(3.0))
        {
            auto n = state.size();
            assert(n == P.rows() && state.size() == P.cols());
            _X.resize(n);
            for (std::size_t i = 0; i < n; i++)
                _X[i]._val = state[i];
            _pdp = p;
            _P = P;
            _pdp.resize(p.size() * 2);
            _theory = theory;
            _eps2 = eps * eps;
        }
        void predict(const Array &p, const Matrix &Q)
        {
            Matrix F;
            auto np = p.size();
            auto n = _X.size();
            assert(np * 2 == _pdp.size() && n == Q.rows() && n == Q.cols());
            for (std::size_t i = 0; i < np; i++)
                _pdp[np + i] = p[i] - _pdp[i];
            // update _X
            for (std::size_t i = 0; i < n; i++)
            {
                _X[i]._dval.resize(n, Cell(0.0));
                _X[i]._dval[i] = Cell(1.0);
            }
            _X = _theory(_X, _pdp);
            assert(_X.size() == n);
            // get F
            F.resize(n, n);
            for (std::size_t i = 0; i < n * n; i++)
                F(i / n, i % n) = _X[i / n]._dval[i % n];
            // update _P
            _P = F * _P * F.transpose() + Q;
        }
        template <typename Measure>
        void update_nolinear(const Array &value, const Matrix &R, const Measure &measure)
        {
            Matrix H, K, dX;
            Matrix dY;
            Cell dxx = Cell(0), dyy = Cell(0);
            auto n = _X.size();
            auto np = _pdp.size() / 2;
            auto nv = value.size();
            std::size_t it = 0;
            assert(nv == R.rows() && nv == R.cols());
            do
            {
                for (std::size_t i = 0; i < n; i++)
                {
                    _X[i]._dval.resize(n, Cell(0.0));
                    _X[i]._dval[i] = Cell(1.0);
                }
                auto Y = measure(_X, _pdp);
                assert(Y.size() == nv);
                dY.resize(nv, 1);
                H.resize(nv, n);
                // get \delta Y
                dyy = Cell(0.0);
                for (std::size_t i = 0; i < nv; i++)
                {
                    dY(i) = value[i] - Y[i]._val;
                    dyy = dyy + dY(i) * dY(i) / R(i, i);
                }
                // get H
                for (std::size_t i = 0; i < nv * n; i++)
                    H(i / n, i % n) = _X[i / n]._dval[i % n];
                // get S[K]
                K = H * _P * H.transpose() + R;
                //modify the error by its error
                for (std::size_t i = 0; i < nv; i++)
                    if (dY(i) * dY(i) / R(i, i) > _eps2 * dyy / nv)
                    {
                        for (std::size_t j = 0; j < nv; j++)
                        {
                            if (i == j)
                            {
                                K(i, i) = K(i, j) + R(i, j) * 99;
                                continue;
                            }
                            K(i, j) = K(i, j) + R(i, j) * 9;
                            K(j, i) = K(j, i) + R(j, i) * 9;
                        }
                    }
                // get K (gain matrix)
                K = K.inverse();
                K = _P * H.transpose() * K;
                // update _X
                dX = K * dY;
                K = _P - K * H * _P;
                dxx = Cell(0.0);
                for (std::size_t i = 0; i < n; i++)
                {
                    dxx = dxx + dX(i) * dX(i) / K(i, i);
                    _X[i]._val = _X[i]._val + dX(i);
                    _X[i]._dval.resize(n, Cell(0.0));
                    _X[i]._dval[i] = Cell(1.0);
                }
//            fprintf(stderr,"%d,%lf,%lf,%lf\n",it,dxx/n,dX(0),K(0,0));
            } while (dxx / n > 1e-4 && it++ < 3);
            // update _P
            _P = K;
            for (std::size_t i = 0; i < np; i++)
                _pdp[i] = _pdp[i] + _pdp[i + np];
        }
        template <typename Measure>
        void update_linear(const Array &value, const Matrix &R, const Measure &measure)
        {
            Matrix H, K, dX;
            Matrix dY;
            auto n = _X.size();
            auto np = _pdp.size() / 2;
            auto nv = value.size();
            assert(nv == R.rows() && nv == R.cols());
            for (std::size_t i = 0; i < n; i++)
            {
                _X[i]._dval.resize(n, Cell(0.0));
                _X[i]._dval[i] = Cell(1.0);
            }
            auto Y = measure(_X, _pdp);
            assert(Y.size() == nv);
            dY.resize(nv, 1);
            H.resize(nv, n);
            // get \delta Y
            for (std::size_t i = 0; i < nv; i++)
                dY(i) = value[i] - Y[i]._val;
            // get H
            for (std::size_t i = 0; i < nv * n; i++)
                H(i / n, i % n) = _X[i / n]._dval[i % n];
            // get S[K]
            K = H * _P * H.transpose() + R;
            // get K (gain matrix)
            K = K.inverse();
            K = _P * H.transpose() * K;
            // update _X
            dX = K * dY;
            // update _X
            for (std::size_t i = 0; i < n; i++)
            {
                _X[i]._val = _X[i]._val + dX(i);
                _X[i]._dval.resize(n, Cell(0.0));
                _X[i]._dval[i] = Cell(1.0);
            }
            // update _P
            _P = _P - K * H * _P;
            for (std::size_t i = 0; i < np; i++)
                _pdp[i] = _pdp[i] + _pdp[i + np];
        }
        Eigen::Matrix<Cell, -1, -1> getState() const
        {
            Eigen::Matrix<Cell, -1, -1> state;
            auto n = _X.size();
            state.resize(n, n + 1);
            for (std::size_t i = 0; i < n; i++)
                for (std::size_t j = 0; j < n + 1; j++)
                    state(i, j) = (j == 0 ? _X[i]._val : _P(i, j - 1));
            return state;
        }
    };

} // namespace exKF