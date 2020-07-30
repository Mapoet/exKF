//
//  Diff.h
//  SphExpress
//
//  Created by Mapoet Niphy on 2020/05/30.
//  Copyright 2020 Mapoet Niphy. All rights reserved.
//

#ifndef __Diff_h__
#define __Diff_h__
#include <valarray>
namespace exKF{
    template<class Cell>
    struct Diff{
        Cell   _val;
        std::valarray<Cell>   _dval;
        Diff(const Cell& val = Cell(), std::valarray<Cell> dval = std::valarray<Cell>()) :_val(val), _dval(dval){}
        Diff<Cell> operator=(const Cell& v){
            return Diff<Cell>(v);
        }
        operator Cell()const{
            return _val;
        }
        friend Diff<Cell> operator +(const Diff<Cell>& a) {
            return a;
        }
        friend Diff<Cell> operator -(const Diff<Cell>& a) {
            Diff<Cell> b = a;
            b._val *= Cell(-1.0);
            b._dval = Cell(-1.0)* b._dval;
            return b;
        }
        Diff<Cell> operator +=(const Diff<Cell>& a) {
            _val = _val + a._val;
            _dval = _dval + a._dval;
            return *this;
        }
        Diff<Cell> operator +=(const Cell& a){
            _val =_val+ a;
            _dval =a._dval;
            return *this;
        }
        Diff<Cell> operator -=(const Diff<Cell>& a){
			_val = _val- a._val;
			_dval = _dval- a._dval;
            return *this;
        }
        Diff<Cell> operator -=(const Cell& a) {
            _val = _val - a;
            _dval = _dval;
            return *this;
        }
        Diff<Cell> operator *=(const Diff<Cell>& a){
            auto  val = _val;
            auto dval = _dval;
            _val = val*a._val;
            _dval = dval*a._val + val*a._dval;
            return *this;
        }
        Diff<Cell> operator *=(const Cell& a) {
            auto  val = _val;
            auto dval = _dval;
            _val = val * a;
            _dval = dval * a;
            return *this;
        }
        Diff<Cell> operator /=(const Diff<Cell>& a){
            auto  val = _val;
            auto dval = _dval;
            _val = val / a._val;
            _dval = (dval*a._val - val*a._dval) / (a._val*a._val);
            return *this;
        }
        Diff<Cell> operator /=(const Cell& a) {
            _val = _val / a;
            _dval = _dval /a;
            return *this;
        }
    }; 
    template<class Cell>
    Diff<Cell> operator+(const Diff<Cell>& a, const Diff<Cell>& b){
        Diff<Cell> v = a;
        v += b;
        return v;
    }
    template<class Cell>
    Diff<Cell> operator-(const Diff<Cell>& a, const Diff<Cell>& b){
        Diff<Cell> v = a;
        v -= b;
        return v;
    }
    template<class Cell>
    Diff<Cell> operator*(const Diff<Cell>& a, const Diff<Cell>& b){
        Diff<Cell> v = a;
        v *= b;
        return v;
    }
    template<class Cell>
    Diff<Cell> operator/(const Diff<Cell>& a, const Diff<Cell>& b){
        Diff<Cell> v = a;
        v /= b;
        return v;
    }
    template<class Cell>
    Diff<Cell> operator+(Cell a, const Diff<Cell>& b){
        Diff<Cell> v =b;
        v += a;
        return v;
    }
    template<class Cell>
    Diff<Cell> operator-(Cell a, const Diff<Cell>& b){
        Diff<Cell> v;
        v._val = a - b._val;
        v._dval = -b._dval;
        return v;
    }
    template<class Cell>
    Diff<Cell> operator*(const Cell& a, const Diff<Cell>& b){
        Diff<Cell> v;
        v._val = a * b._val;
        v._dval = a*b._dval;
        return v;
    }
    template<class Cell>
    Diff<Cell> operator/(const Cell& a, const Diff<Cell>& b){
        Diff<Cell> v;
        v._val = a / b._val;
        v._dval = -a*b._dval / (b._val*b._val);
        return v;
    }
    template<class Cell>
    Diff<Cell> operator+(const Diff<Cell>& a, const Cell& b){
        Diff<Cell> v = a;
		v._val = v._val+ b;
        return v;
    }
    template<class Cell>
    Diff<Cell> operator-(const Diff<Cell>& a, const Cell& b){
        Diff<Cell> v = a;
		v._val = v._val- b;
        return v;
    }
    template<class Cell>
    Diff<Cell> operator*(const Diff<Cell>& a, const Cell& b){
        Diff<Cell> v = a;
		v._val = v._val* b;
		v._dval = v._dval* b;
        return v;
    }
    template<class Cell>
    Diff<Cell> operator/(const Diff<Cell>& a, const Cell& b){
        Diff<Cell> v = a;
		v._val = v._val/ b;
		v._dval = v._dval/ b;
        return v;
	}
	template<class Cell>
	Diff<Cell> pow(const Diff<Cell>& a, const Diff<Cell>& b){
		Diff<Cell> v;
		v._val = std::pow(a._val, b._val);
		v._dval = b._val*std::pow(a._val, b._val - 1)*a._dval + v._val*std::log(a._val)*b._dval;
		return v;
	}
	template<class Cell>
	Diff<Cell> pow(const Cell& a, const Diff<Cell>& b){
		Diff<Cell> v;
		v._val = std::pow(a, b._val);
		v._dval = v._val*std::log(a)*b._dval;
		return v;
	}
	template<class Cell>
	Diff<Cell> pow(const Diff<Cell>& a, const Cell& b){
		Diff<Cell> v;
		v._val = std::pow(a._val, b);
		v._dval = b*std::pow(a._val, b - 1)*a._dval;
		return v;
	}
    template<class Cell>
    Diff<Cell> sin(const Diff<Cell>& a){
        Diff<Cell> v;
        v._val = sin(a._val);
        v._dval = cos(a._val)*a._dval;
        return v;
    }
    template<class Cell>
    Diff<Cell> cos(const Diff<Cell>& a){
        Diff<Cell> v;
        v._val = cos(a._val);
        v._dval = -sin(a._val)*a._dval;
        return v;
    }
    template<class Cell>
    Diff<Cell> tan(const Diff<Cell>& a){
        Diff<Cell> v = sin(a) / cos(a);
        return v;
	}
	template<class Cell>
	Diff<Cell> exp(const Diff<Cell>& a){
		Diff<Cell> v;
		v._val = std::exp(a._val);
		v._dval = std::exp(a._val)*a._dval;
		return v;
	}
	template<class Cell>
	Diff<Cell> asin(const Diff<Cell>& a){
		Diff<Cell> v;
		v._val = asin(a._val);
		v._dval = 1.0/sqrt((1.0-a._val*a._val))*a._dval;
		return v;
	}
	template<class Cell>
	Diff<Cell> acos(const Diff<Cell>& a){
		Diff<Cell> v;
		v._val = acos(a._val);
		v._dval = -1.0/sqrt((1.0-a._val*a._val))*a._dval;
		return v;
	}
	template<class Cell>
	Diff<Cell> atan(const Diff<Cell>& a){
        Diff<Cell> v;
        v._val = atan(a._val);
        v._dval = -1.0/(a._val*a._val+1.0)*a._dval;
		return v;
	}
	template<class Cell>
	Diff<Cell> log(const Diff<Cell>& a){
		Diff<Cell> v;
		v._val = std::log(a._val);
		v._dval = a._dval / (a._val);
		return v;
	}
    template<class Cell>
    Diff<Cell> if_gt(const Diff<Cell>& u, const Diff<Cell>& c, const Diff<Cell>& x, const Diff<Cell>& a, const Diff<Cell>& b) {
        Diff<Cell> cf = Cell(1.0) / (Cell(1.0) + exp(-u * (x - c)));
        return (Cell(1.0) - cf) * a + cf * b;
    }
    template<class Cell>
    Diff<Cell> if_gt(const Cell& u, const Diff<Cell>& c, const Diff<Cell>& x, const Diff<Cell>& a, const Diff<Cell>& b) {
        Diff<Cell> cf = Cell(1.0) / (Cell(1.0) + exp(-u * (x - c)));
        return (Cell(1.0) - cf) * a + cf * b;
    }
    template<class Cell>
    Diff<Cell> if_gt(const Diff<Cell>& u, const Cell& c, const Diff<Cell>& x, const Diff<Cell>& a, const Diff<Cell>& b) {
        Diff<Cell> cf = Cell(1.0) / (Cell(1.0) + exp(-u * (x - c)));
        return (Cell(1.0) - cf) * a + cf * b;
    }
    template<class Cell>
    Diff<Cell> if_gt(const Cell& u, const Cell& c, const Diff<Cell>& x, const Diff<Cell>& a, const Diff<Cell>& b) {
        Diff<Cell> cf = Cell(1.0) / (Cell(1.0) + exp(-u * (x - c)));
        return (Cell(1.0) - cf) * a + cf * b;
    }
    template<class Cell>
    Diff<Cell> if_is(const Diff<Cell>& u, const Diff<Cell>& c, const Diff<Cell>& x, const Diff<Cell>& a, const Diff<Cell>& b) {
        Diff<Cell> cf = Cell(1.0) / (Cell(1.0) + exp(-u*u * (x - c) * (x - c)));
        return (Cell(1.0) - cf) * a + cf * b;
    }
    template<class Cell>
    Diff<Cell> if_is(const Cell& u, const Diff<Cell>& c, const Diff<Cell>& x, const Diff<Cell>& a, const Diff<Cell>& b) {
        Diff<Cell> cf = Cell(1.0) / (Cell(1.0) + exp(-u * u * (x - c) * (x - c)));
        return (Cell(1.0) - cf) * a + cf * b;
    }
    template<class Cell>
    Diff<Cell> if_is(const Diff<Cell>& u, const Cell& c, const Diff<Cell>& x, const Diff<Cell>& a, const Diff<Cell>& b) {
        Diff<Cell> cf = Cell(1.0) / (Cell(1.0) + exp(-u * u * (x - c) * (x - c)));
        return (Cell(1.0) - cf) * a + cf * b;
    }
    template<class Cell>
    Diff<Cell> if_is(const Cell& u, const Cell& c, const Diff<Cell>& x, const Diff<Cell>& a, const Diff<Cell>& b) {
        Diff<Cell> cf = Cell(1.0) / (Cell(1.0) + exp(-u * u * (x - c) * (x - c)));
        return (Cell(1.0) - cf) * a + cf * b;
    }
    template<class Cell>
    Diff<Cell>& abs(const Diff<Cell>& u, const Diff<Cell>& a) {
        return if_gt(u,Cell(0.0),a,a,-a);
    }
    template<class Cell>
    Diff<Cell>& abs(const Cell& u, const Diff<Cell>& a) {
        return if_gt(u, Cell(0.0), a, a, -a);
    }
}
#endif /* __Diff_h__ */
