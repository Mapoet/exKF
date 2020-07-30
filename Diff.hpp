//
//  Diff_tnna.h
//  TNNA
//
//  Created by Mapoet Niphy on 2018/11/2.
//  Copyright © 2018年 Mapoet Niphy. All rights reserved.
//

#ifndef Diff_h
#define Diff_h
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
        Diff<Cell> operator +=(const Diff<Cell>& a){
            _val =_val+ a._val;
            _dval =_dval+ a._dval;
            return *this;
        }
        Diff<Cell> operator -=(const Diff<Cell>& a){
			_val = _val- a._val;
			_dval = _dval- a._dval;
            return *this;
        }
        Diff<Cell> operator *=(const Diff<Cell>& a){
            auto  val = _val;
            auto dval = _dval;
            _val = val*a._val;
            _dval = dval*a._val + val*a._dval;
            return *this;
        }
        Diff<Cell> operator /=(const Diff<Cell>& a){
            auto  val = _val;
            auto dval = _dval;
            _val = val / a._val;
            _dval = (dval*a._val - val*a._dval) / (a._val*a._val);
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
	Diff<Cell>& abs(const Diff<Cell>& a){
		Diff<Cell> v;
		v._val = std::abs(a._val);
		v._dval = (a._val > 0 ? 1.0 : -1.0)*a._dval;
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
    template<typename IOS,typename Cell>
    IOS& operator <<(IOS&ios,const Diff<Cell>&cell){
        ios<<"("<<cell._val<<":";
        for(std::size_t i=0;i<cell._dval.size();i++)
        ios<<cell._dval[i]<<(i==cell._dval.size()-1?')':',');
        return ios;
    }
}
#endif /* Diff_h */
