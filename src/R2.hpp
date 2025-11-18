#ifndef R2_hpp
#define R2_hpp
#include <cmath>
#include <iostream>
#include <cassert> // for assert
#include <cstdlib> // for exit

using namespace std;
typedef double R; 

class R2 {
   public:  
  	R x,y;
  	R2(): x(0.),y(0.){}  
  	R2(R a, R b): x(a),y(b){}
	R2(const R2 &a): x(a.x), y(a.y){}
   	R2(const R2 &a, const R2 &b): x(b.x-a.x), y(b.y-a.y){}
 
 	R2& operator+=(const R2 &P){ x += P.x; y += P.y; return *this; } 
  	R2& operator-=(const R2 &P){ x -= P.x; y -= P.y; return *this; } 
	R2 operator+(const R2 &P)const { return R2(x+P.x, y+P.y); }
  	R2 operator-(const R2 &P)const { return R2(x-P.x, y-P.y); }
  	R  operator,(const R2 &P)const { return x*P.x+y*P.y; }	// produit scalaire
  	R  operator^(const R2 &P)const { return x*P.y-y*P.x; } 	// produit mixte
  	R2 operator*(R c)const {return R2(x*c, y*c);}
  	R2 operator/(R c)const {return R2(x/c, y/c);}
   	R2 operator-()const { return R2(-x, -y); } 
  	const R2& operator+()const { return *this; } //
   	R& operator[](int i) { if(i==0) return x; else if(i==1) return y; else {assert(0);exit(1);} }
    	const R& operator[](int i)const { if(i==0) return x; else if(i==1) return y; else {assert(0);exit(1);} }

  	R2 perp()const { return R2(-y, x); }	// la perpendiculaire
 	R norme()const { return sqrt(x*x+y*y); }

  	// un variable globale dans la classe donc de nom R2::d 
  	// static const int d=2;
};

inline ostream& operator <<(ostream& f, const R2& P){ f << "(" << P.x << "," << P.y << ")" ; return f; }
inline istream& operator >>(istream& f, R2& P){ f >> P.x >> P.y  ; return f; }
inline R2 operator*(R c, const R2 &P){ return P*c; } 
#endif
