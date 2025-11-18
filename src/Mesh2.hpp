#ifndef Mesh2_hpp
#define Mesh2_hpp
#include "R2.hpp"
#include <fstream>
using namespace std;

// class Label:   	lab
// class Vertex2: 	lab  x  y  
// class Triangle: 	lab  *vertices[3]  area
// class Seg:		lab  *vertices[2]  l
// class Mesh2: 	nt nv nbe, area peri, *v *t *be,  

/**************************** class Label *****************************/
class Label{
    public:
	int lab; 
	Label(int r=0):lab(r){}
	int onGamma()const{ return lab; } 
};
inline ostream& operator<<(ostream& f, const Label& r){ f << r.lab; return f; }
inline istream& operator>>(istream& f, Label& r){ f >> r.lab; return f; }


/**************************** class Vertex2 *****************************/
class Vertex2: public R2, public Label{
    public:
	Vertex2(): R2(), Label(){};
	Vertex2(R2 P, int r=0): R2(P), Label(r){}
};
inline ostream& operator <<(ostream& f, const Vertex2& v){ f << (Label&) v << " " << (R2) v; return f; }
inline istream& operator >>(istream& f, Vertex2& v)	 { f >> (R2&) v >> (Label&) v ; return f; }


/*************************** class Triangle *****************************/
class Triangle: public Label{
	Vertex2* v[3];
	R area;  
   public:
	//static const int nv=3; 
	Triangle():area(){}; 
	void init(Vertex2* v0,int i0,int i1,int i2,int r);
	Vertex2& operator[](int i)const{ assert(i>=0 && i <3); return *v[i];} // return the ith vertex
	//R2* pv()const{ return *v;} // return a pointer to the 1st vertex
 	R mesure()const{return area;}  	
};


/*************************** class Seg *****************************/
class Seg: public Label {
	Vertex2 *v[2]; 
	R l; 
    public:
	//static const int nv=2;
	Seg():l(){}; 
	Vertex2& operator[](int i)const{ assert(i>=0 && i <2); return *v[i];} 
	void init(Vertex2 *v0,int *iv,int r);
	R2 operator()(R Phat) const { 
		const R2 &A =*v[0];
		const R2 &B =*v[1];
		return (1-Phat)* A + Phat *B ;
	}
	R mesure() const {return l;}
};


/*************************** class Mesh2 *****************************/
class Mesh2{
    public:
	int nv, nt, nbe;
	R area,peri; 	
	Vertex2 *v; 	
	Triangle *t;	
	Seg *be;	

	Mesh2(const char *filename); 

	int CheckV(int i) const { assert(i>=0 && i < nv); return i; } 
	int CheckT(int i) const { assert(i>=0 && i < nt); return i; }
	int CheckBE(int i) const { assert(i>=0 && i < nbe); return i; }

	// if we have Mesh2 m;
	Triangle& operator[](int i)const{ return t[CheckT(i)]; } // m[0] return the 0th triangle
	Vertex2& operator()(int i)const { return v[CheckV(i)]; } // m(2) return the 2th vertex
	Seg& BE(int i)const{ return be[CheckBE(i)]; }		 // BE(5) return the 5th segment

	// return the number of the input triangle/vertex/segment
	int operator()(const Triangle & tt) const { return CheckT(&tt - t); } 
	int operator()(const Triangle * tt) const { return CheckT( tt - t); } 
	int operator()(const Vertex2 & vv)  const { return CheckV(&vv - v); }
	int operator()(const Vertex2 * vv)  const { return CheckV( vv - v); }
	int operator()(const Seg & s) const { return CheckBE(&s - be); }
	int operator()(const Seg * s) const { return CheckBE( s - be); }

	// return the number of the jth vertex of the ith triangle
	int operator()(int it,int j) const { return (*this)(t[it][j]); }

	~Mesh2(){ delete[] v; delete[] t; delete[] be; }
};

/* extension fct for triangle
R2 Edge(int i)const{ assert(i>=0 && i <3); return R2(*v[(i+1)%3],*v[(i+2)%3]);} //
R  lenEdge(int i)const{ assert(i>=0 && i <3); R2 E=Edge(i);return sqrt((E,E));} // 
R2 H(int i)const{ assert(i>=0 && i <3);	R2 E=Edge(i);return E.perp()*(2*area)/pow(lenEdge(i),2);} // 
void Gradlambda(R2 *GradL)const{ // 
	GradL[1]= H(1);
	GradL[2]= H(2);
	GradL[0]=-GradL[1]-GradL[2];
}
R2 operator()(const R2& Phat) const { //
	const R2 &A =*v[0];
	const R2 &B =*v[1];
	const R2 &C =*v[2];
	return (1-Phat[0]-Phat[1])* A + Phat[0] *B +Phat[1]*C ;
}
*/
#endif
