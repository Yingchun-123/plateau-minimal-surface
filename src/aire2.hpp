#include"R2.hpp"
class R3: public R2 {
    public:
	R z;

	R3(): R2(), z(0.){};
	R3(R x, R y, R z): R2(x,y), z(z){};
	R3(R2 P, R z=0): R2(P), z(z){};
 	R3(const R3 &a, const R3 &b): R2(b.x-a.x, b.y-a.y), z(b.z-a.z){}
	R norme()const { return sqrt(x*x+y*y+z*z); }
	R3 operator^(const R3 &P)const { return R3(y*P.z-z*P.y, z*P.x-x*P.z, x*P.y-y*P.x); }
   	//R& operator[](int i) { if(i==0) return x; else if(i==1) return y; else if(i==2) return z; else {assert(0);exit(1);} }
};

/******************* Q1 **********************/
R aire(const R2& A0, const R2& A1, const R2& A2, R *z){ 
	R3 a(A0,z[0]);
	R3 b(A1,z[1]);
	R3 c(A2,z[2]);
	R3 ab(a,b);
	R3 ac(a,c);
	return 0.5*(ab^ac).norme();
}

/******************* Q2 **********************/
R daire(const R2& A0, const R2& A1, const R2& A2, R *z, R* d){ // d: d0,d1,d2
	R3 a(A0,z[0]);
	R3 b(A1,z[1]);
	R3 c(A2,z[2]);
	//cout << A0 << " " << A1 << " " << A2 << endl;
	R3 ab(a,b);    
	R3 ac(a,c);    
	R3 cross = ab^ac;  
	R s = 0.5*cross.norme();
	
	d[0] = cross.x * (A2.y - A1.y) + cross.y * ((A1.x - A2.x)); 
	d[1] = cross.x * (A0.y - A2.y) + cross.y * ((A2.x - A0.x));
	d[2] = cross.x * (A1.y - A0.y) + cross.y * ((A0.x - A1.x));
	d[0] /= 4*s; d[1] /= 4*s; d[2] /= 4*s;
	return s;
}

/******************* Q3 **********************/
R ddaire(const R2& A0, const R2& A1, const R2& A2, R *z, R* d, R* dd){ // dd: d00 d01 d02, d10 d11 d12, d20 d21 d22
	R3 a(A0,z[0]);
	R3 b(A1,z[1]);
	R3 c(A2,z[2]);

	R3 ab(a,b);    
	R3 ac(a,c);    
	R3 cross = ab^ac; 
 
	R u = cross.x;  // u=(y1-y0)(z2-z0) - (y2-y0)(z1-z0)
	R v = cross.y;  // v=(x2-x0)(z1-z0) - (x1-x0)(z2-z0)
	R w = cross.z;  // w=(x1-x0)(y2-y0) - (y1-y0)(x2-x0)

	R s = 0.5*cross.norme(); // s=0.5*sqrt(u^2+v^2+w^2)

	R u0 = A2.y - A1.y;  // u partial z0 = y2-y1
	R u1 = A0.y - A2.y;  // u partial z1 = y0-y2
	R u2 = A1.y - A0.y;  // u partial z2 = y1-y0
	R v0 = A1.x - A2.x;  // v partial z0 = x1-x2
	R v1 = A2.x - A0.x;  // v partial z1 = x2-x0
	R v2 = A0.x - A1.x;  // v partial z2 = x0-x1

	d[0] = u * u0 + v * v0;   d[0] /= 4*s;
	d[1] = u * u1 + v * v1;   d[1] /= 4*s;
	d[2] = u * u2 + v * v2;   d[2] /= 4*s;
 
	dd[0] = u0*u0 + v0*v0 - 4*d[0]*d[0]; dd[0] /= 4*s; // d00
	dd[1] = u0*u1 + v0*v1 - 4*d[0]*d[1]; dd[1] /= 4*s; // d01
	dd[2] = u0*u2 + v0*v2 - 4*d[0]*d[2]; dd[2] /= 4*s; // d02
	dd[3] = dd[1]; // d10 = d01
	dd[4] = u1*u1 + v1*v1 - 4*d[1]*d[1]; dd[4] /= 4*s; // d11
	dd[5] = u1*u2 + v1*v2 - 4*d[1]*d[2]; dd[5] /= 4*s; // d12
	dd[6] = dd[2]; // d20 = d02
	dd[7] = dd[5]; // d21 = d12
	dd[8] = u2*u2 + v2*v2 - 4*d[2]*d[2]; dd[8] /= 4*s; // d11

	return s;
}

