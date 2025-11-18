#include "Mesh2.hpp"	

void Triangle::init(Vertex2* v0,int i0,int i1,int i2,int r){ 
	v[0]=v0+i0; v[1]=v0+i1; v[2]=v0+i2; 
	R2 AB(*v[0],*v[1]);
	R2 AC(*v[0],*v[2]);
	area = (AB^AC)*0.5;
	lab=r;
	assert(area>0);
}

void Seg::init(Vertex2 *v0,int *iv,int r){
	v[0]=v0+iv[0]; v[1]=v0+iv[1];
	R2 AB(*v[0],*v[1]);
	l= AB.norme();
	lab=r;
	assert(l>0);
}

Mesh2::Mesh2(const char * filename): nv(0), nt(0), nbe(0), area(0), peri(0),
				     v(0),t(0),be(0){ 
	ifstream f(filename);
	if(!f) { cerr << "Failed to open"  << endl; exit(1); }

	// input nv,nt,nbe
	f >> nv >> nt >> nbe ;
	cout << "Nv=" << nv << 
		", Nt=" << nt << 
		", Nbe=" << nbe << endl;
	assert(f.good() && nt && nv ) ;

	t = new Triangle [nt];
	v = new Vertex2[nv];
	be = new Seg[nbe];
	assert(t && v);

	int i,iv[3],ir;

	// input (nv) vertex 
	for (i=0; i<nv; i++){
		f >> v[i];
		assert(f.good());
	}

	// input (nt) triangle 
	for (i=0; i<nt; i++){
		// input the number of the 3 vertex &   the lab of the triangle
 		f >> iv[0] >> iv[1] >> iv[2] >> ir;
		assert(f.good() && iv[0]>0 && iv[0]<=nv && iv[1]>0 && 
			iv[1]<=nv && iv[2]>0 && iv[2]<=nv);

		for (int j=0; j<3; ++j) iv[j]--;
		t[i].init(v,iv[0],iv[1],iv[2],ir); // triangle:init()
		area += t[i].mesure(); 
	}

	// input (nbe) segment
	for (i=0; i<nbe; i++){
		f >> iv[0] >> iv[1] >> ir;
		assert(f.good() && iv[0]>0 && iv[0]<=nv && iv[1]>0 && iv[1]<=nv);
		for (int j=0; j<2; ++j) iv[j]--;
		be[i].init(v,iv,ir);
		peri += be[i].mesure();
	}
	//cout << "Area=" << area << ", Peri=" << peri << endl;
}
