#include"R2.hpp"
#include"Diff.hpp"
typedef Diff<R> DR;
typedef Diff<DR> DDR; 

/******************* Q1 **********************/
R aire(R2 *A, R *z){ 
	R u = (A[1]-A[0])[1]*(z[2]-z[0])-(A[2]-A[0])[1]*(z[1]-z[0]);
	R v = (A[2]-A[0])[0]*(z[1]-z[0])-(A[1]-A[0])[0]*(z[2]-z[0]);
	R w = (A[1]-A[0])[0]*(A[2]-A[0])[1]-(A[1]-A[0])[1]*(A[2]-A[0])[0];
	return 0.5*sqrt(u*u+v*v+w*w);
}

/******************* Q2 **********************/
R daire(R2 *A, R *z, R* d){ // d: d0,d1,d2
	DR z0(z[0],1,0,0);
	DR z1(z[1],0,1,0);
	DR z2(z[2],0,0,1);

	DR u = (A[1]-A[0])[1]*(z2-z0)-(A[2]-A[0])[1]*(z1-z0);
	DR v = (A[2]-A[0])[0]*(z1-z0)-(A[1]-A[0])[0]*(z2-z0);
	R w = (A[1]-A[0])[0]*(A[2]-A[0])[1]-(A[1]-A[0])[1]*(A[2]-A[0])[0];

	DR res = 0.5*sqrt(u*u+v*v+w*w);
	d[0] = res.d[0]; d[1] = res.d[1]; d[2] = res.d[2];
	return res.val;
}

/******************* Q3 **********************/
R ddaire(R2 *A, R *z, R* d, R* dd){ // dd: d00 d01 d02, d10 d11 d12, d20 d21 d22
	DR z0(z[0],1,0,0); 
	DR z1(z[1],0,1,0);
	DR z2(z[2],0,0,1);

	DDR zz0(z0, DR(1,0,0,0), DR(0,0,0,0), DR(0,0,0,0)); 
	DDR zz1(z1, DR(0,0,0,0), DR(1,0,0,0), DR(0,0,0,0)); 
	DDR zz2(z2, DR(0,0,0,0), DR(0,0,0,0), DR(1,0,0,0)); 

	DDR u = (A[1]-A[0])[1]*(zz2-zz0)-(A[2]-A[0])[1]*(zz1-zz0);
	DDR v = (A[2]-A[0])[0]*(zz1-zz0)-(A[1]-A[0])[0]*(zz2-zz0);
	R w = (A[1]-A[0])[0]*(A[2]-A[0])[1]-(A[1]-A[0])[1]*(A[2]-A[0])[0];

	// res.val = (area,'z0,'z1,'z2) 		res.val.d[0]='z0
	// res.d[0] = ('z0, 'z0'z0, 'z0'z1, 'z0'z2)	res.d[0].val='z0
	// res.d[1] = ('z1, 'z1'z0, 'z1'z1, 'z1'z2)
	// res.d[2] = ('z2, 'z2'z0, 'z2'z1, 'z2'z2)
	DDR res = 0.5*sqrt(u*u+v*v+w*w); 
	d[0] = res.val.d[0]; 
	d[1] = res.val.d[1]; 
	d[2] = res.val.d[2];

	dd[0] = res.d[0].d[0]; // d00
	dd[4] = res.d[1].d[1]; // d11 
	dd[8] = res.d[2].d[2]; // d22
	dd[3] = (dd[1] = res.d[0].d[1]); // d10=d01
	dd[6] = (dd[2] = res.d[0].d[2]); // d20=d02
	dd[7] = (dd[5] = res.d[1].d[2]); // d21=d12
	/* pour tester
	for(int i=0; i<3; ++i)
		cout << A[i][0] << " " << A[i][1] << " " << z[i] << endl;
	
	cout << dd[0] << "  " << dd[4] << "  " << dd[8] << endl;
	cout << dd[3] << "  " << dd[6] << "  " << dd[7] << endl;
	cout << endl;
	*/
	return res.val.val;
}


