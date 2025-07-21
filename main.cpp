#include <iostream> 
#include <iomanip>
#include <vector>
#include <ctime>
//#include "aire.hpp"  // include R2, Diff
#include "aire2.hpp" // include R2
#include "Mesh2.hpp" // include R2
//#include "GC.hpp"
using namespace std;

void out(R* data, int n, string name){ // print an array, for test 
	if(n>50) n=50;
	for(int i=0; i<n; ++i){
		cout << name << i << "=" ;
		if(i<10) cout << " ";
		cout << setw(12) << data[i] << "   ";
		if((i+1)%4 == 0) cout << endl;
	}
	cout << endl;
}

void out(const vector<R> data, int n, string name){ // print a vector, for test 
	if(n>50) n=50;
		for(int i=0; i<n; ++i){
			cout << name << i << "=" ;
			if(i<10) cout << " ";
			cout << setw(12) << data[i] << "   ";
			if((i+1)%4 == 0) cout << endl;
		}
	cout << endl;
}

void out(R* data, int n, int m){ // print a matrix, for test
	for(int i=0; i<n; ++i){
	    for(int j=0; j<m; ++j){
		cout << setw(12) << data[i*n+j] ;
		if((j+1)%5 == 0 && (j+1)!=m) cout << endl << "   ";
		if(i>9) cout << " ";
	    }
            cout << endl;
	}
}

/******************* Q5 **********************/
R J(Mesh2* m, const vector<R>& z){
	assert( m->nv == z.size() );
	R area = 0;
	for(int i = 0; i < m->nt; ++i){ // for every triangle m->t[i]
		int n0 = (*m)(i,0); 	// find the number of the vertex for the ith triangle
		int n1 = (*m)(i,1);     // for example: 0-14 in c10.msh 	
		int n2 = (*m)(i,2);      
		R zz[3] = { z[n0], z[n1], z[n2] }; // find the correspond z
		area += aire(m->t[i][0], m->t[i][1], m->t[i][2], zz);     
	}
	return area;
}
/*
/******************* Q6 **********************/
R dJ(Mesh2* m, const vector<R>& z, R* d){ // d: d0,d1,...,dn 
	assert( m->nv == z.size() );
	R area = 0;
	for(int i=0; i<m->nv; ++i ) d[i] = 0.;
	for(int i=0; i<m->nt; ++i){ 	
		int n0 = (*m)(i,0); 
		int n1 = (*m)(i,1);	
		int n2 = (*m)(i,2);      
		R zz[3] = {z[n0], z[n1], z[n2]};
		R tmp[3]; // tmp: d(n0),d(n1),d(n2)	
		area += daire(m->t[i][0], m->t[i][1], m->t[i][2], zz, tmp);
		//cout << &m->t[i][0] << " " << &m->t[i][0]+1 << " " << &m->t[i][1] << endl;
		if( !(m->v+n0)->onGamma() ) d[n0] += tmp[0]; // if (n0) vertex not on gamma
		if( !(m->v+n1)->onGamma() ) d[n1] += tmp[1];
		if( !(m->v+n2)->onGamma() ) d[n2] += tmp[2];
	}	
	return area;
}

/******************* Q7 **********************/
R norme(R* dz, int n){
	R sum = 0.;
	for(int i=0; i<n; ++i)	sum += dz[i]*dz[i];
	return sqrt(sum/n);
}

int descente(Mesh2* m, vector<R>& z, int nbitermax, R eps){

	int n = z.size();
	assert( m->nv == z.size() );
	R rho = 0.2;
	printf("rho=%f, eps=%f \n", rho, eps);
	R d[n];	  dJ(m,z,d);	R Norme = norme(d,n);  printf("initial Norme=%f \n", Norme);
      	if( Norme > eps){ 
	 	for(int iter=0; iter<nbitermax; ++iter) {
			for(int i=0; i<n; ++i){
				z[i] -=  rho*d[i];
			}
			dJ(m,z,d);   Norme = norme(d,n);
			printf("-%d/%d  Norme=%f \n", iter, nbitermax, Norme);
	    		if (Norme < eps ) break;
	   	}
    	}
	return (Norme < eps );
}

/******************* Q8 **********************/
R* ddJ(Mesh2* m, const vector<R>& z, R* u, R* du){ 
	int n = m->nv; 
	assert( n == z.size() );
        for(int k=0; k<n; ++k) du[k]=0.; 
	//R G[n*n];
	//for(int i=0; i<n*n; ++i) G[i] = 0.;
	for(int i=0; i<m->nt; ++i){ 	
		int n0 = (*m)(i,0); 
		int n1 = (*m)(i,1);	
		int n2 = (*m)(i,2); 
		R zz[3]= { z[n0], z[n1], z[n2] };
		R d[3], dd[3*3];	
		ddaire(m->t[i][0], m->t[i][1], m->t[i][2], zz, d, dd);

		/*if(!(m->v+n0)->onGamma()) G[n0*n+n0] += dd[0]; // add d00 into grad2[n0,n0]
		if(!(m->v+n1)->onGamma()) G[n1*n+n1] += dd[4]; // add d11 into [n1,n1]
		if(!(m->v+n2)->onGamma()) G[n2*n+n2] += dd[8]; // add d22 into [n2,n2]
		if(!(m->v+n0)->onGamma() && !(m->v+n1)->onGamma()) G[n1*n+n0] = (G[n0*n+n1] += dd[1]); // add d01 into [n0,n1]
		if(!(m->v+n0)->onGamma() && !(m->v+n2)->onGamma()) G[n2*n+n0] = (G[n0*n+n2] += dd[2]); 
		if(!(m->v+n1)->onGamma() && !(m->v+n2)->onGamma()) G[n2*n+n1] = (G[n1*n+n2] += dd[5]);
	}
	for(int i=0; i<n; ++i) du[i]=0.;
    	for(int i=0; i<n; ++i)
            for(int j=0; j<n; ++j)
             	du[i] += G[i*n+j]*u[j];*/
    
    if(!(m->v+n0)->onGamma()) du[n0]+=dd[0]*u[n0]+dd[1]*u[n1]+dd[2]*u[n2];
    if(!(m->v+n1)->onGamma()) du[n1]+=dd[3]*u[n0]+dd[4]*u[n1]+dd[5]*u[n2];
    if(!(m->v+n2)->onGamma()) du[n2]+=dd[6]*u[n0]+dd[7]*u[n1]+dd[8]*u[n2];

 }

	return du;
}

// the same just ajoute G
R* ddJ(Mesh2* m, const vector<R>& z, R* G, R* u, R* du){ // matrix G = grad2j(z)，for test
	int n = m->nv; 
	assert( n == z.size() );
	for(int i=0; i<n*n; ++i) G[i] = 0.;
	for(int i=0; i<m->nt; ++i){ 	
		int n0 = (*m)(i,0); 
		int n1 = (*m)(i,1);	
		int n2 = (*m)(i,2); 
		R zz[3]= { z[n0], z[n1], z[n2] };
		R d[3], dd[3*3];	
		ddaire(m->t[i][0], m->t[i][1], m->t[i][2], zz, d, dd);

		if(!(m->v+n0)->onGamma()) G[n0*n+n0] += dd[0]; 
		if(!(m->v+n1)->onGamma()) G[n1*n+n1] += dd[4]; 
		if(!(m->v+n2)->onGamma()) G[n2*n+n2] += dd[8]; 
		if(!(m->v+n0)->onGamma() && !(m->v+n1)->onGamma()) G[n1*n+n0] = (G[n0*n+n1] += dd[1]); 
		if(!(m->v+n0)->onGamma() && !(m->v+n2)->onGamma()) G[n2*n+n0] = (G[n0*n+n2] += dd[2]); 
		if(!(m->v+n1)->onGamma() && !(m->v+n2)->onGamma()) G[n2*n+n1] = (G[n1*n+n2] += dd[5]);
	}
	for(int i=0; i<n; ++i) du[i]=0.;
    	for(int i=0; i<n; ++i)
            for(int j=0; j<n; ++j)
             	du[i] += G[i*n+j]*u[j];
	return du;
}

/******************* Q10 **********************/
double* comblineaire(int n, double a, double* A,
                     double b, double* B,double* R ){ // return R = a*A+b*B
    for(int i=0;i< n;++i)
        R[i] = a*A[i] + b*B[i];
    return R;
}

double sdot(int n,double *A, double *B){ // return <A,B>
    double s=0.;
    for(int i=0; i<n; ++i)
        s += A[i]*B[i];
    return s;
}

double* Prod_Id_Vec(void *a, int n, int m, double *x, double *ax) { // return (ax = In * x)
    	assert(n==m);
    	for(int i=0; i<n; ++i)  ax[i]=x[i];
    	return ax;
}
typedef double* (*T_Prod_Mat_Vec)(void *, int ,int , double* x, double* ax); 
// 
int GradienConjugue( // int n,
                    Mesh2* m, const vector<R>& z,//  T_Prod_Mat_Vec PA, void *a, 		// fonction et pointeur data pour A
                    //T_Prod_Mat_Vec PC, void *c,  // fonction et pointeur data pour C
                    double* b, // second membre
                    double* x, // solution qui contient une initialisation
                    int nbitermax, double eps){
    int n = m->nv;
    double * G  = (double*) malloc(sizeof(double)*n);
    double * H  = (double*) malloc(sizeof(double)*n);
    double * CG = (double*) malloc(sizeof(double)*n);
    double * AH = CG;
    double * Ax = CG;
    double Cgg, Cggp, rho, gamma ;
    comblineaire(n,1.,ddJ(m,z,x,Ax),-1,b,G); //  * comblineaire(n,1.,PA(a,n,n,x,Ax),-1,b,G);	// G0  = (Ax=A*x)-b
   // CG0 = C*G0
    for(int i=0; i<n; ++i) CG[i] = G[i];    
	for(int i=0; i<n; ++i) H[i] = - CG[i];   // H0  = -CG0
    Cgg= sdot(n,G,CG);			     // cgg = (G0,G0)c
    if( Cgg > eps){
        for(int iter=0; iter<nbitermax; ++iter){
            ddJ(m,z,H,AH);		     // 改 PA(a,n,n,H,AH);    				// AH  = A*H
            rho = - sdot(n,G,H)/sdot(n,H,AH);// rho = -(Gi,Hi) / (Hi,AHi)
            comblineaire(n,1.,x,rho,H,x);    // xi+1 = xi+rho*Hi	
            comblineaire(n,1.,G,rho,AH,G);   // Gi+1 = Gi+rho*AHi
            Cggp = Cgg;			     // past Cgg
   	 for(int i=0; i<n; ++i) CG[i] = G[i];		     // CGi+1 = C*Gi+1
            Cgg =sdot(n,G,CG);		     // Cgg = (Gi+1,Gi+1)c
	    //assert(Cggp > Cgg);  // 添加
            gamma = Cgg / Cggp;		     // gamma = (Gi+1,Gi+1)c / (Gi,Gi)c
            comblineaire(n,-1.,CG,gamma,H,H);// Hi+1 = -CGi+1
            cout<<"cgg="<<Cgg<<"eps="<<eps<<"gamma="<<gamma<<"rho"<< rho<<endl;
            //printf("  -%d/%d  cgg= %f, eps=%f, gamma=%f, rho=%f\n",iter,nbitermax,Cgg,eps,gamma,rho);
            if( Cgg < eps) break;
        }
    }
    free(G);
    free(H);
    free(CG);
    return Cgg < eps;
}


/******************* Q11 **********************/
int newton(Mesh2* m, vector<R>& z, int nbitermax, R eps){
	int n = z.size();
	assert( m->nv == z.size() );
	 
	R w[n];	
	R d[n];	 dJ(m,z,d);   R Norme = norme(d,n); 
       cout<<"Norme="<<Norme<<"eps="<<eps<<endl; 
       //printf("initial Norme=%d, eps=%d\n", Norme, eps);
	if(Norme > eps){
  		for(int iter=0; iter<nbitermax; ++iter) {
			for(int i=0; i<n; ++i) w[i] =0.;
			assert(GradienConjugue(m,z,d,w,1000,1e-15)); // find w 			
     			for(int i=0; i<n; ++i)
				z[i] -=  w[i]; // wi
			dJ(m,z,d);  Norme = norme(d,n);

			//out(d,n,"b"); cout << endl;
			out(w,n,"w"); cout << endl;
			out(z,n,"z"); cout << endl;
			
			out(d,n,"d"); cout << endl;
                        cout<<iter<<" / "<<nbitermax<<" Norme="<<Norme<<"  eps="<<eps<<endl; 
			//printf("-%d/%d  Norme=%f, eps=%f\n", iter, nbitermax, Norme, eps);
     			if (Norme < eps) break;
   		}
	}
	return (Norme < eps );
}

int main(int argc, char **argv){
	Mesh2 Th(argc>1  ? argv[1] : "c50.msh");
	time_t t1,t2;
	R2 A[3] = { R2(1,0), R2(0,2), R2(0,0)};
	R z[3] = {0,0,3}, d[3], dd[3*3];

	cout << "******* test of Q1 ********" << endl;
	cout << "Aire= " << aire(A[0], A[1], A[2], z) << endl << endl; 
 

	cout << "******* test of Q2 ********" << endl;
	cout << "Aire= " << daire(A[0], A[1], A[2], z,d) << endl; 
 	for(int i=0; i<3; ++i)
		cout << "d" << i << "  = " << d[i] << endl; 
	cout << endl;


	cout << "******* test of Q3 ********" << endl;
	cout << "Aire= " << ddaire(A[0], A[1], A[2], z,d,dd) << endl;
 	for(int i=0; i<3; ++i)
		cout << "d" << i << "  = " << d[i] << endl; 
	out(dd,3,3);
	cout << endl;


	cout << "******* test of Q5 ********" << endl;	
	const int N = Th.nv;
	vector<R> v;
	for(int i = 0; i < N; ++i)
		v.push_back(i*10./N); // v=0.1,0.2,0.3,...
	cout << "Aire= " << J(&Th,v) << endl;
	cout << endl;	


	cout << "******* test of Q6 ********" << endl;
	R dj[N];		
	cout << "Aire= " << dJ(&Th,v,dj) << endl;	
	out(dj, N, "dj");
	cout << endl; 

	
	cout << "******* test of Q7 ********" << endl;
	vector<R> z1=v;
	cout << "initial value" << endl;
	out(z1,N,"z");
	cout << endl;
       // Les fonction test
       /*for(int k=0; k< Th.nv; ++k){

        //z1[k]=cos(Th(k).x)+0.1; 
        //z1[k]=atan(Th(k).y/Th(k).x)+0.1;//Hélicoide
        if(Th(k).lab==0)

        //z1[k]=cos(Th(k).x);//
        //z1[k]=atan(Th(k).y/Th(k).x)//Hélicoide 
         } */
	time(&t1);  assert(descente(&Th,z1,10000,1e-8));  time(&t2);
	cout << endl << "by descente algorithm, result is" << endl;
	out(z1,N,"z");
	/*
        //double err=z1[0]-cos(Th(0).x);
        //double err=z1[0]-atan(Th(k).y/Th(k).x)+0.1;
       
       	for(int k=0; k< Th.nv; ++k){

         //double diff=z1[k]-cos(Th(k).x);
         //double diff=z1[K]-atan(Th(k).y/Th(k).x)+0.1;
         if(err<diff) err=diff;
         } 
        cout << "we spent " << t2-t1 << " second" <<" err=" << err << endl << endl;*/
        cout << "we spent " << t2-t1 << " second" << endl << endl;
	R zero[N];
	dJ(&Th,z1,zero);
	cout << "then, par Descendre grad is" << endl;
	out(zero,N,"G");
	cout << endl;
        
       ofstream gp("Descendre.gp");
        for(int k=0; k< Th.nt; ++k){
       	gp << Th.t[k][0].x << " " << Th.t[k][0].y << " " << z1[Th(k,0)] <<endl;
        gp << Th.t[k][1].x << " " << Th.t[k][1].y << " " << z1[Th(k,1)] <<endl;
 	gp << Th.t[k][2].x << " " << Th.t[k][2].y << " " << z1[Th(k,2)] <<endl;
        gp << endl;
      }

	cout << "******* test of Q8 ********" << endl;
	R du[N];		
	R u[N];  for(int i=0; i<N; ++i) u[i]=i; // u=1,2,3...
	
	cout << "test of fct grad2" << endl;
	R G[N*N]; // conserve la matrix grad2j(z), for test
	ddJ(&Th,v,G,u,du);

	for(int i=0; i<N; ++i){
	    for(int j=0; j<N; ++j)
		cout << G[i*N+j] << " ";
	    cout << endl;
	}
	cout << endl;

	ddJ(&Th,v,u,du);	
	out(du,N,"du");
	cout << endl;


	cout << "******* test of Q10 ********" << endl;
	R b[N], x[N], xx[N];  
   	for(int i=0; i<N; ++i){
		xx[i] = i+1.;  	// xx=1,2,3,...
		x[i]  = 0.;   
	}

	cout << "first, suppose solution xx is" << endl;
	out(xx,N,"xx");
	cout << endl;

	cout << "and, initial x is" << endl;
	out(x,N,"x");
	cout << endl;

   	ddJ(&Th,v,xx,b);    // b=A(m,z)*xx
	cout << "then, b=grad2(Th,v)*xx" << endl;
	out(b,N,"b");
	cout << endl; 

	assert(GradienConjugue(&Th,v,b,x,100,1e-2)); 
	cout << endl << "use fct GC, solution is" << endl;
	out(x,N,"x");
	cout << endl;

  	ddJ(&Th,v,x,b);    // b=A(m,z)*x
	cout << "then b=grad2(Th,v)*x=" << endl;
	out(b,N,"b");
	cout << endl;  


	cout << "******* test of Q11 ********" << endl;
	vector<R> z2=v; 
       //Les fonction test 
       /* for(int k=0; k< Th.nv; ++k){
        z2[k]=cos(Th(k).x)+0.1; 
        if(Th(k).lab==0)
          //z2[k]=cos(Th(k).x);
          //z2[k]=atan(Th(k).y/Th(k).x)+0.1;//Hélicoide
         } */
       
        ofstream gp2("Min1.gp");
        for(int k=0; k< Th.nt; ++k){
       	gp2 << Th.t[k][0].x << " " << Th.t[k][0].y << " " << z2[Th(k,0)] <<endl;
        gp2 << Th.t[k][1].x << " " << Th.t[k][1].y << " " << z2[Th(k,1)] <<endl;
 	gp2 << Th.t[k][2].x << " " << Th.t[k][2].y << " " << z2[Th(k,2)] <<endl;
        gp2 << endl;
      }
	//for(int i=0; i<N; ++i) z2[i] = z1[i]*1.3; // use the descente result as the initial value
	cout << "initial value" << endl;
	out(z2,N,"z");
	cout << endl;

	cout << "by newton algorithm, result is" << endl;
	time(&t1); assert( newton(&Th,z2,80,1e-8));  time(&t2);
	out(z2,N,"z");
        /*//double err=z2[0]-cos(Th(0).x);
        //double err=z2[0]-atan(Th(k).y/Th(k).x)+0.1;//Hélicoide
       	for(int k=0; k< Th.nv; ++k){
         //double diff=z2[k]-cos(Th(k).x);
         //double diff=z2[k]-atan(Th(k).y/Th(k).x)+0.1;//Hélicoide
         if(err<diff) err=diff;
         } 
	cout << "we spent " << t2-t1 << " second" <<" erreur= "<< err<< endl << endl;*/	
        cout << "we spent " << t2-t1 << " second" << endl << endl;
	R zero2[N];
	dJ(&Th,z2,zero2);
	cout << "then, gradJz is" << endl;
	out(zero2,N,"G");
	cout << endl;
        
      
	return 0;
}
