template <class T,int N=3> 
struct Diff{ // val,d[0],d[1],d[2]
  	T val, d[N]; 
  	
  	Diff(): val(T(0)){ for(int i=0; i<N; ++i) d[i]=T(0); }
  	Diff(T a): val(a){ for(int i=0; i<N; ++i) d[i]=T(0); }
  	Diff(T a, T da): val(a){ d[0]=da; for(int i=1; i<N; ++i) d[i]=T(0); }
	Diff(T a, T d0, T d1, T d2): val(a){ d[0]=d0; d[1]=d1; d[2]=d2; for(int i=3; i<N; ++i) d[i]=T(0); }
};


template <class T,int N> // x+y
Diff<T,N> operator+(const Diff<T,N>& x, const Diff<T,N>& y){
  	Diff<T,N> r(x.val+y.val);
  	for(int i=0; i<N; ++i)	     
    		r.d[i] = x.d[i]+y.d[i]; 
  	return r; 
}

template <class T,int N> // x+c
Diff<T,N> operator+(const Diff<T,N>& x, R c){
  	Diff<T,N> r(x.val+c);
  	for(int i=0; i<N; ++i)	     
    		r.d[i] = x.d[i]; 
  	return r; 
}

template <class T,int N> // x-y
Diff<T,N> operator-(const Diff<T,N>& x, const Diff<T,N>& y){
  	Diff<T,N> r(x.val-y.val);
  	for(int i=0; i<N; ++i)	     
    		r.d[i] = x.d[i]-y.d[i]; 
  	return r; 
}

template <class T,int N> // x-c
Diff<T,N> operator-(const Diff<T,N>& x, R c){
   	return x+(-c); 
}

template <class T,int N> // x*y
Diff<T,N> operator*(const Diff<T,N>& x, const Diff<T,N>& y){
  	Diff<T,N> r(x.val*y.val);
  	for(int i=0; i<N; ++i)	     
    		r.d[i] = x.d[i]*y.val + x.val*y.d[i]; 
  	return r; 
}

template <class T,int N> // x*c
Diff<T,N> operator*(const Diff<T,N>& x, R c){ 
  	Diff<T,N> r(x.val*c);
  	for(int i=0; i<N; ++i)	     
    		r.d[i] = x.d[i]*c ; 
  	return r; 
}

template <class T,int N> // c*x
Diff<T,N> operator*(R c, const Diff<T,N>& x){ return x*c; } 

template <class T,int N> // x/y
Diff<T,N> operator/(const Diff<T,N>& x, const Diff<T,N>& y){
  	Diff<T,N> r(x.val/y.val);
  	for(int i=0; i<N; ++i)	     
    		r.d[i] = (x.d[i]*y.val - x.val*y.d[i]) / (y.val*y.val); 
  	return r; 
}

template <class T,int N> // x/c
Diff<T,N> operator/(const Diff<T,N>& x, R c){
  	Diff<T,N> r(x.val/c);
  	for(int i=0; i<N; ++i)	     
    		r.d[i] = x.d[i]/c; 
  	return r; 
}

template <class T,int N> // c/x
Diff<T,N> operator/(R c, const Diff<T,N>& x){
  	Diff<T,N> r(c/x.val);
  	for(int i=0; i<N; ++i)	     
    		r.d[i] = -x.d[i]*c/(x.val*x.val); 
  	return r; 
}
	
template <class T,int N> // sqrt(x)
Diff<T,N> sqrt(const Diff<T,N>& x){
  	Diff<T,N> r(sqrt(x.val));
  	for(int i=0; i<N;++i)	     
    		r.d[i] = x.d[i]*0.5/r.val; 
  	return r; 
}

template <class T,int N> // cout << x
ostream& operator<<(ostream& f, const Diff<T,N>& a){ 
	f << a.val << " ( "; 
	for(int i=0; i<N-1; ++i)
		f << "d" << i << "=" << a.d[i] << ", ";
	f << "d" << N-1 << "=" << a.d[N-1] << " )";
	return f;
} 

