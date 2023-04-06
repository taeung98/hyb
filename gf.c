#include "start.h"
#include <gsl/gsl_integration.h>
#include "inout_data.h"
#include "legendre_transform.h"
#include <sys/stat.h>
#include <sys/types.h>
#define w_n(n, B) ((2*(n)+1)*M_PI/(B))
#define R		1024
#define kmax	1000

double complex w(int, double);
void FT_T(double complex*,int,double complex*,int,double*,double);
void FT_G(double complex*,int,double complex*,int,double*,double);
void IFT_T(double complex*,int,double complex*,int,double*,double);
void IFT_G(double complex*,int,double complex*,int,double*,double*,double);
double E_n(int ,double , double , int,  int, int);

int main() {
    int L,N,amax,D,nei; L=LL; N=NN; amax=AMAX; D=1;
	double a,B,e,dk;	B=BB;
	double ti[L],wi[L],W_n[N],x[L],x_amax[amax],T;
    double complex IGT_G[N],LGT_T[L],LGT_G[L],G_latt[N],DW[N],DT_T[L],GT_T[L],A_GT_T[amax],A_GT_G[amax],A_DT_G[amax],A_DT_T[amax],IGT_T[N],GT_G[L],DT_G[L],IDT_T[N],IDT_G[N],LDT_T[L],LDT_G[L];
//	double complex DD[N],err[L],DOS[R],W[R]; double eat=0.01;
//	double (*E)(int,double, double, int, int, int);
//	E = E_n;
	for(nei=1;nei<5;nei++){
		for(T=1;T<100;T++){
			
	D = 1;a = 1;
	dk = 2*M_PI/a/1000;
	
	//	for(int n=1;n<2;n++){
		char dirname[50];
		sprintf(dirname, "n%d",nei);
		mkdir(dirname, 0700);
//		for(int i=1;i<2;i++){
//			double T = t*i;
	
//	char fn[1024]; sprintf(fn, "%f", t*i);
	
	//sum of k
	
	double TT = T*0.03;
	for(int n=0;n<N;n++){
		double complex g; g=0;
		for(int k=0;k<kmax;k++){    //Q.FBZ-> -pi,pi가 겹치니까 하나 빼줘야 하나?
            e = E_n(D,TT,a,nei,k,1);
            g += 1/( w_n(n,B)*1i - e );
        }
        g *= dk;
        W_n[n] = w_n(n,B);
        G_latt[n] = g;                                  // sum G(k,iw)_latt
        DW[n] = 1i*W_n[n] - 1/g;                        // hyb
    }

	//save_data("G_latt", W_n, G_latt, N, dirname,T);
	//save_data("DW",W_n, DW, N, dirname,T);
	
	//Gaussian Quadrature node & weights
	gsl_integration_glfixed_table *w;
	w = gsl_integration_glfixed_table_alloc(L);
	for(int i=0;i<L;i++){
		x[i] = (B/(L-1))*i;
		gsl_integration_glfixed_point(0,B,i,&ti[i],&wi[i],w);
		//prinf("w_%d = %f, x_%d = %f\n",i,wi[i],i,ti[i]);
	}
	gsl_integration_glfixed_table_free(w);

	//N: before_func_length, L: after_func_length
	FT_T(G_latt,N,GT_T,L,x,B);		//save_data("GT_T",x, GT_T,L,dirname,T);
	FT_G(G_latt,N,GT_G,L,ti,B);		//save_data("GT_G",ti,GT_G,L,dirname,T);
	FT_T(DW,N,DT_T,L,x,B);			//save_data("DT_T",x,DT_T,L,dirname,T);
	FT_G(DW,N,DT_G,L,ti,B);			//save_data("DT_G",ti,DT_G,L,dirname,T);
#if 0
	IFT_T(GT_T,L,IGT_T,N,x,B);		save_data("IGT_T",W_n,IGT_T,N,fn);
	IFT_G(GT_G,L,IGT_G,N,ti,wi,B);	save_data("IGT_G",W_n,IGT_G,N,fn);
	IFT_T(DT_T,L,IDT_T,L,x,B);		save_data("IDT_T",W_n,IDT_T,N,fn);
	IFT_G(DT_G,L,IDT_G,L,ti,wi,B);	save_data("IDT_G",W_n,IDT_G,N,fn);
#endif
	//legendre_transform_traperzoid (b,coeff,af)
	AN_LT_T(GT_T,x_amax,A_GT_T);	save_data("A_GT_T",x_amax,A_GT_T,amax,dirname,TT);
	LT_T(GT_T,x,A_GT_T,LGT_T); 		//save_data("LGT_T",x,LGT_T,L,dirname,T);
	AN_LT_T(DT_T,x_amax,A_DT_T);	save_data("A_DT_T",x_amax,A_DT_T,amax,dirname,TT);
	LT_T(DT_T,x,A_DT_T,LDT_T);		//save_data("LDT_T",x,LGT_T,L,dirname,T);

	//legendre_transform_Gaussian	(b,coeff,node,weihgt,af)
	AN_LT_G(GT_G,x_amax,A_GT_G,ti,wi);	
	save_data("A_GT_G",x_amax,A_GT_G,amax,dirname,TT);
//	LT_G(GT_G,ti,A_GT_G,LGT_G);			save_data("LGT_G",ti,LGT_G,L,fn);
	
	AN_LT_G(DT_G,x_amax,A_DT_G,ti,wi);
	save_data("A_DT_G",x_amax,A_DT_G,amax,dirname,TT);
//	LT_G(DT_G,ti,A_DT_G,LDT_G);			save_data("LDT_G",ti,LDT_G,L,fn);
		}
	}
#if 0
	//Error	
	for(int i=0;i<N;i++){
		err[i] = (G_latt[i]-IGT_G[i]);
	}
	save_err("err",err,N);

	//DOS
	for(int r=0; r<R; r++){
		double complex g2 = 0;
		double complex w_ = (r - 4*0.125*(R-1))*8*pow(R-1,-1) + 1I*ea; 
		for(int k=0; k<kmax; k++){
			e = E(t,a,k,neighbor);
			g2 += dk/( w_ - e);
		W[r] = w_;
		DOS[r] = g2;
		}

	}
	save_data("DOS", W, DOS, R);

	//적분 정밀도 확인	
	//G||DW (iw_n) for au	respectively
	for(int i=0;i<N;i++){
		W_n[i] = w_n(i-511,B);
// 		prinf("%f\n",W_n[i]);
	}
	for(int tau=0;tau<50;tau++){
		for(int i=0;i<N;i++){
			double complex x = 0;
			x = (cos(w_n(i-511,B)*[tau])-1I*sin(w_n(i-511,B)*t[tau]))*DW[i];
			DD[i] = x/B;
		}
		save_d("au",W_n,DD,N,tau);
	}
#endif
	return 0;
}


double complex w(int n, double eta){
    double complex a;
	n = (double)n;
    a= (n-511*0.125*4)*pow(511,-1)*8 + 1I*eta; 
    
    return a;
}

//calculate tight bindign energy
double E_n(int D,double t, double a, int neighbor, int x, int y){
    double e;e=0;
    double dk = 2*M_PI/a/1000; //kmax = 1000

    if(D == 1){
        if(neighbor ==1){
            e = -2*t*cos(-M_PI/a + x*dk);
        }
        else if(neighbor ==2){
            e = -2*t*cos(-M_PI/a + x*dk)-2*(t/2)*cos(-M_PI/a + x*dk);
        }
        else if(neighbor ==3){
            e = -2*t*cos(-M_PI/a + x*dk)-2*(t/2)*cos(2*(-M_PI/a + x*dk))-2*(t/4)*cos(3*(-M_PI/a + x*dk));
        }
		else if(neighbor ==4){
			e = -2*t*cos(-M_PI/a + x*dk)-2*(t/2)*cos(2*(-M_PI/a + x*dk))-2*(t/4)*    cos(3*(-M_PI/a + x*dk)) -2*(t/8)*cos(4*(-M_PI/a + x*dk));
		}
	}
    else if(D == 2){
        if(neighbor == 1){
            e = -2*t*cos(-M_PI/a + x*dk) -2*t*cos(-M_PI/a + y*dk);
        }
        else if(neighbor == 2){
            e = -2*t*cos(-M_PI/a + x*dk)-2*(t/2)*cos(2*(-M_PI/a + x*dk))
            -2*t*cos(-M_PI/a + y*dk)-2*(t/2)*cos(2*(-M_PI/a + y*dk));
        }
        else if(neighbor == 3){
            e = -2*t*cos(-M_PI/a + x*dk)-2*(t/2)*cos(2*(-M_PI/a + x*dk))-2*(t/4)*cos(3*(-M_PI/a + x*dk))
                -2*t*cos(-M_PI/a + y*dk)-2*(t/2)*cos(2*(-M_PI/a + y*dk))-2*(t/4)*cos(3*(-M_PI/a + y*dk));
        }
		else if(neighbor ==4){
			e = -2*t*cos(-M_PI/a + x*dk)-2*(t/2)*cos(2*(-M_PI/a + x*dk))-2*(t/4)*cos(3*(-M_PI/a + x*dk)) - 2*(t/8)*cos(4*(-M_PI/a + x*dk))
                -2*t*cos(-M_PI/a + y*dk)-2*(t/2)*cos(2*(-M_PI/a + y*dk))-2*(t/4)*cos(3*(-M_PI/a + y*dk)) -2*t*(t/8)*cos(4*(-M_PI/a +y*dk));
	
		}
	}
    return e;
}

//Trapozoidal fourier transform
void FT_T(double complex* b,int b_point_num, double complex* af,int af_point_num, double* node,double cycle){
	for(int i=0;i<af_point_num;i++){
		double complex f = 0;
		for(int j=0;j<b_point_num;j++){
			//exp -> Euler formula!!
			f += (cos(w_n(j,cycle)*node[i])-1I*sin(w_n(j,cycle)*node[i]))*b[j];
		}
		af[i] = f/cycle; 
	}
}

//Gaussian fourier ransform
void FT_G(double complex* b, int b_point_num, double complex* af, int af_point_num, double* node, double cycle){
	for(int i=0;i<af_point_num;i++){
		double complex f = 0;
		for(int j=0;j<b_point_num;j++){
			f += (cos(w_n(j,cycle)*node[i])-1I*sin(w_n(j,cycle)*node[i]))*b[j];
		}
		af[i] = f/cycle;
	}
}

//Inverse Fourier Transform_rapozoidal
void IFT_T(double complex* b,int b_point_num, double complex* af,int af_point_num, double* node, double cycle){
	for(int i=0;i<af_point_num;i++){
		double complex x=0;
		for(int j=0;j<b_point_num;j++){
			x += (cos(w_n(i,cycle)*node[j])+1I*sin(w_n(i,cycle)*node[j]))*b[j];
		}
		x -= (cos(w_n(i,cycle)*node[0])+1I*sin(w_n(i,cycle)*node[0]))*b[0];
		af[i] = x*cycle/(b_point_num-1);
	}
}

//Inverse Fourier Transform_Gaussian
void IFT_G(double complex* b,int b_point_num, double complex* af,int af_point_num, double* node,double* weight, double cycle){
	for(int i=0;i<af_point_num;i++){
		double complex x=0;
		for(int j=0;j<b_point_num;j++){
			x += (cos(w_n(i,cycle)*node[j])+1I*sin(w_n(i,cycle)*node[j]))*b[j]*weight[j];
		}
		af[i] = x;
	}
}

