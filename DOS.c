#define N_EIGEN 1000 // Number of eigenenergies
#define N_DOS 1000 // Number of energy interval
#define ETA 0.05 // Broadening rate
#define GREEN(e, eigen_e, fermi) ( ETA / ( pow((e) - ((eigen_e) - (fermi)), 2) + (ETA)*(ETA)) )  // Imaginary Green func

#include <math.h>
#include <stdio.h>
#include <string.h>

void save(char*, double*,double*,int);

double MakeDOS(double *eigen_e, double fermi) {
	double e_min = -8, e_max = 8; // 에너지 범위 지정
	double e[N_DOS],dos[N_DOS]; // DOS

	memset(dos, 0, sizeof(dos)); // DOS array 0으로 초기화

	for(int n=0; n<N_DOS; n++) {
		e[n] = e_min + (e_max - e_min) * n / N_DOS; // DOS 구할 e 계산

		for(int i=0; i<N_EIGEN; i++) {
			dos[n] += GREEN(e[n], eigen_e[i], fermi)/M_PI * (2*M_PI)/N_EIGEN; // DOS 계산
		}
			
//		printf("%16.6f\t%16.6f\n", e, dos[n]); // 출력해보고 싶으면 주석 지우기 (왼쪽이 에너지, 오른쪽이 DOS)
	}
	
	double sum=0; // # of states(# of electorns)
	for(int i=0;i<N_DOS;i++)
		sum += dos[i]*(e_max-e_min)/N_DOS/(2*M_PI);
	//printf("%0.16f\n",sum);
	
	//save("DOS",e,dos,N_DOS);
	
	return sum;
}

double optimizeMakeDOS(double *eigen_e, double *fermi) {
  double threshold = 1e-6; // set threshold for optimization
  double diff = 1.0; // initialize difference
  int max_iter = 100; // set maximum number of iterations

  for (int i = 0; i < max_iter && fabs(diff) > threshold; i++) {
    double dos_sum = MakeDOS(eigen_e, *fermi); // calculate DOS sum using current fermi value
    diff = dos_sum - 1.0; // calculate difference from target value of 1
    *fermi += diff; // adjust fermi value based on difference
  }

  return *fermi;
}


int main() {
	double eigen_e[N_EIGEN], fermi; 
	double t = 0.1;
	double count=0;
	double N; double d = 1.5;
	// 테스트 해보려면 아래 주석을 지우고 고유에너지, 페르미 레벨을 임의로 생성하면 됨
	// 실제 사용할 때는 자기가 구한 고유에너지, 페르미 레벨을 집어넣으면 됨
	for(int i=0; i<N_EIGEN; i++) 
		eigen_e[i] = -2.*t*cos(-M_PI+ 2*M_PI* i / N_EIGEN) -2.*(t/2)*cos(2*(-M_PI+2*M_PI*i/N_EIGEN));
	
	fermi = 3; // set initial fermi value
	double Ni,Nf;
	do {
		printf("%f, fermi=%.16f\n",count,fermi);
		Ni = MakeDOS(eigen_e,fermi);
		printf("Ni= %.16f\n",Ni);
		if ( N > 1 ) {
			fermi += d/pow(2,count);
		}
		else if ( N < 1 ) {
			fermi -= d/pow(2,count);
		}
		Nf = MakeDOS(eigen_e,fermi);
		printf("Nf= %.16f\n",Nf);
		printf("Nf-Ni = %.16f\n",fabs(Nf-Ni));
		count += 1;
	} while(fabs(Nf-1) > 0.001);
	

	
	return 0;	
}

void save(char *head, double* x, double* y, int LN){
    char fname[1024];
    sprintf(fname, "%s.txt", head);
    FILE* fp = fopen(fname,"w");
    for(int i=0; i<LN; i++){
        fprintf(fp, "%lf\t%lf\n", x[i], y[i] );
    }
    fclose(fp);

}

