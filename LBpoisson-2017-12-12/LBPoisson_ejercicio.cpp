#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define N 3
//#define N 7
#define Lx 128
#define tolerance 1.0E-12

using namespace std;

double dt = 1.0;
double cs2 = 1.0/3.0;
//double cs2 = 0.69795332201968308823840905538933;

double potential(int i){
	double g;
	g=0.0;
	return (g);
}
double charge_density(int i){
	return ( 0.0 );
}

double fi[Lx][N], fip[Lx][N]; double taoi, Phi[Lx], Phi_old[Lx], rho[Lx]; int v[N]; double chi_feq[N], chi_S[N]; double w[N];
void CalculoMacro(int t);
double GetPhi(int i){return (Phi[i]);};
double GetRho(int i){return (rho[i]);};
void Inicio(void);
void Evolucion(int t);

void Inicio(void){
	
	taoi=1.0;

	w[0]=2.0/3.0; v[0]=0;
	w[1]=1.0/6.0; v[1]=1;
	w[2]=1.0/6.0; v[2]=-1;

/*
	w[0] = 0.47666988658920744233176499640793;v[0]=0;
	w[1] = 0.23391473782682471337286677188259;v[1]=1;
	w[2] = 0.23391473782682471337286677188259;v[2]=-1;
	w[3] = 0.026938189344825451680983762681790;v[3]=2;
	w[4] = 0.026938189344825451680983762681790;v[4]=-2;
	w[5] = 0.00081212953374611378026696723165760;v[5]=3;
	w[6] = 0.00081212953374611378026696723165760;v[6]=-3;
*/
	for(int n=0;n<N;n++){

		chi_feq[n]=w[n];
		chi_S[n]=w[n]*(1.0 - 0.5*(v[n]*v[n]/cs2 - 1.0) );

/*
		chi_feq[n]=w[n]*(1.0 - (1.0/8.0)*( v[n]*v[n]*v[n]*v[n]/(cs2*cs2) - 6.0*v[n]*v[n]/cs2 + 3.0 ) );
		chi_S[n]=w[n]*(1.0 - 0.5*(v[n]*v[n]/cs2 - 1.0) + (1.0/8.0)*( v[n]*v[n]*v[n]*v[n]/(cs2*cs2) - 6.0*v[n]*v[n]/cs2 + 3.0 ) );
*/
	}
	
	for(int i=0;i<Lx;i++)
		for(int n=0;n<N;n++){
			fi[i][n]=w[n]*potential(i)*chi_feq[n];
		}
	
	
}
void CalculoMacro(int t){
	for(int i=0;i<Lx;i++){
		Phi[i]=0.0;
		for(int n=0;n<N;n++){
			Phi[i]+=fi[i][n];
		}
		rho[i]=charge_density(i);
	}
}
void Evolucion(int t){
	CalculoMacro(t);

	for(int i=0;i<Lx;i++)
		for(int n=0;n<N;n++){
			fip[i][n]=fi[i][n]-(1/taoi)*(fi[i][n]-Phi[i]*chi_feq[n]) - dt*( cs2*(taoi - 0.5)*charge_density(i)*chi_S[n] );
		}

	for(int i=0;i<Lx;i++)
		for(int n=0;n<N;n++){
			int ia=i+v[n];
			
			if(ia>Lx-1 || ia<0)
				ia=(i+v[n]+Lx)%Lx;
			
			fi[ia][n]=fip[i][n];
			
		}
	
}

int main(){
	char filename[3][60];
	Inicio();
	FILE *X[3];
	int t=0;
	double error = 0.0;
	do{
		Evolucion(t);
		if(t%100==0)
			cout << t << endl;
		
		error=0.0;
		for(int i=0;i<Lx;i++){
			error += (Phi_old[i]-Phi[i])*(Phi_old[i]-Phi[i])/Lx;
		}
		error = sqrt(error);
		for(int i=0;i<Lx;i++){
			Phi_old[i]=Phi[i];
		}
		t++;
	}while(error > tolerance || t == 1);
	
	sprintf(filename[0], "Phi.dat");
	sprintf(filename[1], "Rho.dat");
	sprintf(filename[2], "Phi_teo.dat");
	
	for(int a=0;a<3;a++)
		X[a]=fopen(filename[a],"w");
	
	for(int i=0;i<Lx;i++){
		fprintf(X[0],"%.7e \n", GetPhi(i));
		fprintf(X[1],"%.7e \n", GetRho(i));
		fprintf(X[2],"%.7e \n", 0.0 );
	}
	
	for(int a=0;a<2;a++)
		fclose(X[a]);
	
	return 0;

}




