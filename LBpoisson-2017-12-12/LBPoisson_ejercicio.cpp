#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#define N 3  // N: Número de pesos y vectores
//#define N 7
#define Lx 128  // Tamaño - Número de nodos
#define tolerance 1.0E-12  // Parámetro de tolerancia para derivada temporal = 0

using namespace std;

double dt = 1.0;
double cs2 = 1.0/3.0;  // Velocidad térmica/sonido al cuadrado
//double cs2 = 0.69795332201968308823840905538933;

double k = 2*M_PI/Lx;  // Número de onda característico del espacio
double k2 = k*k;
const double A = 0.01;

double potential(int i){  // i por convención es índice de variable espacial
    /*
    Función de condición inicial de potencial.
    */
	double g;
	g=-A*sin(k*i)/k2;
	return (g);
}
double charge_density(int i){
    /*
    Función para definir la densidad de carga
    */
	return ( A*sin(k*i) );
}

double fi[Lx][N], fip[Lx][N];  // Distribuciones (2 copias)
double taoi, Phi[Lx], Phi_old[Lx], rho[Lx];  // Tiempo de relajación, Potencial, densidad de carga
int v[N];  // Vectores de velocidad micro
double chi_feq[N], chi_S[N];
double w[N];  // Pesos
void CalculoMacro(int t);  // Función de cálculo de variables macroscópicas
double GetPhi(int i){return (Phi[i]);};  // Función valor de Phi/Potencial
double GetRho(int i){return (rho[i]);};  // Función que retorna valor de densidad
void Inicio(void);  // FUnción de inicialización de todas las variables
void Evolucion(int t);  // Función evolución con colisión + streaming

void Inicio(void){
    /*
    Función de inicialización de variables
    */
	
	taoi=1.0;

	w[0]=2.0/3.0; v[0]=0;
	w[1]=1.0/6.0; v[1]=1;
	w[2]=1.0/6.0; v[2]=-1;

/*
    // Pesos para expansión de orden 4
	w[0] = 0.47666988658920744233176499640793;v[0]=0;
	w[1] = 0.23391473782682471337286677188259;v[1]=1;
	w[2] = 0.23391473782682471337286677188259;v[2]=-1;
	w[3] = 0.026938189344825451680983762681790;v[3]=2;
	w[4] = 0.026938189344825451680983762681790;v[4]=-2;
	w[5] = 0.00081212953374611378026696723165760;v[5]=3;
	w[6] = 0.00081212953374611378026696723165760;v[6]=-3;
*/
	for(int n=0;n<N;n++){

		chi_feq[n]=w[n];  // Orden 2
		chi_S[n]=w[n]*(1.0 - 0.5*(v[n]*v[n]/cs2 - 1.0) );  // Orden 2

/*
		chi_feq[n]=w[n]*(1.0 - (1.0/8.0)*( v[n]*v[n]*v[n]*v[n]/(cs2*cs2) - 6.0*v[n]*v[n]/cs2 + 3.0 ) );
		chi_S[n]=w[n]*(1.0 - 0.5*(v[n]*v[n]/cs2 - 1.0) + (1.0/8.0)*( v[n]*v[n]*v[n]*v[n]/(cs2*cs2) - 6.0*v[n]*v[n]/cs2 + 3.0 ) );
*/
	}
	
	for(int i=0;i<Lx;i++)
		for(int n=0;n<N;n++){
			fi[i][n] = 0*potential(i)*chi_feq[n];  // Se inicializa f a la función de equilibrio
		}
	
	
}
void CalculoMacro(int t){
	for(int i=0;i<Lx;i++){
		Phi[i]=0.0;
		for(int n=0;n<N;n++){
			Phi[i]+=fi[i][n];  // Valor del pontecial, suma de f
		}
		rho[i]=charge_density(i);  // Cambiando en tiempo valor de carga (de ser necesario)
	}
}
void Evolucion(int t){
    /*
    Función de Evolución.
    */
	CalculoMacro(t);  // Calcular variables macroscópicas

    // Colisión
	for(int i=0;i<Lx;i++)
		for(int n=0;n<N;n++){
			fip[i][n]=fi[i][n]-(1/taoi)*(fi[i][n]-Phi[i]*chi_feq[n]) -
			dt*( cs2*(taoi - 0.5)*charge_density(i)*chi_S[n] );
		}

    // Streaming/Advection
	for(int i=0;i<Lx;i++)
		for(int n=0;n<N;n++){
			int ia=i+v[n];
			
			if(ia>Lx-1 || ia<0)
				ia=(i+v[n]+Lx)%Lx;  // Condiciones de frontera periódicas
			
			fi[ia][n]=fip[i][n];  // Streaming/Advection
			
		}
	
}

double rmsError(void){
    /*
    Function that computes the root-mean-squre error between the analytical and numerical solution
    */
    int i;
    double sum = 0;
    for (i=0; i<Lx; i++){
       sum += (Phi[i] - potential(i))*(Phi[i] - potential(i));
    }
    return sqrt(sum)/Lx;
}

int main(){
	char filename[3][60];  // Tres archivos de salida: densidad carga, potencial, potencial analítico
	Inicio();  // Inicializo lattice
	FILE *X[3];  // Apuntadores para archivos de salida
	int t=0;
	double error = 0.0;
	do{
		Evolucion(t);  // Colisión + streaming/advection
		if(t%100==0)
			cout << t << endl;  // Imprimiendo tiempo para ver progreso
		
		error=0.0;
		for(int i=0;i<Lx;i++){
			error += (Phi_old[i]-Phi[i])*(Phi_old[i]-Phi[i]);  // calcular error/varianza numérico, para chequear tolerancia
		}
		error = sqrt(error)/Lx;
		for(int i=0;i<Lx;i++){
			Phi_old[i]=Phi[i];
		}
		t++;
	}while(error > tolerance || t == 1);

	// Nombres de archivos de salida
	sprintf(filename[0], "Phi.dat");
	sprintf(filename[1], "Rho.dat");
	sprintf(filename[2], "Phi_teo.dat");

	// Abrir archivos para escritura
	for(int a=0;a<3;a++)
		X[a]=fopen(filename[a],"w");

	// Imprimir valores a archivos
	for(int i=0;i<Lx;i++){
		fprintf(X[0],"%.7e \n", GetPhi(i));
		fprintf(X[1],"%.7e \n", GetRho(i));
		fprintf(X[2],"%.7e \n", potential(i) );
	}

	// Cierre de archivos
	for(int a=0;a<2;a++)
		fclose(X[a]);

    // Printing rms error
    cout << "RMS error is: " << rmsError() << endl;

	return 0;
}
