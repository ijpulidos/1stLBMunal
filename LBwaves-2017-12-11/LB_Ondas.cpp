//Mi primer programa en C++
#include <iostream> 
#include <fstream> 
#include <cmath>
using namespace std;

const int Lx=100;
const int Ly=100;

const int Q=5;
const double W0=1.0/3;

const double C=0.5; // C<0.707 celdas/click
const double TresC2=3*C*C;
const double AUX0=1-TresC2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

class LatticeBoltzmann{
private:
  double w[Q];
  int v[2][Q];//v_{ix}=v[0][i] ,v_{iy}=v[1][i]
  double f[Lx][Ly][Q],fnew[Lx][Ly][Q];//f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UsarNew);
  double Jx(int ix,int iy,bool UsarNew);
  double Jy(int ix,int iy,bool UsarNew);
  double feq(int i,double rho0,double Jx0,double Jy0);
  void Inicie(double rho0,double Jx0,double Jy0);
  void ImponerCampos(int ix,int iy,double & rho0,double & Jx0,double & Jy0,int t);
  void Colisione(int t);
  void Adveccione(void);
  void Imprimase(const char * NombreArchivo,int t);
};

LatticeBoltzmann::LatticeBoltzmann(void){
  //los pesos
  w[0]=1.0/3; w[1]=w[2]=w[3]=w[4]=1.0/6;
  //los vectores
  v[0][0]=0;
  v[1][0]=0;
  
  v[0][1]=1;  v[0][2]=0;  v[0][3]=-1;  v[0][4]=0;
  v[1][1]=0;  v[1][2]=1;  v[1][3]=0;   v[1][4]=-1;
}

double LatticeBoltzmann::rho(int ix,int iy,bool UsarNew){
  int i; double suma=0;
  for(i=0;i<Q;i++)
    if(UsarNew)
      suma+=fnew[ix][iy][i];
    else
      suma+=f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jx(int ix,int iy,bool UsarNew){
  int i; double suma=0;
  for(i=0;i<Q;i++)
    if(UsarNew)
      suma+=v[0][i]*fnew[ix][iy][i];
    else
      suma+=v[0][i]*f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jy(int ix,int iy,bool UsarNew){
  int i; double suma=0;
  for(i=0;i<Q;i++)
    if(UsarNew)
      suma+=v[1][i]*fnew[ix][iy][i];
    else
      suma+=v[1][i]*f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::feq(int i,double rho0,double Jx0,double Jy0){
  if(i==0)
    return AUX0*rho0;
  else
    return w[i]*(TresC2*rho0 + 3* (v[0][i]*Jx0+v[1][i]*Jy0) );
}
void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++)
	f[ix][iy][i]=fnew[ix][iy][i]=feq(i,rho0,Jx0,Jy0);
}
void LatticeBoltzmann::ImponerCampos(int ix,int iy,double & rho0,double & Jx0,double & Jy0,int t){
  double A=10, lambda=10, omega=2*M_PI*C/lambda;
  if(ix==Lx/2 && iy==Ly/2)
    rho0=A*sin(omega*t);
}
void LatticeBoltzmann::Colisione(int t){
  int ix,iy,i; double rho0,Jx0,Jy0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){//Para cada celda
      //Calcular las cantidades macroscopicas
      rho0=rho(ix,iy,false);  Jx0=Jx(ix,iy,false);  Jy0=Jy(ix,iy,false);
      ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
      //Calcular las fnew en cada direccion
      for(i=0;i<Q;i++) 
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(i,rho0,Jx0,Jy0);
    }
}
void LatticeBoltzmann::Adveccione(void){
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++) 
	f[(ix+v[0][i]+Lx)%Lx][(iy+v[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}
void LatticeBoltzmann::Imprimase(const char * NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0;
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);   Jx0=Jx(ix,iy,true);  Jy0=Jy(ix,iy,true);
      ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
      MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=100;
  
  Ondas.Inicie(0,0,0);
  for(t=0;t<tmax;t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
  }
  Ondas.Imprimase("Ondas.dat",t);

  return 0;
}
