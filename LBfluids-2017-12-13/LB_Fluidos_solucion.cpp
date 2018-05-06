//Mi primer programa en C++
#include <iostream> 
#include <fstream> 
#include <cmath>
using namespace std;

const int Lx=512;
const int Ly=128;

const int Q=9;

const double tau=0.53;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

const double RHO0=1.0, UX0=0.1, UY0=0;

enum TipoCelda{aire,obstaculo,ventilador};

class LatticeBoltzmann{
private:
  double w[Q];
  //int v[2][Q];//v_{ix}=v[0][i] ,v_{iy}=v[1][i]
  int** v;
  //double f[Lx][Ly][Q],fnew[Lx][Ly][Q];//f[ix][iy][i]
  double*** f;
  double*** fnew;
  //TipoCelda Celda[Lx][Ly];
  TipoCelda** Celda;
public:
  LatticeBoltzmann(void);
  void ConstruyeTuGeometria(void);
  double rho(int ix,int iy,bool UsarNew);
  double Ux(int ix,int iy,bool UsarNew);
  double Uy(int ix,int iy,bool UsarNew);
  double feq(int i,double rho0,double Ux0,double Uy0);
  void Inicie(void);
  void Colisione(int t);
  void Adveccione(void);
  void Imprimase(const char * NombreArchivo,int t);
};
void LatticeBoltzmann::ConstruyeTuGeometria(void){
  int ix,iy;
  //Aire
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      Celda[ix][iy]=aire;
  //Ventilador
  for(iy=0;iy<Ly;iy++)
    Celda[0][iy]=ventilador;
  //Obstaculo
  int xc=Lx/8.0, yc=Ly/2.0, R=Ly/5.0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      if((ix-xc)*(ix-xc)+(iy-yc)*(iy-yc)<=R*R)
	Celda[ix][iy]=obstaculo;
  Celda[xc][yc+R+1]=obstaculo;
}
LatticeBoltzmann::LatticeBoltzmann(void){
  //los pesos
  w[0]=4.0/9; 
  w[1]=w[2]=w[3]=w[4]=1.0/9;
  w[5]=w[6]=w[7]=w[8]=1.0/36;
  // Allocating dynamic memory arrays
  // Cell Types
  Celda = new TipoCelda*[Lx];
  for (int i=0; i<Lx; ++i)
    Celda[i] = new TipoCelda[Ly];
  // Velocities
  v = new int*[2];
  for (int i=0; i<2; ++i)
    v[i] = new int[Q];
  v[0][0]=0;
  v[1][0]=0;
  
  v[0][1]=1;  v[0][2]=0;  v[0][3]=-1;  v[0][4]=0;
  v[1][1]=0;  v[1][2]=1;  v[1][3]=0;   v[1][4]=-1;
  
  v[0][5]=1;  v[0][6]=-1; v[0][7]=-1;  v[0][8]=1;
  v[1][5]=1;  v[1][6]=1;  v[1][7]=-1;  v[1][8]=-1;
  
  // Equilibrium distributions
  f =  new double**[Lx];
  for (int i=0; i<Lx; ++i)
    f[i] = new double*[Ly];
  for (int i=0; i<Lx; ++i)
    for (int j=0; j<Ly; ++j)
      f[i][j] = new double[Q];
  
  fnew =  new double**[Lx];
  for (int i=0; i<Lx; ++i)
    fnew[i] = new double*[Ly];
  for (int i=0; i<Lx; ++i)
    for (int j=0; j<Ly; ++j)
      fnew[i][j] = new double[Q];
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
double LatticeBoltzmann::Ux(int ix,int iy,bool UsarNew){
  int i; double suma=0;
  for(i=0;i<Q;i++)
    if(UsarNew)
      suma+=v[0][i]*fnew[ix][iy][i];
    else
      suma+=v[0][i]*f[ix][iy][i];
  return suma/rho(ix,iy,UsarNew);
}
double LatticeBoltzmann::Uy(int ix,int iy,bool UsarNew){
  int i; double suma=0;
  for(i=0;i<Q;i++)
    if(UsarNew)
      suma+=v[1][i]*fnew[ix][iy][i];
    else
      suma+=v[1][i]*f[ix][iy][i];
  return suma/rho(ix,iy,UsarNew);
}
double LatticeBoltzmann::feq(int i,double rho0,double Ux0,double Uy0){
  double UdotVi=Ux0*v[0][i]+Uy0*v[1][i];
  double U2=Ux0*Ux0+Uy0*Uy0;

  return w[i]*rho0*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}
void LatticeBoltzmann::Inicie(void){
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++)
	if(Celda[ix][iy]==obstaculo)
	  f[ix][iy][i]=fnew[ix][iy][i]=feq(i,RHO0,0,0);
	else
	  f[ix][iy][i]=fnew[ix][iy][i]=feq(i,RHO0,UX0,UY0);
}
void LatticeBoltzmann::Colisione(int t){
  int ix,iy,i; double rho0,Ux0,Uy0; TipoCelda Celda0;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){//Para cada celda
      //Calcular las cantidades macroscopicas
      rho0=rho(ix,iy,false); Ux0=Ux(ix,iy,false); Uy0=Uy(ix,iy,false);
      Celda0=Celda[ix][iy];
      //Calcular las fnew en cada direccion
      for(i=0;i<Q;i++)
	if(Celda0==ventilador) //ventilador
	  fnew[ix][iy][i]=feq(i,RHO0,UX0,UY0);
	else if(Celda0==obstaculo) //obstaculo
	  fnew[ix][iy][i]=feq(i,RHO0,0,0);
	else //aire
	  fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(i,rho0,Ux0,Uy0);
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
  ofstream MiArchivo(NombreArchivo);
  for(int ix=0;ix<Lx;ix+=4){
    for(int iy=0;iy<Ly;iy+=4)
      MiArchivo<<ix<<" "<<iy<<" "
	       <<4.0/UX0*(Ux(ix,iy,true)-UX0)<<" "
	       <<4.0/UX0*Uy(ix,iy,true)<<endl;
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

int main(void){
  LatticeBoltzmann Fluido;
  int t,tmax=8000;
  
  Fluido.ConstruyeTuGeometria();
  Fluido.Inicie();
  for(t=0;t<tmax;t++){
    Fluido.Colisione(t);
    Fluido.Adveccione();
  }
  Fluido.Imprimase("fluid.dat",t);

  return 0;
}
