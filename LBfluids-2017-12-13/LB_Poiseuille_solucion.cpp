//Mi primer programa en C++
#include <iostream> 
#include <fstream> 
#include <cmath>
using namespace std;

const int Lx=1;
const int Ly=256;

const int Q=9;

const double tau=12.85;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
const double UmU2tau = 1 - 0.5*Utau;
const double Deltat = 1.0;

const double RHO0=1.0, UX0=0.06, UY0=0;
const double g = 0.001;

enum TipoCelda{aire,obstaculo,ventilador};

class LatticeBoltzmann{
private:
  double w[Q];
  int v[2][Q];//v_{ix}=v[0][i] ,v_{iy}=v[1][i]
  double f[Lx][Ly][Q],fnew[Lx][Ly][Q];//f[ix][iy][i]
  TipoCelda Celda[Lx][Ly];
public:
  LatticeBoltzmann(void);
  void ConstruyeTuGeometria(void);
  double rho(int ix,int iy,bool UsarNew);
  double Ux(int ix,int iy,bool UsarNew);
  double Uy(int ix,int iy,bool UsarNew);
  double feq(int i,double rho0,double Ux0,double Uy0);
  double Fi(double rho0, double Ux0, double Uy0, int i);  // Función de forzamiento para la feq por forzamiento Guo
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
  //Obstaculo
  Celda[0][0]=obstaculo;
  Celda[0][Ly-1]=obstaculo;
}
LatticeBoltzmann::LatticeBoltzmann(void){
  //los pesos
  w[0]=4.0/9; 
  w[1]=w[2]=w[3]=w[4]=1.0/9;
  w[5]=w[6]=w[7]=w[8]=1.0/36;
  //los vectores
  v[0][0]=0;
  v[1][0]=0;
  
  v[0][1]=1;  v[0][2]=0;  v[0][3]=-1;  v[0][4]=0;
  v[1][1]=0;  v[1][2]=1;  v[1][3]=0;   v[1][4]=-1;
  
  v[0][5]=1;  v[0][6]=-1; v[0][7]=-1;  v[0][8]=1;
  v[1][5]=1;  v[1][6]=1;  v[1][7]=-1;  v[1][8]=-1;
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
  double FextX;
  double rho0 = rho(ix,iy,UsarNew);
  FextX = g;
  for(i=0;i<Q;i++)
    if(UsarNew)
      suma+=v[0][i]*fnew[ix][iy][i];
    else
      suma+=v[0][i]*f[ix][iy][i];
  suma += 0.5*Deltat*FextX; // Adiciono fuerza externa
  return suma/rho0;
}
double LatticeBoltzmann::Uy(int ix,int iy,bool UsarNew){
  int i; double suma=0;
  double FextY;
  FextY = 0.0;
  double rho0 = rho(ix,iy,UsarNew);
  for(i=0;i<Q;i++)
    if(UsarNew)
      suma+=v[1][i]*fnew[ix][iy][i];
    else
      suma+=v[1][i]*f[ix][iy][i];
  suma += 0.5*Deltat*FextY;  // Adiciono fuerza externa en y
  return suma/rho0;
}
double LatticeBoltzmann::feq(int i,double rho0,double Ux0,double Uy0){
  double UdotVi=Ux0*v[0][i]+Uy0*v[1][i];
  double U2=Ux0*Ux0+Uy0*Uy0;

  return w[i]*rho0*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}

double LatticeBoltzmann::Fi(double rho0, double Ux0, double Uy0, int i){
    double Fext[2]; Fext[0] = rho0*g; Fext[1] = 0;
    double VidotFext = v[0][i]*Fext[0] + v[1][i]*Fext[1];
    double UdotFext = Ux0*Fext[0] + Uy0*Fext[1];
    double VidotU = v[0][i]*Ux0 + v[1][i]*Uy0;
    return UmU2tau*w[i]*(3*VidotFext - 3*UdotFext + 9*VidotU*VidotFext);
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
      for(i=0;i<Q;i++){
        if(Celda0==obstaculo){ //obstaculo
          fnew[ix][iy][0]=f[ix][iy][0];
          fnew[ix][iy][1]=f[ix][iy][3];
          fnew[ix][iy][2]=f[ix][iy][4];
          fnew[ix][iy][3]=f[ix][iy][1];
          fnew[ix][iy][4]=f[ix][iy][2];
          fnew[ix][iy][5]=f[ix][iy][7];
          fnew[ix][iy][6]=f[ix][iy][8];
          fnew[ix][iy][7]=f[ix][iy][5];
          fnew[ix][iy][8]=f[ix][iy][6];
        }
        else //aire
          fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(i,rho0,Ux0,Uy0) +  Deltat*Fi(rho0, Ux0, Uy0, i);
      }
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
    for(int iy=0;iy<Ly;iy+=4)
      MiArchivo<<iy<<" "<<Ux(0,iy,true)<<endl;
  MiArchivo.close();
}

int main(void){
  LatticeBoltzmann Fluido;
  int t,tmax=100000;
  
  Fluido.ConstruyeTuGeometria();
  Fluido.Inicie();
  for(t=0;t<tmax;t++){
    Fluido.Colisione(t);
    Fluido.Adveccione();
  }
  Fluido.Imprimase("Pouseuille2.dat",t);

  double nu = (tau - 0.5)/3;
  cout << "Viscosidad cinematica = " << nu << endl;
  cout << "Reynolds = " << Ly*Fluido.Ux(0, Ly/2, true)/nu << endl;  // d = ly unidades celdas/click

  return 0;
}
