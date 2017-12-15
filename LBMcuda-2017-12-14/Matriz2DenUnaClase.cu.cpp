// Programa que hace c=a+b para a,b,c vectores en CUDA
#include <iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;

#define Lx 16
#define Ly 8
#define Nx 8
#define Ny 8
const int Mx=(Lx+Nx-1)/Nx;
const int My=(Ly+Ny-1)/Nx;

//--------------------KERNELS----------------
__global__ void IncrementarMatriz(float *d_a,size_t pitcha){
  int ix,iy; float *a;
  ix=blockIdx.x*blockDim.x+threadIdx.x;  iy=blockIdx.y*blockDim.y+threadIdx.y;

  a=d_a+(ix*pitcha)/sizeof(float)+iy;
  
  (*a)++;
}
//------------------- CLASES ----------------
class LatticeBoltzmann{
private:
  float h_a[Lx][Ly]; float*d_a; size_t pitcha; 
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  void Inicie(void);
  void Incremente(void);
  void Imprimase(void);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Construir las matrices virtuales en el Device
  cudaMallocPitch((void**) &d_a,&pitcha,Ly*sizeof(float),Lx);
}
LatticeBoltzmann::~LatticeBoltzmann(void){
  cudaFree(d_a);
}
void LatticeBoltzmann::Inicie(void){
  int ix,iy;
  //Cargar valores en el Host
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      h_a[ix][iy]=ix*Ly+iy;
  //Llevar al Device
  cudaMemcpy2D(d_a,pitcha,h_a,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
}
void LatticeBoltzmann::Imprimase(void){
  int ix,iy;
  //Devolver los datos al Host
  cudaMemcpy2D(h_a,Ly*sizeof(float),d_a,pitcha,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  //Mostrar
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;
}
void LatticeBoltzmann::Incremente(void){
  //Procesar en el Device
  dim3 ThreadsPerBlock(Nx,Ny,1);
  dim3 BlocksPerGrid(Mx,My,1);
  IncrementarMatriz<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,pitcha);
}
// ----------------- FUNCIONES GLOBALES ----------
int main(void){
  LatticeBoltzmann Ondas;

  Ondas.Inicie();
  Ondas.Imprimase();
  Ondas.Incremente();
  Ondas.Imprimase();

  return 0;
}
