// Suma de dos vectores c=a+b en CUDA
#include<iostream>
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
const int My=(Ly+Ny-1)/Ny;

//--------------- KERNELS ----------------
__constant__ float d_w[5];
__constant__ int d_Vx[5];
__constant__ int d_Vy[5];

__device__ float SumeleUno(float x){
  return x+1;
}
__global__ void SumarConstante(float * d_a,size_t pitcha){
  int ix,iy; float *aux;
  ix=blockIdx.x*blockDim.x+threadIdx.x;  
  iy=blockIdx.y*blockDim.y+threadIdx.y;
  
  aux=d_a+(ix*pitcha)/sizeof(float)+iy; // aux es &(d_a[ix][iy])
  
  (*aux)++; //  (*aux) es d_a[ix][iy]
}

int main(){
  float h_w[5]; // w[i]
  int h_Vx[5],h_Vy[5]; // Vx[i],Vy[i]
  int ix,iy;

  //DECLARAR LAS MATRICES
  //CONSTANTES
  //---Cargar las constantes en el Host-----------------
  //Cargar los pesos
  h_w[0]=1.0/3;    h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/6;
  //Cargar los vectores
  h_Vx[0]=0;  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;  
  h_Vy[0]=0;  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;
  //------Enviarlas al Device-----------------
  cudaMemcpyToSymbol(d_w,h_w,5*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,5*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,5*sizeof(int),0,cudaMemcpyHostToDevice);

  //DECLARAR
  //Declarar matrices en el Host
  float h_a[Lx][Ly];
  //Declarar matrices en el Device
  float*d_a; size_t pitcha; 
  cudaMallocPitch((void**) &d_a,&pitcha,Ly*sizeof(float),Lx);

  //INICIALIZAR LOS DATOS
  //Cargar los datos en el Host
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      h_a[ix][iy]=Ly*ix+iy;
  //Mostrar
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;

  
  //PROCESAR EN LA TARJETA GRAFICA
  //Enviarlos al Device
  cudaMemcpy2D(d_a,pitcha,h_a,Ly*sizeof(float),Ly*sizeof(float),Lx,cudaMemcpyHostToDevice);
  dim3 ThreadsPerBlock(Nx,Ny,1);
  dim3 BlocksPerGrid(Mx,My,1);
  SumarConstante<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,pitcha);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy2D(h_a,Ly*sizeof(float),d_a,pitcha,Ly*sizeof(float),Lx,cudaMemcpyDeviceToHost);
  //Imprimirlos
  //Mostrar
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++)
      cout<<h_a[ix][iy]<<" ";
    cout<<endl;
  }
  cout<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);

  return 0;
}
