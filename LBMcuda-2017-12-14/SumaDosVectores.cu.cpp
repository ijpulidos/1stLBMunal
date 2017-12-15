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
#define Nx 8
const int Mx=(Lx+Nx-1)/Nx;

//--------------------KERNELS----------------
__global__ void SumaVectores(float *d_a,float *d_b,float *d_c){
  int ix;
  ix=blockIdx.x*blockDim.x+threadIdx.x;
  
  d_c[ix]=d_a[ix]+d_b[ix];
}

int main(void){
  int ix;
  //DECLARAR
  //Declarar matrices en el Host
  float h_a[Lx],h_b[Lx],h_c[Lx]; 
  //Declarar matrices en el Device
  float*d_a; cudaMalloc((void**) &d_a,Lx*sizeof(float));
  float*d_b; cudaMalloc((void**) &d_b,Lx*sizeof(float));
  float*d_c; cudaMalloc((void**) &d_c,Lx*sizeof(float));

  //CARGAR
  //Cargar valores en el Host
  for(ix=0;ix<Lx;ix++){
    h_a[ix]=ix; h_b[ix]=2*ix;
  } 
  //Llevar al Device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(d_b,h_b,Lx*sizeof(float),cudaMemcpyHostToDevice);
  
  //PROCESAR
  //Procesar en el Device
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  SumaVectores<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,d_b,d_c);

  //MOSTRAR RESULTADOS
  //Devolver al Host
  cudaMemcpy(h_c,d_c,Lx*sizeof(float),cudaMemcpyDeviceToHost);
  for(ix=0;ix<Lx;ix++) cout<<ix<<" "<<h_c[ix]<<endl;

  cudaFree(d_a);  cudaFree(d_b);  cudaFree(d_c);

  return 0;
}
