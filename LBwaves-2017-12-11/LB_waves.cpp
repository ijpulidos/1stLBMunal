/*
First Workshop on Lattice-Boltzmann Methods at Universidad Nacional de Colombia
First day of practice session 2017-12-11
Notes from Iván Pulido
*/
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double W0 = 1/3.0;  // Constante para el peso cero
const double C = 0.5;  // Velocidad "física", debe ser menor que la velocidad de la información sqrt(2) en celdas/click
const double ThreeC2 = 3*C*C;
const double AUX0 = 1-ThreeC2*(1-W0);

const double tau = 0.5;
const double Utau = 1.0/tau;
const double UmUtau = 1.0-Utau;

const int Lx = 512;  // Tamaño dominio en x
const int Ly = 512;  // Tamaño dominio en y
const int Q = 5;  // Número de pesos/vectores en cada celda


class LatticeBoltzmann{
    /*
    Generic class for simple Lattice-Boltzmann method for waves with D2Q5 scheme
    */
    // TODO: Tratar de hacer las cosas más eficientes al usar punteros en lugar de arreglos bidimensionales.
    private:
        double weights[Q];  // Pesos del esquema de LB
//        int velocities[2][Q];  // Velocidades de la celda (micro). Convención: v_{ix}=velocities[0][i]
                                  // Dejar la i de último índice hace que sea más eficiente por memoria contigua
                                  // Tienen que ser int/enteros para que se puedan sumar dentro de la advección
        int** velocities;
//        double f[Lx][Ly][Q], f_new[Lx][Ly][Q];  // Probability densities arrays f[ix][iy][i]
        double*** f;
        double*** f_new;
    public:
        LatticeBoltzmann(void);
        double rho(int ix, int iy, bool isNew);
        double Jx(int ix, int iy, bool isNew);
        double Jy(int ix, int iy, bool isNew);
        double feq(int i, double rho0, double Jx0, double Jy0);
        void initialize(double rho0, double Jx0, double Jy0);
        void collision(int t);
        void streaming(void);
        void imposeFields(int ix, int iy, double & rho0, double & Jx0, double & Jy0, int t);
        void print_data(const char * filename, int t);
};

LatticeBoltzmann::LatticeBoltzmann(void){
    /*
    Constructor de la clase LatticeBoltzmann. Asigna valores a los pesos y vectores v.
    */
    // Weights
    weights[0] = W0;
    weights[1] = weights[2] = weights[3] = weights[4] = 1.0/6;
    // Allocation of dynamic arrays
    velocities = new int*[2];
    for (int i=0; i<2; ++i)
            velocities[i] = new int[Q];
    f = new double**[Lx];
    for (int i=0; i<Lx; ++i)
        f[i] = new double*[Ly];
    for (int i=0; i<Lx; ++i){
        for (int j=0; j<Ly; ++j)
            f[i][j] = new double[Q];
    }
    f_new = new double**[Lx];
    for (int i=0; i<Lx; ++i)
        f_new[i] = new double*[Ly];
    for (int i=0; i<Lx; ++i){
        for (int j=0; j<Ly; ++j)
            f_new[i][j] = new double[Q];
    }
    // Velocities D2Q5
    velocities[0][0] = 0;  // v_0
    velocities[1][0] = 0;
    velocities[0][1] = 1;  // v_1
    velocities[1][1] = 0;
    velocities[0][2] = 0;  // v_2
    velocities[1][2] = 1;
    velocities[0][3] = -1;  // v_3
    velocities[1][3] = 0;
    velocities[0][4] = 0;  // v_4
    velocities[1][4] = -1;
}

double LatticeBoltzmann::rho(int ix, int iy, bool isNew){
    /*
    Method that computes the macroscopic density rho. Adding densities f together in each cell.
    */
    int i;
    double _sum = 0;
    if(isNew)
        for (i=0; i<Q; i++)
            _sum += f_new[ix][iy][i];
    else
        for (i=0; i<Q; i++)
            _sum += f[ix][iy][i];
    return _sum;
}

double LatticeBoltzmann::Jx(int ix, int iy, bool isNew){
    /*
    Method that computes the x component of the macroscopic 2nd moment (linear momentum). Adding velocities times f
    densities.
    */
    int i;
    double _sum = 0;
    if(isNew)
        for (i=0; i<Q; i++)
            _sum += velocities[0][i]*f_new[ix][iy][i];
    else
        for (i=0; i<Q; i++)
            _sum += velocities[0][i]*f[ix][iy][i];
    return _sum;
}

double LatticeBoltzmann::Jy(int ix, int iy, bool isNew){
    /*
    Method that computes the x component of the macroscopic 2nd moment (linear momentum). Adding velocities times f
    densities.
    */
    int i;
    double _sum = 0;
    if(isNew)
        for (i=0; i<Q; i++)
            _sum += velocities[1][i]*f_new[ix][iy][i];
    else
        for (i=0; i<Q; i++)
            _sum += velocities[1][i]*f[ix][iy][i];
    return _sum;
}

double LatticeBoltzmann::feq(int i, double rho0, double Jx0, double Jy0){
    /*
    Method to compute the equilibrium function for the simple Lattice Boltzmann method for waves.
    */
    if(i==0)
        return AUX0*rho0;
    else
        return weights[i]*(ThreeC2*rho0 + 3*(velocities[0][i]*Jx0 + velocities[1][i]*Jy0));
}

void LatticeBoltzmann::initialize(double rho0, double Jx0, double Jy0){
    /*
    Method that initializes the Lattice-Boltzmann model, filling up the density functions f_i and f_new_i
    */
    for (int ix=0; ix<Lx; ix++)
        for (int iy=0; iy<Ly; iy++)
            for (int i=0; i<Q; i++)
                f[ix][iy][i] = f_new[ix][iy][i] = feq(i, rho0, Jx0, Jy0);
}

void LatticeBoltzmann::collision(int t){
    /*
    Method that implements the collision step for the Lattice-Boltzmann method for waves.
    */
    int ix, iy, i;
    double rho0, Jx0, Jy0;
    for (ix=0; ix<Lx; ix++)
        for (iy=0; iy<Ly; iy++){  // For each cell
            rho0 = rho(ix, iy, false);
            Jx0 = Jx(ix, iy, false);
            Jy0 = Jy(ix, iy, false);
            imposeFields(ix, iy, rho0, Jx0, Jy0, t);
            for (i=0; i<Q; i++)  // Compute f_new in each direction
                f_new[ix][iy][i] = UmUtau*f[ix][iy][i] + Utau*feq(i, rho0, Jx0, Jy0);
        }
}

void LatticeBoltzmann::streaming(void){
    /*
    Method that implements the streaming/advection step for the Lattice-Boltzmann method for waves.
    */
    // TODO: Mirar manera más eficiente de hacer las condiciones de frontera periódicas, tal vez con ciclos separados.
    for (int ix=0; ix<Lx; ix++)
        for (int iy=0; iy<Ly; iy++)
            for (int i=0; i<Q; i++)
                f[(ix+velocities[0][i]+Lx)%Lx][(iy+velocities[1][i]+Ly)%Ly][i] = f_new[ix][iy][i];
}

void LatticeBoltzmann::imposeFields(int ix, int iy, double & rho0, double & Jx0, double & Jy0, int t){
    /*
    This method forces the fields to be something specific. This basically perturbates the domain with a sine function.
    */
    double A = 10, lambda = 10, omega = 2*M_PI*C/lambda;
    if (ix == Lx/2 && iy == Ly/2)
        rho0 = A * sin(omega*t);
}

void LatticeBoltzmann::print_data(const char * filename, int t){
    /*
    Method that prints data as a 2D array that GNUplot can handle, meshgrid type of print.
    */
    ofstream myFile(filename); double rho0, Jx0, Jy0;
    for (int ix=0; ix<Lx; ix++){
        for (int iy=0; iy<Ly; iy++){
            rho0 = rho(ix, iy, true);
            Jx0 = Jx(ix, iy, true);
            Jy0 = Jy(ix, iy, true);
            imposeFields(ix, iy, rho0, Jx0, Jy0, t);
            myFile << ix << " " << iy << " " << rho0 << endl;
        }
        myFile << endl;
    }
    myFile.close();
}

double fsinc(double x){
    /*
    Función Sinc de ejemplo para saber cómo definir funciones
    */
    return sin(x)/x;
}

int main() {
    /*
    Función principal para el método de Lattice-Boltzmann para ondas
    */
    LatticeBoltzmann ondas;  // Instancio objeto ondas de clase LatticeBoltzmann
    int time;  // variable for storing instant of time
    int tmax=100;

    ondas.initialize(0, 0, 0);

    for (time=0; time<tmax; time++){
        ondas.collision(time);
        ondas.streaming();
    }

    ondas.print_data("waves.dat", time);

    return 0;
}