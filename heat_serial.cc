#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
using namespace std;

/*
void PrintState(double n, double t, const double *x) {
  printf("%15.8f", t);
  for (int i = 0; i < n; ++i) {
    printf("%15.8f", x[i]);
  }
  printf("\n");
}*/

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("USAGE: %s <grid size> \n", argv[0]);
        exit(1);
    }
    //Input argument
    const int nx = atoi(argv[1]);

    //Constants and definitions
    const double pi = 4.0*atan(1);
    const double domainLength = pi; //square domain, domain length is the same in x and y
    const double dx = domainLength/nx; //dy is the same
    const double kappa = 1.0;
    const double endtime = 0.5*pi*pi/kappa; //Defined in problem PDF
    const double dt = (0.5)*(dx*dx)/(4.0*kappa); //Half the maximum stable time step
  
    //Storage for the temperature
    double T[nx][nx];

    printf("nx: %15.8f \npi: %15.8f \ndomainLength: %15.8f \ndx: %15.8f \nkappa: %15.8f \nendtime: %15.8f \ndt: %15.8f \n", nx, pi, domainLength, dx, kappa, endtime, dt);    
    
    return 0;
}
