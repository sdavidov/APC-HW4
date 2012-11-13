#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
//#include <fstream>
using namespace std;

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
    const double dx = domainLength/(nx-1.0); //dy is the same
    const double kappa = 1.0;
    const double endtime = 0.5*pi*pi/kappa; //Defined in problem PDF
    const double dt = (0.5)*(dx*dx)/(4.0*kappa); //Half the maximum stable time step
    const int numsteps = 4*(nx-1)*(nx-1); // = endtime/dt but done analytically as one way to ensure integer value
    const double mult = dt*kappa/pow(dx,2); //mutliplicative factor for time advance

    printf("nx: %i \npi: %15.8f \ndomainLength: %15.8f \ndx: %15.8f \nkappa: %15.8f \nendtime: %15.8f \ndt: %15.8f \n numsteps: %i \n", nx, pi, domainLength, dx, kappa, endtime, dt,numsteps);
  
    //Storage for the temperature
    double T[nx][nx]; //T[x][y]
    double Tprev[nx][nx];
    /*Note: In order for our domain to truly have a length of pi in the x (periodic) direction, we must equate T[0][j] = T[nx-1][j] and then skip one of these values appropriately in evaluating the finite difference (i.e., we should NOT just do a*T[nx-2][j]+b*T[nx-1][j]+c*T[0][j]), or we effectively lengthen the domain by dx. If we want to do this, we should have T[nx-1][nx], but this makes things a bit uglier when we want to domain decompose with even numbers of nx as assigned*/
    
    //Set the boundary temperature for the y direction
    for (int i=0;i < nx;++i){
        Tprev[i][0] = pow(cos(i*dx),2);
        T[i][0] = Tprev[i][0];
        Tprev[i][nx-1] = pow(sin(i*dx),2);
        T[i][nx-1] = Tprev[i][nx-1];
    }
    //Set the initial condition for the rest of the grid
    for (int i=0;i < nx;++i){
        for (int j=1;j < (nx-1);++j){
            Tprev[i][j] = 0.0;
            T[i][j] = 0.0;
        }
    }
    //Advance solution
    for (int t=0;t < numsteps;++t){
        //Domain interior
        for(int i=1;i < (nx-1);++i){
            for(int j=1;j < (nx-1);++j){
                T[i][j] = Tprev[i][j]+mult*(Tprev[i+1][j]+Tprev[i-1][j]+Tprev[i][j+1]+Tprev[i][j-1]-4.0*T[i][j]);
            }
        }
        //Boundaries. Only need to take care of periodic boundary, since y boundary values are const
        for(int j=1;j < (nx-1);++j){
            T[0][j] = Tprev[0][j]+mult*(Tprev[1][j]+Tprev[nx-2][j]+Tprev[0][j+1]+Tprev[0][j-1]-4.0*T[0][j]);
            T[nx-1][j] = T[0][j];
        }
        
        //Update Tprev to prepare for the next step
        for(int i=0;i < nx;++i){
            for(int j=0;j < nx;++j){
                Tprev[i][j] = T[i][j];
            }
        }
    }

    //Calculate volume averaged temperature
    double volAve = 0;    
    for(int i=0;i < nx;++i){
        for(int j=0;j < nx;++j){
            volAve = volAve + T[i][j];
        }
    }
    volAve = volAve/(pow(nx-1,2)); //All cells the same volume, so we can just do this 

    cout << "Volume average: " << volAve << "\n\n";

    //Dump matrix
    for(int j=0;j < nx;++j){
        for(int i=0;i < nx;++i){
            cout << Tprev[i][nx-j-1] << ",";
        }
        cout << "\n\n";
    }
      
    return 0;
}
