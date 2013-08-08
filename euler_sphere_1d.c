#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#define PI    3.14159265
#define G     .1
#define GAMMA 1.4
#define X 10000
#define XMIN 0.0
#define XMAX 1.0
#define THETA 2.0
#define tmax 1.0


void Grid(double *gridX) {
  *gridX = (XMAX - XMIN) / (X - 1);
}

double sgn(double x) { 				// Gives the sign of x
  double S;
  if (x>0.) S = 1.;
  else S = -1.;
  return S;
}

double minmod(double x, double y, double z) {		// The minmod function, described in the paper
  double min, M;
  min = fabs(x);
  if (fabs(y) < min) min = fabs(y);
  if (fabs(z) < min) min = fabs(z);
  M = .25 * fabs(sgn(x) + sgn(y)) * (sgn(x) + sgn(z)) * min;
  return M;
}

// Given the density, speed and pressure, we calculate the flux F:
void FluxCalcPX(double *FX, double *phys, int N) {
  int i=0;
  for (i=0;i<N;i++) {
      FX[3*i+0] = phys[3*i+0] * phys[3*i+1];
      FX[3*i+1] = phys[3*i+0] * phys[3*i+1] * phys[3*i+1] + phys[3*i+2];
      FX[3*i+2] = phys[3*i+1] *(phys[3*i+0] *(phys[3*i+1] * phys[3*i+1]) * .5 + phys[3*i+2] * GAMMA / (GAMMA - 1));

  }
}

// Given the density, speed and pressure, we calculate the conserved variable U:
void Ucalc(double *U, double *physicalvar, int N) {     // physicalvar[] = rho, vx, vy, vz, p
  int i=0;
  for (i=0;i<N;i++) {
    U[3*i+0] = physicalvar[3*i+0];
    U[3*i+1] = physicalvar[3*i+0] * physicalvar[3*i+1];
    U[3*i+2] = physicalvar[3*i+2] / (GAMMA - 1) + 0.5 * physicalvar[3*i+0] * (physicalvar[3*i+1] * physicalvar[3*i+1]);
  }
}

// And the inverse:
void Ucalcinv(double *physicalvar, double *U, int N) {     // physicalvar[] = rho, v, p
  int i=0;
  for (i=0;i<N;i++) {
    physicalvar[3*i+0] = U[3*i+0];
    physicalvar[3*i+1] = U[3*i+1] / U[3*i+0];
    physicalvar[3*i+2] = (U[3*i+2] - 0.5 * (U[3*i+1] * U[3*i+1])/ U[3*i+0]) * (GAMMA - 1);
  }
}



///SNR initial condition

void sphere(double *physical) {
  int i=0, N;
  double   dx;
  Grid (&dx);
  double x;

  for (i=0;i<X;i++) {
    
    x = XMIN + dx * i;
    N = 3 * i;
    if (x*x<0.01) {
      physical[N+0] = 1.0;
      physical[N+1] = 0.0;
      physical[N+2] = 10.0;
    }
    else {
      physical[N+0] = 0.1/(x*x);
      physical[N+1] = 0.0;
      physical[N+2] = 0.1;
    }
  }
}

// Calculates the flux at a half integer coordinate using the riemann method.
void riemansolverX(double *F_mid, double *U, double *max) {
  double *phys  = (double*) malloc (3*(X)*sizeof (double));
  double *phys_temp= (double*) malloc (3*(X+4)*sizeof (double));
  double *physL = (double*) malloc (3*(X+1)*sizeof (double));
  double *physR = (double*) malloc (3*(X+1)*sizeof (double));
  double *FLX   = (double*) malloc (3*(X+1)*sizeof (double));
  double *FRX   = (double*) malloc (3*(X+1)*sizeof (double));
  double *UL    = (double*) malloc (3*(X+1)*sizeof (double));
  double *UR    = (double*) malloc (3*(X+1)*sizeof (double));
  double SoundSpeedL, SoundSpeedR;											//L, SoundSpeedR, SoundSpeedLL, SoundSpeedRR;
  double AlphaPlus=0, AlphaMinus=0;
  int i, l, N;
  int Sx = 3;

  Ucalcinv (phys,U,(X));

  *max=0;

/**********periodic***************/
  for (i=0;i<X+4;i++) {
    for (l=0;l<3;l++) {
      N = Sx*i+l;
      if (i<2) {
      phys_temp[N] = phys[Sx*(i+X-2)+l];
      }
      else if (i>X+1) {
        phys_temp[N] = phys[Sx*(i-X-2)+l];
      }
      else phys_temp[N] = phys[N-2*Sx];
    }
  }
/************************************/

/***************outflow*****************
  for (i=0;i<X+4;i++) {
    for (j=0;j<Y;j++) {
      for (k=0;k<Z;k++) {
        for (l=0;l<5;l++) {
          N = Sx*i+Sy*j+Sz*k+l;
          if      (i<2) {
            phys_temp[N] = phys[Sx*(0)+Sy*j+Sz*k+l];
          }
          else if (i>X+1) {
            phys_temp[N] = phys[Sx*(X-1)+Sy*j+Sz*k+l];
          }
          else phys_temp[N] = phys[N-2*Sx];
        }
      }
    }
  }
/***************************************/

  for (i=0;i<X+1;i++) {
        for (l=0; l<3; l++) {
          N = Sx*i+l;
          physL[N] = phys_temp[N+Sx]   + 0.5 * minmod (THETA*(phys_temp[N+Sx]   - phys_temp[N]),    0.5*(phys_temp[N+2*Sx] - phys_temp[N]),    THETA*(phys_temp[N+2*Sx] - phys_temp[N+Sx]));
          physR[N] = phys_temp[N+2*Sx] - 0.5 * minmod (THETA*(phys_temp[N+2*Sx] - phys_temp[N+Sx]), 0.5*(phys_temp[N+3*Sx] - phys_temp[N+Sx]), THETA*(phys_temp[N+3*Sx] - phys_temp[N+2*Sx]));
        }
      
  }
  
  FluxCalcPX (FLX, physL,(X+1));
  FluxCalcPX (FRX, physR,(X+1));
  Ucalc     (UL,           physL,(X+1));
  Ucalc     (UR,           physR,(X+1));

  for (i=0;i<X+1;i++) {
    N = Sx*i;
    AlphaPlus=0.;
    AlphaMinus=0.;
    SoundSpeedL = sqrt (GAMMA * physL[N+2] / physL[N+0]);
    SoundSpeedR = sqrt (GAMMA * physR[N+2] / physR[N+0]);
    if (AlphaPlus  <  physL[N+1] + SoundSpeedL ) AlphaPlus  =  physL[N+1] + SoundSpeedL;
    if (AlphaMinus < -physL[N+1] + SoundSpeedL ) AlphaMinus = -physL[N+1] + SoundSpeedL;
    if (AlphaPlus  <  physR[N+1] + SoundSpeedR ) AlphaPlus  =  physR[N+1] + SoundSpeedR;
    if (AlphaMinus < -physR[N+1] + SoundSpeedR ) AlphaMinus = -physR[N+1] + SoundSpeedR;

    for (l=0; l<3; l++) {
      N = Sx*i+l;
      F_mid[N] = (AlphaPlus*FLX[N]+AlphaMinus*FRX[N]-AlphaMinus*AlphaPlus*(UR[N]-UL[N]))/(AlphaPlus+AlphaMinus);
    }

    if (*max < AlphaPlus ) *max = AlphaPlus;
    if (*max < AlphaMinus) *max = AlphaMinus;
  }

  free (phys);
  free (phys_temp);
  free (physL);
  free (physR);
  free (FLX);
  free (FRX);
  free (UL);
  free (UR);
}


void Advance(double *U_new,double *U_old, double *dt) {
  int i=0, l=0, N;								// counters
  double *phys   = (double*) malloc (3*(X)*sizeof (double));
  double *FmidX  = (double*) malloc (3*(X+1)*sizeof (double));
  double *LU     = (double*) malloc (3*(X)*sizeof (double));
  double *U1     = (double*) malloc (3*(X)*sizeof (double));
  double *U2     = (double*) malloc (3*(X)*sizeof (double));
//  double *Source = (double*) malloc (3*(X)*sizeof (double));
  double max1, maxX=0;
  int Sx  = 3;
  int Sxx = 3;
  int N1, Nx;

  double  GridRatioX;
  Grid(&GridRatioX);
  GridRatioX = *dt/GridRatioX; 

  riemansolverX(FmidX,U_old,&max1);
  if (max1> maxX) maxX=max1;
//  Fluxsource(Source, U_old, dt);

  for (i=0;i<X;i++) {
    
    for (l=0; l<3; l++) {
      N  = i*Sx + l;
      N1 = i*Sx;
      Nx = i*Sxx + l;
      LU[N] = - GridRatioX * (FmidX[Nx+Sxx] - FmidX[Nx]);
      U1[N] = U_old[N] + LU[N];
    }
  }
  
  riemansolverX(FmidX,U1,&max1);
  if (max1> maxX) maxX=max1;
  //Fluxsource(Source, U1, dt);

  for (i=0;i<X;i++) {
    for (l=0; l<3; l++) {
    N  = i*Sx +  l;
    N1 = i*Sx ;
    Nx = i*Sxx + l;
    LU[N] = - GridRatioX * (FmidX[Nx+Sxx] - FmidX[Nx]);
    U2[N] = 0.75 * U_old[N] + 0.25 * U1[N] + 0.25 * LU[N];
    }
  }
  
  riemansolverX(FmidX,U2,&max1);
  if (max1> maxX) maxX=max1;
 // Fluxsource(Source, U2, dt);

  for (i=0;i<X;i++) {
    for (l=0; l<3; l++) {
    N  = i*Sx + l;
    N1 = i*Sx;
    Nx = i*Sxx + l;
    LU[N] = - GridRatioX * (FmidX[Nx+Sxx] - FmidX[Nx]);
    U_new[N] = (1./3.) * U_old[N] + (2./3.) * U2[N] + (2./3.) * LU[N]; //+ Source[N];
    }
  }
  
  Grid(&GridRatioX);
  double dtX = 0.3 * GridRatioX / maxX;
  *dt=dtX;

  free (phys);
  free (FmidX);
  free (LU);
  free (U1);
  free (U2);
  //free (Source);
}

void output_file(double *U,double t){

  int i=0, N;
  double xstep;
  Grid  (&xstep);
  int Sx = 3;
  
  char command[160];
  char name[160];
  double *phys = (double*) malloc (3*X*sizeof(double));
  Ucalcinv(phys, U, (X));
  

  sprintf(command, "rm density_%d_2d.dat",X);
  system(command);
  sprintf(name,"density_%d_2d.dat",X);
  FILE *dens = fopen(name, "a+");
    
  sprintf(command, "rm velocity_%d_2d.dat",X);
  system(command);
  sprintf(name, "velocity_%d_2d.dat",X);

  FILE *vel = fopen(name,"a+");
    
  sprintf(command, "rm pressure_%d_2d.dat",X);
  system(command);
  sprintf(name, "pressure_%d_2d.dat",X);

  FILE *press = fopen(name,"a+");

    


  for (i=0; i<X; i++) {
    N = Sx*i;
    fprintf(dens , "%f\t",phys[N]);
    fprintf(vel  , "%f\t",phys[N+1]);
    fprintf(press, "%f\t",phys[N+2]);
  }  
  
   
  fclose(dens);
  fclose(vel);
  fclose(press);

  sprintf(command, "cp density_%d_2d.dat %d/density/%f.dat",X,X,t);
  system(command);
  sprintf(command, "cp velocity_%d_2d.dat %d/velocity/%f.dat",X,X,t);
  system(command);
  sprintf(command, "cp pressure_%d_2d.dat %d/pressure/%f.dat",X,X,t);
  system(command);
  free (phys);

}

void Creat_folder(){
  char command[120];
  sprintf(command, "mkdir %d",X);
  system(command);
  sprintf(command, "mkdir %d/velocity",X);
  system(command);
  sprintf(command, "mkdir %d/density",X);
  system(command);
  sprintf(command, "mkdir %d/pressure",X);
  system(command);
}


int main (int argc, char **argv) {
  double *phys = (double*) malloc(3*X*sizeof(double));
  double *U    = (double*) malloc(3*X*sizeof(double));
  double *U_adv = (double*) malloc(3*X*sizeof(double));
  double dt=0.00001;
  double t=0;
  int i, k=0, l=0;

  Creat_folder();//Creat folder to store data;

   
  // Initial conditions:
  sphere(phys);
  Ucalc(U,phys, (X));
  printf ("tmax = %f\n",tmax);

#if 1

  while (t<tmax) {

    t+=dt;
    Advance(U_adv, U, &dt);
    for ( i=0;i<3*X;i++) {
      U[i]=U_adv[i];
    }

    if (t>0.05 * k) {
      output_file(U, t);
      k++;
    }
    printf ("%d| t = %f\n",l,t);
    l++;
  }

  free (phys);
  free (U);
  free (U_adv);

#endif
  return 0;
}
