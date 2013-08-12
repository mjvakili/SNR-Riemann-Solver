#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


#define PI    3.14159265
#define G     .1
#define GAMMA 1.4
#define X 1001
#define XMIN 0.1
#define XMAX 1.1
#define THETA 2.0
///#define tmax 0.1
#define R 0.5
#define Rc 0.05
void Grid(double *gridX) {
  *gridX = (XMAX - XMIN) / (X - 1);
}

double sgn(double x) {					// Gives the sign of x
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
      FX[3*i+2] = phys[3*i+1] *(phys[3*i+0] *phys[3*i+1] * phys[3*i+1] * .5 + phys[3*i+2]*GAMMA / (GAMMA - 1));

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

// Surfaces


void CellSurface(double *SX) {
  int i=0;
  for (i=0;i<X;i++) {
    double dx;
    Grid(&dx);
    double x;
    x = XMIN +dx*i;
    SX[i]=4*PI*x*x;
  }
}


void CellVolume(double *VX){

 int i=0;
 for (i=0;i<X;i++) {
   double dx;
   Grid(&dx);
   double x;
   x = XMIN + dx*(i+1);
   VX[i] = 4*PI*x*x*dx;
 }
}



   
///SNR initial condition

void constant(double *physical) {
  int i=0, N;
  double   dx;
  Grid (&dx);
  double x;

  for (i=0;i<X;i++) {
    
    x = XMIN + dx * i;
    N = 3 * i;
    if (x<R) {
      physical[N+0] = 1.0;
      physical[N+1] = 0.0;
      physical[N+2] = 1.0;
    }
    else {
      physical[N+0] = 0.1;
      physical[N+1] = 0.0;
      physical[N+2] = 0.125;
    }
  }
}

void sphere(double *physical) {
  int i=0, N;
  double   dx;
  Grid (&dx);
  double x;

  for (i=0;i<X;i++) {
    
    x = XMIN + dx * i;
    N = 3 * i;
    if (x<Rc) {
      physical[N+0] = 1.0;
      physical[N+1] = 0.0;
      physical[N+2] = 1.0;
    }
    else {
      physical[N+0] = (Rc*Rc)/(x*x);
      physical[N+1] = 0.0;
      physical[N+2] = 0.125;
    }
  }
}

double MAX2 (double a, double b){
  double z;
  if (a>b)
    z=a;
  else
    z=b;
  return z;
}


double MAX3 (double a, double b, double c){
    double z;
    if (a>=b && a>=c)
        z=a;
    else if (b>=a && b>=c)
        z=b;
    else if (c>=a && c>=b)
        z=c;
    return z;
    
}


// Calculates the flux at a half integer coordinate using the riemann method.
void riemansolverX(double *F_mid, double *U, double *S_mid, double *max) {
  double *phys  = (double*) malloc (3*(X)*sizeof (double));
  double *phys_temp= (double*) malloc (3*(X+4)*sizeof (double));
  double *physL = (double*) malloc (3*(X+1)*sizeof (double));
  double *physR = (double*) malloc (3*(X+1)*sizeof (double));
  double *FLX   = (double*) malloc (3*(X+1)*sizeof (double));
  double *FRX   = (double*) malloc (3*(X+1)*sizeof (double));
  double *UL    = (double*) malloc (3*(X+1)*sizeof (double));
  double *UR    = (double*) malloc (3*(X+1)*sizeof (double));
  double *SX    = (double*) malloc ((X)*sizeof (double));
  double *SLX   = (double*) malloc ((X+1)*sizeof (double));
  double *SRX   = (double*) malloc ((X+1)*sizeof (double));
  double *S_temp = (double*) malloc ((X+4)*sizeof (double));
  double SoundSpeedL, SoundSpeedR;											//
//L, SoundSpeedR, SoundSpeedLL, SoundSpeedRR;
  double AlphaPlus=0., AlphaMinus=0.;
  int i, l, N;  
  Ucalcinv (phys,U,(X));

  *max=0;
  int Sx =3;


/**********reflective at r=0 and outflow at rmax **** /                                    
  for (i=0;i<X+4;i++) {
    for (l=0;l<3;l++) {
      N = Sx*i+l;
      if      (i==0) {
        if    (l==1)  phys_temp[N] = -phys[Sx*(1)+l];
        else          phys_temp[N] = phys[Sx*(0)+l];
      }
      else if (i==1) {
        if    (l==1)  phys_temp[N] = -phys[Sx*(0)+l];
        else          phys_temp[N] = phys[Sx*(0)+l];
      }
      else if (i==X+2) {
        if    (l==1)  phys_temp[N] = phys[Sx*(X-1)+l];
        else          phys_temp[N] = phys[Sx*(X-1)+l];
      }
      else if (i==X+3) {
        if    (l==1) phys_temp[N] = phys[Sx*(X-2)+l];
        else         phys_temp[N] = phys[Sx*(X-1)+l];
      }
      else            phys_temp[N] = phys[Sx*(i-2)+l];
    }
  }
/**********reflective on both sides ***************/ 
  for (i=0;i<X+4;i++) {
    for (l=0;l<3;l++) {
      N = Sx*i+l;
      if      (i==0) {
        if    (l==1)  phys_temp[N] = -phys[Sx*(1)+l];
        else          phys_temp[N] = phys[Sx*(0)+l];
      }
      else if (i==1) {
        if    (l==1)  phys_temp[N] = -phys[Sx*(0)+l];
        else          phys_temp[N] = phys[Sx*(0)+l];
      }
      else if (i==X+2) {
        if    (l==1)  phys_temp[N] = -phys[Sx*(X-1)+l];
        else          phys_temp[N] = phys[Sx*(X-1)+l];
      }
      else if (i==X+3) {
        if    (l==1) phys_temp[N] = -phys[Sx*(X-2)+l];
        else         phys_temp[N] = phys[Sx*(X-1)+l];
      }
      else            phys_temp[N] = phys[Sx*(i-2)+l];
    }
  }
 
/************************************/
  for (i=0;i<X+4;i++) {
    double dx;
    Grid(&dx);
    double x;
    x = XMIN + dx*i;
    
  
    if (i==0){
      S_temp[i] = 0;
    }
    else if (i==1){
      S_temp[i] = 0;
    }
    else if (i==X+2){
      S_temp[i] = 4*PI*x*x;
    }
    else if (i==X+3){
      S_temp[i] = 4*PI*x*x;
    }
    else{
      S_temp[i] = SX[i-2];
    }
  }

  for (i=0;i<X+1;i++) {
	
      SLX[i] = S_temp[i+1] + 0.5 * minmod (THETA*(S_temp[i+1] - S_temp[i]), 0.5*(S_temp[i+2] -S_temp[i]), THETA*(S_temp[i+2] - S_temp[i+1]));
      SRX[i] = S_temp[i+2] - 0.5 * minmod (THETA*(S_temp[i+2] - S_temp[i+1]), 0.5*(S_temp[i+3] - S_temp[i+1]), THETA*(S_temp[i+3] - S_temp[i+2]));
  }
      
  

  for (i=0;i<X+1;i++) {
    for (l=0; l<3; l++) {
      N = Sx*i+l;
	
      physL[N] = phys_temp[N+Sx] + 0.5 * minmod (THETA*(phys_temp[N+Sx] - phys_temp[N]), 0.5*(phys_temp[N+2*Sx] -phys_temp[N]), THETA*(phys_temp[N+2*Sx] - phys_temp[N+Sx]));
      physR[N] = phys_temp[N+2*Sx] - 0.5 * minmod (THETA*(phys_temp[N+2*Sx] - phys_temp[N+Sx]), 0.5*(phys_temp[N+3*Sx] - phys_temp[N+Sx]), THETA*(phys_temp[N+3*Sx] - phys_temp[N+2*Sx]));
    }
      
  }
  

  FluxCalcPX (FLX, physL, (X+1));
  FluxCalcPX (FRX, physR, (X+1));
  Ucalc (UL, physL, (X+1));
  Ucalc (UR, physR, (X+1));

  for (i=0;i<X+1;i++) {
    N = Sx*i;
    AlphaPlus=0.;
    AlphaMinus=0.;
    SoundSpeedL = sqrt (GAMMA * physL[N+2] / physL[N+0]);
    SoundSpeedR = sqrt (GAMMA * physR[N+2] / physR[N+0]);
    AlphaMinus = MAX3(-physL[N+1] + SoundSpeedL , -physR[N+1] + SoundSpeedR , 0);
    AlphaPlus = MAX3(physR[N+1] + SoundSpeedR , physL[N+1] + SoundSpeedL , 0);
    S_mid[i] = (AlphaPlus*SLX[i]+AlphaMinus*SRX[i]-AlphaMinus*AlphaPlus*(UR[N]-UL[N]))/(AlphaPlus+AlphaMinus);

    for (l=0; l<3; l++) {
      N = Sx*i+l;
      F_mid[N] = (AlphaPlus*FLX[N]+AlphaMinus*FRX[N]-AlphaMinus*AlphaPlus*(UR[N]-UL[N]))/(AlphaPlus+AlphaMinus);
    }

    *max = MAX2 (AlphaPlus , AlphaMinus);
  }

  free (phys);
  free (phys_temp);
  free (physL);
  free (physR);
  free (FLX);
  free (FRX);
  free (UL);
  free (UR);
  free (S_mid);
  free (SX);
  free (SLX);   
  free (SRX);   
  free (S_temp); 
}




void Fluxsource (double *Source, double *U, double *dt) {
  
  
  double *S_mid= (double*) malloc ((X+1)*sizeof (double));
  double *phys = (double*) malloc (3*(X)*sizeof(double));
  double *VX   = (double*) malloc ((X)*sizeof (double));
  double *UL    = (double*) malloc (3*(X+1)*sizeof (double));
  double *UR    = (double*) malloc (3*(X+1)*sizeof (double));
  double *SX    = (double*) malloc ((X)*sizeof (double));
  double *SLX   = (double*) malloc ((X+1)*sizeof (double));
  double *SRX   = (double*) malloc ((X+1)*sizeof (double));
  double *S_temp = (double*) malloc ((X+4)*sizeof (double));
  double *phys_temp= (double*) malloc (3*(X+4)*sizeof (double));
  double *physL = (double*) malloc (3*(X+1)*sizeof (double));
  double *physR = (double*) malloc (3*(X+1)*sizeof (double));
  double SoundSpeedL, SoundSpeedR;
  double AlphaPlus=0., AlphaMinus=0.;
  int i, l, N;  
  int Sx =3; 							
  

  Ucalcinv(phys, U, (X));
  CellVolume(VX);
  CellSurface(SX);
  for (i=0;i<X+4;i++) {
    for (l=0;l<3;l++) {
      N = Sx*i+l;
      if      (i==0) {
        if    (l==1)  phys_temp[N] = -phys[Sx*(1)+l];
        else          phys_temp[N] = phys[Sx*(0)+l];
      }
      else if (i==1) {
        if    (l==1)  phys_temp[N] = -phys[Sx*(0)+l];
        else          phys_temp[N] = phys[Sx*(0)+l];
      }
      else if (i==X+2) {
        if    (l==1)  phys_temp[N] = -phys[Sx*(X-1)+l];
        else          phys_temp[N] = phys[Sx*(X-1)+l];
      }
      else if (i==X+3) {
        if    (l==1) phys_temp[N] = -phys[Sx*(X-2)+l];
        else         phys_temp[N] = phys[Sx*(X-1)+l];
      }
      else            phys_temp[N] = phys[Sx*(i-2)+l];
    }
  }
 
/************************************/
  for (i=0;i<X+4;i++) {
    double dx;
    Grid(&dx);
    double x;
    x = XMIN + dx*i;
    
  
    if (i==0){
      S_temp[i] = 0;
    }
    else if (i==1){
      S_temp[i] = 0;
    }
    else if (i==X+2){
      S_temp[i] = 4*PI*x*x;
    }
    else if (i==X+3){
      S_temp[i] = 4*PI*x*x;
    }
    else{
      S_temp[i] = SX[i-2];
    }
  }

  for (i=0;i<X+1;i++) {
	
      SLX[i] = S_temp[i+1] + 0.5 * minmod (THETA*(S_temp[i+1] - S_temp[i]), 0.5*(S_temp[i+2] -S_temp[i]), THETA*(S_temp[i+2] - S_temp[i+1]));
      SRX[i] = S_temp[i+2] - 0.5 * minmod (THETA*(S_temp[i+2] - S_temp[i+1]), 0.5*(S_temp[i+3] - S_temp[i+1]), THETA*(S_temp[i+3] - S_temp[i+2]));
  }
      
  

  for (i=0;i<X+1;i++) {
    for (l=0; l<3; l++) {
      N = Sx*i+l;
	
      physL[N] = phys_temp[N+Sx] + 0.5 * minmod (THETA*(phys_temp[N+Sx] - phys_temp[N]), 0.5*(phys_temp[N+2*Sx] -phys_temp[N]), THETA*(phys_temp[N+2*Sx] - phys_temp[N+Sx]));
      physR[N] = phys_temp[N+2*Sx] - 0.5 * minmod (THETA*(phys_temp[N+2*Sx] - phys_temp[N+Sx]), 0.5*(phys_temp[N+3*Sx] - phys_temp[N+Sx]), THETA*(phys_temp[N+3*Sx] - phys_temp[N+2*Sx]));
    }
      
  }
   
  
  for (i=0;i<X+1;i++) {
    N = Sx*i;
    AlphaPlus=0.;
    AlphaMinus=0.;
    SoundSpeedL = sqrt (GAMMA * physL[N+2] / physL[N+0]);
    SoundSpeedR = sqrt (GAMMA * physR[N+2] / physR[N+0]);
    AlphaMinus = MAX3(-physL[N+1] + SoundSpeedL , -physR[N+1] + SoundSpeedR , 0);
    AlphaPlus = MAX3(physR[N+1] + SoundSpeedR , physL[N+1] + SoundSpeedL , 0);
    S_mid[i] = (AlphaPlus*SLX[i]+AlphaMinus*SRX[i]-AlphaMinus*AlphaPlus*(UR[N]-UL[N]))/(AlphaPlus+AlphaMinus);

  }
  


  for (i=0;i<X;i++) {
	int Nf;
        Nf = 3*i;
        Source[Nf+0] = 0.0;
        Source[Nf+1] = phys[Nf+2]*(S_mid[i+1]-S_mid[i])*dt/VX[i];
        Source[Nf+2] = 0.0; 
  }

  free (phys);
  free (physL);
  free (physR);
  free (phys_temp);
  free (S_mid);
  free (phys);
  free (VX);
  free (UL);
  free (UR);
  free (SX);
  free (SLX);   
  free (SRX);   
  free (S_temp); 
  
}



void Advance(double *U_new,double *U_old, double *dt) {
  int i=0, l=0, N;								// counters
  double *phys   = (double*) malloc (3*(X)*sizeof (double));
  double *VX   = (double*) malloc ((X)*sizeof (double));
  double *FmidX  = (double*) malloc (3*(X+1)*sizeof (double));
  double *LU     = (double*) malloc (3*(X)*sizeof (double));
  double *U1     = (double*) malloc (3*(X)*sizeof (double));
  double *U2     = (double*) malloc (3*(X)*sizeof (double));
  double *Source = (double*) malloc (3*(X)*sizeof (double));
  double max1, maxX=0;
  int Sx  = 3;
  
  //double  GridRatioX;
  //Grid(&GridRatioX);
  //GridRatioX = *dt/GridRatioX;
  CellVolume(VX);

  riemansolverX(FmidX,U_old,Smid,&max1);
  if (max1> maxX) maxX=max1;
  Fluxsource(U_old, dt);

  for (i=0;i<X;i++) {
    for (l=0; l<3; l++) {
      N  = i*Sx + l;
      LU[N] = - dt * (FmidX[N+Sx]*Smid[i+1] - FmidX[N]*Smid[i])/VX[i];
      U1[N] = U_old[N] + LU[N];
    }
   
  
  }
  
  riemansolverX(FmidX,U1,Smid,&max1);
  if (max1> maxX) maxX=max1;
  Fluxsource(U1, dt);

  for (i=0;i<X;i++) { 
    for (l=0; l<3; l++) {
      N  = i*Sx +  l;
      LU[N] = -1.* dt * (FmidX[N+Sx]*Smid[i+1] - FmidX[N]*Smid[i])/VX[i];
      U2[N] = 0.75 * U_old[N] + 0.25 * U1[N] + 0.25 * LU[N];
    }
  }
  
  riemansolverX(FmidX,U2,Smid,&max1);
  if (max1> maxX) maxX=max1;
  Fluxsource(U2, dt);

  for (i=0;i<X;i++) {  
    for (l=0; l<3; l++) {
      N  = i*Sx + l;
      LU[N] = -1.* dt * (FmidX[N+Sx]*Smid[i+1] - FmidX[N]*Smid[i])/VX[i];
      U_new[N] = (1./3.) * U_old[N] + (2./3.) * U2[N] + (2./3.) * LU[N] + Source[N];
    }
  }
  double GridRatioX;
  Grid(&GridRatioX);
  double dtX = 0.7 * GridRatioX / maxX;
  *dt=dtX;

  free (phys);
  free (FmidX);
  free (LU);
  free (U1);
  free (U2);
  free (Source);
}

// Output:
void output_file1(double *U, double t, int stat) {
  int i=0, l, N;
  double xstep;
  Grid  (&xstep);
  int Sx = 3;
  double xproj;

  double *phys = (double*) malloc (3*(X)*sizeof(double));
  Ucalcinv(phys, U, (X));

  if (stat==0) system("rm numerical.dat");
  FILE *full = fopen("numerical.dat", "a+");
  for (i=0; i<X; i++) {
    N = Sx*i;
    fprintf(full, "%f\t%f\t%f\t%f\t%f\t\n",XMIN + i*xstep,t,phys[N+0],phys[N+1],phys[N+2]);
  }
  fprintf(full, "\n\n");
  fclose(full);
  free (phys);
}

void output_file(double *U,double t){

  int i=0, N;
  double xstep;
  Grid  (&xstep);
  int Sx = 3;
  
  char command[160];
  char name[160];
  double *phys = (double*) malloc (3*(X)*sizeof(double));
  Ucalcinv(phys, U, X);
  

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
  double *Fmid = (double*) malloc(3*(X+1)*sizeof(double));
  double dt=0.00001, dx;
  char command[40];
  double t=0, tmax=0.1;
  int i,j=0,k=0;

  Creat_folder();//Creat folder to store data;

   
  // Initial conditions:
  sphere(phys);
  Ucalc(U,phys,(X));
  double p;
  //riemansolverX(Fmid, U, &p);
  //for (i=0;i<X+1;i++) {
    //int N;
    //N = 3*i;
    //printf ("flux = %f\n",Fmid[N]);
  //}

  printf ("tmax = %f\n",tmax);
  
  while (t<tmax) {
    output_file1(U, t, j);
    t+=dt;
    Advance(U_adv, U, &dt);
    for (i=0;i<3*X;i++) {
      U[i]=U_adv[i];
    }
    j++;
    if (t>0.5 * k) {
      sprintf(command, "cp yproj.dat %d/%f.dat",X,t);
      system (command);
    k++;
    }
    printf ("%d| t = %f\n",k,t);
    output_file(U,t);
  }
 

  free (phys);
  free (U);
  free (U_adv);

  return 0;
}

