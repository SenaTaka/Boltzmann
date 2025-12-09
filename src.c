//*******************************************************************************************
//
// Non-dimentional boltzmann equation BGK model 2D3V Numerical simulation
// <shocktube>
//reference: Mieussens, L. and Struchtrup, H., “Numerical comparison of Bhatnagar-Gross-Krook models with proper Prandtl number,” Physics of Fluids, Vol. 16, No. 8, 2004, pp. 2797–2813.
//
//gcc -fopenmp -O2 src.c -o src -lm
//*******************************************************************************************
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>
//*******************************************************************************************
//definition constant value
//*******************************************************************************************
//mesh(2D)
//y-direction
#define xmin -4.0e-2
#define nx 55
#define Dx 2.0e-2
//y-direction
#define ymin -5.0e-2
#define ny 45
#define Dy 1.25e-2
//time (iteration)
//#define dt 1.0e-4
#define nstep 30
#define tstop 10.0
//micro velocity
#define cmin -5.0
#define ncj 100
#define nck 100
#define ncl 60
#define dcj 0.125
#define dck 0.125
#define dcl 0.2
//initial condition (macro volumes)
//left side
#define pu 1.000
#define Tu 1.000 
#define nu 1.000
#define v1u 4.564 
#define v2u 0.000
#define T1 1.000
//right side
//#define pd 0.100
//#define Td 0.800 
//#define nd 0.125
//#define v1d 0.000
//#define v2d 0.000
//local velocity distribution(newton method)
#define MAX_ITER_NEWTON 100    // ニュートン法の最大反復回数
#define TOLERANCE_NEWTON 1e-8  // ニュートン法の収束判定のための許容誤差
//initial condition (micro volume(gas:algon))
#define M 6.6339e-26
#define r 1.85e-10
//constant
#define pi 3.1415926535897932
#define BC 1.380658e-23
#define R BC/M
#define sigma 4*pi*r*r
#define Tref 273
#define myuref 2.117e-5
#define omega 0.81
#define nref 1.894e-6    //[kg/m^3]
#define pref 1.076e-1 
#define L 0.9     //reference length
#define l0 pow(pi*R*Tref/2.0,0.5)*myuref/pref
#define tref pow(2*R*Tref,0.5)/L
#define alpha 1.0
//*******************************************************************************************
//definition <global variables>
//*******************************************************************************************
//initial condition
double t,dt;
double *x,*y,*c1,*c2,*c3;
double **n,**v1,**v2,**p,**T;
double *****f,*****feq;
double **micVel_ave,**freq_coll,**mfp,**myu;
double ****nb1,***nb2;
//flux
double *****F, *****G;
double *****phi_x, *****phi_y;
//time step
double *****fn,*****fb;
//integral
double ****sum1_n,***sum2_n;
double ****sum1_v1,***sum2_v1;
double ****sum1_v2,***sum2_v2;
double ****sum1_p,***sum2_p;
//local velocity distribution
typedef struct {
    double a;     // a_i
    double Gamma; // Γ_i
    double gamma_x;
    double gamma_y; 
} Coefficients;
//*******************************************************************************************
//definition <function>
//*******************************************************************************************
void initialCondition(void);
void boundaryCondition(void);
void flux_x(void);
void flux_y(void);
void timeStep(void);
void replaceDistribution(void);
void integral(void);
void vol_re_coll(void);
static inline double calculate_f_gamma_approx(int j, int k, int l, const Coefficients* coeffs, double vi, double vi2);
static void solve_newton_for_coeffs(Coefficients* coeffs, int i, int i2);
void localf_eq_Newton(void); // 元のlocalf_eqと置き換える関数
double det3x3(
    double a11, double a12, double a13,
    double a21, double a22, double a23,
    double a31, double a32, double a33
);
void output_csv(int);
void output_csv_f(int);
double *****alloc5D(int nx1, int nc1, int nc2, int nc3, int nc4);
void free5D(double *****arr, int nx1, int nc1, int nc2, int nc3);
double ****alloc4D(int nx1, int nc1, int nc2, int nc3);
void free4D(double ****arr, int nx1, int nc1, int nc2);
double ***alloc3D(int nx1, int nc1, int nc2);
void free3D(double ***arr, int nx1, int nc1);
double **alloc2D(int nx1, int nc1);
void free2D(double **arr, int nx1);
static inline double max3(double a, double b, double c);
static inline double min3(double a, double b, double c);
static inline double minmod(double a, double b, double c);
void fai_x(void);
void fai_y(void);
void cfl_check(void);
void timeStep_RK2_1st(void);
void timeStep_RK2_2nd(void);
//*******************************************************************************************
//main program
//*******************************************************************************************
int main(void){
    int i,j,m;
//initial Condition
    x   = (double *)malloc(nx * sizeof(double));
    y   = (double *)malloc(ny * sizeof(double));
    c1  = (double *)malloc(ncj * sizeof(double));
    c2  = (double *)malloc(nck * sizeof(double));
    c3  = (double *)malloc(ncl * sizeof(double));
    n   = alloc2D(nx, ny);
    T   = alloc2D(nx, ny);
    p   = alloc2D(nx, ny);
    v1  = alloc2D(nx, ny);
    v2  = alloc2D(nx, ny);
    f   = alloc5D(nx, ny, ncj, nck, ncl);
    feq  = alloc5D(nx, ny, ncj, nck, ncl);
    micVel_ave = alloc2D(nx, ny);
    freq_coll  = alloc2D(nx, ny);
    mfp        = alloc2D(nx, ny);
    myu        = alloc2D(nx, ny);
    nb1        = alloc4D(nx, ny, ncj,ncl);
    nb2        = alloc3D(nx, ny, ncl);
//flux
    F   = alloc5D(nx, ny, ncj, nck, ncl);
    G   = alloc5D(nx, ny, ncj, nck, ncl);
//phi (limiter)
    phi_x = alloc5D(nx, ny, ncj, nck, ncl);
    phi_y = alloc5D(nx, ny, ncj, nck, ncl);
//time step
    fb  = alloc5D(nx, ny, ncj, nck, ncl);
    fn  = alloc5D(nx, ny, ncj, nck, ncl);
//integral
//density
    sum1_n   = alloc4D(nx, ny, nck, ncl);
    sum2_n   = alloc3D(nx, ny, ncl);
//velosity(x-direction)
    sum1_v1   = alloc4D(nx, ny, nck, ncl);
    sum2_v1   = alloc3D(nx, ny, ncl);
//velosity(y-direction)
    sum1_v2   = alloc4D(nx, ny, ncj, ncl);
    sum2_v2   = alloc3D(nx, ny, ncl);
//pressure
    sum1_p   = alloc4D(nx, ny, nck, ncl);
    sum2_p   = alloc3D(nx, ny, ncl);
printf("malloc OK!\n");
    initialCondition();
printf("initial OK!\n");
    omp_set_num_threads(20);
    t=0.0;
    m=1;
    for(m=1;m<=nstep;m++){//while(t<tstop){
        printf("===step : %d\n",m);
        boundaryCondition();
        printf("boundary OK!\n");
        localf_eq_Newton();
        printf("local f OK!\n");
        fai_x();
        flux_x();
        fai_y();
        flux_y();
        printf("flux OK!\n");
        cfl_check();
        timeStep_RK2_1st();
        printf("timestep RK1 OK!\n");
        integral();
        printf("integral OK!\n");
        vol_re_coll();

        boundaryCondition();
        printf("boundary OK!\n");
        localf_eq_Newton();
        printf("local f OK!\n");
        fai_x();
        flux_x();
        fai_y();
        flux_y();
        printf("flux OK!\n");
        cfl_check();
        timeStep_RK2_2nd();
        printf("timestep RK2 OK!\n");
        replaceDistribution();
        printf("replaceD OK!\n");
        integral();
        printf("integral OK!\n");
        vol_re_coll();

        if(m%1000==0){
            output_csv(m);  
        }
        if(m%5000==0){ 
            output_csv_f(m);
        }
        t=t+dt;
        m=m+1;
    }
    output_csv(m);  
    //output_csv_f(m);

//initial condition
    free(x); 
    free(c1); 
    free(c2); 
    free(c3);
    free(n); 
    free(T); 
    free(p); 
    free(v1);
    free5D(f, nx, ny, ncj, nck);
    free5D(feq, nx, ny, ncj, nck);
    free4D(nb1, nx, ny, ncj);
    free3D(nb2, nx, ny);
    free2D(micVel_ave, nx); 
    free2D(freq_coll, nx); 
    free2D(mfp, nx); 
    free2D(myu,nx);
//flux
    free5D(F, nx, ny, ncj, nck);
    free5D(G, nx, ny, ncj, nck);
//phi (limiter)
    free5D(phi_x, nx, ny, ncj, nck);
    free5D(phi_y, nx, ny, ncj, nck);
//time step
    free5D(fb, nx, ny, ncj, nck);
    free5D(fn, nx, ny, ncj, nck);
//integral
//density
    free4D(sum1_n, nx, ny, nck);
    free3D(sum2_n, nx, ny);
//velocity(x-direction)
    free4D(sum1_v1, nx, ny, nck);
    free3D(sum2_v1, nx, ny);
//velocity(y-direction)
    free4D(sum1_v2, nx, ny, ncj);
    free3D(sum2_v2, nx, ny);
//pressure
    free4D(sum1_p, nx, ny, nck);
    free3D(sum2_p, nx, ny);
    return 0;
}
//*******************************************************************************************
//function
//*******************************************************************************************
//****************************************************************************
//initial Condition
//shock tube (2D3V)
//****************************************************************************
void initialCondition(void){
    int i,i2,j,k,l;
//region(real space (2-dimentional))
//x-direction
    #pragma omp parallel for 
    for(i=0;i<nx;i++){
        x[i]=xmin+i*Dx;
    }
//y-direction
    #pragma omp parallel for 
    for(i2=0;i2<ny;i2++){
        y[i2]=ymin+i2*Dy;
    }
//micro velocity
//x-direction
    #pragma omp parallel for 
    for(j=0;j<ncj;j++){
        c1[j]=cmin+j*dcj; 
       }
//y-direction
       #pragma omp parallel for 
       for(k=0;k<nck;k++){
        c2[k]=cmin+k*dck;
       }
//z-direction
       #pragma omp parallel for 
       for(l=0;l<ncl;l++){
        c3[l]=cmin+l*dcl;
       }
//macro physical volumes
       #pragma omp parallel for collapse(2)
       for(i=0;i<nx;i++){
            for(i2=0;i2<ny;i2++){
                if(i2 == 0){
                    T[i][0]=T1;
                    n[i][0]=0.0;

                    for(j=0;j<ncj;j++){
                        for(l=0;l<ncl;l++){
                            nb1[i][0][j][l]=0.0;
                            nb2[i][0][l]=0.0;
                        }
                    }
                    
                    for(j=0;j<ncj;j++){
                        for(k=0;k<0.5*nck-1;k++){
                            for(l=0;l<ncl;l++){
                                nb1[i][0][j][l]+=0.5*(c2[k]*f[i][2][j][k][l]+c2[k+1]*f[i][2][j][k+1][l])*dck;
                            }
                        }
                    }
                    for(j=0;j<ncj-1;j++){
                        for(l=0;l<ncl;l++){
                            nb2[i][0][l]+=0.5*(nb1[i][0][j][l]+nb1[i][0][j+1][l])*dcj;
                        }
                    }
                    for(l=0;l<ncl-1;l++){
                        n[i][0]+=(-1.0)*pow(pi/T1,0.5)*2.0*0.5*(nb2[i][0][l]+nb2[i][0][l+1])*dcl;
                    }
                    p[i][0]=n[i][0]*T[i][0];
                }
                //if(i<0.5*nx){
                else{
                    n[i][i2]=nu;
                    v1[i][i2]=v1u;
                    v2[i][i2]=v2u;
                    p[i][i2]=pu;
                    T[i][i2]=Tu;
                }
                //}
                /*else{
                    n[i][i2]=nd;
                    v1[i][i2]=v1d;
                    v2[i][i2]=v2d;
                    p[i][i2]=pd;
                    T[i][i2]=Td;
                }*/
            }
        }
//velocity distribution
     #pragma omp parallel for collapse(2)
       for(i=0;i<nx;i++){
            for(i2=0;i2<ny;i2++){
                for(j=0;j<ncj;j++){
                    for(k=0;k<nck;k++){
                        for(l=0;l<ncl;l++){
                            f[i][i2][j][k][l]=n[i][i2]*pow(pi*T[i][i2],-1.5)*exp(-((c1[j]-v1[i][i2])*(c1[j]-v1[i][i2])+(c2[k]-v2[i][i2])*(c2[k]-v2[i][i2])+(c3[l]*c3[l]))/T[i][i2]);
                        }
                    }
                }
            }
        }
//local velocity distribution
        #pragma omp parallel for collapse(2)
        for(i=0;i<nx;i++){
            for(i2=0;i2<ny;i2++){
                for(j=0;j<ncj;j++){
                    for(k=0;k<nck;k++){
                        for(l=0;l<ncl;l++){
                            feq[i][i2][j][k][l]=n[i][i2]*pow(pi*T[i][i2],-1.5)*exp(-((c1[j]-v1[i][i2])*(c1[j]-v1[i][i2])+(c2[k]-v2[i][i2])*(c2[k]-v2[i][i2])+(c3[l]*c3[l]))/T[i][i2]);
                        }
                    }
                }
            }
        }
//average micro velocity
#pragma omp parallel for private(i, i2) collapse(2)
       for(i=0;i<nx;i++){
            for(i2=0;i2<ny;i2++){
                micVel_ave[i][i2]=pow(8.0*R*Tref*T[i][i2]/pi,0.5);
            }
       }
//viscosity : myu
 #pragma omp parallel for private(i, i2) collapse(2)
       for(i=0;i<nx;i++){
            for(i2=0;i2<ny;i2++){
                myu[i][i2]=myuref*pow(T[i][i2],omega);
            }
        }
//mean free path
#pragma omp parallel for private(i, i2) collapse(2)
       for(i=0;i<nx;i++){
            for(i2=0;i2<ny;i2++){
                mfp[i][i2]=pow(pi*R*Tref*T[i][i2]/2.0,0.5)*myu[i][i2]/(pref*p[i][i2]);
            }
       }
//collsion frequency
#pragma omp parallel for private(i, i2) collapse(2)
        for(i=0;i<nx;i++){
            for(i2=0;i2<ny;i2++){
                freq_coll[i][i2]=micVel_ave[i][i2]/mfp[i][i2];
            }
        }
}
//****************************************************************************
//Boundary Condition
//non-gradient(i=0,nx-1)
//****************************************************************************
void boundaryCondition(void){
    int i,i2,j,k,l;
    double f_diff,f_spec;
//wall density
    #pragma omp parallel for private(j, k, l)
    for(i=0;i<nx;i++){
        T[i][0]=T1;
        n[i][0]=0.0;

        for(j=0;j<ncj;j++){
            for(l=0;l<ncl;l++){
                nb1[i][0][j][l]=0.0;
                nb2[i][0][l]=0.0;
            }
        }

        for(j=0;j<ncj;j++){
            for(k=0;k<0.5*nck-1;k++){
                for(l=0;l<ncl;l++){
                    nb1[i][0][j][l]+=0.5*(c2[k]*f[i][2][j][k][l]+c2[k+1]*f[i][2][j][k+1][l])*dck;
                }
            }
        }
        for(j=0;j<ncj-1;j++){
            for(l=0;l<ncl;l++){
                nb2[i][0][l]+=0.5*(nb1[i][0][j][l]+nb1[i][0][j+1][l])*dcj;
            }
        }
        for(l=0;l<ncl-1;l++){
            n[i][0]+=(-1.0)*2.0*pow(pi/T1,0.5)*0.5*(nb2[i][0][l]+nb2[i][0][l+1])*dcl;
        }
    }

//left side
    #pragma omp parallel for private(j, k, l)
    for(i2=0;i2<ny;i2++){
        for(j=0;j<ncj;j++){
            for(k=0;k<nck;k++){
                for(l=0;l<ncl;l++){
                    double val = nu*pow(pi*Tu,-1.5)*exp(-((c1[j]-v1u)*(c1[j]-v1u)+(c2[k]-v2u)*(c2[k]-v2u)+(c3[l]*c3[l]))/Tu);
                    f[0][i2][j][k][l]= val;
                    f[1][i2][j][k][l]= val;
                }   
            }
        }
    }
//right side
     #pragma omp parallel for private(j, k, l)
    for(i2=0;i2<ny;i2++){
        for(j=0;j<ncj;j++){
            for(k=0;k<nck;k++){
                for(l=0;l<ncl;l++){
                    f[nx-1][i2][j][k][l]=f[nx-3][i2][j][k][l];
                    f[nx-2][i2][j][k][l]=f[nx-3][i2][j][k][l];
                }
            }
        } 
    }
//upside
    #pragma omp parallel for private(j, k, l)
    for(i=0;i<nx;i++){
        for(j=0;j<ncj;j++){
            for(k=0;k<nck;k++){
                for(l=0;l<ncl;l++){
                    f[i][ny-1][j][k][l]=f[i][ny-3][j][k][l];
                    f[i][ny-2][j][k][l]=f[i][ny-3][j][k][l];
                }
            }
        } 
    }
//downside
    #pragma omp parallel for private(j, k, l, f_spec, f_diff)
    for(i=0;i<nx;i++){
        for(j=0;j<ncj;j++){
            for(k=0;k<nck;k++){
                for(l=0;l<ncl;l++){
                    if(x[i]<=0.1){
                        if(c2[k]<0.0){
                            f[i][1][j][k][l]=f[i][2][j][k][l];
                        }
                        else{
                            f[i][1][j][k][l]=f[i][2][j][nck-k-1][l];
                        }
                        f[i][0][j][k][l]=f[i][1][j][k][l];
                    }
                    else{
                        if(c2[k]<0.0){
                            f[i][1][j][k][l]=f[i][2][j][k][l];
                            f[i][0][j][k][l]=f[i][1][j][k][l];
                        }
                        else{
                            f_spec = f[i][2][j][nck-k-1][l];
                            f_diff = n[i][0]*pow(pi*T[i][0],-1.5)*exp(-(c1[j]*c1[j]+c2[k]*c2[k]+c3[l]*c3[l])/T[i][0]);
                            f[i][1][j][k][l] = (1.0 - alpha)*f_spec + alpha * f_diff;
                            f[i][0][j][k][l] = f[i][1][j][k][l];
                        }
                    }
                }
            }
        } 
    }
}
//****************************************************************************
//flux (x-direction)
//one-order upwind difference scheme
//****************************************************************************
void flux_x(void){
    int i,i2,j,k,l;
    #pragma omp parallel for private(i2,j, k, l)
    for(i=0;i<nx-1;i++){
        for(i2=0;i2<ny;i2++){
            for(j=0;j<ncj;j++){
                for(k=0;k<nck;k++){
                    for(l=0;l<ncl;l++){
                        F[i][i2][j][k][l]=0.5*c1[j]*(f[i+1][i2][j][k][l]+f[i][i2][j][k][l])-0.5*fabs(c1[j])*(f[i+1][i2][j][k][l]-f[i][i2][j][k][l]-phi_x[i][i2][j][k][l]);
                    }
                }
            }
        }
    }
}
//****************************************************************************
//flux (y-direction)
//one-order upwind difference scheme
//****************************************************************************
void flux_y(void){
    int i,i2,j,k,l;
    #pragma omp parallel for private(i2,j, k, l)
    for(i=0;i<nx;i++){
        for(i2=0;i2<ny-1;i2++){
            for(j=0;j<ncj;j++){
                for(k=0;k<nck;k++){
                    for(l=0;l<ncl;l++){
                        G[i][i2][j][k][l]=0.5*c2[k]*(f[i][i2+1][j][k][l]+f[i][i2][j][k][l])-0.5*fabs(c2[k])*(f[i][i2+1][j][k][l]-f[i][i2][j][k][l]-phi_y[i][i2][j][k][l]);
                    }
                }
            }
        }
    }
}
//****************************************************************************
//time step
//2nd order explicit scheme
//****************************************************************************
// 1段階目: 予測ステップ (Predictor)
// f^n から f^* (仮の次の値) を計算して fn に格納
void timeStep_RK2_1st(void){
    int i, i2, j, k, l;

    #pragma omp parallel for private(i2, j, k, l)
    for(i = 2; i < nx-1; i++){
        for(i2 = 2; i2 < ny-1; i2++){
            for(j = 0; j < ncj; j++){
                for(k = 0; k < nck; k++){
                    for(l = 0; l < ncl; l++){
                        // RHS(f^n) の計算
                        double rhs = -1.0/Dx * (F[i][i2][j][k][l] - F[i-1][i2][j][k][l])
                                     -1.0/Dy * (G[i][i2][j][k][l] - G[i][i2-1][j][k][l])
                                     -freq_coll[i][i2]/tref * (f[i][i2][j][k][l] - feq[i][i2][j][k][l]);
                        
                        fb[i][i2][j][k][l]=f[i][i2][j][k][l];
                        // f^* = f^n + dt * RHS(f^n)
                        f[i][i2][j][k][l] = fb[i][i2][j][k][l] + dt * rhs;
                    }
                }
            }
        }
    }
}

// 2段階目: 修正ステップ (Corrector)
// f^n と f^* (現在のfn) を用いて f^{n+1} を計算し fn に上書き
// ※注意: この関数を呼ぶ前に、Step1で求めたfnを使って F, G, feq を更新しておく必要があります
void timeStep_RK2_2nd(void){
    int i, i2, j, k, l;

    #pragma omp parallel for private(i2, j, k, l)
    for(i = 2; i < nx-1; i++){
        for(i2 = 2; i2 < ny-1; i2++){
            for(j = 0; j < ncj; j++){
                for(k = 0; k < nck; k++){
                    for(l = 0; l < ncl; l++){
                        // RHS(f^*) の計算
                        // ここでの F, G, feq は「Step1の結果(fn)」に基づいて更新されたものを使う
                        double rhs_star = -1.0/Dx * (F[i][i2][j][k][l] - F[i-1][i2][j][k][l])
                                          -1.0/Dy * (G[i][i2][j][k][l] - G[i][i2-1][j][k][l])
                                          -freq_coll[i][i2]/tref * (f[i][i2][j][k][l] - feq[i][i2][j][k][l]);
                        
                        // f^{n+1} = 0.5 * (f^n + f^* + dt * RHS(f^*))
                        // 数式的には f^{n+1} = f^n + 0.5*dt*(RHS(f^n) + RHS(f^*)) と等価です
                        fn[i][i2][j][k][l] = 0.5 * (fb[i][i2][j][k][l] + f[i][i2][j][k][l] + dt * rhs_star);
                    }
                }
            }
        }
    }
}
//****************************************************************************
//cfl check
//****************************************************************************
void cfl_check(void) {
    int i, i2, j, k;
    double max_nu = 0.0;
    double max_cx = 0.0;
    double max_cy = 0.0;
    double cfl = 0.25; // 安全係数 (通常 0.8 ~ 0.95)
    double denominator;

    // 1. 最大衝突周波数 (max v_i,j) の探索
    // freq_coll が無次元化されている場合は /tref が必要かどうか確認してください
    // ここでは freq_coll[i][i2]/tref が実際の周波数(1/s)の次元を持つと仮定します
    
    // OpenMP 3.1以降なら reduction(max:...) が使えます
    #pragma omp parallel for reduction(max: max_nu) private(i2)
    for(i = 0; i < nx; i++) {
        for(i2 = 0; i2 < ny; i2++) {
            double nyu = freq_coll[i][i2] / tref; 
            if(nyu > max_nu) {
                max_nu = nyu;
            }
        }
    }

    // 2. 最大粒子速度 (max |cx|, |cy|) の探索
    // x方向
    for(j = 0; j < ncj; j++) {
        if(fabs(c1[j]) > max_cx) {
            max_cx = fabs(c1[j]);
        }
    }
    // y方向
    for(k = 0; k < nck; k++) {
        if(fabs(c2[k]) > max_cy) {
            max_cy = fabs(c2[k]);
        }
    }

    // 3. dt の計算
    // 分母 = (衝突項の最大値) + (x方向の移流項の最大値) + (y方向の移流項の最大値)
    denominator = max_nu + (max_cx / Dx) + (max_cy / Dy);

    if(denominator > 0.0) {
        dt = cfl / denominator;
    } else {
        // 万が一分母が0の場合（初期状態など）、安全な小さい値を入れるかエラー処理
        dt = 1.0e-6; 
    }

    // 4. 計算終了時刻を超えないよう調整
    if ((t + dt > tstop) && (t < tstop)) {
        dt = tstop - t;
    }

    // デバッグ用（必要ならコメントアウトを外す）
    printf("dt determined: %e (max_nu: %e, Adv_x: %e, Adv_y: %e)\n", 
            dt, max_nu, max_cx/Dx, max_cy/Dy);
}
//****************************************************************************
//replace distrubution
//****************************************************************************
void replaceDistribution(void){
    int i,i2,j,k,l;
    #pragma omp parallel for private(i2,j, k, l)
    for(i=2;i<nx-2;i++){
        for(i2=2;i2<ny-2;i2++){
            for(j=0;j<ncj;j++){
                for(k=0;k<nck;k++){
                    for(l=0;l<ncl;l++){
                        f[i][i2][j][k][l]=fn[i][i2][j][k][l];
                    }
                }
            }
        }
    }
}
//****************************************************************************
//integral macro volume(number density,temprature,pressure)
//numerical integration : trapezoidal rule
//****************************************************************************
void integral(void){
    int i,i2,j,k,l;
//initialize macro volume
#pragma omp parallel private(i,i2,j, k, l)
{
    #pragma omp for
    for(i=0;i<nx;i++){
        for(i2=0;i2<ny;i2++){
            T[i][i2]=0.0;
            n[i][i2]=0.0;
            v1[i][i2]=0.0;
            v2[i][i2]=0.0;
            p[i][i2]=0.0;
                for(l=0;l<ncl;l++){
                    sum2_n[i][i2][l]=0.0;
                    sum2_v1[i][i2][l]=0.0;
                    sum2_v2[i][i2][l]=0.0;
                    sum2_p[i][i2][l]=0.0;
                        for(k=0;k<nck;k++){
                            sum1_n[i][i2][k][l]=0.0;
                            sum1_v1[i][i2][k][l]=0.0;
                            sum1_p[i][i2][k][l]=0.0;
                        }
                }
                for(j=0;j<ncj;j++){
                        for(l=0;l<ncl;l++){
                            sum1_v2[i][i2][j][l]=0.0;
                    }
                }
        }
    }
//integral calculation
//number density
    #pragma omp for
    for(i=0;i<nx;i++){
        for(i2=0;i2<ny;i2++){
            for(j=0;j<ncj-1;j++){
                for(k=0;k<nck;k++){
                    for(l=0;l<ncl;l++){
                        sum1_n[i][i2][k][l]+=0.5*(f[i][i2][j][k][l]+f[i][i2][j+1][k][l])*dcj;
                    }
                }
            }
            for(k=0;k<nck-1;k++){
                for(l=0;l<ncl;l++){
                    sum2_n[i][i2][l]+=0.5*(sum1_n[i][i2][k][l]+sum1_n[i][i2][k+1][l])*dck;
                }
            }
            for(l=0;l<ncl-1;l++){
                n[i][i2]+=0.5*(sum2_n[i][i2][l]+sum2_n[i][i2][l+1])*dcl;
            }
        }
    }
//velocity(x-direction)
    #pragma omp for 
    for(i=0;i<nx;i++){
        for(i2=0;i2<ny;i2++){
            for(j=0;j<ncj-1;j++){
                for(k=0;k<nck;k++){
                    for(l=0;l<ncl;l++){
                        sum1_v1[i][i2][k][l]+=0.5*(c1[j]*f[i][i2][j][k][l]+c1[j+1]*f[i][i2][j+1][k][l])*dcj;
                    }
                }
            }
            for(k=0;k<nck-1;k++){
                for(l=0;l<ncl;l++){
                    sum2_v1[i][i2][l]+=0.5*(sum1_v1[i][i2][k][l]+sum1_v1[i][i2][k+1][l])*dck;
                }
            }
            for(l=0;l<ncl-1;l++){
                v1[i][i2]+=0.5*(sum2_v1[i][i2][l]+sum2_v1[i][i2][l+1])*dcl/n[i][i2];
            }
        }
    }
//velocity(y-direction)
    #pragma omp for 
    for(i=0;i<nx;i++){
        for(i2=0;i2<ny;i2++){
            for(j=0;j<ncj;j++){
                for(k=0;k<nck-1;k++){
                    for(l=0;l<ncl;l++){
                        sum1_v2[i][i2][k][l]+=0.5*(c2[k]*f[i][i2][j][k][l]+c2[k+1]*f[i][i2][j][k+1][l])*dck;
                    }
                }
            }
            for(j=0;j<ncj-1;j++){
                for(l=0;l<ncl;l++){
                    sum2_v2[i][i2][l]+=0.5*(sum1_v2[i][i2][j][l]+sum1_v2[i][i2][j+1][l])*dcj;
                }
            }
            for(l=0;l<ncl-1;l++){
                v2[i][i2]+=0.5*(sum2_v2[i][i2][l]+sum2_v2[i][i2][l+1])*dcl/n[i][i2];
            }
        }
    }
//pressure
    #pragma omp for 
    for(i=0;i<nx;i++){
        for(i2=0;i2<ny;i2++){
            for(j=0;j<ncj-1;j++){
                for(k=0;k<nck;k++){
                    for(l=0;l<ncl;l++){
                        sum1_p[i][i2][k][l]+=0.5*(pow((c1[j]-v1[i][i2]),2)*f[i][i2][j][k][l]+pow((c1[j+1]-v1[i][i2]),2)*f[i][i2][j+1][k][l])*dcj;
                    }
                }
            }
            for(k=0;k<nck-1;k++){
                for(l=0;l<ncl;l++){
                    sum2_p[i][i2][l]+=0.5*(sum1_p[i][i2][k][l]+sum1_p[i][i2][k+1][l])*dck;
                }
            }
            for(l=0;l<ncl-1;l++){
                p[i][i2]+=0.5*(sum2_p[i][i2][l]+sum2_p[i][i2][l+1])*dcl*2.0;
            }
        }
    }   
    //temparature
    for(i=0;i<nx;i++){
        for(i2=0;i2<ny;i2++){
            T[i][i2]=p[i][i2]/n[i][i2];
        }
    }
}
}
//****************************************************************************
//local velocity distribution (Newton method version)
//****************************************************************************
void localf_eq_Newton(void){
    int i, i2, j, k, l;

    // --- 空間ループ: 各空間点 (i,i2) について計算を実行 ---
    #pragma omp parallel for collapse(2) private(j, k, l) schedule(dynamic)
    for (i = 0; i < nx; i++) {
        for(i2=0;i2<ny;i2++){
        // 1. ニュートン法で解くべき係数 (a_i, Γ_i, γ_i) の初期値を設定
        // マクスウェル分布に近い値から始めることで収束を速める
        Coefficients coeffs;
        coeffs.a = n[i][i2]*pow(pi * T[i][i2], -1.5);
        coeffs.Gamma = 1.0 / T[i][i2];
        coeffs.gamma_x = 0.0;
        coeffs.gamma_y = 0.0;
        // 2. ニュートン法を実行して、現在の空間点iの係数を収束させる
        solve_newton_for_coeffs(&coeffs,i,i2);
        // 3. 収束した係数を使い、全ての速度点で f_gamma を計算して feq に格納
            for (j = 0; j < ncj; j++) {
                for (k = 0; k < nck; k++) {
                    for (l = 0; l < ncl; l++) {
                        feq[i][i2][j][k][l] = calculate_f_gamma_approx(j, k, l, &coeffs, v1[i][i2],v2[i][i2]);
                    }
                }
            }
        }
    }
}
/**
 * @brief (補助関数) 近似分布関数 f_γ を計算する
 * @param j,k,l 速度グリッドのインデックス
 * @param coeffs 計算済みの係数 (a_i, Γ_i, γ_i,γ_i2)
 * @param vi,vi2 空間点iにおける平均速度 v^n_i,v^n_i2
 * @return 計算された f_γ の値
 */
static inline double calculate_f_gamma_approx(int j, int k, int l, const Coefficients* coeffs, double vi,double vi2) {
    // グローバル変数から速度座標値を取得
    double cx = c1[j];
    double cy = c2[k];
    double cz = c3[l];

    // C^2_{j} = (c_{x} - v_i)^2 + c_{y}^2 + c_{z}^2
    double C2 = pow(cx - vi, 2) + pow(cy-vi2, 2) + pow(cz, 2);
    
    // 指数部分: -Γ_i * C^2 + γ_i * (c_x - v_i)
    double exponent = -coeffs->Gamma * C2 + coeffs->gamma_x * (cx - vi)+ coeffs->gamma_y * (cy - vi2);
    
    // f_γ = a_i * exp(...)
    return coeffs->a * exp(exponent);
}

/**
 * @brief (補助関数) ニュートン法を用いて非線形連立方程式を解き、係数を求める
 * @param[in,out] coeffs 係数の初期値を与え、計算結果がこの中に格納される
 * @param i 現在計算中の空間グリッドのインデックス
 */
static void solve_newton_for_coeffs(Coefficients* coeffs, int i,int i2) {
    // 速度空間の離散化体積
    int iter;
    const double dV = dcj * dck * dcl;

    //printf("\n");
    //printf("calculating Neewton method.....\n");
    //printf("\n");

    for (iter = 0; iter < MAX_ITER_NEWTON; ++iter) {
        double F[4] = {0.0};
        double J[4][4] = {{0.0}};
         //#pragma omp parallel for reduction(+:F[:4], J[:4][:4])
        // --- 速度グリッド全体で和を計算 (Σ_{j,k,l}) ---
        for (int j = 0; j < ncj; j++) {
            for (int k = 0; k < nck; k++) {
                for (int l = 0; l < ncl; l++) {
                    
                    double cur_cx = c1[j];
                    double cur_cy = c2[k];
                    double cur_cz = c3[l];

                    // 現在の係数でf_gammaを計算
                    double f_gamma = calculate_f_gamma_approx(j, k, l, coeffs, v1[i][i2],v2[i][i2]);
                    //printf("fgamma OK!\n");
                    // 差分 (f^n_{i,j} - f^n_{γ,i,j})
                    double diff = f[i][i2][j][k][l] - f_gamma;

                    // 各項に共通する係数部分 ν^n_{i,j} * m * Δc
                    // userコードではνが速度に依存しない freq_coll[i] のため、それを使用
                    double common_prefix = freq_coll[i][i2]  * dV;
                    
                    double c_sq = cur_cx*cur_cx + cur_cy*cur_cy + cur_cz*cur_cz;

                    // --- 残差ベクトル F の計算 ---
                    F[0] += common_prefix * diff;
                    F[1] += common_prefix * cur_cx * diff;
                    F[2] += common_prefix * cur_cy * diff;
                    F[3] += common_prefix * (c_sq / 2.0) * diff;
                    //printf("F OK!\n");
                    // --- ヤコビアン行列 J の計算 ---
                    double C2 = pow(cur_cx - v1[i][i2], 2) + pow(cur_cy-v2[i][i2], 2) + pow(cur_cz, 2);
                    double df_da = f_gamma / coeffs->a;
                    double df_dGamma = -C2 * f_gamma;
                    double df_dgamma_x = (cur_cx - v1[i][i2]) * f_gamma;
                    double df_dgamma_y = (cur_cy - v2[i][i2]) * f_gamma;

                    J[0][0] += common_prefix * df_da;
                    J[0][1] += common_prefix * df_dGamma;
                    J[0][2] += common_prefix * df_dgamma_x;
                    J[0][3] += common_prefix * df_dgamma_y;
//-------------------------------------------------
                    J[1][0] += common_prefix * cur_cx * df_da;
                    J[1][1] += common_prefix * cur_cx * df_dGamma;
                    J[1][2] += common_prefix * cur_cx * df_dgamma_x;
                    J[1][3] += common_prefix * cur_cx * df_dgamma_y;
//-------------------------------------------------
                    J[2][0] += common_prefix * cur_cy * df_da;
                    J[2][1] += common_prefix * cur_cy * df_dGamma;
                    J[2][2] += common_prefix * cur_cy * df_dgamma_x;
                    J[2][3] += common_prefix * cur_cy * df_dgamma_y;
//-------------------------------------------------
                    J[3][0] += common_prefix * (c_sq / 2.0) * df_da;
                    J[3][1] += common_prefix * (c_sq / 2.0) * df_dGamma;
                    J[3][2] += common_prefix * (c_sq / 2.0) * df_dgamma_x;
                    J[3][3] += common_prefix * (c_sq / 2.0) * df_dgamma_y;
                    //printf("jacobian OK!\n");
                }
            }
        }
         //printf("jacobian OK!\n");
// --- 線形連立方程式 J * dx = F を解く ---
//4x4の行列式を余因子展開で計算
//3x3の行列式を計算
        double detJ1 = det3x3(
            J[1][1], J[1][2], J[1][3],
            J[2][1], J[2][2], J[2][3],
            J[3][1], J[3][2], J[3][3]
        );
        double detJ2 = det3x3(
            J[0][1], J[0][2], J[0][3],
            J[2][1], J[2][2], J[2][3],
            J[3][1], J[3][2], J[3][3]
        );
        double detJ3 = det3x3(
            J[0][1], J[0][2], J[0][3],
            J[1][1], J[1][2], J[1][3],
            J[3][1], J[3][2], J[3][3]
        );
        double detJ4 = det3x3(
            J[0][1], J[0][2], J[0][3],
            J[1][1], J[1][2], J[1][3],
            J[2][1], J[2][2], J[2][3]
        );
        double detJ = 
            J[0][0] * detJ1
          - J[1][0] * detJ2
          + J[2][0] * detJ3
          - J[3][0] * detJ4;
          //printf("detJ Ok!\n");
        if (fabs(detJ) < 1e-20) {
            //行列式が0に近い場合、更新せずループを抜ける
            return;
        }
       // printf("detJ Ok!\n");

        double inv_detJ = 1.0 / detJ;
        double invJ[4][4],detM[4][4];
        detM[0][0] = det3x3(
            J[1][1], J[1][2], J[1][3],
            J[2][1], J[2][2], J[2][3],
            J[3][1], J[3][2], J[3][3]
        );
        detM[0][1] = -det3x3(
            J[1][0], J[1][2], J[1][3],
            J[2][0], J[2][2], J[2][3],
            J[3][0], J[3][2], J[3][3]
        );
        detM[0][2] = det3x3(
            J[1][0], J[1][1], J[1][3],
            J[2][0], J[2][1], J[2][3],
            J[3][0], J[3][1], J[3][3]
        );
        detM[0][3] = -det3x3(
            J[1][0], J[1][1], J[1][2],
            J[2][0], J[2][1], J[2][2],
            J[3][0], J[3][1], J[3][2]
        );
//-------------------------------------------------
        detM[1][0] = -det3x3(
            J[0][1], J[0][2], J[0][3],
            J[2][1], J[2][2], J[2][3],
            J[3][1], J[3][2], J[3][3]
        );
        detM[1][1] = det3x3(
            J[0][0], J[0][2], J[0][3],
            J[2][0], J[2][2], J[2][3],
            J[3][0], J[3][2], J[3][3]
        );
        detM[1][2] = -det3x3(
            J[0][0], J[0][1], J[0][3],
            J[2][0], J[2][1], J[2][3],
            J[3][0], J[3][1], J[3][3]
        );
        detM[1][3] = det3x3(
            J[0][0], J[0][1], J[0][2],
            J[2][0], J[2][1], J[2][2],
            J[3][0], J[3][1], J[3][2]
        );
//-------------------------------------------------
        detM[2][0] = det3x3(
            J[0][1], J[0][2], J[0][3],
            J[1][1], J[1][2], J[1][3],
            J[3][1], J[3][2], J[3][3]
        );
        detM[2][1] = -det3x3(
            J[0][0], J[0][2], J[0][3],
            J[1][0], J[1][2], J[1][3],
            J[3][0], J[3][2], J[3][3]
        );
        detM[2][2] = det3x3(
            J[0][0], J[0][1], J[0][3],
            J[1][0], J[1][1], J[1][3],
            J[3][0], J[3][1], J[3][3]
        );
        detM[2][3] = -det3x3(
            J[0][0], J[0][1], J[0][2],
            J[1][0], J[1][1], J[1][2],
            J[3][0], J[3][1], J[3][2]
        );
//-------------------------------------------------
        detM[3][0] = -det3x3(
            J[0][1], J[0][2], J[0][3],
            J[1][1], J[1][2], J[1][3],
            J[2][1], J[2][2], J[2][3]
        );
        detM[3][1] = det3x3(
            J[0][0], J[0][2], J[0][3],
            J[1][0], J[1][2], J[1][3],
            J[2][0], J[2][2], J[2][3]
        );
        detM[3][2] = -det3x3(
            J[0][0], J[0][1], J[0][3],
            J[1][0], J[1][1], J[1][3],
            J[2][0], J[2][1], J[2][3]
        );
        detM[3][3] = det3x3(
            J[0][0], J[0][1], J[0][2],
            J[1][0], J[1][1], J[1][2],
            J[2][0], J[2][1], J[2][2]
        );
//-------------------------------------------------
// 逆行列 = 余因子行列の転置 / 行列式
        invJ[0][0] = detM[0][0] * inv_detJ;
        invJ[1][0] = detM[0][1] * inv_detJ;
        invJ[2][0] = detM[0][2] * inv_detJ;
        invJ[3][0] = detM[0][3] * inv_detJ;
        invJ[0][1] = detM[1][0] * inv_detJ;
        invJ[1][1] = detM[1][1] * inv_detJ;
        invJ[2][1] = detM[1][2] * inv_detJ;
        invJ[3][1] = detM[1][3] * inv_detJ;
        invJ[0][2] = detM[2][0] * inv_detJ;
        invJ[1][2] = detM[2][1] * inv_detJ;
        invJ[2][2] = detM[2][2] * inv_detJ;
        invJ[3][2] = detM[2][3] * inv_detJ;
        invJ[0][3] = detM[3][0] * inv_detJ;
        invJ[1][3] = detM[3][1] * inv_detJ;
        invJ[2][3] = detM[3][2] * inv_detJ;
        invJ[3][3] = detM[3][3] * inv_detJ;
        
            //printf("gyaku OK!\n");

        double dx[4]; // 更新量 (δa, δΓ, δγ_x,δγ_y)
        dx[0] = invJ[0][0]*F[0] + invJ[0][1]*F[1] + invJ[0][2]*F[2] + invJ[0][3]*F[3];
        dx[1] = invJ[1][0]*F[0] + invJ[1][1]*F[1] + invJ[1][2]*F[2] + invJ[1][3]*F[3];
        dx[2] = invJ[2][0]*F[0] + invJ[2][1]*F[1] + invJ[2][2]*F[2] + invJ[2][3]*F[3];
        dx[3] = invJ[3][0]*F[0] + invJ[3][1]*F[1] + invJ[3][2]*F[2] + invJ[3][3]*F[3];

        // --- 係数を更新 ---
        coeffs->a     += dx[0];
        coeffs->Gamma += dx[1];
        coeffs->gamma_x += dx[2];
        coeffs->gamma_y += dx[3];

        // --- 収束判定 ---
        //printf("\n");
        //printf("calculating Neewton method.....\n")
        //printf("\n");
        double norm_dx = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + dx[3]*dx[3]);
        if (norm_dx < TOLERANCE_NEWTON) {
            break; // 収束したら関数を抜ける
        }
    }
    //printf("\n");
    //printf("convergence tolerance newton = %d\n",iter);
    //printf("\n");
}
//****************************************************************************
//3*3 matrix determinant
//****************************************************************************
// 3x3行列の行列式を計算する関数（各要素を直接引数で入力）
double det3x3(
    double a11, double a12, double a13,
    double a21, double a22, double a23,
    double a31, double a32, double a33
) {
    return
        a11 * (a22 * a33 - a23 * a32)
      - a12 * (a21 * a33 - a23 * a31)
      + a13 * (a21 * a32 - a22 * a31);
}
//****************************************************************************
//average micro velocity, collision frequency, mean free path,
//viscosity
//****************************************************************************
void vol_re_coll(void){
    int i,i2;
//initialize average micro velocity, viscosity, collision frequency, mean free path,
  #pragma omp parallel for collapse(2)
  for(i=0;i<nx;i++){
    for(i2=0;i2<ny;i2++){
        micVel_ave[i][i2]=0.0;
        myu[i][i2]=0.0;
        freq_coll[i][i2]=0.0;
        mfp[i][i2]=0.0;
    }
 }
//average micro velocity, viscosity, collision frequency, mean free path,
    for(i=0;i<nx;i++){
        for(i2=0;i2<ny;i2++){
            micVel_ave[i][i2]=pow(8.0*R*Tref*T[i][i2]/pi,0.5);
            myu[i][i2]=myuref*pow(T[i][i2],omega);
            mfp[i][i2]=pow(pi*R*Tref*T[i][i2]/2.0,0.5)*myu[i][i2]/(pref*p[i][i2]);
            freq_coll[i][i2]=micVel_ave[i][i2]/mfp[i][i2];
        }
    }
}
//*****************************************************************************************
//high order scheme (x-direction)
//*****************************************************************************************
void fai_x(void){
    int i,i2,j,k,l;
    //double df_minus, df_mid, df_plus;
    #pragma omp parallel for private(i2,j, k, l)
    for(i=1;i<nx-2;i++){
        for(i2=0;i2<ny;i2++){
            for(j=0;j<ncj;j++){
                for(k=0;k<nck;k++){
                    for(l=0;l<ncl;l++){
                        double df_minus=f[i][i2][j][k][l]-f[i-1][i2][j][k][l];
                        double df_mid=f[i+1][i2][j][k][l]-f[i][i2][j][k][l];
                        double df_plus=f[i+2][i2][j][k][l]-f[i+1][i2][j][k][l];
                        phi_x[i][i2][j][k][l]=minmod(df_minus,df_mid,df_plus);
                    }
                }
            }
        }
    }
}
//*****************************************************************************************
//high order scheme (y-direction)
//*****************************************************************************************
void fai_y(void){
    int i,i2,j,k,l;
    //double df_minus, df_mid, df_plus;
    #pragma omp parallel for private(i2,j, k, l)
    for(i=0;i<nx;i++){
        for(i2=1;i2<ny-2;i2++){
            for(j=0;j<ncj;j++){
                for(k=0;k<nck;k++){
                    for(l=0;l<ncl;l++){
                        double df_minus=f[i][i2][j][k][l]-f[i][i2-1][j][k][l];
                        double df_mid=f[i][i2+1][j][k][l]-f[i][i2][j][k][l];
                        double df_plus=f[i][i2+2][j][k][l]-f[i][i2+1][j][k][l];
                        phi_y[i][i2][j][k][l]=minmod(df_minus,df_mid,df_plus);
                    }
                }
            }
        }
    }
}
//*****************************************************************************************
//max
//*****************************************************************************************
static inline double max3(double a, double b, double c) {
    double m = a;
    if (b > m) m = b;
    if (c > m) m = c;
    return m;
}
//*****************************************************************************************
//min
//*****************************************************************************************
static inline double min3(double a, double b, double c) {
    double m = a;
    if (b < m) m = b;
    if (c < m) m = c;
    return m;
}
//*****************************************************************************************
//minmod
//*****************************************************************************************
static inline double minmod(double a, double b, double c) {
    if (a > 0 && b > 0 && c > 0) {
        return min3(a, b, c);
    } else if (a < 0 && b < 0 && c < 0) {
        return max3(a, b, c);
    } else {
        return 0.0;
    }
}
//*****************************************************************************************
//output CSV (macro volumes etc.(time series))
//*****************************************************************************************
void output_csv(int step){
    int i,i2;
    char fname[256];
    sprintf(fname, "2d3v_maxwellRefrection_flatPlate_%07d_%f.csv", step,t); 
    FILE *ofcsv=fopen(fname,"w");
    if(!ofcsv){
        perror("You cannot open this file.");
        exit(1);
    }
    fprintf(ofcsv,"i,j,x,y,n,v1,v2,p,T,mean free path\n");
    for(i=0;i<nx;i++){
        for(i2=0;i2<ny;i2++){
        fprintf(ofcsv,"%d,%d,%f,%f,%f,%f,%f,%f,%f,%f\n", i,i2,x[i],y[i2],n[i][i2],v1[i][i2],v2[i][i2],p[i][i2],T[i][i2],mfp[i][i2]);
        }
    }
    fclose(ofcsv);
}
//*****************************************************************************************
//output CSV (velocity distribution(Time series))
//*****************************************************************************************
void output_csv_f(int step){
    int i,i2,j,k,l;
    char dirname[256], fname[1024];
    sprintf(dirname, "f_output_%07d_%f", step,t);
    #ifdef _WIN32
    _mkdir(dirname);
    #else
    mkdir(dirname, 0777);
    #endif
    for(i=0;i<nx;i++){
        for(i2=0;i2<ny;i2++){
            snprintf(fname, sizeof(fname),"%s/f_%d,%d.csv", dirname, i,i2);
            FILE *fp=fopen(fname,"w");
            if(!fp){
                perror("You cannot open this file.");
                exit(1);
            }
            fprintf(fp,"j,cx,cy,cz,f\n");
            for(j=0;j<ncj;j++){
                for(k=0;k<nck;k++){
                    for(l=0;l<ncl;l++){
                        fprintf(fp,"%d,%d,%d,%f,%f,%f,%e\n", j,k,l,c1[j],c2[k],c3[l],f[i][i2][j][k][l]);
                    }
                }
            }
            fclose(fp);
        }
    }
}
//*****************************************************************************************
//Dynamic Array Allocation
//*****************************************************************************************
//5次元配列の動的確保関数
double *****alloc5D(int nx1, int nc1, int nc2, int nc3, int nc4) {
    double *****arr;
    arr = (double *****)malloc(nx1 * sizeof(double ****));
    for (int i = 0; i < nx1; i++) {
        arr[i] = (double ****)malloc(nc1 * sizeof(double ***));
        for (int j = 0; j < nc1; j++) {
            arr[i][j] = (double ***)malloc(nc2 * sizeof(double **));
            for (int k = 0; k < nc2; k++) {
                arr[i][j][k] = (double **)malloc(nc3 * sizeof(double *));
                for (int l = 0; l < nc3; l++) {
                    arr[i][j][k][l] = (double *)malloc(nc4 * sizeof(double));
                    for (int m = 0; m < nc4; m++) arr[i][j][k][l][m] = 0.0;
                }
            }
        }
    }
    return arr;
}

// 5次元配列の解放関数
void free5D(double *****arr, int nx1, int nc1, int nc2, int nc3) {
    for (int i = 0; i < nx1; i++) {
        for (int j = 0; j < nc1; j++) {
            for (int k = 0; k < nc2; k++) {
                for (int l = 0; l < nc3; l++) {
                    free(arr[i][j][k][l]);
                }
                free(arr[i][j][k]);
            }
            free(arr[i][j]);
        }
        free(arr[i]);
    }
    free(arr);
}

// 4次元配列動的確保関数
double ****alloc4D(int nx1, int nc1,int nc2,int nc3){
    double ****arr;
    arr = (double ****)malloc(nx1 * sizeof(double ***));
    for(int i=0;i<nx1;i++){
        arr[i] = (double ***)malloc(nc1 * sizeof(double **));
        for(int j=0;j<nc1;j++){
            arr[i][j] = (double **)malloc(nc2 * sizeof(double *));
            for(int k=0;k<nc2;k++){
                arr[i][j][k] = (double *)malloc(nc3 * sizeof(double));
                for(int l=0;l<nc3;l++) arr[i][j][k][l] = 0.0;
            }
        }
    }
    return arr;
}

// 4次元配列解放
void free4D(double ****arr, int nx1, int nc1,int nc2){
    for(int i=0;i<nx1;i++){
        for(int j=0;j<nc1;j++){
            for(int k=0;k<nc2;k++){
                free(arr[i][j][k]);
            }
            free(arr[i][j]);
        }
        free(arr[i]);
    }
    free(arr);
}

// 3次元配列動的確保関数
double ***alloc3D(int nx1, int nc1, int nc2){
    double ***arr;
    arr = (double ***)malloc(nx1 * sizeof(double **));
    for(int i=0;i<nx1;i++){
        arr[i] = (double **)malloc(nc1 * sizeof(double *));
        for(int j=0;j<nc1;j++){
            arr[i][j] = (double *)malloc(nc2 * sizeof(double));
            for(int k=0;k<nc2;k++) arr[i][j][k] = 0.0;
        }
    }
    return arr;
}

// 3次元配列解放
void free3D(double ***arr, int nx1, int nc1){
    for(int i=0;i<nx1;i++){
        for(int j=0;j<nc1;j++) free(arr[i][j]);
        free(arr[i]);
    }
    free(arr);
}

// 2次元配列動的確保関数
double **alloc2D(int nx1, int nc1){
    double **arr;
    arr = (double **)malloc(nx1 * sizeof(double *));
    for(int i=0;i<nx1;i++){
        arr[i] = (double *)malloc(nc1 * sizeof(double));
        for(int j=0;j<nc1;j++) arr[i][j] = 0.0;
    }
    return arr;
}

// 2次元配列解放
void free2D(double **arr, int nx1){
    for(int i=0;i<nx1;i++) free(arr[i]);
    free(arr);
}
