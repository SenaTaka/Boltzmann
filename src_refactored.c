//*******************************************************************************************
//
// Non-dimensional Boltzmann Equation BGK Model 2D3V Numerical Simulation
// Reference: Mieussens, L. and Struchtrup, H., Physics of Fluids, Vol. 16, No. 8, 2004
//
// Compile: gcc -fopenmp -O3 -march=native -ffast-math -funroll-loops src_refactored.c -o sim -lm
//*******************************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>
#include <stdbool.h>

//*******************************************************************************************
// Physical Constants
//*******************************************************************************************
typedef struct {
    double boltzmann;       // Boltzmann constant
    double mass;            // Molecular mass
    double radius;          // Molecular radius
    double pi;              // Pi
    double gas_constant;    // R = BC/M
    double cross_section;   // sigma = 4*pi*r*r
    double T_ref;           // Reference temperature
    double mu_ref;          // Reference viscosity
    double omega;           // Viscosity index
    double n_ref;           // Reference density
    double p_ref;           // Reference pressure
    double L_ref;           // Reference length
    double l0;              // Mean free path reference
    double t_ref;           // Time reference
    double alpha;           // Accommodation coefficient
} PhysicalConstants;

//*******************************************************************************************
// Grid Configuration
//*******************************************************************************************
typedef struct {
    // Spatial grid
    int nx, ny;
    double xmin, ymin;
    double dx, dy;
    
    // Velocity grid
    int ncj, nck, ncl;
    double cmin;
    double dcj, dck, dcl;
    
    // Time
    int nstep;
    double tstop;
} GridConfig;

//*******************************************************************************************
// Initial Conditions
//*******************************************************************************************
typedef struct {
    double n_upstream;
    double v1_upstream;
    double v2_upstream;
    double T_upstream;
    double p_upstream;
    double T_wall;
} InitialConditions;

//*******************************************************************************************
// Newton Method Configuration
//*******************************************************************************************
typedef struct {
    int max_iter;
    double tolerance;
} NewtonConfig;

//*******************************************************************************************
// Simulation State
//*******************************************************************************************
typedef struct {
    // Time
    double t;
    double dt;
    
    // Spatial coordinates
    double *x;
    double *y;
    
    // Velocity coordinates
    double *c1;
    double *c2;
    double *c3;
    
    // Macroscopic quantities
    double **n;         // Number density
    double **v1;        // Velocity x-direction
    double **v2;        // Velocity y-direction
    double **p;         // Pressure
    double **T;         // Temperature
    
    // Distribution functions
    double *****f;      // Distribution function
    double *****feq;    // Equilibrium distribution
    
    // Transport properties
    double **mic_vel_ave;   // Average microscopic velocity
    double **freq_coll;     // Collision frequency
    double **mfp;           // Mean free path
    double **mu;            // Viscosity
    
    // Boundary arrays
    double ****nb1;
    double ***nb2;
    
    // Flux arrays
    double *****F;
    double *****G;
    
    // Time stepping arrays
    double *****fb;
    
    // Integration arrays
    double ****sum1_n, ***sum2_n;
    double ****sum1_v1, ***sum2_v1;
    double ****sum1_v2, ***sum2_v2;
    double ****sum1_p, ***sum2_p;
} SimulationState;

//*******************************************************************************************
// Coefficients for Local Distribution
//*******************************************************************************************
typedef struct {
    double a;
    double Gamma;
    double gamma_x;
    double gamma_y;
} Coefficients;

//*******************************************************************************************
// Function Declarations
//*******************************************************************************************

// Initialization
PhysicalConstants init_physical_constants(void);
GridConfig init_grid_config(void);
InitialConditions init_initial_conditions(void);
NewtonConfig init_newton_config(void);
SimulationState* create_simulation_state(const GridConfig *config);
void destroy_simulation_state(SimulationState *state, const GridConfig *config);
void setup_initial_condition(SimulationState *state, const GridConfig *config, 
                            const InitialConditions *ic, const PhysicalConstants *pc);

// Memory management
double **alloc2D(int n1, int n2);
double ***alloc3D(int n1, int n2, int n3);
double ****alloc4D(int n1, int n2, int n3, int n4);
double *****alloc5D(int n1, int n2, int n3, int n4, int n5);
void free2D(double **arr, int n1);
void free3D(double ***arr, int n1, int n2);
void free4D(double ****arr, int n1, int n2, int n3);
void free5D(double *****arr, int n1, int n2, int n3, int n4);

// Boundary conditions
void apply_boundary_conditions(SimulationState *state, const GridConfig *config,
                              const InitialConditions *ic, const PhysicalConstants *pc);

// Distribution functions
void compute_local_equilibrium_newton(SimulationState *state, const GridConfig *config,
                                     const NewtonConfig *nc, const PhysicalConstants *pc);
void solve_newton_coefficients(Coefficients *coeffs, const SimulationState *state,
                              const GridConfig *config, const NewtonConfig *nc,
                              const PhysicalConstants *pc, int i, int i2);
double calculate_f_gamma(int j, int k, int l, const Coefficients *coeffs,
                        const SimulationState *state, double vi, double vi2);

// Flux computation
void compute_flux_x(SimulationState *state, const GridConfig *config);
void compute_flux_y(SimulationState *state, const GridConfig *config);

// Time integration
void determine_timestep(SimulationState *state, const GridConfig *config,
                       const PhysicalConstants *pc);
void runge_kutta_step1(SimulationState *state, const GridConfig *config,
                      const PhysicalConstants *pc);
void runge_kutta_step2(SimulationState *state, const GridConfig *config,
                      const PhysicalConstants *pc);

// Integration
void integrate_macroscopic_quantities(SimulationState *state, const GridConfig *config);
void update_transport_properties(SimulationState *state, const GridConfig *config,
                                const PhysicalConstants *pc);

// Utility functions
double det3x3(double a11, double a12, double a13,
             double a21, double a22, double a23,
             double a31, double a32, double a33);
double minmod(double a, double b, double c);
double max3(double a, double b, double c);
double min3(double a, double b, double c);

// Output
void output_macroscopic_csv(const SimulationState *state, const GridConfig *config,
                           int step);
void output_distribution_csv(const SimulationState *state, const GridConfig *config,
                            int step);

//*******************************************************************************************
// Main Program
//*******************************************************************************************
int main(void) {
    // Initialize configurations
    PhysicalConstants pc = init_physical_constants();
    GridConfig config = init_grid_config();
    InitialConditions ic = init_initial_conditions();
    NewtonConfig nc = init_newton_config();
    
    // Create simulation state
    SimulationState *state = create_simulation_state(&config);
    if (!state) {
        fprintf(stderr, "Error: Failed to allocate simulation state\n");
        return EXIT_FAILURE;
    }
    
    printf("Memory allocation successful!\n");
    
    // Setup initial conditions
    setup_initial_condition(state, &config, &ic, &pc);
    printf("Initial conditions set!\n");
    
    // Set OpenMP threads
    omp_set_num_threads(20);
    
    // Time integration loop
    state->t = 0.0;
    for (int step = 1; step <= config.nstep; step++) {
        printf("=== Step: %d ===\n", step);
        
        // First RK2 stage
        apply_boundary_conditions(state, &config, &ic, &pc);
        compute_local_equilibrium_newton(state, &config, &nc, &pc);
        compute_flux_x(state, &config);
        compute_flux_y(state, &config);
        determine_timestep(state, &config, &pc);
        runge_kutta_step1(state, &config, &pc);
        integrate_macroscopic_quantities(state, &config);
        update_transport_properties(state, &config, &pc);
        
        // Second RK2 stage
        apply_boundary_conditions(state, &config, &ic, &pc);
        compute_local_equilibrium_newton(state, &config, &nc, &pc);
        compute_flux_x(state, &config);
        compute_flux_y(state, &config);
        determine_timestep(state, &config, &pc);
        runge_kutta_step2(state, &config, &pc);
        integrate_macroscopic_quantities(state, &config);
        update_transport_properties(state, &config, &pc);
        
        // Output results
        if (step % 1000 == 0) {
            output_macroscopic_csv(state, &config, step);
        }
        if (step % 5000 == 0) {
            output_distribution_csv(state, &config, step);
        }
        
        state->t += state->dt;
    }
    
    // Final output
    output_macroscopic_csv(state, &config, config.nstep);
    
    // Cleanup
    destroy_simulation_state(state, &config);
    
    return EXIT_SUCCESS;
}

//*******************************************************************************************
// Initialization Functions
//*******************************************************************************************

PhysicalConstants init_physical_constants(void) {
    PhysicalConstants pc;
    pc.boltzmann = 1.380658e-23;
    pc.mass = 6.6339e-26;
    pc.radius = 1.85e-10;
    pc.pi = 3.1415926535897932;
    pc.gas_constant = pc.boltzmann / pc.mass;
    pc.cross_section = 4.0 * pc.pi * pc.radius * pc.radius;
    pc.T_ref = 273.0;
    pc.mu_ref = 2.117e-5;
    pc.omega = 0.81;
    pc.n_ref = 1.894e-6;
    pc.p_ref = 1.076e-1;
    pc.L_ref = 0.9;
    pc.l0 = sqrt(pc.pi * pc.gas_constant * pc.T_ref / 2.0) * pc.mu_ref / pc.p_ref;
    pc.t_ref = sqrt(2.0 * pc.gas_constant * pc.T_ref) / pc.L_ref;
    pc.alpha = 1.0;
    return pc;
}

GridConfig init_grid_config(void) {
    GridConfig config;
    config.nx = 55;
    config.ny = 45;
    config.xmin = -4.0e-2;
    config.ymin = -5.0e-2;
    config.dx = 2.0e-2;
    config.dy = 1.25e-2;
    config.ncj = 100;
    config.nck = 100;
    config.ncl = 60;
    config.cmin = -5.0;
    config.dcj = 0.125;
    config.dck = 0.125;
    config.dcl = 0.2;
    config.nstep = 30;
    config.tstop = 10.0;
    return config;
}

InitialConditions init_initial_conditions(void) {
    InitialConditions ic;
    ic.n_upstream = 1.000;
    ic.v1_upstream = 4.564;
    ic.v2_upstream = 0.000;
    ic.T_upstream = 1.000;
    ic.p_upstream = 1.000;
    ic.T_wall = 1.000;
    return ic;
}

NewtonConfig init_newton_config(void) {
    NewtonConfig nc;
    nc.max_iter = 100;
    nc.tolerance = 1e-8;
    return nc;
}

//*******************************************************************************************
// Memory Management Functions
//*******************************************************************************************

double **alloc2D(int n1, int n2) {
    double **arr = (double **)malloc(n1 * sizeof(double *));
    if (!arr) return NULL;
    
    for (int i = 0; i < n1; i++) {
        arr[i] = (double *)calloc(n2, sizeof(double));
        if (!arr[i]) {
            for (int j = 0; j < i; j++) free(arr[j]);
            free(arr);
            return NULL;
        }
    }
    return arr;
}

void free2D(double **arr, int n1) {
    if (!arr) return;
    for (int i = 0; i < n1; i++) {
        free(arr[i]);
    }
    free(arr);
}

double ***alloc3D(int n1, int n2, int n3) {
    double ***arr = (double ***)malloc(n1 * sizeof(double **));
    if (!arr) return NULL;
    
    for (int i = 0; i < n1; i++) {
        arr[i] = (double **)malloc(n2 * sizeof(double *));
        if (!arr[i]) {
            for (int j = 0; j < i; j++) {
                for (int k = 0; k < n2; k++) free(arr[j][k]);
                free(arr[j]);
            }
            free(arr);
            return NULL;
        }
        for (int j = 0; j < n2; j++) {
            arr[i][j] = (double *)calloc(n3, sizeof(double));
            if (!arr[i][j]) {
                for (int k = 0; k < j; k++) free(arr[i][k]);
                free(arr[i]);
                for (int k = 0; k < i; k++) {
                    for (int l = 0; l < n2; l++) free(arr[k][l]);
                    free(arr[k]);
                }
                free(arr);
                return NULL;
            }
        }
    }
    return arr;
}

void free3D(double ***arr, int n1, int n2) {
    if (!arr) return;
    for (int i = 0; i < n1; i++) {
        if (arr[i]) {
            for (int j = 0; j < n2; j++) {
                free(arr[i][j]);
            }
            free(arr[i]);
        }
    }
    free(arr);
}

double ****alloc4D(int n1, int n2, int n3, int n4) {
    double ****arr = (double ****)malloc(n1 * sizeof(double ***));
    if (!arr) return NULL;
    
    for (int i = 0; i < n1; i++) {
        arr[i] = alloc3D(n2, n3, n4);
        if (!arr[i]) {
            for (int j = 0; j < i; j++) {
                free3D(arr[j], n2, n3);
            }
            free(arr);
            return NULL;
        }
    }
    return arr;
}

void free4D(double ****arr, int n1, int n2, int n3) {
    if (!arr) return;
    for (int i = 0; i < n1; i++) {
        free3D(arr[i], n2, n3);
    }
    free(arr);
}

double *****alloc5D(int n1, int n2, int n3, int n4, int n5) {
    double *****arr = (double *****)malloc(n1 * sizeof(double ****));
    if (!arr) return NULL;
    
    for (int i = 0; i < n1; i++) {
        arr[i] = alloc4D(n2, n3, n4, n5);
        if (!arr[i]) {
            for (int j = 0; j < i; j++) {
                free4D(arr[j], n2, n3, n4);
            }
            free(arr);
            return NULL;
        }
    }
    return arr;
}

void free5D(double *****arr, int n1, int n2, int n3, int n4) {
    if (!arr) return;
    for (int i = 0; i < n1; i++) {
        free4D(arr[i], n2, n3, n4);
    }
    free(arr);
}

SimulationState* create_simulation_state(const GridConfig *config) {
    SimulationState *state = (SimulationState *)malloc(sizeof(SimulationState));
    if (!state) return NULL;
    
    // Allocate coordinate arrays
    state->x = (double *)calloc(config->nx, sizeof(double));
    state->y = (double *)calloc(config->ny, sizeof(double));
    state->c1 = (double *)calloc(config->ncj, sizeof(double));
    state->c2 = (double *)calloc(config->nck, sizeof(double));
    state->c3 = (double *)calloc(config->ncl, sizeof(double));
    
    // Allocate macroscopic arrays
    state->n = alloc2D(config->nx, config->ny);
    state->v1 = alloc2D(config->nx, config->ny);
    state->v2 = alloc2D(config->nx, config->ny);
    state->p = alloc2D(config->nx, config->ny);
    state->T = alloc2D(config->nx, config->ny);
    
    // Allocate distribution functions
    state->f = alloc5D(config->nx, config->ny, config->ncj, config->nck, config->ncl);
    state->feq = alloc5D(config->nx, config->ny, config->ncj, config->nck, config->ncl);
    
    // Allocate transport properties
    state->mic_vel_ave = alloc2D(config->nx, config->ny);
    state->freq_coll = alloc2D(config->nx, config->ny);
    state->mfp = alloc2D(config->nx, config->ny);
    state->mu = alloc2D(config->nx, config->ny);
    
    // Allocate boundary arrays
    state->nb1 = alloc4D(config->nx, config->ny, config->ncj, config->ncl);
    state->nb2 = alloc3D(config->nx, config->ny, config->ncl);
    
    // Allocate flux arrays
    state->F = alloc5D(config->nx, config->ny, config->ncj, config->nck, config->ncl);
    state->G = alloc5D(config->nx, config->ny, config->ncj, config->nck, config->ncl);
    
    // Allocate time stepping arrays
    state->fb = alloc5D(config->nx, config->ny, config->ncj, config->nck, config->ncl);
    
    // Allocate integration arrays
    state->sum1_n = alloc4D(config->nx, config->ny, config->nck, config->ncl);
    state->sum2_n = alloc3D(config->nx, config->ny, config->ncl);
    state->sum1_v1 = alloc4D(config->nx, config->ny, config->nck, config->ncl);
    state->sum2_v1 = alloc3D(config->nx, config->ny, config->ncl);
    state->sum1_v2 = alloc4D(config->nx, config->ny, config->ncj, config->ncl);
    state->sum2_v2 = alloc3D(config->nx, config->ny, config->ncl);
    state->sum1_p = alloc4D(config->nx, config->ny, config->nck, config->ncl);
    state->sum2_p = alloc3D(config->nx, config->ny, config->ncl);
    
    // Check if any allocation failed
    if (!state->x || !state->y || !state->c1 || !state->c2 || !state->c3 ||
        !state->n || !state->v1 || !state->v2 || !state->p || !state->T ||
        !state->f || !state->feq || !state->mic_vel_ave || !state->freq_coll ||
        !state->mfp || !state->mu || !state->nb1 || !state->nb2 ||
        !state->F || !state->G ||
        !state->fb || !state->sum1_n || !state->sum2_n ||
        !state->sum1_v1 || !state->sum2_v1 || !state->sum1_v2 || !state->sum2_v2 ||
        !state->sum1_p || !state->sum2_p) {
        destroy_simulation_state(state, config);
        return NULL;
    }
    
    return state;
}

void destroy_simulation_state(SimulationState *state, const GridConfig *config) {
    if (!state) return;
    
    free(state->x);
    free(state->y);
    free(state->c1);
    free(state->c2);
    free(state->c3);
    
    free2D(state->n, config->nx);
    free2D(state->v1, config->nx);
    free2D(state->v2, config->nx);
    free2D(state->p, config->nx);
    free2D(state->T, config->nx);
    
    free5D(state->f, config->nx, config->ny, config->ncj, config->nck);
    free5D(state->feq, config->nx, config->ny, config->ncj, config->nck);
    
    free2D(state->mic_vel_ave, config->nx);
    free2D(state->freq_coll, config->nx);
    free2D(state->mfp, config->nx);
    free2D(state->mu, config->nx);
    
    free4D(state->nb1, config->nx, config->ny, config->ncj);
    free3D(state->nb2, config->nx, config->ny);
    
    free5D(state->F, config->nx, config->ny, config->ncj, config->nck);
    free5D(state->G, config->nx, config->ny, config->ncj, config->nck);
    
    free5D(state->fb, config->nx, config->ny, config->ncj, config->nck);
    
    free4D(state->sum1_n, config->nx, config->ny, config->nck);
    free3D(state->sum2_n, config->nx, config->ny);
    free4D(state->sum1_v1, config->nx, config->ny, config->nck);
    free3D(state->sum2_v1, config->nx, config->ny);
    free4D(state->sum1_v2, config->nx, config->ny, config->ncj);
    free3D(state->sum2_v2, config->nx, config->ny);
    free4D(state->sum1_p, config->nx, config->ny, config->nck);
    free3D(state->sum2_p, config->nx, config->ny);
    
    free(state);
}

//*******************************************************************************************
// Initial Condition Setup
//*******************************************************************************************

void setup_initial_condition(SimulationState *state, const GridConfig *config,
                            const InitialConditions *ic, const PhysicalConstants *pc) {
    // Setup spatial coordinates
    #pragma omp parallel for
    for (int i = 0; i < config->nx; i++) {
        state->x[i] = config->xmin + i * config->dx;
    }
    
    #pragma omp parallel for
    for (int i2 = 0; i2 < config->ny; i2++) {
        state->y[i2] = config->ymin + i2 * config->dy;
    }
    
    // Setup velocity coordinates
    #pragma omp parallel for
    for (int j = 0; j < config->ncj; j++) {
        state->c1[j] = config->cmin + j * config->dcj;
    }
    
    #pragma omp parallel for
    for (int k = 0; k < config->nck; k++) {
        state->c2[k] = config->cmin + k * config->dck;
    }
    
    #pragma omp parallel for
    for (int l = 0; l < config->ncl; l++) {
        state->c3[l] = config->cmin + l * config->dcl;
    }
    
    // Initialize macroscopic quantities
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < config->nx; i++) {
        for (int i2 = 0; i2 < config->ny; i2++) {
            if (i2 == 0) {
                state->T[i][0] = ic->T_wall;
                state->n[i][0] = 0.0;
            } else {
                state->n[i][i2] = ic->n_upstream;
                state->v1[i][i2] = ic->v1_upstream;
                state->v2[i][i2] = ic->v2_upstream;
                state->p[i][i2] = ic->p_upstream;
                state->T[i][i2] = ic->T_upstream;
            }
        }
    }
    
    // Initialize distribution function (Maxwellian)
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < config->nx; i++) {
        for (int i2 = 0; i2 < config->ny; i2++) {
            double n = state->n[i][i2];
            double v1 = state->v1[i][i2];
            double v2 = state->v2[i][i2];
            double T = state->T[i][i2];
            double coeff = n * pow(pc->pi * T, -1.5);
            
            for (int j = 0; j < config->ncj; j++) {
                for (int k = 0; k < config->nck; k++) {
                    for (int l = 0; l < config->ncl; l++) {
                        double c1 = state->c1[j];
                        double c2 = state->c2[k];
                        double c3 = state->c3[l];
                        double exponent = -((c1 - v1) * (c1 - v1) +
                                          (c2 - v2) * (c2 - v2) +
                                          c3 * c3) / T;
                        state->f[i][i2][j][k][l] = coeff * exp(exponent);
                        state->feq[i][i2][j][k][l] = state->f[i][i2][j][k][l];
                    }
                }
            }
        }
    }
    
    // Initialize transport properties
    update_transport_properties(state, config, pc);
}

//*******************************************************************************************
// Transport Properties Update
//*******************************************************************************************

void update_transport_properties(SimulationState *state, const GridConfig *config,
                                const PhysicalConstants *pc) {
    #pragma omp parallel
    {
        #pragma omp for collapse(2) nowait
        for (int i = 0; i < config->nx; i++) {
            for (int i2 = 0; i2 < config->ny; i2++) {
                double T = state->T[i][i2];
                double p = state->p[i][i2];
                
                state->mic_vel_ave[i][i2] = sqrt(8.0 * pc->gas_constant * pc->T_ref * T / pc->pi);
                state->mu[i][i2] = pc->mu_ref * pow(T, pc->omega);
                state->mfp[i][i2] = sqrt(pc->pi * pc->gas_constant * pc->T_ref * T / 2.0) *
                                   state->mu[i][i2] / (pc->p_ref * p);
                state->freq_coll[i][i2] = state->mic_vel_ave[i][i2] / state->mfp[i][i2];
            }
        }
    }
}

//*******************************************************************************************
// Boundary Conditions (Continued in next part due to length)
//*******************************************************************************************

void apply_boundary_conditions(SimulationState *state, const GridConfig *config,
                              const InitialConditions *ic, const PhysicalConstants *pc) {
    const double pi_sqrt = sqrt(pc->pi / ic->T_wall);
    
    // Wall density calculation
    #pragma omp parallel for
    for (int i = 0; i < config->nx; i++) {
        state->T[i][0] = ic->T_wall;
        state->n[i][0] = 0.0;
        
        for (int j = 0; j < config->ncj; j++) {
            for (int l = 0; l < config->ncl; l++) {
                state->nb1[i][0][j][l] = 0.0;
            }
        }
        
        for (int j = 0; j < config->ncj; j++) {
            for (int k = 0; k < config->nck / 2 - 1; k++) {
                for (int l = 0; l < config->ncl; l++) {
                    state->nb1[i][0][j][l] += 0.5 * (state->c2[k] * state->f[i][2][j][k][l] +
                                                     state->c2[k+1] * state->f[i][2][j][k+1][l]) *
                                             config->dck;
                }
            }
        }
        
        for (int l = 0; l < config->ncl; l++) {
            state->nb2[i][0][l] = 0.0;
        }
        
        for (int j = 0; j < config->ncj - 1; j++) {
            for (int l = 0; l < config->ncl; l++) {
                state->nb2[i][0][l] += 0.5 * (state->nb1[i][0][j][l] +
                                             state->nb1[i][0][j+1][l]) * config->dcj;
            }
        }
        
        for (int l = 0; l < config->ncl - 1; l++) {
            state->n[i][0] += (-2.0) * pi_sqrt * 0.5 * (state->nb2[i][0][l] +
                                                         state->nb2[i][0][l+1]) * config->dcl;
        }
    }
    
    // Left side (inflow)
    #pragma omp parallel for collapse(4)
    for (int i2 = 0; i2 < config->ny; i2++) {
        for (int j = 0; j < config->ncj; j++) {
            for (int k = 0; k < config->nck; k++) {
                for (int l = 0; l < config->ncl; l++) {
                    double c1 = state->c1[j];
                    double c2 = state->c2[k];
                    double c3 = state->c3[l];
                    double exponent = -((c1 - ic->v1_upstream) * (c1 - ic->v1_upstream) +
                                       (c2 - ic->v2_upstream) * (c2 - ic->v2_upstream) +
                                       c3 * c3) / ic->T_upstream;
                    double val = ic->n_upstream * pow(pc->pi * ic->T_upstream, -1.5) * exp(exponent);
                    state->f[0][i2][j][k][l] = val;
                    state->f[1][i2][j][k][l] = val;
                }
            }
        }
    }
    
    // Right side (outflow)
    #pragma omp parallel for collapse(4)
    for (int i2 = 0; i2 < config->ny; i2++) {
        for (int j = 0; j < config->ncj; j++) {
            for (int k = 0; k < config->nck; k++) {
                for (int l = 0; l < config->ncl; l++) {
                    state->f[config->nx-1][i2][j][k][l] = state->f[config->nx-3][i2][j][k][l];
                    state->f[config->nx-2][i2][j][k][l] = state->f[config->nx-3][i2][j][k][l];
                }
            }
        }
    }
    
    // Upper side
    #pragma omp parallel for collapse(4)
    for (int i = 0; i < config->nx; i++) {
        for (int j = 0; j < config->ncj; j++) {
            for (int k = 0; k < config->nck; k++) {
                for (int l = 0; l < config->ncl; l++) {
                    state->f[i][config->ny-1][j][k][l] = state->f[i][config->ny-3][j][k][l];
                    state->f[i][config->ny-2][j][k][l] = state->f[i][config->ny-3][j][k][l];
                }
            }
        }
    }
    
    // Wall (lower side) - Maxwell reflection
    #pragma omp parallel for
    for (int i = 0; i < config->nx; i++) {
        for (int j = 0; j < config->ncj; j++) {
            for (int k = 0; k < config->nck; k++) {
                for (int l = 0; l < config->ncl; l++) {
                    double c1 = state->c1[j];
                    double c2 = state->c2[k];
                    double c3 = state->c3[l];
                    
                    if (state->x[i] <= 0.1) {
                        if (c2 < 0.0) {
                            state->f[i][1][j][k][l] = state->f[i][2][j][k][l];
                        } else {
                            state->f[i][1][j][k][l] = state->f[i][2][j][config->nck-k-1][l];
                        }
                        state->f[i][0][j][k][l] = state->f[i][1][j][k][l];
                    } else {
                        if (c2 < 0.0) {
                            state->f[i][1][j][k][l] = state->f[i][2][j][k][l];
                            state->f[i][0][j][k][l] = state->f[i][1][j][k][l];
                        } else {
                            double f_spec = state->f[i][2][j][config->nck-k-1][l];
                            double T = state->T[i][0];
                            double exponent = -(c1*c1 + c2*c2 + c3*c3) / T;
                            double f_diff = state->n[i][0] * pow(pc->pi * T, -1.5) * exp(exponent);
                            state->f[i][1][j][k][l] = (1.0 - pc->alpha) * f_spec + pc->alpha * f_diff;
                            state->f[i][0][j][k][l] = state->f[i][1][j][k][l];
                        }
                    }
                }
            }
        }
    }
}

//*******************************************************************************************
// Local Equilibrium Distribution (Newton Method)
//*******************************************************************************************

double calculate_f_gamma(int j, int k, int l, const Coefficients *coeffs,
                        const SimulationState *state, double vi, double vi2) {
    double cx = state->c1[j];
    double cy = state->c2[k];
    double cz = state->c3[l];
    
    double C2 = (cx - vi) * (cx - vi) + (cy - vi2) * (cy - vi2) + cz * cz;
    double exponent = -coeffs->Gamma * C2 + coeffs->gamma_x * (cx - vi) +
                     coeffs->gamma_y * (cy - vi2);
    
    return coeffs->a * exp(exponent);
}

void solve_newton_coefficients(Coefficients *coeffs, const SimulationState *state,
                              const GridConfig *config, const NewtonConfig *nc,
                              const PhysicalConstants *pc, int i, int i2) {
    const double dV = config->dcj * config->dck * config->dcl;
    
    for (int iter = 0; iter < nc->max_iter; iter++) {
        double F[4] = {0.0};
        double J[4][4] = {{0.0}};
        
        for (int j = 0; j < config->ncj; j++) {
            for (int k = 0; k < config->nck; k++) {
                for (int l = 0; l < config->ncl; l++) {
                    double cx = state->c1[j];
                    double cy = state->c2[k];
                    double cz = state->c3[l];
                    
                    double f_gamma = calculate_f_gamma(j, k, l, coeffs, state,
                                                      state->v1[i][i2], state->v2[i][i2]);
                    double diff = state->f[i][i2][j][k][l] - f_gamma;
                    double prefix = state->freq_coll[i][i2] * dV;
                    double c_sq = cx*cx + cy*cy + cz*cz;
                    
                    // Residual vector
                    F[0] += prefix * diff;
                    F[1] += prefix * cx * diff;
                    F[2] += prefix * cy * diff;
                    F[3] += prefix * (c_sq / 2.0) * diff;
                    
                    // Jacobian matrix
                    double C2 = (cx - state->v1[i][i2]) * (cx - state->v1[i][i2]) +
                               (cy - state->v2[i][i2]) * (cy - state->v2[i][i2]) + cz * cz;
                    double df_da = f_gamma / coeffs->a;
                    double df_dGamma = -C2 * f_gamma;
                    double df_dgamma_x = (cx - state->v1[i][i2]) * f_gamma;
                    double df_dgamma_y = (cy - state->v2[i][i2]) * f_gamma;
                    
                    J[0][0] += prefix * df_da;
                    J[0][1] += prefix * df_dGamma;
                    J[0][2] += prefix * df_dgamma_x;
                    J[0][3] += prefix * df_dgamma_y;
                    
                    J[1][0] += prefix * cx * df_da;
                    J[1][1] += prefix * cx * df_dGamma;
                    J[1][2] += prefix * cx * df_dgamma_x;
                    J[1][3] += prefix * cx * df_dgamma_y;
                    
                    J[2][0] += prefix * cy * df_da;
                    J[2][1] += prefix * cy * df_dGamma;
                    J[2][2] += prefix * cy * df_dgamma_x;
                    J[2][3] += prefix * cy * df_dgamma_y;
                    
                    J[3][0] += prefix * (c_sq / 2.0) * df_da;
                    J[3][1] += prefix * (c_sq / 2.0) * df_dGamma;
                    J[3][2] += prefix * (c_sq / 2.0) * df_dgamma_x;
                    J[3][3] += prefix * (c_sq / 2.0) * df_dgamma_y;
                }
            }
        }
        
        // Compute determinant
        double detJ1 = det3x3(J[1][1], J[1][2], J[1][3],
                             J[2][1], J[2][2], J[2][3],
                             J[3][1], J[3][2], J[3][3]);
        double detJ2 = det3x3(J[0][1], J[0][2], J[0][3],
                             J[2][1], J[2][2], J[2][3],
                             J[3][1], J[3][2], J[3][3]);
        double detJ3 = det3x3(J[0][1], J[0][2], J[0][3],
                             J[1][1], J[1][2], J[1][3],
                             J[3][1], J[3][2], J[3][3]);
        double detJ4 = det3x3(J[0][1], J[0][2], J[0][3],
                             J[1][1], J[1][2], J[1][3],
                             J[2][1], J[2][2], J[2][3]);
        double detJ = J[0][0] * detJ1 - J[1][0] * detJ2 + J[2][0] * detJ3 - J[3][0] * detJ4;
        
        if (fabs(detJ) < 1e-20) {
            return;
        }
        
        // Compute inverse Jacobian using cofactor method
        double inv_detJ = 1.0 / detJ;
        double invJ[4][4];
        
        // First row of inverse
        invJ[0][0] = detJ1 * inv_detJ;
        invJ[0][1] = -detJ2 * inv_detJ;
        invJ[0][2] = detJ3 * inv_detJ;
        invJ[0][3] = -detJ4 * inv_detJ;
        
        // Second row (compute cofactors)
        invJ[1][0] = -det3x3(J[1][0], J[1][2], J[1][3],
                            J[2][0], J[2][2], J[2][3],
                            J[3][0], J[3][2], J[3][3]) * inv_detJ;
        invJ[1][1] = det3x3(J[0][0], J[0][2], J[0][3],
                           J[2][0], J[2][2], J[2][3],
                           J[3][0], J[3][2], J[3][3]) * inv_detJ;
        invJ[1][2] = -det3x3(J[0][0], J[0][1], J[0][3],
                            J[2][0], J[2][1], J[2][3],
                            J[3][0], J[3][1], J[3][3]) * inv_detJ;
        invJ[1][3] = det3x3(J[0][0], J[0][1], J[0][2],
                           J[2][0], J[2][1], J[2][2],
                           J[3][0], J[3][1], J[3][2]) * inv_detJ;
        
        // Third row
        invJ[2][0] = det3x3(J[1][0], J[1][1], J[1][3],
                           J[2][0], J[2][1], J[2][3],
                           J[3][0], J[3][1], J[3][3]) * inv_detJ;
        invJ[2][1] = -det3x3(J[0][0], J[0][1], J[0][3],
                            J[1][0], J[1][1], J[1][3],
                            J[3][0], J[3][1], J[3][3]) * inv_detJ;
        invJ[2][2] = det3x3(J[0][0], J[0][1], J[0][2],
                           J[1][0], J[1][1], J[1][2],
                           J[3][0], J[3][1], J[3][2]) * inv_detJ;
        invJ[2][3] = -det3x3(J[0][0], J[0][1], J[0][2],
                            J[1][0], J[1][1], J[1][2],
                            J[2][0], J[2][1], J[2][2]) * inv_detJ;
        
        // Fourth row
        invJ[3][0] = -det3x3(J[1][0], J[1][1], J[1][2],
                            J[2][0], J[2][1], J[2][2],
                            J[3][0], J[3][1], J[3][2]) * inv_detJ;
        invJ[3][1] = det3x3(J[0][0], J[0][1], J[0][2],
                           J[1][0], J[1][1], J[1][2],
                           J[2][0], J[2][1], J[2][2]) * inv_detJ;
        invJ[3][2] = -det3x3(J[0][0], J[0][1], J[0][2],
                            J[1][0], J[1][1], J[1][2],
                            J[3][0], J[3][1], J[3][2]) * inv_detJ;
        invJ[3][3] = det3x3(J[0][0], J[0][1], J[0][2],
                           J[1][0], J[1][1], J[1][2],
                           J[2][0], J[2][1], J[2][2]) * inv_detJ;
        
        // Compute update
        double dx[4];
        dx[0] = invJ[0][0]*F[0] + invJ[0][1]*F[1] + invJ[0][2]*F[2] + invJ[0][3]*F[3];
        dx[1] = invJ[1][0]*F[0] + invJ[1][1]*F[1] + invJ[1][2]*F[2] + invJ[1][3]*F[3];
        dx[2] = invJ[2][0]*F[0] + invJ[2][1]*F[1] + invJ[2][2]*F[2] + invJ[2][3]*F[3];
        dx[3] = invJ[3][0]*F[0] + invJ[3][1]*F[1] + invJ[3][2]*F[2] + invJ[3][3]*F[3];
        
        // Update coefficients
        coeffs->a += dx[0];
        coeffs->Gamma += dx[1];
        coeffs->gamma_x += dx[2];
        coeffs->gamma_y += dx[3];
        
        // Check convergence
        double norm_dx = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] + dx[3]*dx[3]);
        if (norm_dx < nc->tolerance) {
            break;
        }
    }
}

void compute_local_equilibrium_newton(SimulationState *state, const GridConfig *config,
                                     const NewtonConfig *nc, const PhysicalConstants *pc) {
    #pragma omp parallel for collapse(2) schedule(static)
    for (int i = 0; i < config->nx; i++) {
        for (int i2 = 0; i2 < config->ny; i2++) {
            Coefficients coeffs;
            coeffs.a = state->n[i][i2] * pow(pc->pi * state->T[i][i2], -1.5);
            coeffs.Gamma = 1.0 / state->T[i][i2];
            coeffs.gamma_x = 0.0;
            coeffs.gamma_y = 0.0;
            
            solve_newton_coefficients(&coeffs, state, config, nc, pc, i, i2);
            
            for (int j = 0; j < config->ncj; j++) {
                for (int k = 0; k < config->nck; k++) {
                    for (int l = 0; l < config->ncl; l++) {
                        state->feq[i][i2][j][k][l] = calculate_f_gamma(j, k, l, &coeffs, state,
                                                                       state->v1[i][i2],
                                                                       state->v2[i][i2]);
                    }
                }
            }
        }
    }
}

//*******************************************************************************************
// Flux Computation
//*******************************************************************************************

void compute_flux_x(SimulationState *state, const GridConfig *config) {
    #pragma omp parallel
    {
        #pragma omp for collapse(2) nowait
        for (int i = 0; i < config->nx - 1; i++) {
            for (int i2 = 0; i2 < config->ny; i2++) {
                for (int j = 0; j < config->ncj; j++) {
                    double c1 = state->c1[j];
                    for (int k = 0; k < config->nck; k++) {
                        for (int l = 0; l < config->ncl; l++) {
                            double f_left = state->f[i][i2][j][k][l];
                            double f_right = state->f[i+1][i2][j][k][l];
                            
                            // Compute phi_x inline (flux limiter)
                            double phi;
                            if (i >= 1 && i <= config->nx - 2) {
                                double df_minus = state->f[i][i2][j][k][l] - state->f[i-1][i2][j][k][l];
                                double df_mid = state->f[i+1][i2][j][k][l] - state->f[i][i2][j][k][l];
                                double df_plus = state->f[i+2][i2][j][k][l] - state->f[i+1][i2][j][k][l];
                                phi = minmod(df_minus, df_mid, df_plus);
                            } else {
                                phi = 0.0;
                            }
                            
                            state->F[i][i2][j][k][l] = 0.5 * c1 * (f_right + f_left) -
                                                      0.5 * fabs(c1) * (f_right - f_left - phi);
                        }
                    }
                }
            }
        }
    }
}

void compute_flux_y(SimulationState *state, const GridConfig *config) {
    #pragma omp parallel
    {
        #pragma omp for collapse(2) nowait
        for (int i = 0; i < config->nx; i++) {
            for (int i2 = 0; i2 < config->ny - 1; i2++) {
                for (int j = 0; j < config->ncj; j++) {
                    for (int k = 0; k < config->nck; k++) {
                        double c2 = state->c2[k];
                        for (int l = 0; l < config->ncl; l++) {
                            double f_lower = state->f[i][i2][j][k][l];
                            double f_upper = state->f[i][i2+1][j][k][l];
                            
                            // Compute phi_y inline (flux limiter)
                            double phi;
                            if (i2 >= 1 && i2 <= config->ny - 2) {
                                double df_minus = state->f[i][i2][j][k][l] - state->f[i][i2-1][j][k][l];
                                double df_mid = state->f[i][i2+1][j][k][l] - state->f[i][i2][j][k][l];
                                double df_plus = state->f[i][i2+2][j][k][l] - state->f[i][i2+1][j][k][l];
                                phi = minmod(df_minus, df_mid, df_plus);
                            } else {
                                phi = 0.0;
                            }
                            
                            state->G[i][i2][j][k][l] = 0.5 * c2 * (f_upper + f_lower) -
                                                      0.5 * fabs(c2) * (f_upper - f_lower - phi);
                        }
                    }
                }
            }
        }
    }
}

//*******************************************************************************************
// Time Integration
//*******************************************************************************************

void determine_timestep(SimulationState *state, const GridConfig *config,
                       const PhysicalConstants *pc) {
    double max_nu = 0.0;
    double max_cx = 0.0;
    double max_cy = 0.0;
    const double cfl = 0.25;
    
    #pragma omp parallel for reduction(max: max_nu) collapse(2)
    for (int i = 0; i < config->nx; i++) {
        for (int i2 = 0; i2 < config->ny; i2++) {
            double nu = state->freq_coll[i][i2] / pc->t_ref;
            if (nu > max_nu) {
                max_nu = nu;
            }
        }
    }
    
    for (int j = 0; j < config->ncj; j++) {
        if (fabs(state->c1[j]) > max_cx) {
            max_cx = fabs(state->c1[j]);
        }
    }
    
    for (int k = 0; k < config->nck; k++) {
        if (fabs(state->c2[k]) > max_cy) {
            max_cy = fabs(state->c2[k]);
        }
    }
    
    double denominator = max_nu + (max_cx / config->dx) + (max_cy / config->dy);
    
    if (denominator > 0.0) {
        state->dt = cfl / denominator;
    } else {
        state->dt = 1.0e-6;
    }
    
    if ((state->t + state->dt > config->tstop) && (state->t < config->tstop)) {
        state->dt = config->tstop - state->t;
    }
    
    printf("dt: %e (max_nu: %e, Adv_x: %e, Adv_y: %e)\n",
           state->dt, max_nu, max_cx/config->dx, max_cy/config->dy);
}

void runge_kutta_step1(SimulationState *state, const GridConfig *config,
                      const PhysicalConstants *pc) {
    #pragma omp parallel
    {
        #pragma omp for collapse(2) nowait
        for (int i = 2; i < config->nx - 1; i++) {
            for (int i2 = 2; i2 < config->ny - 1; i2++) {
                for (int j = 0; j < config->ncj; j++) {
                    for (int k = 0; k < config->nck; k++) {
                        for (int l = 0; l < config->ncl; l++) {
                            double rhs = -1.0 / config->dx * (state->F[i][i2][j][k][l] -
                                                              state->F[i-1][i2][j][k][l]) -
                                        1.0 / config->dy * (state->G[i][i2][j][k][l] -
                                                           state->G[i][i2-1][j][k][l]) -
                                        state->freq_coll[i][i2] / pc->t_ref *
                                        (state->f[i][i2][j][k][l] - state->feq[i][i2][j][k][l]);
                            
                            state->fb[i][i2][j][k][l] = state->f[i][i2][j][k][l];
                            state->f[i][i2][j][k][l] = state->fb[i][i2][j][k][l] + state->dt * rhs;
                        }
                    }
                }
            }
        }
    }
}

void runge_kutta_step2(SimulationState *state, const GridConfig *config,
                      const PhysicalConstants *pc) {
    #pragma omp parallel
    {
        #pragma omp for collapse(2) nowait
        for (int i = 2; i < config->nx - 1; i++) {
            for (int i2 = 2; i2 < config->ny - 1; i2++) {
                for (int j = 0; j < config->ncj; j++) {
                    for (int k = 0; k < config->nck; k++) {
                        for (int l = 0; l < config->ncl; l++) {
                            double f_current = state->f[i][i2][j][k][l];
                            double rhs_star = -1.0 / config->dx * (state->F[i][i2][j][k][l] -
                                                                   state->F[i-1][i2][j][k][l]) -
                                             1.0 / config->dy * (state->G[i][i2][j][k][l] -
                                                                state->G[i][i2-1][j][k][l]) -
                                             state->freq_coll[i][i2] / pc->t_ref *
                                             (f_current - state->feq[i][i2][j][k][l]);
                            
                            // Write directly to f, eliminating the need for fn array
                            state->f[i][i2][j][k][l] = 0.5 * (state->fb[i][i2][j][k][l] +
                                                              f_current +
                                                              state->dt * rhs_star);
                        }
                    }
                }
            }
        }
    }
}

//*******************************************************************************************
// Integration of Macroscopic Quantities
//*******************************************************************************************

void integrate_macroscopic_quantities(SimulationState *state, const GridConfig *config) {
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for (int i = 0; i < config->nx; i++) {
            for (int i2 = 0; i2 < config->ny; i2++) {
                state->n[i][i2] = 0.0;
                state->v1[i][i2] = 0.0;
                state->v2[i][i2] = 0.0;
                state->p[i][i2] = 0.0;
                
                for (int k = 0; k < config->nck; k++) {
                    for (int l = 0; l < config->ncl; l++) {
                        state->sum1_n[i][i2][k][l] = 0.0;
                        state->sum1_v1[i][i2][k][l] = 0.0;
                        state->sum1_p[i][i2][k][l] = 0.0;
                    }
                }
                
                for (int j = 0; j < config->ncj; j++) {
                    for (int l = 0; l < config->ncl; l++) {
                        state->sum1_v2[i][i2][j][l] = 0.0;
                    }
                }
                
                for (int l = 0; l < config->ncl; l++) {
                    state->sum2_n[i][i2][l] = 0.0;
                    state->sum2_v1[i][i2][l] = 0.0;
                    state->sum2_v2[i][i2][l] = 0.0;
                    state->sum2_p[i][i2][l] = 0.0;
                }
            }
        }
        
        // Number density
        #pragma omp for collapse(2)
        for (int i = 0; i < config->nx; i++) {
            for (int i2 = 0; i2 < config->ny; i2++) {
                for (int j = 0; j < config->ncj - 1; j++) {
                    for (int k = 0; k < config->nck; k++) {
                        for (int l = 0; l < config->ncl; l++) {
                            state->sum1_n[i][i2][k][l] += 0.5 * (state->f[i][i2][j][k][l] +
                                                                state->f[i][i2][j+1][k][l]) *
                                                         config->dcj;
                        }
                    }
                }
                for (int k = 0; k < config->nck - 1; k++) {
                    for (int l = 0; l < config->ncl; l++) {
                        state->sum2_n[i][i2][l] += 0.5 * (state->sum1_n[i][i2][k][l] +
                                                         state->sum1_n[i][i2][k+1][l]) *
                                                  config->dck;
                    }
                }
                for (int l = 0; l < config->ncl - 1; l++) {
                    state->n[i][i2] += 0.5 * (state->sum2_n[i][i2][l] +
                                             state->sum2_n[i][i2][l+1]) * config->dcl;
                }
            }
        }
        
        // Velocity x-direction
        #pragma omp for collapse(2)
        for (int i = 0; i < config->nx; i++) {
            for (int i2 = 0; i2 < config->ny; i2++) {
                for (int j = 0; j < config->ncj - 1; j++) {
                    for (int k = 0; k < config->nck; k++) {
                        for (int l = 0; l < config->ncl; l++) {
                            state->sum1_v1[i][i2][k][l] += 0.5 * (state->c1[j] * state->f[i][i2][j][k][l] +
                                                                 state->c1[j+1] * state->f[i][i2][j+1][k][l]) *
                                                          config->dcj;
                        }
                    }
                }
                for (int k = 0; k < config->nck - 1; k++) {
                    for (int l = 0; l < config->ncl; l++) {
                        state->sum2_v1[i][i2][l] += 0.5 * (state->sum1_v1[i][i2][k][l] +
                                                          state->sum1_v1[i][i2][k+1][l]) *
                                                   config->dck;
                    }
                }
                for (int l = 0; l < config->ncl - 1; l++) {
                    state->v1[i][i2] += 0.5 * (state->sum2_v1[i][i2][l] +
                                              state->sum2_v1[i][i2][l+1]) * config->dcl;
                }
                if (state->n[i][i2] > 0.0) {
                    state->v1[i][i2] /= state->n[i][i2];
                }
            }
        }
        
        // Velocity y-direction
        #pragma omp for collapse(2)
        for (int i = 0; i < config->nx; i++) {
            for (int i2 = 0; i2 < config->ny; i2++) {
                for (int j = 0; j < config->ncj; j++) {
                    for (int k = 0; k < config->nck - 1; k++) {
                        for (int l = 0; l < config->ncl; l++) {
                            state->sum1_v2[i][i2][j][l] += 0.5 * (state->c2[k] * state->f[i][i2][j][k][l] +
                                                                 state->c2[k+1] * state->f[i][i2][j][k+1][l]) *
                                                          config->dck;
                        }
                    }
                }
                for (int j = 0; j < config->ncj - 1; j++) {
                    for (int l = 0; l < config->ncl; l++) {
                        state->sum2_v2[i][i2][l] += 0.5 * (state->sum1_v2[i][i2][j][l] +
                                                          state->sum1_v2[i][i2][j+1][l]) *
                                                   config->dcj;
                    }
                }
                for (int l = 0; l < config->ncl - 1; l++) {
                    state->v2[i][i2] += 0.5 * (state->sum2_v2[i][i2][l] +
                                              state->sum2_v2[i][i2][l+1]) * config->dcl;
                }
                if (state->n[i][i2] > 0.0) {
                    state->v2[i][i2] /= state->n[i][i2];
                }
            }
        }
        
        // Pressure
        #pragma omp for collapse(2)
        for (int i = 0; i < config->nx; i++) {
            for (int i2 = 0; i2 < config->ny; i2++) {
                for (int j = 0; j < config->ncj - 1; j++) {
                    for (int k = 0; k < config->nck; k++) {
                        for (int l = 0; l < config->ncl; l++) {
                            double diff1 = state->c1[j] - state->v1[i][i2];
                            double diff2 = state->c1[j+1] - state->v1[i][i2];
                            state->sum1_p[i][i2][k][l] += 0.5 * (diff1 * diff1 * state->f[i][i2][j][k][l] +
                                                                diff2 * diff2 * state->f[i][i2][j+1][k][l]) *
                                                         config->dcj;
                        }
                    }
                }
                for (int k = 0; k < config->nck - 1; k++) {
                    for (int l = 0; l < config->ncl; l++) {
                        state->sum2_p[i][i2][l] += 0.5 * (state->sum1_p[i][i2][k][l] +
                                                         state->sum1_p[i][i2][k+1][l]) *
                                                  config->dck;
                    }
                }
                for (int l = 0; l < config->ncl - 1; l++) {
                    state->p[i][i2] += 0.5 * (state->sum2_p[i][i2][l] +
                                             state->sum2_p[i][i2][l+1]) * config->dcl * 2.0;
                }
            }
        }
        
        // Temperature
        #pragma omp for collapse(2)
        for (int i = 0; i < config->nx; i++) {
            for (int i2 = 0; i2 < config->ny; i2++) {
                if (state->n[i][i2] > 0.0) {
                    state->T[i][i2] = state->p[i][i2] / state->n[i][i2];
                }
            }
        }
    }
}

//*******************************************************************************************
// Utility Functions
//*******************************************************************************************

double det3x3(double a11, double a12, double a13,
             double a21, double a22, double a23,
             double a31, double a32, double a33) {
    return a11 * (a22 * a33 - a23 * a32) -
           a12 * (a21 * a33 - a23 * a31) +
           a13 * (a21 * a32 - a22 * a31);
}

double max3(double a, double b, double c) {
    double m = a;
    if (b > m) m = b;
    if (c > m) m = c;
    return m;
}

double min3(double a, double b, double c) {
    double m = a;
    if (b < m) m = b;
    if (c < m) m = c;
    return m;
}

double minmod(double a, double b, double c) {
    if (a > 0 && b > 0 && c > 0) {
        return min3(a, b, c);
    } else if (a < 0 && b < 0 && c < 0) {
        return max3(a, b, c);
    } else {
        return 0.0;
    }
}

//*******************************************************************************************
// Output Functions
//*******************************************************************************************

void output_macroscopic_csv(const SimulationState *state, const GridConfig *config, int step) {
    char fname[256];
    snprintf(fname, sizeof(fname), "2d3v_refactored_%07d_%.6f.csv", step, state->t);
    
    FILE *fp = fopen(fname, "w");
    if (!fp) {
        perror("Cannot open output file");
        return;
    }
    
    fprintf(fp, "i,j,x,y,n,v1,v2,p,T,mfp\n");
    for (int i = 0; i < config->nx; i++) {
        for (int i2 = 0; i2 < config->ny; i2++) {
            fprintf(fp, "%d,%d,%f,%f,%f,%f,%f,%f,%f,%f\n",
                   i, i2, state->x[i], state->y[i2],
                   state->n[i][i2], state->v1[i][i2], state->v2[i][i2],
                   state->p[i][i2], state->T[i][i2], state->mfp[i][i2]);
        }
    }
    
    fclose(fp);
}

void output_distribution_csv(const SimulationState *state, const GridConfig *config, int step) {
    char dirname[256];
    snprintf(dirname, sizeof(dirname), "f_output_%07d_%.6f", step, state->t);
    
    #ifdef _WIN32
    _mkdir(dirname);
    #else
    mkdir(dirname, 0777);
    #endif
    
    for (int i = 0; i < config->nx; i++) {
        for (int i2 = 0; i2 < config->ny; i2++) {
            char fname[512];
            snprintf(fname, sizeof(fname), "%s/f_%d_%d.csv", dirname, i, i2);
            
            FILE *fp = fopen(fname, "w");
            if (!fp) {
                perror("Cannot open distribution file");
                continue;
            }
            
            fprintf(fp, "j,k,l,cx,cy,cz,f\n");
            for (int j = 0; j < config->ncj; j++) {
                for (int k = 0; k < config->nck; k++) {
                    for (int l = 0; l < config->ncl; l++) {
                        fprintf(fp, "%d,%d,%d,%f,%f,%f,%e\n",
                               j, k, l, state->c1[j], state->c2[k], state->c3[l],
                               state->f[i][i2][j][k][l]);
                    }
                }
            }
            
            fclose(fp);
        }
    }
}
