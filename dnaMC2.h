// Include all the libraries we need across different files
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cblas.h>
#include <string.h>
#include <ctype.h>

// Defining some constants
#define PI 3.14159265358979

// System size needs to be hardcoded here
#define N 147

// The current code does not do tension, but
// I'll leave this here for the future.
#define f1pN 0.02414  // Pulling force of 1pN/(kT_r) * (m/Angstrom);
                      // 1pN of force, defined to give E/(kT) for room temperature if distance units are in Angstrom
                      


// Global random number generator
gsl_rng * r;  

// DNAState structure
typedef struct {
    double pos[3];
    double R[3][3];
} DNAState;

// Base variable type
typedef enum Bases {
	A, T, C, G, X
} Base;


// Some variables that are defined in the C files but
// need to be available across files.
extern double Eq[4][4][6];                     
extern double Kf[4][4][6][6];
extern Base *Sequence;
extern char *BaseString[5];
extern int boundbp[29];


// Function prototypes

void fillK();

// Move functions
double SBPMove(int base, double stepsize_t, double stepsize_a, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1]);
double SBPMoveFMP(int base, int complement, double stepsize_t, double stepsize_a, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1], DNAState midplanes[28]);
double Mutation(DNAState state[N], Base* seq, int a, Base* b, double BPSEnergies[N-1], double BPSEnergiesTrial[N-1]);

// Energy functions
double BPSEnergy(DNAState BP1, DNAState BP2, Base B1, Base B2);
double Energy(DNAState state[N], Base sequence[N]);


// Auxiliary functions

//// Related to bound base pairs in the nucleosome
int isBound(int base);
int BaseToBond(int base);

//// Related to midplanes
void midplane(DNAState BP1, DNAState BP2, DNAState* out);
void midplanestate(DNAState state[N], DNAState midstate[N-1]);

//// Matrix manipulations
void mmult2(double A[3][3], double B[3][3], double C[3][3]);
void LmmultT2(double A[3][3], double B[3][3], double C[3][3]);
void RmmultT2(double A[3][3], double B[3][3], double C[3][3]);
void OrthogonalizeR(double R[3][3]);

//// Numerical calculations
void randUnitVectorB(double n[3]);
int LinRandBP(int bp, double a);
void RfromAA(double angle, double axis[3], double R[3][3]);

//// Outputting and copying data structures
void PrintBPState(FILE* output, DNAState state);
void PrintDNAState(FILE* output, DNAState state[N]);
void CopyDNAState(DNAState state1[N], DNAState state2[N], int start, int end);
void ReadDNAState(FILE* file, DNAState state[N]);
void CopyBPSEnergies(double Esrc[N-1], double Etrg[N-1], int start, int end);
void PrintSequence(FILE* file, Base* sequence);
Base CharToBase(char c);
