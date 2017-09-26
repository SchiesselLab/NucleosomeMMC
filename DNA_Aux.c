#include "dnaMC2.h"

// Print Base Pair State Configuration (PrintBPState)
//
// Prints a single bp DNAState to the specified output 
// filestream (can also be stdout/stderr). The format is x, y, 
// z, R00, R01, R02, etc., tab-separate, no trailing tab. The 
// printing commands are separated out purely for readability.
//
// Return value: 
// none
//
// Arguments: 
// output - output filestream
// state  - single DNAState instance representing one base pair
void PrintBPState(FILE* output, DNAState state) {
    // Print position vector components
    fprintf(output, "%e\t%e\t%e\t", state.pos[0], state.pos[1], state.pos[2]);
    // Print orientation matrix components
    fprintf(output, "%e\t%e\t%e\t", state.R[0][0], state.R[0][1], state.R[0][2]);
    fprintf(output, "%e\t%e\t%e\t", state.R[1][0], state.R[1][1], state.R[1][2]);
    fprintf(output, "%e\t%e\t%e", state.R[2][0], state.R[2][1], state.R[2][2]);
}

// Print DNA Chain Configuration (PrintDNAState)
//
// Calls PrintBPState for every base pair in a DNAState
// array, with a tab character in between.
//
// Return value:
// none
//
// Arguments:
// output - output filestream
// state  - array of DNAStates representing a DNA molecule
void PrintDNAState(FILE* output, DNAState state[N]) {
    int n;
    // For every base pair...
    for ( n = 0; n < N; n++ ) {
        // ...print its state.
        PrintBPState(output, state[n]);
        // Separate by tab characters, but don't add
        // a trailing tab.
        if ( n == N-1 ) { }
        else { fprintf(output, "\t"); }
    }
    fprintf(output, "\n");
}

// Copy DNA Chain Configuration (CopyDNAState)
//
// Copies (a subset of) the base pair configurations in
// a DNAState array to another array. Allows copying only
// a part of the configuration. This allows the function
// to be more efficient when copying trial states to current 
// states in the MC simulation when the MC move only affects
// a small part of the chain.
//
// Return value:
// none
//
// Arguments:
// state1 - source array
// state2 - target array
// start  - base pair index from which to start copying
// stop   - base pair index at which to stop copying
void CopyDNAState(DNAState state1[N], DNAState state2[N], int start, int end) {
    int n;
    // These checks can be optimized away if you're sure you're
    // calling the function correctly.
    
    // In case the start and end arguments are accidentally reversed.
    if ( end < start ) { int x = end; end = start; start = x; }
    
    // In case the start and end are not in the valid range.
    if ( start < 0 ) { start = 0; }
    if ( end > N-1 ) { end = N-1; }
    
    // Copy.
    for ( n = start; n <= end; n++ ) {
        state2[n] = state1[n];
    }
}

// Copy Base Pair Step Elastic Energies (CopyBPSEnergies)
//
// Copies (a subset of) the base pair step energies in
// an array to another array. Its use is similar to that of
// CopyDNAState, but it operates on the BPSEnergies(Trial)
// bookkeeping arrays.
//
// Return value:
// none
//
// Arguments:
// Esrc  - source array
// Etrg  - target array
// start - base pair index from which to start copying
// stop  - base pair index at which to stop copying
void CopyBPSEnergies(double Esrc[N-1], double Etrg[N-1], int start, int end) {
    int n;
    // These checks can be optimized away if you're sure you're
    // calling the function correctly.
    
    // In case the start and end arguments are accidentally reversed.
    if ( end < start ) { int x = end; end = start; start = x; }
    
    // In case the start and end are not in the valid range.
    if ( start < 0 ) { start = 0; }
    if ( end > N-2 ) { end = N-2; }
    
    // Copy.
    for ( n = start; n <= end; n++ ) {
        Etrg[n] = Esrc[n];
    }
}

// Read DNAState From File (ReadDNAState)
//
// Reads in a DNAState from the given file. It expects a 
// specific format: Nx12 tab-separated floating point numbers
// representing in order, x, y, z, R00, R01, R02, etc. for
// N base pairs consecutively. (Note: this is the format in
// which PrintDNAState outputs a DNAState array.) N must be
// equal to the system size set in the header file.
//
// Return value:
// none
//
// Arguments:
// file  - input filestream
// state - DNAState array into which to read the data
void ReadDNAState(FILE* file, DNAState state[N]) {
    int n = 0;
    // Using fscanf to read things in. The syntax is a bit
    // ugly, but that's how it works.
    while ( fscanf(file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf", &state[n].pos[0], &state[n].pos[1], &state[n].pos[2], &state[n].R[0][0], &state[n].R[0][1], &state[n].R[0][2], &state[n].R[1][0], &state[n].R[1][1], &state[n].R[1][2], &state[n].R[2][0], &state[n].R[2][1], &state[n].R[2][2]) != EOF ) {
        n++;
    }
    // Check that we read in the right amount of data.
    if ( n != N ) { fprintf(stderr, "Error reading data from file: length mismatch, n = %d.\n", n); exit(1); }
}

// Print DNA Sequence (PrintSequence)
//
// Prints a Base array to the specified filestream.
//
// Return value:
// none
//
// Arguments:
// file     - output filestream
// sequence - sequence array to be output
void PrintSequence(FILE* file, Base* sequence) {
    int i;
    char *SeqString;
    // Allocate memory for the string version of the sequence
	SeqString = malloc((N+1)*sizeof(char));
	SeqString[0] = '\0';
	// Convert each character in the sequence to a string
	// and save in the allocated memory
    for ( i = 0; i < N; i++ ) { strcat(SeqString, BaseString[sequence[i]]); }
    // Print the resulting string.
    fprintf(file, "%s\n", SeqString);
}


// Simplified wrappers for some CBLAS matrix multiplication
// functions.

// Basic Matrix Multiplication:
// Multiplies A with B and stores result in C.
void mmult2(double A[3][3], double B[3][3], double C[3][3])  {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, *A, 3, *B, 3, 0.0, *C, 3);
}

// Matrix Multiplication with left matrix transposed:
// Multiplies A^T with B and stores result in C.
void LmmultT2(double A[3][3], double B[3][3], double C[3][3])  {
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, 3, 3, 1.0, *A, 3, *B, 3, 0.0, *C, 3);
}

// Matrix Multiplication with right matrix transposed:
// Multiplies A with B^T and stores result in C.
void RmmultT2(double A[3][3], double B[3][3], double C[3][3])  {
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, 3, 3, 1.0, *A, 3, *B, 3, 0.0, *C, 3);
}


// Calculate BPS Midplane (midplane)
//
// Calculates the midplane position and orientation
// between two base pairs.
//
// Return value:
// none
//
// Arguments:
// BP1 - first base pair DNAState
// BP2 - second base pair DNAState
// out - DNAState object into which to write the result
void midplane(DNAState BP1, DNAState BP2, DNAState *out) {
    double R[3][3], Rm[3][3];
    double a[3];
    double rot_ang, frac, t;
    int i;
    
    // Get the rotation matrix between BP1 and BP2
    RmmultT2(BP2.R, BP1.R, R);
    
    // Get the rotation angle
    rot_ang = acos(0.5*(R[0][0] + R[1][1] + R[2][2] - 1.0));
    
    // Get the rotation axis
    frac = 0.5*rot_ang / sin(rot_ang);
    a[0] = frac*(R[2][1] - R[1][2]);
    a[1] = frac*(R[0][2] - R[2][0]);
    a[2] = frac*(R[1][0] - R[0][1]);
	t = 1.0/sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	for ( i = 0; i < 3; i++ ) {
	    a[i] *= t;
	}
	
	// Generate rotation matrix from the calculated
	// rotation axis, and half the rotation angle.
	// This rotation matrix rotates from BP1 to the
	// midframe.
	RfromAA(0.5*rot_ang, a, Rm);
	
	// Get the midframe orientation
	mmult2(Rm, BP1.R, (*out).R);
	
	// Get the midframe position
	for ( i = 0; i < 3; i++ ) {
	    (*out).pos[i] = 0.5*(BP1.pos[i] + BP2.pos[i]);
	}
}

// Calculate Midplane State (midplanestate)
// 
// Calculates the midplanes between all successive base pairs
// in a DNA chain state.
//
// Return value:
// none
//
// Arguments:
// 
// state    - DNA molecule state
// midstate - DNAState array into which to write the results
void midplanestate(DNAState state[N], DNAState midstate[N-1]) {
    int n;
    for ( n = 0; n < N-1; n++ ) {
        midplane(state[n], state[n+1], &midstate[n]);
    }
}


// Pick random base pair from lineair distribution around bp
// i.e. bp has the highest probability and away from bp the
// probability scales as 1 - a*|n - bp|/N.
// (This function not in use in the Nucleosome code, but can
// come in handy for crankshaft moves and such. Leaving it in
// in case.)
//
// Return value:
// m  - selected base pair index
//
// Arguments:
//
// bp - base pair index around which to center the distribution
// a  - slope of the linear decay
int LinRandBP(int bp, double a) {
    double random, func;
    int m, accept = 0;
    while ( accept == 0 ) {
        // Random float between 0 and 1
        random = gsl_rng_uniform(r);
        
        // Random integer between 0 and N-1
        m = gsl_rng_uniform_int(r, N-1);
        
        // Linear decay function
        func = 1.0 - a*fabs((double) (m-bp))/N;
        
        // Accept the trial with probability
        // given by func
        if ( random < func ) { accept = 1; }
    }
    return m;
}


// Generate a random unit vector, such that the
// probability density on the surface of the
// unit sphere is uniform. See 
// mathworld.wolfram.com/SpherePointPicking.html
//
// Return value:
// None
//
// Arguments:
// n - 3-array into which to store the result
void randUnitVectorB(double n[3]) {

    double a, x, b;
   
    // Random azimuth
    a = gsl_rng_uniform(r)*2*PI;
    
    // Random float between 0 and 1
    x = gsl_rng_uniform(r);
    
    // Random elevation
    b = acos(1-2*x);
    
    // Generate unit vector
    n[0] = sin(b)*cos(a);
    n[1] = sin(b)*sin(a);
    n[2] = cos(b);
}

// Rotation Matrix from Axis-Angle (RfromAA)
//
// Generates a rotation matrix from a given
// unit vector and angle.
//
// Return value:
// none
//
// Arguments:
// angle - rotation angle
// axis  - rotation axis
// R     - 3x3 array into which to store the result
void RfromAA(double angle, double axis[3], double R[3][3]) {
    double c = cos(angle);
    double s = sin(angle);
    double c1 = 1.0-c;
    double s2 = 2.0*s;
    
    R[0][0] = c + axis[0]*axis[0]*c1;
    
    R[1][0] = axis[0]*axis[1]*c1 + axis[2]*s;
    R[2][0] = axis[0]*axis[2]*c1 - axis[1]*s;
    R[0][1] = R[1][0] - axis[2]*s2;
    
    R[1][1] = c + axis[1]*axis[1]*c1;
    
    R[2][1] = axis[1]*axis[2]*c1 + axis[0]*s;
    R[0][2] = R[2][0] + axis[1]*s2;
    R[1][2] = R[2][1] - axis[0]*s2;
    
    R[2][2] = c + axis[2]*axis[2]*c1;
}

// Convert characters to Base variables (CharToBase)
//
// Given a character (one-letter string), give the
// corresponding Base variable. Character can be
// upper or lower case. Any character that is not
// A/a, T/t, C/c or G/g is converted into X (i.e.
// homogeneous DNA).
Base CharToBase(char c) {
	Base B;
	if ( c == 'A' || c == 'a' ) { B = A; }
	else if ( c == 'T' || c == 't' ) { B = T; }
	else if ( c == 'C' || c == 'c' ) { B = C; }
	else if ( c == 'G' || c == 'g' ) { B = G; }
	else { B = X; }
	return B;
}
