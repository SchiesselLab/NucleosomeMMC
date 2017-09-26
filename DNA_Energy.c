#include "dnaMC2.h"

// Sets the stretching force on the system. The numerical factor is to be given in units of pN.
double f = 0.0*f1pN;

// Averaged (over Kf) stiffness matrix             
double KAVG[6][6] = {{2.706075e+00, 2.768619e-15, 1.633416e-14, -2.282971e+00, -3.892719e-15, -2.578493e-14},
					{2.768619e-15, 4.581385e+00, 2.576186e+00, -1.101584e-13, -2.271799e+00, -1.558278e+01},
					{1.633416e-14, 2.576186e+00, 1.277129e+01, -2.063349e-13, -3.061377e-01, -1.277594e+01},
					{-2.282971e+00, -1.101584e-13, -2.063349e-13, 1.958269e+02, -6.039058e-13, -1.458195e-12},
					{-3.892719e-15, -2.271799e+00, -3.061377e-01, -6.039058e-13, 1.200616e+02, 4.039089e+01},
					{-2.578493e-14, -1.558278e+01, -1.277594e+01, -1.458195e-12, 4.039089e+01, 2.234921e+02}};

// Naive equilibrium values	
double EqDiff[6] = {0.0, 0.0, 3.4, 0.0, 0.0, 0.6};
           
// Tensor containing the 16 (first two indices) stiffness matrices (last two indices).          
double Kf[4][4][6][6];

// Tensor containing the 16 (first two indices) equilibrium vectors (last index).
double Eq[4][4][6];

// Array that converts between the enumerated values A, T, C, G and X, and their corresponding string values,
// in the sense that BaseString[A] = "A", etc.
char *BaseString[5] = {"A", "T", "C", "G", "X"};

// Load in the parameterization files
void fillK() {
    // Loop counters
	int i, j, k, l, m;
	
	// Buffer to contain the filename to open
	char filename[100];
	
	// Here we specify the parameterization files to input. The full set
	// (stiffness and equilibrium) is included for both the crystallography
	// and molecular dynamics parameterizations. From these, the hybrids
	// can also be mixed.
	
	// Base of the filenames for the stiffnesses
    char* basefilename = "Parameterization/Crystallography/Stiffness-C-";
    
    // Base of the filenames for the equilibrium values
	char* baseeqfilename = "Parameterization/MolecularDynamics/Equilibrium-MD-";
	
	// File object
	FILE* file;
	
	// Loop over all dinucleotides
	for ( i = A; i <= G; i++ ) { for ( j = A; j <= G; j++ ) {
	    // Load the base file name for stiffness into filename
		strcpy(filename, basefilename);
		// Add the first nucleotide as a letter
		strcat(filename, BaseString[i]);
		// Add the second nucleotide as a letter
		strcat(filename, BaseString[j]);
		// Append ".txt" to complete the filename
		strcat(filename, ".txt");
		// Open the file
		file = fopen(filename, "r");
		
		// The files have the rotational degrees of freedom listed first, but we want the translational ones first, therefore
		// we have a somewhat odd loop construction. But in essence, this just reads the files into Kf.
		k = 3;
		while (!feof(file)) {
       		if (fscanf(file,"%lf %lf %lf %lf %lf %lf", &Kf[i][j][k][3], &Kf[i][j][k][4], &Kf[i][j][k][5], &Kf[i][j][k][0], &Kf[i][j][k][1], &Kf[i][j][k][2]) != 6) {
         		continue;
       		}
       		k++;
       		if ( k > 5 ) { k = 0; }
    	}
    	fclose(file);
    	
    	// The files use nm as their unit of length. The program uses Angstrom, so we convert the stiffnesses appropriately.
    	for ( l = 0; l < 6; l++ ) { 
    		for ( m = 0; m < 6; m++ ) {
    			if ( l < 3 ) { if ( m < 3 ) { Kf[i][j][l][m] /= 100.0; } else { Kf[i][j][l][m] /= 10.0; } }
    			else { if ( m < 3 ) { Kf[i][j][l][m] /= 10.0; } }
    		}
		}
		
		// We now follow a similar procedure to load in the equilibrium values.
		strcpy(filename, baseeqfilename);
		strcat(filename, BaseString[i]);
		strcat(filename, BaseString[j]);
		strcat(filename, ".txt");
		file = fopen(filename, "r");
		k = 3;
		while (!feof(file)) {
       		if (fscanf(file,"%lf", &Eq[i][j][k]) != 1) {
         		continue;
       		}
       		k++;
       		if ( k > 5 ) { k = 0; }
    	}
   		for ( l = 0; l < 3; l++ ) { Eq[i][j][l] *= 10.0; }
	}}
}

// Base Pair Step Energy (BSPEnergy)
//
// This function calculates the elastic energy between two base pairs, given their
// locations and orientations, and their identities. Note:
//  1) All calculations need to take place in the midframe between the two base pairs,
//     since this is how the stiffness matrices and equilibrium values are defined.
//  2) The calculations require conversion between rotations given as rotation matrices
//     and as axis-angle pairs. The full details are not documented here.
//  3) The calculations require various transformations between different frames of
//     reference. The full details are not documented here.
//
// Return value:
// e    - calculated energy
//
// Arguments:
// BP1  - position and orientation of the first base pair, encoded as a DNAState variable
// BP2  - as BP1, for the second base pair
// B1   - identity of the first base pair
// B2   - identity of the second base pair

double BPSEnergy(DNAState BP1, DNAState BP2, Base B1, Base B2) {
    // Defining variables
    double e;
    double R[3][3];
    double Rtemp[3][3]; 
    double x[6], xrel[6], temp[6];
    int i, j;
    double a[3], b[3];
    double t;
    double rot_ang, frac;
    
    // We need the relative distance vector between the two base pairs in the midframe.
    // We also need the relative rotation between the two base pairs, in axis-angle
    // representation, in the midframe. The following process gives us both.
    
    // Calculate relative distance vector in the coordinates of the first base pair.
    // Note the use of the transpose/inverse of the rotation matrix.
    for ( i = 0; i < 3; i++ ) { ;
        x[i] = 0.0;
        for ( j = 0; j < 3; j++ ) {
            x[i] += BP1.R[j][i]*(BP2.pos[j] - BP1.pos[j]); 
        }
    }
    
    // Calculate the relative rotation in the coordinate frame of the first base pair
    // and store it in R.
    LmmultT2(BP1.R, BP2.R, R); 
    
    // Extract the rotation angle and unit vector from R. We store the angle in rot_ang,
    // the normalized axis in a and the non-normalized axis rot_ang*a (i.e. the degrees
    // of freedom we need) in the last three components of x.
    // (Note: we don't need to rotate this vector to the midframe, because the rotation
    // from the first base pair to the midframe is a rotation around a.)
    rot_ang = acos(0.5*(R[0][0] + R[1][1] + R[2][2] - 1.0));
    frac = 0.5*rot_ang / sin(rot_ang);
	x[3] = frac*(R[2][1] - R[1][2]);
	x[4] = frac*(R[0][2] - R[2][0]);
	x[5] = frac*(R[1][0] - R[0][1]);
	t = 1.0/sqrt(x[3]*x[3] + x[4]*x[4] + x[5]*x[5]);
	for ( i = 0; i < 3; i++ ) { a[i] = x[i+3]*t; }
	
	// Calculating the rotation matrix from axis a and angle rot_ang/2, we get the rotation
	// between the first base pair and the midframe.
	RfromAA(0.5*rot_ang, a, Rtemp);
	
	// Using that rotation, we can get the difference vector (which was stored in the first
	// 3 components of x) in the midframe coordinates.
	cblas_dgemv(CblasRowMajor, CblasTrans, 3, 3, 1.0, *Rtemp, 3, x, 1, 0.0, b, 1);
	for ( i = 0; i < 3; i++ ) { x[i] = b[i]; }
	
	// We now have all our degrees of freedom, and we can calculate the energy. First in the
	// case of A/T/C/G dinucleotides...
	if ( B1 < 4 && B2 < 4 ) {
	    // ...we compare our degrees of freedom with their equilibrium values...
	    for ( i = 0; i < 6; i++ ) { xrel[i] = x[i] - Eq[B1][B2][i]; }
	    // ...and calculate the quadratic form in two steps:
	    cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, *Kf[B1][B2], 6, xrel, 1, 0.0, temp, 1);
	    e = cblas_ddot(6, xrel, 1, temp, 1); }
    else { 
        // If the base pair step includes one or more non-A/T/C/G dinucleotides, assumed to
        // represent homogeneous DNA, we do the same calculation, but with the naive equilibrium
        // values, and averaged stiffness.
	    for ( i = 0; i < 6; i++ ) { xrel[i] = x[i] - EqDiff[i]; }
	    cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, *KAVG, 6, xrel, 1, 0.0, temp, 1);
	    e = cblas_ddot(6, xrel, 1, temp, 1);
    }
    
    // Finally we apply the factor 1/2 and return the energy.
    return 0.5*e;
}

// Total Elastic Energy (Energy)
//
// This function simply loops over all base pair steps and returns the total
// elastic energy stored in them.
//
// Return value:
// e        - total elastic energy
//
// Arguments:
// state    - DNA conformation for which to calculate energy
// sequence - DNA sequence for which to calculate energy

double Energy(DNAState state[N], Base sequence[N]) {
    double e = 0.0;
    int n;
       
    for ( n = 1; n < N; n++ ) {
        e += BPSEnergy(state[n-1], state[n], sequence[n-1], sequence[n]);
    }
    
    return e;
}
