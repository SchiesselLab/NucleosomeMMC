#include "dnaMC2.h"

// Single Base Pair Move (SBPMove)
//
// This function performs a random move (translation and rotation)
// on a single base pair.
//
// Return value:
// dE               - Energy difference between trial state and current state
//
// Arguments:
// base             - (zero-indexed) number of the base pair to be moved
// stepsize_t       - maximum step allowed in translational DOF
// stepsize_a       - maximum step allowed in rotational DOF
// currentState     - current accepted state of the system
// newState         - DNAState structure in which to store the new state after the move
// BPSEnergies      - Bookkeeping array containing the current elastic energies between successive base pairs
// BPSEnergiesTrial - Bookkeeping array in which to store the modified elastic energies between successive base pairs

double SBPMove(int base, double stepsize_t, double stepsize_a, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1]) {
    // Declaring variables
    int i;
    double dE = 0.0;
    double step;
    
    // For each of the translational degrees of freedom...
    for ( i = 0; i < 3; i++ ) {
        // ...generate a random value between -stepsize_t and +stepsize_t...
        step = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_t;
        // ...and add this to the current position of the base pair, storing the result in newState.
        newState[base].pos[i] = currentState[base].pos[i] + step;
    }    
    
    // Declaring variables for random rotation:
    
    // Random rotation axis
    double n[3];
    randUnitVectorB(n);
                   
    // Random rotation angle
    double theta = gsl_rng_uniform(r)*stepsize_a;         
    
    // Obtain rotation matrix            
    double R[3][3];
    RfromAA(theta, n, R);           
    
    // Rotate base pair, storing the result in newState
    mmult2(R, currentState[base].R, newState[base].R);      
    
    // Calculate the energy difference on either side of the base pair
    
    // Unless the very first base pair was moved...
    if ( base != 0 )    {     
        // ...calculate the new energy between the moved base pair and the previous one...
        BPSEnergiesTrial[base-1] = BPSEnergy(newState[base-1], newState[base], Sequence[base-1], Sequence[base]);
        // ...and store the difference.
        dE += BPSEnergiesTrial[base-1] - BPSEnergies[base-1];    
    }
    
    // Unless the very last base pair was moved...
    if ( base != N-1 )  {     
        // ...do the same for the moved base pair and the next one.
        BPSEnergiesTrial[base] = BPSEnergy(newState[base], newState[base+1], Sequence[base], Sequence[base+1]);    
        dE += BPSEnergiesTrial[base] - BPSEnergies[base];     
    }
    
    // Return the difference in energy between the trial state and the current state.
    return dE;
} 

// Single Base Pair Move, with a Fixed MidPlane (SBPMoveFMP)
//
// This function performs a random move (translation and rotation)
// on a single base pair, but with the constraint that the base pair
// has a complementary base pair that needs to be moved in the opposite
// direction in order to preserve the midplane between the two.
//
// Return value:
// dE               - Energy difference between trial state and current state
//
// Arguments:
// base             - (zero-indexed) number of the base pair to be moved
// complement       - (zero-indexed) number of the complementary base pair
// stepsize_t       - maximum step allowed in translational DOF
// stepsize_a       - maximum step allowed in rotational DOF
// currentState     - current accepted state of the system
// newState         - DNAState structure in which to store the new state after the move
// BPSEnergies      - Bookkeeping array containing the current elastic energies between successive base pairs
// BPSEnergiesTrial - Bookkeeping array in which to store the modified elastic energies between successive base pairs
// midplanes        - Array of the orientations and locations of the 28 midplanes that are fixed in the nucleosome

double SBPMoveFMP(int base, int complement, double stepsize_t, double stepsize_a, DNAState currentState[N], DNAState newState[N], double BPSEnergies[N-1], double BPSEnergiesTrial[N-1], DNAState midplanes[28]) {
    // Declaring variables
    int i;
    double dE = 0.0;
    double step;
    
    // Translational move identical to that in SBPMove(...), but the complement is
    // moved by the same amount in the opposite direction.
    for ( i = 0; i < 3; i++ ) {
        step = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_t;
        newState[base].pos[i] = currentState[base].pos[i] + step;
        newState[complement].pos[i] = currentState[complement].pos[i] - step;
    }

    // The implementation of the random rotation is somewhat non-trivial. One could just generate
    // a random rotation matrix and apply it to one base pair, while applying its inverse to the
    // complementary base pair. However, this approach could lead to a numerically unstable midplane.
    // To ensure that the midplane is absolutely fixed, we instead do the following:
    //  * Generate random rotation R
    //  * Calculate the rotation matrix Rmid between the orientation of the base pair and the midplane
    //  * Calculate the combined rotation matrix R.Rmid, i.e. the rotation between base pair and midplane plus a random rotation
    //  * Apply this to the midplane orientation to generate a randomly modified orientation for the base pair
    //  * Apply the inverse rotation to the complementary base pair to generate the appropriate orientation for
    //    the complementary base pair, regardless of its original orientation.
    // This scheme will rigidly enforce that the simulation always respects the original midframes.
    
    // Random rotation matrix
    double n[3];
    double R[3][3];
    double theta = (2.0*gsl_rng_uniform(r)-1.0)*stepsize_a; 
    randUnitVectorB(n);
    RfromAA(theta, n, R);

    // Here we find out which bond we are dealing with, in order to get the right midplane.
    // Note that we assume that the function has been passed a correct set of base pairs.
    // ToDo: Error Checking on this point.
    int bond = BaseToBond(base);
    
    // If the base pair is not in the list of bound base pairs, we assume the other one must be.
    if ( bond == -1 ) { bond = BaseToBond(complement); }
        
    // We obtain the rotation matrix Rmid between the base pair and the midplane...
    double Rmid[3][3], Rmove[3][3];
    RmmultT2(currentState[base].R, midplanes[bond].R, Rmid);
    
    // ...reorthogonalize it (because all the matrix multiplications
    // tend to be numerical unstable)...
    OrthogonalizeR(Rmid);

    // ...and combine it with our random matrix to obtain our move relative to the midframe
    mmult2(R, Rmid, Rmove);
    
    // Rotate one base pair by Rmove... 
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, *Rmove, 3, *midplanes[bond].R, 3, 0.0, *newState[base].R, 3);
    // ...and the other by the inverse of Rmove.
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, 3, 3, 1.0, *Rmove, 3, *midplanes[bond].R, 3, 0.0, *newState[complement].R, 3);
    
    // Calculate the new energies and the energy difference
    
    BPSEnergiesTrial[base-1] = BPSEnergy(newState[base-1], newState[base], Sequence[base-1], Sequence[base]);    
    dE += BPSEnergiesTrial[base-1] - BPSEnergies[base-1];

    BPSEnergiesTrial[complement] = BPSEnergy(newState[complement], newState[complement+1], Sequence[complement], Sequence[complement+1]);
    dE += BPSEnergiesTrial[complement] - BPSEnergies[complement]; 
                 
    BPSEnergiesTrial[base] = BPSEnergy(newState[base], newState[complement], Sequence[base], Sequence[complement]);       
    dE += BPSEnergiesTrial[base] - BPSEnergies[base];
    
    return dE;
} 

// Mutation Move
//
// This function performs a random mutation on a single base pair.
//
// Return value:
// dE               - Energy difference between trial state and current state
//
// Arguments:
// state            - current state of the simulation
// seq              - current sequence of the simulation
// a                - (zero-indexed) number of the base to mutate
// b                - Base variable in which to store the trial mutation
// BPSEnergies      - Bookkeeping array containing the current elastic energies between successive base pairs
// BPSEnergiesTrial - Bookkeeping array in which to store the modified elastic energies between successive base pairs

double Mutation(DNAState state[N], Base* seq, int a, Base* b, double BPSEnergies[N-1], double BPSEnergiesTrial[N-1]) {
    // Declare variables
	double dE = 0.0;
	double s;
	
	// Randomly pick a new base identity. The mutation must constitute an actual alteration, hence the while loop.
	*b = seq[a];
    while ( *b == seq[a] ) {
        s = gsl_rng_uniform(r);
        if ( s < 0.25 ) { *b = C; }
        else if ( s < 0.50 ) { *b = T; }
        else if ( s < 0.75 ) { *b = A; }
        else { *b = G; }
    }
	
	// Calculate the energy change
    if ( a > 0 ) {      BPSEnergiesTrial[a-1] = BPSEnergy(state[a-1], state[a], seq[a-1], b[0]);
                        dE += BPSEnergiesTrial[a-1] - BPSEnergies[a-1];  }
    if ( a < N-1 ) {    BPSEnergiesTrial[a] = BPSEnergy(state[a], state[a+1], b[0], seq[a+1]);
                        dE += BPSEnergiesTrial[a] - BPSEnergies[a];  }
        
    return dE;
}
