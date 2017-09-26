// Header file
#include "dnaMC2.h"

// Includes
#include <signal.h>
#include <unistd.h>

// MC Sweep size
const int ENSSTEP = 147;



// START OF SIGNAL HANDLING CODE
//
// consider as black box
int sigint_received = 0;
int sigterm_received = 0;
int sigquit_received = 0;

void handle_sigint(int sig)  { sigint_received = 1;  }
void handle_sigterm(int sig) { sigterm_received = 1; }
void handle_sigquit(int sig) { sigquit_received = 1; }

static void setup_signal_handler(int sig, void (*handler)(  )) {
    #if _POSIX_VERSION > 198800L
    struct sigaction action;
    
    action.sa_handler = handler;
    sigemptyset(&(action.sa_mask));
    sigaddset(&(action.sa_mask), sig);
    action.sa_flags = 0;
    sigaction(sig, &action, 0);
    #else
    signal(sig, handler);
    #endif
}

static int signal_was_caught(void)
{
    if (sigint_received) fprintf(stderr, "SIGINT received!\n");
    if (sigterm_received) fprintf(stderr, "SIGTERM received!\n");
    if (sigquit_received) fprintf(stderr, "SIGQUIT received!\n");
    return (sigint_received || sigterm_received || sigquit_received);
}
// END OF SIGNAL HANDLING CODE



// Globally accessible sequence arrays
Base *Sequence;
Base *PullSequence;

// Array that keeps track of the right-hand base pairs of the 28 base pair steps that have
// fixed midplanes in the RBP nucleosome model.
int boundbp[29] = {3, 7, 15, 18, 25, 30, 35, 39, 46, 50, 56, 60, 66, 70, 77, 81, 87, 91, 97, 101, 108, 112, 117, 122, 129, 132, 140, 144, 5000};


// Initialization function. Currently only sets up the 
// random number generator.
void init(int seed) {
    
    gsl_rng_env_setup();
    r = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (r, seed);
    
}


// Return 1 if base is bound to base+1, -1 if to base-1, 0 if free.
// The function determines this based on the boundbp array defined above.
int isBound(int base) {    
    int base2 = base + 1;
    int n;
    
    // Copy base into the last (unused) element of boundbp
    boundbp[28] = base;
    
    // Loop over boundbp until we find the index at which
    // base is located
    for ( n = 0; boundbp[n] != base; n++ );
    
    // If base was not found until at the end of the array,
    // it's not present in the array.
    if ( n == 28 ) {
        // But boundbp only contains one bp for every set
        // of bound base pairs, so we also have to check
        // for base+1. We use the same procedure.
        boundbp[28] = base2;
        for ( n = 0; boundbp[n] != base2; n++ );
        
        // If the potential complement was also not found,
        // the base pair is not bound. Return 0.
        if ( n == 28 ) { return 0; }
        
        // If it was found, the base pair is bound to 
        // base+1. Return 1.
        else { return 1; }
    }
    
    // If the outer if-statement above triggered, the 
    // following line will never be reached. If we do
    // get here, it means base was found in boundbp.
    // Return -1.
    return -1;
    
}

// Convert Base number to Bond number (BaseToBond)
//
// This function takes the number of the base along a nucleosome and checks
// if that base is one of the ones listed in the boundbp array.
int BaseToBond(int base) {
    int n;
    for ( n = 0; n < 28; n++ ) {
        if ( boundbp[n] == base ) {
            return n;
        }
    }
    return -1;
}

// Procedure to make sure the given matrix is orthogonal.
// It does so by converting to axis-angle representation
// and back to a rotation matrix. If the rotation matrix 
// is not perfectly orthogonal, the axis and angle will
// also come out imperfectly. The procedure forcibly
// renormalizes the axis and assumes the angle is a good
// enough approximation. Should only be used on matrices
// that are very approximately orthogonal.
void OrthogonalizeR(double R[3][3]) {
    double a[3];
    int i;
    
    // Get rotation angle
    double rot_ang = acos(0.5*(R[0][0] + R[1][1] + R[2][2] - 1.0));
    
    // Get rotation axis. 
    a[0] = R[2][1] - R[1][2];
    a[1] = R[0][2] - R[2][0];
    a[2] = R[1][0] - R[0][1];
    
    // Normalize axis.
    double t = 1.0/sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    for ( i = 0; i < 3; i++ ) {
        a[i] *= t;
    }   
    
    // Rebuild rotation matrix. Store in same array.
    RfromAA(rot_ang, a, R);
}

// Main function, starting point of the program. The arguments
// encode the number and values of the command line arguments
// passed when calling the program.
int main(int argc, char *argv[])
{
    // Start signal handling
    setup_signal_handler(SIGINT, handle_sigint);
    setup_signal_handler(SIGTERM, handle_sigterm);
    setup_signal_handler(SIGQUIT, handle_sigquit);
    
    // Loop counters
    int i, n, n1;
    
    // Auxiliary
    int accept, tries;
    
    // Temporary base pair identity storage
    Base stor;
    
    // Step size along long sequence
    int dec = 1;
    
    // MC loop counter and length
    long int c;
    long int stop = 1e3;
    
    // Inverse temperature
    double beta = 1.0;
    
    // Equipartition energy, currently not used
    // double Eth = 3.0*(147.0-28.0)/beta;
    
    // DNAState arrays
    DNAState state[N];
    DNAState start[N];
    DNAState trial[N];
    
    // Midplane arrays
    DNAState midstate[N-1];
    DNAState midplanes[28];
    
    // BPS energy bookkeeping arrays. These arrays record the
    // local elastic energies in all base pair steps in the system
    // for the current state of the simulation, as well as for the
    // trial state under consideration during a move attempt.
    // The reason for these arrays' existence is that in order to
    // calculate the energy difference of a trial move with the current
    // state we need the energy of both the trial and current states.
    // However, we can reduce the number of calculations required by
    // realizing the following:
    // 1) When we accept a move, the new "current" energy becomes the
    //    trial energy. We can save calculations in the next move by 
    //    remembering this energy.
    // 2) We perform many localized moves that only affect the energy
    //    of a small part of the system. By breaking up our recollection
    //    of the energy down to the BPS level, we can recalculate and
    //    update only the base pair steps affected by a move.
    // The above two optimizations are addressed by introducing these
    // two bookkeeping arrays. The updating of the Trial array is done
    // within the various Move functions, copying over of the appropriate
    // values is done below in the MC loop.
    double BPSEnergies[N-1], BPSEnergiesTrial[N-1];
    
    // Stepsize for random MC moves
    double stepsize_t = 1.0/(sqrt(beta)*15.0); 
    double stepsize_a = PI/(sqrt(beta)*150.0);
    
    // Variables for the MC loop:
    // e  - energy of current state
    // d  - random floating point container
    // w  - Boltzmann weight container
    // dE - energy difference for trial move
    // s  - random floating point for move type selection
    double e, d, w, dE, s;
    
    // Mutation move probability
    double Pmut = 0.2;
    
    // Variables related to bound base pairs
    int bb, bplo, bphi;
    
    // Random seed for the random number generator.
    // It is set here to an arbitrary number. If a
    // seed value is provided on the command line,
    // this value is overwritten. Not providing a
    // value will lead to this same value always
    // being used, with reproducible results.
    int seed = 10089;
    
    // Start position and number of base pairs
    // to indicate region along long sequence
    // to analyse
    int startpos, nstep;
    
    // Auxiliary variables for reading in sequences
    int numCharacters = 0, cnt = 0;
    
    // Read in sequence start position and number of steps 
    // to take from the command line. See README for full
    // definition of the command line arguments. See comments
    // describing the function of the two main nested loops
    // below for a description of how these variables are used.
    if (argc > 1) {
        startpos = atoi(argv[1]);
        nstep = atoi(argv[2]);
        if ( argc > 3 ) {
            seed = atoi(argv[3]);
        }
    }
    else { // Default is to start at the beginning and take only the first position
        startpos = 0;
        nstep = 1;
    }
    
    // Set the nucleosome-bound sequence to be the 147-bp sequence starting at 
    // the starting position
    Sequence += startpos;
    
    // Initialization routines
    init(seed);
    fillK();
    
    // Read in the sequence once to determine the length
    fprintf(stderr, "Reading sequence length.\n");
    FILE* seqfile = fopen("./State/Ctract.seq", "r");
    char nextChar = getc(seqfile);
    while (nextChar != EOF) {
        if ( isspace(nextChar) ) { nextChar = getc(seqfile); }
        else {
            numCharacters++;
            nextChar = getc(seqfile);
        }
    }
    fclose(seqfile);
    
    // Allocate memory to hold the sequence
    fprintf(stderr, "Allocating space for sequence.\n");
    PullSequence = (Base*) malloc(numCharacters*sizeof(int));
    Sequence = PullSequence;
    
    // Read in the sequence again and store in memory
    fprintf(stderr, "Reading sequence.\n");
    seqfile = fopen("./State/Ctract.seq", "r");
    nextChar = getc(seqfile);
    while (nextChar != EOF) {
        if ( isspace(nextChar) ) { nextChar = getc(seqfile); }
        else {
            Sequence[cnt] = CharToBase(nextChar);
            cnt++;
            nextChar = getc(seqfile);
        }
    }
    fclose(seqfile);
    fprintf(stderr, "Sequence read.\n");
    
    // Read in nucleosome configuration file
    fprintf(stderr, "Reading state.\n");
    FILE* statefile = fopen("State/Nucleosome.state", "r");
    ReadDNAState(statefile, state);
    fclose(statefile);
    fprintf(stderr, "State read.\n");
    
    // Calculate midplanes for the state just read in
    midplanestate(state, midstate);
    
    // Store the 28 relevant midframes in a new array
    for ( i = 0; i < N; i++ ) {
        if ( BaseToBond(i) >= 0 ) {
            //fprintf(stderr, "%d\n", i);
            midplanes[BaseToBond(i)] = midstate[i-1];
        }
    }
    
    // Set the starting configuration for the MC loops
    // to the state just read in
    CopyDNAState(state, start, 0, N-1);
    
    // Starting energy with the given configuration
    // and sequence.
    e = Energy(state, Sequence);
    fprintf(stderr, "AVG Energy: %e\n", e);
    
    // Done initializing, let's start
    fprintf(stderr, "Starting MC.\n");
    
    // The program has two nested loops:
    //
    // 1) Loop with counter i. Starting at position startpos along a
    // sequence of length at least 147 + startpos + nstep, it
    // loops over all positions (startpos + n*dec) < startpos + nstep
    // along the given sequence and runs MC for the 147-bp subsequence
    // at each position
    //.
    // 2) The actual MC loop with counter c. It runs a number of sweeps
    // defined by the variable stop. In each sweep, it makes ENSSTEP moves.
    // Only accepted moves are counted.
    
    for ( i = startpos; i < startpos+nstep; i+=dec ) {
        // Initialize variables
        
        // Always start in the same state
        CopyDNAState(start, state, 0, N-1);
        
        // Clean trial state
        CopyDNAState(state, trial, 0, N-1);
        
        // Start energy tracking
        e = Energy(state, Sequence);
        
        for ( n = 0; n < N-1; n++ ) {
            // For safety, orthogonalize the orientations before starting
            OrthogonalizeR(state[n].R);
            
            // Set up bookkeeping arrays
            BPSEnergies[n] = BPSEnergy(state[n], state[n+1], Sequence[n], Sequence[n+1]);
            BPSEnergiesTrial[n] = BPSEnergies[n];
        }
        
        // Reset number of attempted MC moves
        tries = 0;
        
        // Start MC loop
        // We make a number of sweeps defined by stop
        // Each sweep consists of ENSSTEP moves
        for ( c = 0; c < (ENSSTEP*stop); c++ ) {
            
            // Loop until accepted move is found
            accept = 0;
            while ( accept == 0 ) {
                // Increase attempt number
                tries += 1;
                
                // Select random base pair
                n1 = gsl_rng_uniform_int(r, N);     
                
                // Check if base pair is bound
                bb = isBound(n1);                   
                
                // If the base pair is bound, bb is either +1 or -1,
                // indicating whether it is bound to its left or right
                // neighbour. The following ternary expressions lead to
                // the following results, depending on this bound state:
                //
                // Case 1: n1 is not bound -> bplo = bphi = n1
                // Case 2: n1 is bound to the left neighbour -> bplo = n1-1, bphi = n1
                // Case 3: n1 is bound to the right neighbour -> bplo = n1, bphi = n1+1
                //
                // The bphi and bplo variables thus encode which base pairs
                // need to be moved, and thus which energies neeed to be
                // recalculated.
                //
                // Explanation of the ternary expressions:
                // 
                // bplo : if n1 < n1 + bb, then the base pair is bound to its right neighbour,
                // and bplo needs to be n1. If not, then there are two possibilities:
                // 1) bb = 0  : base pair not bound, bplo should be equal to n1, but we can safely add bb as it is 0.
                // 2) bb = -1 : base pair bound on the left, bplo needs to be n1 - 1 = n1 + bb.
                // Thus, in both cases bplo = n1 + bb.
                bplo = ( n1 < n1+bb ? n1 : n1+bb ); 
                // bphi : analogous to the above. If n1 < n1 + bb, we are bound on the right, so we
                // need to add bb to n1. In the other two cases, bphi = n1.
                bphi = ( n1 < n1+bb ? n1+bb : n1 );
                
                // Select random move type
                s = gsl_rng_uniform(r);     
                
                // Store the sequence identity of the selected base pair
                // (only needed in case of mutation).    
                stor = Sequence[n1];
                
                // Select mutation with probability Pmut
                if ( s < Pmut ) {
                    dE = Mutation(state, Sequence, n1, &stor, BPSEnergies, BPSEnergiesTrial);
                }
                // Else, we do a configurational move. If the base pair 
                // is bound, keep midplane fixed.
                else if ( bb ) { 
                    dE = SBPMoveFMP(bplo, bphi, stepsize_a, stepsize_t, state, trial, BPSEnergies, BPSEnergiesTrial, midplanes);
                }
                // Otherwise, move freely.
                else {      
                    dE = SBPMove(n1, 2.0*stepsize_a, 2.0*stepsize_t, state, trial, BPSEnergies, BPSEnergiesTrial); 
                }
                
                // Acceptance conditions.
                // If the move lowers the energy, we always accept it.
                if ( dE < 0 ) { 
                    // Setting accept = 1 will end the while loop after this iteration
                    accept = 1;
                    
                    // Copy the trial state to the current state
                    CopyDNAState(trial, state, bplo, bphi);
                    
                    // Copy the trial energies as well
                    CopyBPSEnergies(BPSEnergiesTrial, BPSEnergies, bplo-1, bphi);
                    
                    // Copy the mutation (does nothing if we didn't mutate)
                    Sequence[n1] = stor;
                    
                    // Update the energy
                    e += dE;
                }
                // If the move doesn't lower the energy, we accept it with
                // probability exp(-beta*dE).
                else { 
                    // Random float between 0 and 1.
                    d = gsl_rng_uniform(r);
                    
                    // Boltzmann weight of the move.
                    w = exp(-beta*dE); 
                    
                    // With probability given by w, accept the move anyway.
                    if ( d < w ) {
                        // Same acceptance operations as above.
                        accept = 1;
                        CopyDNAState(trial, state, bplo, bphi);
                        CopyBPSEnergies(BPSEnergiesTrial, BPSEnergies, bplo-1, bphi);
                        e += dE;
                        Sequence[n1] = stor;
                    }
                    // With probability 1-w, reject the move.
                    else {
                        // Clean up trial state
                        CopyDNAState(state, trial, bplo, bphi);
                        
                        // Clean up bookkeeping
                        CopyBPSEnergies(BPSEnergies, BPSEnergiesTrial, bplo-1, bphi);
                    }
                }
            }
            
            // We usually want to output something during the simulation.
            // Here are some examples. The current statement will print
            // the sweep number, the tracked energy, and a freshly calculated
            // energy. The latter two should be identical, otherwise there
            // is a leak somewhere.
            // The lines that are commented out would instead print the
            // full DNA configuration or the sequence.
            if ( c % (1*ENSSTEP) == 0 ) {
                fprintf(stdout, "%ld\t%e\t%e\n", c/ENSSTEP, e, Energy(state, Sequence));
                //PrintDNAState(stdout, state);
                //PrintSequence(stdout, Sequence);
                fflush(stdout);
            }
            
            // If any signal was caught, break out of the loop
            if (signal_was_caught(  )) break;
        }
        
        // If any signal was caught, break out of the loop
        if (signal_was_caught(  )) break;
    }
    
    // Clean up the random number generator.
    gsl_rng_free (r);
    
    // Cleanly exit, returning exit status 0 (success).
    return 0;
}
