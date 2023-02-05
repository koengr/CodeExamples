import cirq
from cirq.contrib import circuit_to_latex_using_qcircuit

import numpy as np
import matplotlib.pyplot as plt

from pprint import pprint # for pretty printing of lists



##################  Useful functions ########################
print_circuit = lambda circuit : "  " + (str(circuit).replace('\n','\n  ') if len(circuit) > 0 else "<<This circuit contains no gates.>>")

def show( circuit ):

    # Sorted qubits are required to 1) get the right order and 2) get rid of Frozenset
    qubits = sorted( circuit.all_qubits() )

    print("The circuit you implemented is:")
    print()
    print(print_circuit(circuit))
    print()

    print("The state probabilities are:")
    show_state_ampl( circuit.apply_unitary_effect_to_state(0) )


def tex( circuit ):
    print( circuit_to_latex_using_qcircuit(circuit) )


def show_state_prob( vector, max_len = 16 ):
    v = np.abs(vector)**2
    formatter = '{0:0' + str(int( np.ceil( np.log( len(v) ) / np.log(2) ))) + 'b}'
    if isinstance(max_len, int) and len(v) > max_len :
        v = v[:12]

    for j in range( len(v) ):
        print( formatter.format(j), v[j] )


def show_state_ampl( vector, max_len = 16 ):
    v = vector * 1
    formatter = '{0:0' + str(int( np.ceil( np.log( len(v) ) / np.log(2) ))) + 'b}'
    if isinstance(max_len, int) and len(v) > max_len :
        v = v[:12]

    for j in range( len(v) ):
        print( formatter.format(j), v[j] )

# Given a quantum state of the form |coin> |position>, return the position probabilities
#TODO add some "asserts" on the length of the input vector
def quantum_walk_position_prob( vector ):
    halflength =int( len(vector) / 2 )
    probs = np.abs( vector )**2
    posprobs = probs[:halflength] + probs[halflength:]
    return(posprobs)

# Given the quantum state vector, and the number of X and Y qubits nx, ny, return
# the probabilities to be at positions (x,y), as a matrix.
def quantum_walk_position_prob_2D( vector, nx, ny ):
    probs = np.abs( vector ) ** 2
    posmat = np.zeros( (2**nx, 2**ny) )
    coinvalue = 2**(nx+ny)
    xvalue = 2**ny
    for x in range(0, 2**nx) :
        posmat[x] = probs[x * xvalue : x * xvalue + 2**ny] + \
                    probs[x * xvalue + coinvalue : coinvalue + x * xvalue + 2**ny]
    return(posmat)

def rotate_list( vector, n ):
    return vector[n:] + vector[:n]

# Map the position probabilities of a "compact" walk (going only forward) into a "standard" walk (going both back and forth, and having twice as many sites)
# Index j goes to 2j-M, where M is the number of step operations performed in the walk.
def positions_compact_to_standard( vector, M ):
    n = len(vector)
    newvec = np.zeros( 2 * n )
    for j in range(0,n):
        newvec[ np.mod(2*j - M, 2*n) ] = vector[j]
    return newvec

# Map the position probabilities of a "compact" walk (going only forward) into a "standard" walk (going both back and forth, and having twice as many sites)
# Index j goes to 2j-M both for the X and Y coordinate, with M the number of steps.
# Note that this mapping is ONLY valid for a compact walker that goes both RIGHT and UP each step.
def positions_compact_to_standard_2D( matrix, M ):
    nx, ny = np.shape(matrix)
    newmat = np.zeros( (2*nx, 2*ny) )
    for x in range(0,nx):
        for y in range(0,ny):
            newmat[ np.mod(2*x-M,2*nx), np.mod(2*y-M,2*ny) ] = matrix[x,y]
    return newmat

# Force any operation to become a circuit - whether it's a circuit, gate or matrix.
def circuitform( c, qubits ):

    if isinstance( c, cirq.circuits.circuit.Circuit ):
        return c
    if isinstance( c, cirq.ops.gate_operation.GateOperation ):
        return( cirq.circuit( [c1(qubits)] ) )
    if isinstance( c, numpy.ndarray ):
        gate = cirq.SingleQubitMatrixGate( c )
        return( cirq.circuit( [gate(qubits)] ) )




# Returns the inverse of a circuit, by running the circuit backwards and inverting all gates.
def Conjugate(circuit):
    newcircuit = cirq.Circuit()
    newcircuit.append( reversed(
        [ cirq.inverse(op) for op in circuit.all_operations() ]
    ) )
    return newcircuit




#################### Basic operations #############################

''' A controlled-phase gate as used in fourier transfations.
    It implements the controlled version of Z^(+2^-exponent), e.g. exp=0 is the normal Z gate, exp=1 the normal S gate with phases (1,i).
    The exponent differs by 1 compared to Draper's definitions, where k=1 implements the Z gate.   '''
def fractional_cphase( qubits, exponent ):
    c = cirq.Circuit()
    c.append( cirq.CZ( qubits[0], qubits[1] ) ** ( 2**(-exponent) ) )
    return c


'''Convention:  Rotate as if we have a Hamiltonian of the form (0,0,0,-1):
        A phase of Pi should be a rotation with eigenvalues (1,1,1,-1)
        A phase of Pi/2 should be (1,1,1,i)        '''
def MS_CPhase( qubits, phase=np.pi ):
    exponent = phase/np.pi/2
    c = cirq.Circuit()
    c.append( [ cirq.H(qubits[0]), cirq.H(qubits[1]) ])
    c += MS_controlledX( qubits, phase )
    c.append( [ cirq.H(qubits[0]), cirq.H(qubits[1]) ])

    return c


''' The Controlled-X is basically the CPhase, but in the X-basis.
    Convention:  Rotate as if we have a Hamiltonian with eigenvalues (0,0,0,-1):
    A phase of Pi should be a rotation with eigenvalues (1,1,1,-1)
    A phase of Pi/2 should be (1,1,1,i)        '''
def MS_controlledX( qubits, phase=np.pi ):
    exponent = phase/np.pi/2
    c = cirq.Circuit()

    c.append( cirq.XXPowGate( exponent=-exponent, global_shift=0 )( qubits[0], qubits[1]) )
    c.append( [cirq.XPowGate( exponent=+exponent, global_shift=0 )( qubits[0] ),
                    cirq.XPowGate( exponent=+exponent, global_shift=0 )( qubits[1] ) ] )

    return c

''' Builds a unitary identical to fractional_cphase( qubits,exponent), but constructed from the MS interaction. '''
def MS_fractional_cphase( qubits, exponent ):
    phase = np.pi * 2**(-exponent)

    c = cirq.Circuit()
    c += Hadamards( qubits )
    c += MS_controlledX( qubits, phase )
    c += Hadamards( qubits )

    return c


''' Builds a unitary identical to MS_fractional_cphase(q1,q2,exponent),
    but without the Hadamards (to be sandwiched by the user manually) '''
def MS_fractional_cx( qubits, exponent ):
    phase = np.pi * 2**(-exponent)

    c = cirq.Circuit()
    c += MS_controlledX( qubits, phase )

    return c


#################### Functions for initialization ###################

''' Apply a Hadamard gate to all qubits, in parallel '''
def Hadamards(qubits):
    c = cirq.Circuit()
    c.append( [ cirq.H(qubits[j]) for j in range(0, len(qubits))] )
    return c

''' Apply the regular Fourier Transform to the given qubits (including an all qubits swap at the end) '''
def QFT(qubits):
    c = cirq.Circuit()
    n = len(qubits)

    for j in range( 0, n ):
        c.append( cirq.H( qubits[j]) )
        for i in range( j+1, n ) :
            k = 1 + (i-j)
            c.append( cirq.CZ( qubits[j], qubits[i] ) ** (2**(1-k) )  )

    for j in range( 0, len(qubits) // 2  ):
        c.append( cirq.SWAP( qubits[j], qubits[n-j-1] ) )

    return c

# Inverse Fourier Transform (inverses of other functions can be manually produced using Conjugate( circ )).
def QFTinverse(qubits):
    return Conjugate( QFT(qubits) )



# Apply the flipped QFT (flips qubit ordering)
def QFT_Flipped(qubits):

    c = cirq.Circuit()

    for j in range( 0, len(qubits) ):
        c.append( cirq.H(qubits[j]) )
        for i in range( j+1, len(qubits) ) :
            k = 1 + (i-j)
            c.append( cirq.CZ( qubits[j], qubits[i] ) ** (2**(1-k) )  )

    return c

# Apply the flipped QFT (flips qubit ordering) where each CZ is replaced by an MS gate
def QFT_Flipped_MS(qubits):
    c = cirq.Circuit()
    n = len(qubits)

    for j in range( 0, n ):
        c.append( cirq.H(qubits[j]) )

        # If it's not the last qubit, then perform controlled phases
        if j < n-1 :
            c += Hadamards(qubits[j:])
            for i in range( j+1, n ) :
                exponent = i-j
                c += MS_fractional_cx( [qubits[j], qubits[i]], exponent)
            c += Hadamards(qubits[j:])
    return c

# Apply the flipped QFT (flips qubit ordering) asssuming the input is in the x-basis (i.e. was acted on by Hadamards)
def QFT_Flipped_MS_Xbasis(qubits):
    c = cirq.Circuit()
    n = len(qubits)

    for j in range( 0, n ):
        c.append( cirq.H(qubits[j]) )

        # If it's not the last qubit, then perform controlled phases
        if j < n-1 :
            for i in range( j+1, n ) :
                exponent = i-j
                c += MS_fractional_cx( [qubits[j], qubits[i]], exponent)
    return c




################## Step operators ######################

# Assume qubit 0 is the coin, qubits 1...n-1 are fourier-transformed position register.
# Add the number 1 to the position.
def ConditionalAddOne( qubits ):
    c = cirq.Circuit()
    for j in range( 1, len(qubits) ): #note: start counting at 1
        # Here the *first* qubit is the most significant bit, rotating the least
        c += fractional_cphase( [qubits[0], qubits[j] ], j-1 )
    return c

def ConditionalAddOne_Flipped( qubits ):
    n = len(qubits)
    c = cirq.Circuit()
    for j in range( 1, n ): #note: start counting at 1
        # Here the *last* qubit is the most significant bit, rotating the least
        c += fractional_cphase( [qubits[0], qubits[j]], n-j-1 )
    return c

def ConditionalAddOne_Flipped_Xbasis( qubits ):
    n = len(qubits)
    c = cirq.Circuit()
    for j in range( 1, n ): #note: start counting at 1
        # Here the *last* qubit is the most significant bit, rotating the least
        c += MS_fractional_cx( [qubits[0], qubits[j]], n-j-1 )
    return c


################### Full Quantum Walk on the line ###########

def QuantumWalk1DFourier( qubits, num_steps ):

    c = cirq.Circuit()

    # initialization: Fourier transform
    c += QFT( qubits[1:] )

    # perform the steps
    for j in range( 1, num_steps+1 ):
        c.append( cirq.H(qubits[0]) )
        c += ConditionalAddOne( qubits )

    # post processing
    c+= Conjugate( QFT( qubits[1:] ) )

    return c


def QuantumWalk1DFourier_Flipped( qubits, num_steps ):

    c = cirq.Circuit()

    # initialization: Fourier transform
    c += QFT_Flipped( qubits[1:])

    # perform the steps
    for j in range( 1, num_steps+1 ):
        c.append( cirq.H(qubits[0]) )
        c+=ConditionalAddOne_Flipped(qubits)

    # post processing
    c += Conjugate( QFT_Flipped(qubits[1:]) )

    return c

# For the 2D walk (Compact Quantum Walk, moving right and moving up with some probability),
# we require as input: the coin-bit, the qubits representing the X and Y register,
# two coins c1, c2, and a number of steps.
def QuantumWalk2DFourier_Flipped( coinbit, qx, qy, c1, c2, num_steps ):

    c = cirq.Circuit()

    # initialization: Fourier transform
    c += QFT_Flipped( qx )
    c += QFT_Flipped( qy )

    # perform the steps
    for j in range( 1, num_steps+1 ):
        c.append( [c1(coinbit)] )
        c += ConditionalAddOne_Flipped( [coinbit] + qx )
        c.append( [c2(coinbit)] )
        c += ConditionalAddOne_Flipped( [coinbit] + qy )

    # post processing
    c += Conjugate( QFT_Flipped(qx) )
    c += Conjugate( QFT_Flipped(qy) )

    return c


# Note that the basic 1D DTQW circuit consists only of controlled-phase gates and Hadamards.
# Hence, we may sandwich the circuit between Hadamards at the beginning and end,
# and perform only MS_controlledX operations.
def TrappedIon1DWalk( qubits, num_steps ):

    c = cirq.Circuit()

    # special initialization: move to X basis
    c += Hadamards( qubits )

    # initialization: Fourier transform
    c += QFT_Flipped_MS_Xbasis( qubits[1:] )

    # perform the steps
    for j in range( 1, num_steps + 1 ):
        c.append( cirq.H(qubits[0]))
        c+=ConditionalAddOne_Flipped_Xbasis(qubits)

    # post-processing
    c += Conjugate( QFT_Flipped_MS_Xbasis( qubits[1:] ) )
    c += Hadamards( qubits )
    return c




#######################



if __name__ == "__main__":


    ## Force pretty printing of numpy arrays
    np.set_printoptions(suppress=True, precision=3)
    np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

    ## SETUP
    num_steps = 3
    n = 3
    qubits = [cirq.GridQubit(i,0) for i in range(n)]


    c1 = QuantumWalk1DFourier_Flipped( qubits, num_steps )
    state = c1.apply_unitary_effect_to_state(0)

    print("Quantum walk position probabilities:")
    print( quantum_walk_position_prob( state ) )


    c = QFT_Flipped_MS(qubits)
    print(c)
    u = c.to_unitary_matrix()

    c2 = QFT_Flipped(qubits)
    u2 = c2.to_unitary_matrix()



    # plt.bar( range(0,2**n), s1 )
    # plt.show()

    # for m in range(0,2):
    #
    #     circuit = QuantumWalk1DFourier_Flipped( qubits, m )
    #     finalstate = np.abs( circuit.apply_unitary_effect_to_state( 0 ) )
    #     print(m, finalstate)
