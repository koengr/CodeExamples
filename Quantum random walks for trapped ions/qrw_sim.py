from qrw import *
from qrw_tests import *

# This function pretends that we 'double' the number of sites (compared to the actual number
# of sites availlable in the n-1 qubit register).
def simulate_1D_walk(n, num_steps, initstate):
    ## SETUP
    qubits = [cirq.GridQubit(i,0) for i in range(n)]

    # map initstate to an actual position... not that initstate must be less than the number of positions!
    initpos = np.zeros(2**n)
    initpos[ initstate * 2 ] = 1.

    # Store a vector with probabilities per position, for each step.
    poss = [ initpos ]
    for s in range(1,num_steps+1):
        c1 = QuantumWalk1DFourier_Flipped( qubits, s )
        poss.append(  positions_compact_to_standard(
            quantum_walk_position_prob(  c1.apply_unitary_effect_to_state(2)  )
            , s
        ) )

    return poss

def plot_positions_per_step( poss ):
    ax = plt.subplot(111)
    im = ax.imshow( np.abs( poss ))
    plt.colorbar(im)
    plt.show()



# This function pretends that we 'double' the number of sites (compared to the actual number
# of sites availlable in the n-1 qubit register).
def simulate_2D_walk( nx, ny, num_steps ):
    coinbit = cirq.GridQubit(0,0)
    qx = [cirq.GridQubit(1,i) for i in range(nx)]
    qy = [cirq.GridQubit(2,i) for i in range(ny)]
    c1 = cirq.H
    c2 = cirq.H

    initstate = np.zeros(2**(nx+ny+1))
    initstate[0] = 1.
    amps = [initstate]
    probs = [ positions_compact_to_standard_2D( quantum_walk_position_prob_2D( initstate, nx, ny ), 0 ) ]

    for m in range(1,num_steps):
        walk = QuantumWalk2DFourier_Flipped( coinbit, qx, qy, c1, c2, m)
        amp = walk.apply_unitary_effect_to_state(0)
        amps.append( amp )
        probs.append( positions_compact_to_standard_2D( quantum_walk_position_prob_2D( amp, nx, ny ), m ) )

    return( probs )

def plot_2D_per_step( prob_matrices ):

    (num_steps, nx, ny) = np.shape( prob_matrices )

    # Attempt to make a square array of plots
    figwidth = int( np.ceil( np.sqrt( num_steps ) ) )
    fig, axes = plt.subplots(nrows=figwidth, ncols=figwidth )

    # In each plot, show the position matrix (with values in range (0,1))
    for j in range(0,num_steps):
        im = axes[j // figwidth, j % figwidth ].imshow( prob_matrices[j], vmin=0, vmax=1 )

    # Make space for the colorbar, taking the rightmost 10%
    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.91, 0.25, 0.03, 0.5])
    fig.colorbar(im, cax=cbar_ax)

    plt.show()


def example_2D_simulation():
    nx = 2
    ny = 2
    num_steps = 25
    posmats = np.sqrt( simulate_2D_walk( nx, ny, num_steps ) )
    plot_2D_per_step( posmats )

if __name__ == "__main__":
    ## Force pretty printing of numpy arrays
    np.set_printoptions(suppress=True, precision=3)
    np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

    example_2D_simulation()
