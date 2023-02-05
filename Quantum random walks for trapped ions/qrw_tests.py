import cirq
from qrw import *



def test_qrw_functions():
    num_qubits = 4
    qubits = [cirq.GridQubit(i,0) for i in range(num_qubits)]

    #
    # print("QFT_Flipped_MS compared to QFT_Flipped (using CZ):")
    # u1 = QFT_Flipped_MS(qubits).to_unitary_matrix()
    # u2 = QFT_Flipped(qubits).to_unitary_matrix()
    # print( np.linalg.norm(u1-u2) )


    print("fractional_cphase versus MS_fractional_cphase:")
    exponents = [0, 1, 4, 2.5]
    for exp in exponents :
        u1 = MS_fractional_cphase( qubits[0:2], exp).to_unitary_matrix()
        u2 = fractional_cphase( qubits[0:2], exp).to_unitary_matrix()
        print( np.linalg.norm( u1-u2) )


    print("QFT_Flipped_MS_Xbasis versus QFT_Flipped:")
    c1 = Hadamards(qubits) + QFT_Flipped_MS_Xbasis(qubits) + Hadamards(qubits)
    u1 = c1.to_unitary_matrix()
    u2 = QFT_Flipped(qubits).to_unitary_matrix()
    print( np.linalg.norm(u1-u2) )



    print("ConditionalAddOne_Flipped_Xbasis, versus ConditionalAddOne_Flipped:")
    c1 = Hadamards(qubits) + ConditionalAddOne_Flipped_Xbasis(qubits) + Hadamards(qubits)
    u1 = c1.to_unitary_matrix()
    u2 = ConditionalAddOne_Flipped(qubits).to_unitary_matrix()
    print( np.linalg.norm( u1-u2 ))

    print("The MS-based 1D QW, versus the standard gate 1D QW:")
    num_stepss = [ 1, 4]
    for num_steps in num_stepss :
        u1 = QuantumWalk1DFourier( qubits, num_steps ).to_unitary_matrix()
        u2 = TrappedIon1DWalk( qubits, num_steps ).to_unitary_matrix()
        print( np.linalg.norm( u1-u2) )





def tex_1D_walk() :
        ## SETUP
        num_steps = 1
        n = 4
        qubits = [cirq.GridQubit(i,0) for i in range(n)]

        c1 = QuantumWalk1DFourier_Flipped( qubits, num_steps )

        tex(c1)

## check me! This function was not tested yet.
## Print the TeX circuit for the 1+2+2 lazy split-step walk on a square grid.
## The coin operators are assumed to be Hadamard gates.
def tex_2D_walk() :
        ## SETUP
        num_steps = 1
        qubits = [cirq.GridQubit(i,0) for i in range(5)]
        c1 = QuantumWalk2DFourier_Flipped( qubits[0], qubits[1:3], qubits[3:5], cirq.H, cirq.H, num_steps )
        tex(c1)



def tex_ms_CZ() :
    n=2
    qubits = [cirq.GridQubit(i,0) for i in range(n)]
    c = MS_controlledZ(qubits, phase=np.pi/8, prettyphase = False )
    tex(c)


def tex_MS_QFT() :
    n = 3
    qubits = [cirq.GridQubit(i,0) for i in range(n)]
    c = Hadamards(qubits) + QFT_Flipped_MS_Xbasis(qubits)
    tex(c)

def tex_MS_invQFT():
    n = 3
    qubits = [cirq.GridQubit(i,0) for i in range(n)]
    c = Hadamards(qubits) + QFT_Flipped_MS_Xbasis(qubits)
    tex(Conjugate(c))

def tex_MS_step() :
    n = 4
    qubits = [cirq.GridQubit(i,0) for i in range(n)]
    c = ConditionalAddOne_Flipped_Xbasis(qubits)
    tex(c)


def quirk_walk(num_qubits, num_steps):

    # Get the right string to be placed behind the Chance quirk object
    if num_qubits > 2 :
        numstring = str(num_qubits-1)
    else :
        numstring = ''

    # Single step string
    onestep = '["H","Chance' + numstring + '"],["•","inc' + str(num_qubits-1) + \
        '"],["◦","dec' + str(num_qubits-1) + '"]'

    # Concatenate single steps
    allsteps = ''
    for j in range(0,num_steps):
        allsteps += onestep
        if j<num_steps-1 : allsteps += ','
    else : allsteps += ',[1,"Chance' + numstring + '"]'
    print('https://algassert.com/quirk#circuit={"cols":[' + allsteps + ']}'  )




if __name__ == "__main__":

    quirk_walk(3,3)
    n = 3
    qubits = [cirq.GridQubit(i,0) for i in range(n)]
    test_qrw_functions()

    c = QFT(qubits[1:])
    c.append(cirq.X(qubits[0]))
    for i in range(4):
        c += ConditionalAddOne( qubits )

    c.append(cirq.X(qubits[0]))
    c += QFTinverse(qubits[1:])

    print(c)
    print(c.apply_unitary_effect_to_state(0))
