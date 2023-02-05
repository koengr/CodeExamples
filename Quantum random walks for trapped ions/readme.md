Most experimental implementations of quantum random walks use physical "walkers" that can be in superposition. From a computational perspective, that's silly: using a quantum computer, one can encode the position of a walker in Log(N) qubits, and have full flexibility to use any possible coin or graph. 

These scripts attempt to efficiently compile quantum walks for specific types of hardware, and show some intersting experiments (like measuring topological winding numbers) that could be performed on NISQ devices. 
