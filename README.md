# Simulating long-range 1-dimensional Ising model with tensor networks

A script that computes the ground state of the quantum Ising model one a transverse field and next-to-nearest neighbour interactions, using DMRG.

## Hamiltonian

The transverse field Ising Hamiltonian with nearest and next-to-nearest neighbour interactions is given by:

$$
H = J_1 \sum_{i=1}^{N-1} \sigma_i^x \sigma_{i+1}^x + J_2 \sum_{i=1}^{N-2} \sigma_i^x \sigma_{i+2}^x - g \sum_{i=1}^N \sigma_i^z
$$

where $J_1$ and $J_2$ describe how strong the nearest neighbour and next-to-nearest neigbour interactions are, respectively. The term $g$ describes the strength of the external magnetic field. Note that the above Hamiltonian approaches the nearest-neighbour Ising model for $J_2 \rightarrow 0$.

## Matrix Product Operator

In order to find the ground state of a system, the Hamiltonian can be described by a matrix product operator (MPO). The Hamiltonian defined above can be rewritten as a function of terms defined on the left and right side of the $k$-th term:

$$
\hat{H} = (I_{1..k-1}) \otimes J_1 \sigma_k^x \otimes (\sigma_{k+1}^x \otimes I_{k+2..N}) + (I_{1..k-2} \otimes J_1 \sigma_{k-1}^x) \otimes \sigma_k^x \otimes (I_{k+1..N}) +
$$

$$
( \sum_{i=1}^{k-2} I_{1..i-1} \otimes J_1 \sigma_i^x \otimes \sigma_{i+1}^x \otimes I_{i+2..k-1}) \otimes I_k \otimes (I_{k+1..N}) +
$$

$$
(I_{1..k-1}) \otimes I_k \otimes (\sum_{i=k+1}^{N-1} I_{k+1..i-1} \otimes J_1 \sigma_i^x \otimes \sigma_{i+1}^x \otimes I_{i+2..N}) +
$$

$$
(I_{1..k-1}) \otimes J_2 \sigma_k^x \otimes (I_{k+1} \otimes \sigma_{k+2}^x \otimes I_{k+3..N} ) + (I_{1..k-2} \otimes J_2 \sigma_{k-1}^x) \otimes I_k \otimes (\sigma_{k+1}^x \otimes I_{k+2..N}) +
$$

$$
(I_{1..k-3} \otimes J_2 \sigma_{k-2} \otimes I_{k-1}) \otimes \sigma_k^x \otimes (I_{k+1..N}) +
$$

$$
( \sum_{i=1}^{k-3} I_{1..i-1} \otimes J_2 \sigma_i^x \otimes I \otimes \sigma_{i+2}^x \otimes I_{i+3..k-1}) \otimes I_k \otimes (I_{k+1..N}) +
$$

$$
(I_{1..k-1}) \otimes I_k \otimes (\sum_{i=k+1}^{N-1} I_{k+1..i-1} \otimes J_2 \sigma_i^x \otimes I \otimes \sigma_{i+2}^x \otimes I_{i+3..N}) +
$$

$$
...
$$

where the $g$-terms were not written for visibility and $1_{i..j} = \otimes_{l=i}^k I_l$.
This can be brought into the form of an MPO by introducing matrices $\hat{W}_{[N]}$ whose product reproduces the expression above:

$$
\hat{H} = \hat{W}_{[1]} \otimes \hat{W}_{[2]} \otimes \dots \otimes \hat{W}_{[N-1]} \otimes \hat{W}_{[N]},
$$

where each tensor acts on one site in the chain.

The transverse field Ising Hamiltonian with nearest and next-to-nearest neighbour interactions can now be rewritten as:

$$
\hat{W}_{[l]} =
\begin{pmatrix}
I & \sigma_x & 0 & -g \sigma_z\\
0 & 0 & I & J_1 \sigma_x \\
0 & 0 & 0 & J_2 \sigma_x \\
0 & 0 & 0 & I
\end{pmatrix}, \quad
\hat{W}_{[1]} =
\begin{pmatrix}
I & \sigma_x & 0 & 0
\end{pmatrix}, \quad
\hat{W}_{[N]} =
\begin{pmatrix}
-h \sigma_z \\
J_1 \sigma_x \\
J_2 \sigma_x \\
I
\end{pmatrix}.
$$

The matrix $\hat{W}_{[l]}$ is implemented in the function `get_Ising_MPO(L, J1, J2, g)` which takes the length of the chain `L`, the coupling terms `J1` and `J2` and the field strength `g` and returns the Hamiltonian in MPO form, so that it can be used in the DMRG algorithm.

## DMRG Algorithm

The Density Matrix Renormalization Group (DMRG) algorithm is implemented in the function `run_dmrg(model_pams)` which takes the parameters of the model as the input and returns the ground state energy of the system as well as the magnetization in the $x$- and $z$-direction.

First, the wavefunction is initialized from the model where we assume that all spins in the chain are up. The algorithm variationally optimizes neighbouring tensors while minimizing the ground state energy until the wavefunction converges to the ground state.

## Examples

This algorithm can be used to find the ground state energies and magnetic susceptibility of the long-range Ising model given by the Hamiltonian above.

### Ground state energy as a function of the external field

The ground state energy was plotted as a function of the external field for both the finite and infinite version of the model.

#### Finite Case

For this case, a chain of length 4 was chosen. In order to test the accuracy of the DMRG algorithm, also a manually calculated result is shown in the same graph. This result was calculated by diagonalizing the Hamiltonian using the `eig` method which is available in numpy and is implemented in function `Ising_diag_4sites`. Open boundary conditions were used for the lattice.

![energy_next](https://user-images.githubusercontent.com/49079733/221142987-7e3539e1-f6ef-4906-9b10-e3d2c25b7018.png)

#### Infinite Case

For this case, periodic boundary conditions were used for the lattice. Due to the computational cost, less values were calculated than in the limited case. The ground state energy is considerably higher than in the finite case, but the decay w.r.t. the external field is preserved.

![GroundEvsFieldInfinite](https://user-images.githubusercontent.com/49079733/221143139-beac395a-d696-48a6-b9f4-63ac7d62a347.png)

### Magnetization as a function of the external field

The magnetization in the x-direction of the Ising model can be found as $M = < \sum_i s_i^x >$.

![Magnetization](https://user-images.githubusercontent.com/49079733/221143260-7ee66849-94f5-4c21-85d2-f8c86d9db1a8.png)
