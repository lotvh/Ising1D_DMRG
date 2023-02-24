import numpy as np
import tenpy.linalg.np_conserved as npc
import matplotlib.pyplot as plt
from numpy import linalg as LA

from tenpy.networks.mpo import MPO
from tenpy.networks.mps import MPS
from tenpy.networks.site import SpinHalfSite
from tenpy.models.model import MPOModel
from tenpy.models.lattice import Chain
from tenpy.algorithms import dmrg
from tenpy.tools.params import asConfig


def get_Ising_MPO(L, J1, J2, g):
    r"""Returns MPO Hamiltonian for Transverse Field Ising Model With Nearest and
    Next-to-nearest Interactions

    The Hamiltonian reads:

    .. math ::

        H = J_1 \sum_{i=1}^{N-1} \sigma_i^x \sigma_{i+1}^x + J_2 \sum_{i=1}^{N-2}
        \sigma_i^x \sigma_{i+2}^x - g \sum_{i=1}^N \sigma_i^z

    All parameters are collected in a single dictionary `model_params`.

    Parameters
    ----------
    L : int
        Length of the chain.
    J1, J2, g: float
        Coupling parameters as defined for the Hamiltonian above.
    bc_MPS : {'finite' | 'infinte'}
        MPS boundary conditions.
    bc: {'periodic' | 'open'}
        Chain boundary conditions
    """
    # Define the sites and the lattice, which in this case is a simple uniform chain
    # of spin 1/2 sites
    site = SpinHalfSite(conserve=None)
    lat = Chain(L, site, bc_MPS="finite", bc="open")

    # Define the operators that appear in the Hamiltonian.
    sx, sz, Id = site.Sigmax, site.Sigmaz, site.Id

    # The grid (list of lists) that defines the MPO
    grid = [
        [Id, J2 * sx, J1 * sx, -g * sz],
        [None, None, Id, None],
        [None, None, None, sx],
        [None, None, None, Id],
    ]
    grids = [grid] * L

    # Generate the MPO Hamiltonian from the grid.
    H = MPO.from_grids(lat.mps_sites(), grids, bc="finite", IdL=0, IdR=-1)
    return H, lat


def run_dmrg(model_pams):
    r"""Runs the DMRG algorithm and returns the ground state energy, as well as the
        magnetization in the x- and z-direction

    Parameters
    ----------
    model_pams: dict
        dictionary of the model parameters L, J1, J2 and g
    """
    H, lat = get_Ising_MPO(
        model_pams["L"], model_pams["J1"], model_pams["J2"], model_pams["g"]
    )

    # Get the model from the MPO Hamiltonian and lattice
    model = MPOModel(lat, H)

    # Initialize the wavefunction from the model, using spin up for all lattice sites
    psi = MPS.from_product_state(
        model.lat.mps_sites(), ["up"] * model_pams["L"], "finite"
    )

    # Define the parameters to run the DMRG optimization algorithm
    dmrg_pams = {
        "mixer": True,
        "chi_list": {0: 100},
        "trunc_params": {"svd_min": 1.0e-10},
        "verbose": 1,
    }

    # Run the DMRG algorithm
    results = dmrg.run(psi, model, dmrg_pams)
    mag_x = np.sum(psi.expectation_value("Sigmax"))
    mag_z = np.sum(psi.expectation_value("Sigmaz"))

    return results["E"], mag_x, mag_z


def Ising_diag_4sites(J1, J2, g):
    r"""Finds the ground state of the 4-site TFIM by diagonalizing the Hamiltonian and
        finding the smallest eigenvalue.

    Parameters
    ----------
    J1, J2, g: float
        Coupling parameters as defined for the Hamiltonian in 'get_Ising_MPO'
    """
    s1, s3, iD = (
        np.matrix([[0, 1], [1, 0]]),
        np.matrix([[1, 0], [0, -1]]),
        np.identity(2),
    )

    H1 = (
        np.kron(s1, np.kron(s1, np.kron(iD, iD)))
        + np.kron(iD, np.kron(s1, np.kron(s1, iD)))
        + np.kron(iD, np.kron(iD, np.kron(s1, s1)))
    )
    H2 = np.kron(s1, np.kron(iD, np.kron(s1, iD))) + np.kron(
        iD, np.kron(s1, np.kron(iD, s1))
    )
    H3 = (
        np.kron(s3, np.kron(iD, np.kron(iD, iD)))
        + np.kron(iD, np.kron(s3, np.kron(iD, iD)))
        + np.kron(iD, np.kron(iD, np.kron(s3, iD)))
        + np.kron(iD, np.kron(iD, np.kron(iD, s3)))
    )

    w, v = LA.eig(J1 * H1 + J2 * H2 - g * H3)

    return min(w)


def plot_ground_energy(
    ground_energy_range, ground_energy_range_manual, g_range, filename
):
    plt.plot(g_range, ground_energy_range, ".", label="result DMRG")
    plt.plot(g_range, ground_energy_range_manual, label="manually computed result")
    plt.title("Ground State Energy of infinite 1D-TFIM as a Function of Ext Field")
    plt.ylabel("Energy")
    plt.ylim(min(ground_energy_range), max(ground_energy_range))
    plt.xlabel("Strength of Ext Field")
    plt.xlim(min(g_range), max(g_range))
    plt.legend(loc="upper right")
    plt.savefig(filename)


def plot_magnetization(magx, g_range, filename):
    plt.plot(g_range, magx)
    plt.title("Magnetization of 1D-TFIM as a Function of Ext Field")
    plt.ylabel("Magnetization in x-direction")
    plt.xlabel("Strength of Ext Field")
    plt.xlim(min(g_range), max(g_range))
    plt.ylim(min(magx), max(magx))
    plt.savefig(filename)


if __name__ == "__main__":

    g_range = np.linspace(0, 2.0, 21)
    ground_energy_range = []
    ground_energy_range_manual = []
    magx = []
    magz = []
    for g in g_range:
        # set the model parameters
        L, J1, J2 = 4, -1.0, -1.0
        model_pams = {"L": L, "J1": J1, "J2": J2, "g": g}

        # run the DMRG optimization algorithm
        ground_energy, mag_x, mag_z = run_dmrg(model_pams)

        # Find the ground state energy to compare to the DMRG result
        ground_energy_manual = Ising_diag_4sites(J1, J2, g)

        # Save the results
        ground_energy_range.append(ground_energy)
        ground_energy_range_manual.append(ground_energy_manual)
        magx.append(mag_x)
        magz.append(mag_z)

    # plot the ground energy as a function of the external field
    filename = "GroundEvsField.pdf"
    plot_ground_energy(
        ground_energy_range, ground_energy_range_manual, g_range, filename
    )

    # plot the magnetization in the x-direction as a function of the external field
    filename = "MagnetizationXvsField.pdf"
    plot_magnetization(magx, g_range, filename)
