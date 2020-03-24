import math
import Bio

from Bio.PDB import PPBuilder
from Bio.PDB.Structure import Structure
from typing import List, Tuple

import matplotlib.pyplot as plt

# custom type

def phi_psi_angles(pdb_struct: Structure) -> Tuple[List[float], List[float]]:
    """
    Returns tortions angles (φ and ψ) for the proteins from Structure object[1],
    which represents PDB or MMCIF files in biopython package.

    Torsion angles are dihedral angles, which are defined by 4 points in space. 
    In proteins the two torsion angles φ and ψ describe the rotation of the 
    polypeptide chain around the two bonds on both sides of the Cα atom (N-Cα, C-Cα)[2].

    References
    ----------
    [1]: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc159
    [2]: https://proteinstructures.com/Structure/Structure/Ramachandran-plot.html

    Parameters
    ----------
    pdb_struct: Bio.PBD.Structure.Structure
        Structure object from biopython package       

    Returns
    -------
    phi_values: List[float]
        amino acids phi angles in degrees
    psi_values: List[float]
        amino acids psi angles in degrees

    The order of angles is maintained, so zip(phi_values, psi_values) gives phi, psi angle for each amino acid.
    """    
    phis, psis = [], []
    peptides = PPBuilder().build_peptides(pdb_struct)
    for peptide in peptides:
        for phi, psi in peptide.get_phi_psi_list():
            if phi and psi:
                phi, psi = map(math.degrees, (phi, psi))
                phis.append(phi)
                psis.append(psi)
    
    return phis, psis


def ram_plot(phi: List[float], psi: List[float], name: str):
    """
    Plots The Ramachandran plot from Phi & Psi values (in degrees).

    The Ramachandran plot provides a convenient way to view the distribution of torsion angles
    in a protein structure. It also provides an overview of excluded regions that show which 
    rotations of the polypeptide are not allowed due to steric hindrance (collisions between atoms)[1].

    References
    ----------
    [1]: https://proteinstructures.com/Structure/Structure/Ramachandran-plot.html

    Parameters
    ----------
    phi: List[float]
        array with values of amino acids phi angles (in degrees)
    psi: List[float]
        array with values of amino acids psi angles (in degrees)
    name: str
        PDB Id, which will be shown in plot title

    """
    plt.figure(figsize=(9, 9))
    plt.scatter(phi, psi)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.plot([0, 0], [-180, 180], color="black")
    plt.plot([-180, 180], [0, 0], color="black")
    plt.xlabel("$\phi$", size=14)
    plt.ylabel("$\psi$", size=14)
    plt.title(f"Ramachandran Plot. PDB ID: {name}", size=16)
    plt.grid()