import math
import Bio

from Bio.PDB import PPBuilder
from typing import List, Tuple

import matplotlib.pyplot as plt

# custom types
AngleList = List[float]
PhiPsiTuple = Tuple[AngleList, AngleList]


def phi_psi_angles(pdb_struct: str) -> PhiPsiTuple:
    """
    Set docstring here.

    Parameters
    ----------
    pdb_struct:str: 

    Returns
    -------

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


def ram_plot(phi: AngleList, psi: AngleList, name: str):
    """
    Set docstring here.

    Parameters
    ----------
    phi:AngleList: 
    psi:AngleList: 
    name:str: 

    Returns
    -------

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