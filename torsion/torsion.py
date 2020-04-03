import os
import requests
import matplotlib

import numpy as np
import matplotlib.pyplot as plt

from math import degrees
from typing import List, Tuple
from scipy.stats import gaussian_kde
from .utils import density_estimation

from Bio.PDB import PPBuilder
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import MMCIFParser


def load_pdb(pdb_path: str) -> Structure:
    """
    Load PDB file in .cif format. 

    Parameters
    ----------
    pdb_path: str
        path to .cif file. 

    Returns
    -------
    pdb_struct: Bio.PBD.Structure.Structure
        Structure object from biopython package for specified PDB .cif file.
    
    """
    pdb_struct = MMCIFParser().get_structure("tmp", pdb_path)
    
    return pdb_struct


def fetch_pdb(pdb_id: str, remove: bool=True) -> Structure:
    """
    Downloads a protein structure in .cif format from Protein Data Bank 
    by PDB ID and returns its as a Structure object from biopython package. 
    The downloaded file will be saved in current directory and then deleted
    if remove is True. PDB ID — the 4-character unique identifier of 
    every entry in the Protein Data Bank[1]. All IDs can be found on PDB site[2].

    References
    ----------
    [1]: https://www.rcsb.org/pages/help/advancedsearch/pdbIDs
    [2]: https://www.rcsb.org/pages/general/summaries

    Parameters
    ----------
    pdb_id: str
        PDB ID associated with the needed structure.
    remove: bool, default=True
        Remove or save downloaded file.

    Returns
    -------
    pdb_struct: Bio.PBD.Structure.Structure
        Structure object from biopython package for specified PDB ID.

    """
    if len(pdb_id) != 4:
        raise ValueError("Wrong PDB ID format. PDB ID should have 4-character length.")
    
    filename = f"{pdb_id}.cif"
    pdb = requests.get(f"http://files.rcsb.org/download/{filename}")
    if pdb.status_code != 200:
        raise ConnectionError(f"Retrieval failed. Status: {pdb.status_code}")

    with open(filename, "wb") as file:
        file.write(pdb.content) 

    pdb_struct = MMCIFParser().get_structure(pdb_id, filename)
    if remove: os.remove(filename)

    return pdb_struct    


def phi_psi_angles(pdb_struct: Structure) -> Tuple[List[float], List[float]]:
    """
    Returns torsions angles (φ and ψ) for the proteins from Structure object[1],
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
        Structure object from biopython package.

    Returns
    -------
    phi_values: List[float]
        amino acids phi angles in degrees.
    psi_values: List[float]
        amino acids psi angles in degrees.

    The order of angles is maintained, so zip(phi_values, psi_values) gives phi, psi angle for each amino acid.
    """
    if not isinstance(pdb_struct, Structure):
        raise TypeError("Incorrect type, use Bio.PDB.Structure.Structure from biopython.")

    phis, psis = [], []
    peptides = PPBuilder().build_peptides(pdb_struct)
    for peptide in peptides:
        for phi, psi in peptide.get_phi_psi_list():
            if phi and psi:
                phis.append(degrees(phi))
                psis.append(degrees(psi))
    
    return phis, psis


def ram_plot(phi: List[float], psi: List[float], pdb_id: str, 
            density: bool=False, contour: bool=False, secondary: bool=False, 
            dpi: int=128, figsize: Tuple[int, int]=(5, 5)):
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
        array with values of amino acids phi angles (in degrees).
    psi: List[float]
        array with values of amino acids psi angles (in degrees).
    pdb_id: str
        PDB Id, which will be shown in plot title.
    density: bool, default=False
        density estimation using Gaussian kernels.
    contour: bool, defualt=False
        draw contour lines using Gaussian kernels.
    secondary: bool, defualt=False
        plot labels for common secondary structures (a-helices, b-sheets, la-helices).
    figsize: int, defualt=(5, 5)
        plot size, same as in matplotlib (actually, exactly the same).
    dpi: int, defualt=128
        plot dpi.
    """
    matplotlib.rcParams['axes.linewidth'] = 0.4

    plt.figure(figsize=figsize, dpi=dpi)
    plt.xlabel("$\phi$", size=6)
    plt.ylabel("$\psi$", size=6)
    plt.title(f"Ramachandran Plot. PDB ID: {pdb_id}", size=10, pad=10)

    # angles for common secondary structures
    if secondary:
        plt.text(-57, -47, r"$\alpha$", fontsize=20, zorder=12)
        plt.text(-119, 113, r"$\beta$", fontsize=18, zorder=12)
        plt.text(57, 47, r"$L\alpha$", fontsize=18, zorder=12)
        plt.text(75, -65, r"$\gamma$", fontsize=18, zorder=12)

    plt.axvline(0, c="black", linewidth=0.4, zorder=4)
    plt.axhline(0, c="black", linewidth=0.4, zorder=4)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)

    if density:
        xy = np.vstack([phi, psi])
        z = gaussian_kde(xy)(xy)

        plt.scatter(phi, psi, c=z, s=10, zorder=6)
    else:
        plt.scatter(phi, psi, s=10, zorder=6)

    if contour:
        x, y, z = density_estimation(np.array(phi), np.array(psi))

        plt.contour(x, y, z, colors="black", linewidths=0.4, zorder=8, levels=12)
        # well, if you really want, you can uncomment this line
        # plt.imshow(np.rot90(z), cmap="Greens", extent=[-180, 180, -180, 180], alpha=0.6)

    plt.grid(b=None, which='major', axis='both', color='#d1d1d1', alpha=0.5)
    plt.xticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], fontsize=6)
    plt.yticks([-180, -135, -90, -45, 0, 45, 90, 135, 180], fontsize=6)
