# Ramachandran Plots with python

Simple module to plot Ramachandran plots with python. It can extract torsion angles from pdb/cif files, plot Ram plots,
or fetch and save pdb files from PDB by ID. The rest of the functionality is described in `torsion.py` docstrings. It's simple!

- `torsion.py` - main module with all needed functions
- `pdb_plots.ipynb` - example plots & description of Ramplot

# Quick example

```python
from torsion import fetch_pdb, ram_plot
from torsion import phi_psi_angles

# load file from PDB and return as Structure object from biopython
pdb_struct = fetch_pdb("6PCY")  

# extract torsion angles
phi_angles, psi_angles = phi_psi_angles(pdb_struct)

# plot!
ram_plot(phi_angles, psi_angles, density=True, secondary=True)
```
![ramplot](figs/example.png)


# Limitations

In the current state it is possible to add contour on a plot with gaussian kde `contour=True`, but this is not the same as "core/not-core" or “core/allowed“generously allowed” regions.