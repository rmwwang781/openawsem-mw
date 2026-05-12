try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except ModuleNotFoundError:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *

import numpy as np
from Bio.PDB.PDBParser import PDBParser

def pullAway_term(oa, k=100, R0=0, forceGroup=19, center=(0,0,1.3)):
    K_pull = k * 4.184 * oa.k_awsem  
    
    pull_Away = CustomCVForce("0.5 * K_pull * (d - R0)^2")  # Harmonic potential
    pull_Away.addGlobalParameter("K_pull", K_pull)
    pull_Away.addGlobalParameter("R0", R0)

    com_dist = CustomCentroidBondForce(1, "sqrt((x1 - x0)^2 + (y1 - y0)^2 + (z1 - z0)^2)")
    
    com_dist.addGlobalParameter("x0", center[0])
    com_dist.addGlobalParameter("y0", center[1])
    com_dist.addGlobalParameter("z0", center[2])

    # Add all atoms to compute COM
    com_dist.addGroup([i for i in range(oa.natoms)])
    com_dist.addBond([0], [])  

    # Attach the collective variable
    pull_Away.addCollectiveVariable("d", com_dist)
    pull_Away.setForceGroup(forceGroup)

    return pull_Away

def pullAwayNt_term(oa, k=100, R0=0, forceGroup=19, center=(0, 0, 0)):
    K_pull = k * 4.184 * oa.k_awsem  

    pull_Away = CustomCVForce("0.5 * K_pull * (d - R0)^2")
    pull_Away.addGlobalParameter("K_pull", K_pull)
    pull_Away.addGlobalParameter("R0", R0)

    com_dist = CustomCentroidBondForce(1, "sqrt((x1 - x0)^2 + (y1 - y0)^2 + (z1 - z0)^2)")
    com_dist.addGlobalParameter("x0", center[0])
    com_dist.addGlobalParameter("y0", center[1])
    com_dist.addGlobalParameter("z0", center[2])

    # Get atom indices of the first 17 residues
    n17_atom_indices = []
    for res in oa.residues[:17]:
        for atom in res.atoms():
            n17_atom_indices.append(atom.index)

    com_dist.addGroup(n17_atom_indices)
    com_dist.addBond([0], [])

    pull_Away.addCollectiveVariable("d", com_dist)
    pull_Away.setForceGroup(forceGroup)

    return pull_Away


def compute_COM_distance(oa, center=(0,0,1.3), forceGroup=3):

    com_x = CustomCentroidBondForce(1, "x1")
    com_x.addGroup([i for i in range(oa.natoms)])
    com_x.addBond([0], [])

    com_y = CustomCentroidBondForce(1, "y1")
    com_y.addGroup([i for i in range(oa.natoms)])
    com_y.addBond([0], [])

    com_z = CustomCentroidBondForce(1, "z1")
    com_z.addGroup([i for i in range(oa.natoms)])
    com_z.addBond([0], [])

    com_distance = CustomCVForce("sqrt((x1 - x0)^2 + (y1 - y0)^2 + (z1 - z0)^2)")

    com_distance.addGlobalParameter("x0", center[0])
    com_distance.addGlobalParameter("y0", center[1])
    com_distance.addGlobalParameter("z0", center[2])

    com_distance.addCollectiveVariable("x1", com_x)
    com_distance.addCollectiveVariable("y1", com_y)
    com_distance.addCollectiveVariable("z1", com_z)

    com_distance.setForceGroup(forceGroup)

    return com_distance

def CA9_Center_dist(oa, center=(0, 0, 1.3), forceGroup=4):
# def compute_CA9_COM_distance(oa, center=(0, 0, 1.3), forceGroup=3):
    # Get CA atom index from residue 9 (index 8 if 0-based)
    residue9 = oa.residues[8]  # Python is 0-indexed, so index 8 is residue 9
    ca9_index = None
    for atom in residue9.atoms():
        if atom.name == 'CA':
            ca9_index = atom.index
            break
    if ca9_index is None:
        raise ValueError("Cα atom not found in residue 9.")

    # Create centroid forces for x, y, and z
    com_x = CustomCentroidBondForce(1, "x1")
    com_x.addGroup([ca9_index])
    com_x.addBond([0], [])

    com_y = CustomCentroidBondForce(1, "y1")
    com_y.addGroup([ca9_index])
    com_y.addBond([0], [])

    com_z = CustomCentroidBondForce(1, "z1")
    com_z.addGroup([ca9_index])
    com_z.addBond([0], [])

    # Define the CV force to measure distance to vesicle center
    com_distance = CustomCVForce("sqrt((x1 - x0)^2 + (y1 - y0)^2 + (z1 - z0)^2)")
    com_distance.addGlobalParameter("x0", center[0])
    com_distance.addGlobalParameter("y0", center[1])
    com_distance.addGlobalParameter("z0", center[2])

    com_distance.addCollectiveVariable("x1", com_x)
    com_distance.addCollectiveVariable("y1", com_y)
    com_distance.addCollectiveVariable("z1", com_z)

    com_distance.setForceGroup(forceGroup)

    return com_distance

def compute_ntCOM_distance(oa, center=(0, 0, 0), forceGroup=3):
    # Get atom indices of the first 17 residues
    n17_atom_indices = []
    for res in oa.residues[:17]:
        for atom in res.atoms():
            n17_atom_indices.append(atom.index)

    # Create COM components for the first 17 residues only
    com_x = CustomCentroidBondForce(1, "x1")
    com_x.addGroup(n17_atom_indices)
    com_x.addBond([0], [])

    com_y = CustomCentroidBondForce(1, "y1")
    com_y.addGroup(n17_atom_indices)
    com_y.addBond([0], [])

    com_z = CustomCentroidBondForce(1, "z1")
    com_z.addGroup(n17_atom_indices)
    com_z.addBond([0], [])

    com_distance = CustomCVForce("sqrt((x1 - x0)^2 + (y1 - y0)^2 + (z1 - z0)^2)")
    com_distance.addGlobalParameter("x0", center[0])
    com_distance.addGlobalParameter("y0", center[1])
    com_distance.addGlobalParameter("z0", center[2])

    com_distance.addCollectiveVariable("x1", com_x)
    com_distance.addCollectiveVariable("y1", com_y)
    com_distance.addCollectiveVariable("z1", com_z)

    com_distance.setForceGroup(forceGroup)

    return com_distance

