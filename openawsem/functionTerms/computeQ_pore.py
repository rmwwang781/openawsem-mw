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

def read_reference_structure_for_q_calculation_pore(oa, pdb_file, reference_chain_name="ALL", min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers, Qflag=0, a=0.1, removeDNAchains=True):
    # default use all chains in pdb file.
    # this change use the canonical Qw/Qo calculation for reference Q
    # for Qw calculation is 0; Qo is 1;
    structure_interactions = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', pdb_file)
    model = structure[0]
    chain_start = 0
    count = 0
    proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR', 'ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
    proteinResidues += ["NGP", "IGL", "IPR"]

    for chain in model.get_chains():
        chain_start += count
        count = 0
        for i, residue_i in enumerate(chain.get_residues()):
            count +=1
            for j, residue_j in enumerate(chain.get_residues()):
                if abs(i-j) >= min_seq_sep and abs(i-j) <= max_seq_sep:  # taking the signed value to avoid double counting
                    ca_i = residue_i['CA']
                    ca_j = residue_j['CA']

                    # Get Z-coordinates (index 2 of the coordinate array)
                    z_i = ca_i.get_coord()[2]
                    z_j = ca_j.get_coord()[2]

                    # Only proceed if both residues are within the membrane boundaries
                    if not (-15 < z_i < 15 and -15 < z_j < 15):
                        continue

                    r_ijN = abs(ca_i - ca_j)/10.0*nanometers # convert to nm
                    if Qflag ==1 and r_ijN >= contact_threshold: continue
                    sigma_ij = a*(abs(i-j)**0.15)  # 0.1 nm = 1 A
                    gamma_ij = 1.0

                    if reference_chain_name != "ALL" and (chain.id not in reference_chain_name):
                        continue
                    i_index = oa.ca[i+chain_start]
                    j_index = oa.ca[j+chain_start]
                    structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                    structure_interactions.append(structure_interaction)
    return structure_interactions

def q_value_pore(oa, reference_pdb_file, reference_chain_name="ALL", min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers, forceGroup=3):
    ### Modified by Mingchen to compute canonical QW/QO

    # create bonds
    # structure_interactions = oa.read_reference_structure_for_q_calculation(reference_pdb_file, reference_chain_name, min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold)
    structure_interactions = read_reference_structure_for_q_calculation_pore(oa, reference_pdb_file, reference_chain_name=reference_chain_name,
        min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold, Qflag=0)
    # print(len(structure_interactions))
    # print(structure_interactions)
    # create bond force for q calculation
    normalization = len(structure_interactions)
    # print(normalization)
    qvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    if oa.periodic:
        qvalue.setUsesPeriodicBoundaryConditions(True)
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")

    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)

    is_periodic=qvalue.usesPeriodicBoundaryConditions()
    print("\nqvalue is in PBC",is_periodic)
    qvalue.setForceGroup(forceGroup)
    return qvalue