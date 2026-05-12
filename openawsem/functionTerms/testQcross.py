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

def read_reference_structure_for_qinter_calculation(oa, pdb_file, reference_chain_name="ALL",min_seq_sep=3, max_seq_sep=np.inf,contact_threshold=0.95*nanometers, Qflag=0,a=0.1):

    structure_interactions = []
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', pdb_file)
    model = structure[0]
    chain_start = 0
    # out = open("ijpairs.dat", 'w')

    for i, res_i in enumerate(structure.get_residues()):
        chain_i = res_i.get_parent().id
        for j, res_j in enumerate(structure.get_residues()):
            chain_j = res_j.get_parent().id
            if j-i >= min_seq_sep and chain_i != chain_j:
                # out.write(f"{ci},{i},{cj},{j}")
                # out.write("\n")
                ca_i = res_i['CA']
                ca_j = res_j['CA']

                r_ijN = abs(ca_i - ca_j)/10.0*nanometers # convert to nm
                if Qflag ==1 and r_ijN >= contact_threshold: continue
                sigma_ij = a*(abs(i-j)**0.15)  # 0.1 nm = 1 A
                gamma_ij = 1.0

                i_index = oa.ca[i+chain_start]
                j_index = oa.ca[j+chain_start]
                structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                structure_interactions.append(structure_interaction)
    # out.close()
    return structure_interactions


def q_inter_value(oa, reference_pdb_file, min_seq_sep=3, a=0.1, forceGroup=4):
    structure_interactions = read_reference_structure_for_qinter_calculation(oa, reference_pdb_file, min_seq_sep=min_seq_sep, a=a)
    normalization = len(structure_interactions)
    qvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")
    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)
    qvalue.setForceGroup(forceGroup)
    return qvalues