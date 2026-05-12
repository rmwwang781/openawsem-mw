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

def read_reference_structure_for_q_calculation_3(oa, pdb_file, reference_chain_name="ALL", min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers, Qflag=0, a=0.1, removeDNAchains=True):
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
    rnaResidues = ['A', 'G', 'C', 'U', 'I']
    dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']
    # out = open("ijpairs_q.dat", 'w')

    for chain in model.get_chains():
        chain_start += count
        count = 0
        if removeDNAchains and np.alltrue([a.get_resname().strip() in dnaResidues for a in chain.get_residues()]):
            print(f"chain {chain.id} is a DNA chain. will be ignored for Q evaluation")
            continue
        elif removeDNAchains and np.alltrue([a.get_resname().strip() not in proteinResidues for a in chain.get_residues()]):
            print(f"chain {chain.id} is a ligand chain. will be ignored for Q evaluation")
            continue
        # print(chain)
        for i, residue_i in enumerate(chain.get_residues()):
            #  print(i, residue_i)
            count +=1
            for j, residue_j in enumerate(chain.get_residues()):
                # out.write(f"{i},{j}")
                # out.write("\n")
                if abs(i-j) >= min_seq_sep and abs(i-j) <= max_seq_sep:  # taking the signed value to avoid double counting
                    ca_i = residue_i['CA']

                    ca_j = residue_j['CA']

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
    # print("Done reading")
    # print(structure_interactions)
    return structure_interactions

def read_reference_structure_for_qc_calculation(oa, pdb_file, min_seq_sep=3,contact_threshold=0.95*nanometers,  Qflag=1,a=0.1, startResidueIndex=0, endResidueIndex=-1, residueIndexGroup=None):
    # default use all chains in pdb file.
    # this change use the canonical Qw/Qo calculation for reference Q
    # for Qw calculation is 0; Qo is 1;
    structure_interactions = []
    parser = PDBParser()

    # out = open("ijpairs_qc.dat", 'w')

    structure = parser.get_structure('X', pdb_file)
    if endResidueIndex == -1:
        endResidueIndex = len(list(structure.get_residues()))
    if residueIndexGroup is None:
        # residueIndexGroup is used for non-continuous residues that used for Q computation.
        residueIndexGroup = list(range(len(list(structure.get_residues()))))
    for i, res_i in enumerate(structure.get_residues()):
        chain_i = res_i.get_parent().id
        if i < startResidueIndex:
            continue
        if i not in residueIndexGroup:
            continue
        for j, res_j in enumerate(structure.get_residues()):
            # out.write(f"{i},{j}")
            # out.write("\n")
            if j > endResidueIndex:
                continue
            if j not in residueIndexGroup:
                continue
            chain_j = res_j.get_parent().id
            if j-i >= min_seq_sep and chain_i == chain_j:
                    ca_i = res_i['CA']
                    ca_j = res_j['CA']

                    r_ijN = abs(ca_i - ca_j)/10.0  # convert to nm
                    if Qflag ==1 and r_ijN >= contact_threshold: continue
                    # if Qflag ==1 and r_ijN >= 0.95: continue
                    sigma_ij = a*(abs(i-j)**0.15)  # 0.1 nm = 1 A
                    gamma_ij = 1.0
                    i_index = oa.ca[i]
                    j_index = oa.ca[j]
                    structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                    structure_interactions.append(structure_interaction)
    return structure_interactions

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
            # if j-i >= min_seq_sep and chain_i != chain_j:
            if chain_i != chain_j:
                # out.write(f"{ci},{i},{cj},{j}")
                # out.write("\n")
                ca_i = res_i['CA']
                ca_j = res_j['CA']

                r_ijN = abs(ca_i - ca_j)/10.0*nanometers # convert to nm
                # if Qflag ==1 and r_ijN >= contact_threshold: continue
                if r_ijN >= contact_threshold: continue
                # print(res_i,res_j)

                sigma_ij = a*(abs(i-j)**0.15)  # 0.1 nm = 1 A
                gamma_ij = 1.0

                i_index = oa.ca[i+chain_start]
                j_index = oa.ca[j+chain_start]
                structure_interaction = [i_index, j_index, [gamma_ij, r_ijN, sigma_ij]]
                structure_interactions.append(structure_interaction)
    # out.close()
    return structure_interactions


def q_value(oa, reference_pdb_file, reference_chain_name="ALL", min_seq_sep=3, max_seq_sep=np.inf, contact_threshold=0.95*nanometers, forceGroup=1):
    ### Modified by Mingchen to compute canonical QW/QO

    # create bonds
    # structure_interactions = oa.read_reference_structure_for_q_calculation(reference_pdb_file, reference_chain_name, min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold)
    structure_interactions = read_reference_structure_for_q_calculation_3(oa, reference_pdb_file, reference_chain_name=reference_chain_name,
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


def qc_value(oa, reference_pdb_file, min_seq_sep=10, contact_threshold=0.95, a=0.2, forceGroup=3):
    # create bonds
    # structure_interactions = oa.read_reference_structure_for_q_calculation(reference_pdb_file, reference_chain_name, min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold)
    structure_interactions = read_reference_structure_for_qc_calculation(oa, reference_pdb_file, min_seq_sep=min_seq_sep,contact_threshold=contact_threshold, a=a)
    # print(len(structure_interactions))
    # print(structure_interactions)

    normalization = len(structure_interactions)
    qvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")
    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)
    qvalue.setForceGroup(forceGroup)
    return qvalue

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
    return qvalue

def partial_q_value(oa, reference_pdb_file, min_seq_sep=3, a=0.1, startResidueIndex=0, endResidueIndex=-1, residueIndexGroup=None, contact_threshold=0.95, forceGroup=4):
    print(f"Including partial q value computation, start residue index: {startResidueIndex}, end residue index: {endResidueIndex}, residueIndexGroup: {residueIndexGroup}")
    # create bonds
    # structure_interactions = oa.read_reference_structure_for_q_calculation(reference_pdb_file, reference_chain_name, min_seq_sep=min_seq_sep, max_seq_sep=max_seq_sep, contact_threshold=contact_threshold)
    structure_interactions = read_reference_structure_for_qc_calculation(oa, reference_pdb_file, min_seq_sep=min_seq_sep, a=a, startResidueIndex=startResidueIndex, endResidueIndex=endResidueIndex, contact_threshold=contact_threshold, residueIndexGroup=residueIndexGroup)
    # print(len(structure_interactions))
    # print(structure_interactions)
    # if len(structure_interactions) == 0:
    #     print("No atom found, Please check your startResidueIndex and endResidueIndex.")
    #     exit()
    normalization = len(structure_interactions)
    qvalue = CustomBondForce(f"(1/{normalization})*gamma_ij*exp(-(r-r_ijN)^2/(2*sigma_ij^2))")
    qvalue.addPerBondParameter("gamma_ij")
    qvalue.addPerBondParameter("r_ijN")
    qvalue.addPerBondParameter("sigma_ij")
    for structure_interaction in structure_interactions:
        qvalue.addBond(*structure_interaction)
    qvalue.setForceGroup(forceGroup)
    return qvalue

def qbias_term(oa,reference_pdb_file, q0, reference_chain_name="ALL", k_qbias=200*kilocalorie_per_mole, qbias_min_seq_sep=3, qbias_max_seq_sep=np.inf, qbias_contact_threshold=0.8*nanometers, forceGroup=4):
    k_qbias = k_qbias.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    # print(k_qbias)
    qbias = CustomCVForce(f"0.5*{k_qbias}*(q-{q0})^2")
    # qbias = CustomCVForce(f"0.5*{k_qbias}*(q-q0)^2")
    q = q_value(oa, reference_pdb_file, reference_chain_name, min_seq_sep=qbias_min_seq_sep, max_seq_sep=qbias_max_seq_sep, contact_threshold=qbias_contact_threshold)
    if oa.periodic:
        q.setUsesPeriodicBoundaryConditions(True)
    qbias.addCollectiveVariable("q", q)
    # qbias.addGlobalParameter("k_qbias", k_qbias)
    # qbias.addGlobalParameter("q0", q0)
    is_periodic=qbias.usesPeriodicBoundaryConditions()
    print("\nqbias_term is in PBC",is_periodic)
    print("\nqbias is",q0)
    qbias.setForceGroup(forceGroup)
    return qbias

def create_dist(oa,fileA="groupA.dat",fileB="groupB.dat",forceGroup=4):
    #groupA = list(range(68))
    #groupB = list(range(196, oa.natoms))
    #print(cA)
    #print(cB)
    groupA = np.array(np.loadtxt(fileA,dtype=int)).tolist()
    #print (groupA)
    #groupA = [0,1]
    groupB = np.array(np.loadtxt(fileB,dtype=int)).tolist()
    #groupA, groupB = get_contact_atoms('crystal_structure-openmmawsem.pdb', chainA=cA, chainB=cB)
    #pull_d = CustomCentroidBondForce(2, 'distance(g1,g2)-R0') # 为什么这里给了R0，实际distance变成了2倍？
    #pull_d.addGlobalParameter("R0", 0.0*angstroms)
    pull_d = CustomCentroidBondForce(2, 'distance(g1,g2)')
    pull_d.addGroup(groupA)
    pull_d.addGroup(groupB) # addGroup(groupB)
    pull_d.addBond([0, 1])
    pull_d.setForceGroup(forceGroup)
    return pull_d

def create_centroid_system(oa, fileA="groupA.dat",fileB="groupB.dat",k=100,R0=0,forceGroup=26):
    K_pull = k*4.184 * oa.k_awsem
    #R0 = R0*nm
    pull_force = CustomCVForce("0.5*K_pull*(d-R0)^2")
    d = create_dist(oa)#    cA = list(range(68)),
    pull_force.addCollectiveVariable("d", d)
    pull_force.addGlobalParameter("K_pull", K_pull)
    pull_force.addGlobalParameter("R0", R0)
    pull_force.setForceGroup(forceGroup)
    return pull_force

def create_dist_vector(oa,fileA="groupA.dat",fileB="groupB.dat",fileC="groupC.dat",forceGroup=4):
    #groupA = list(range(68))
    #groupB = list(range(196, oa.natoms))
    #print(cA)
    #print(cB)
    groupA = np.array(np.loadtxt(fileA,dtype=int)).tolist()
    groupB = np.array(np.loadtxt(fileB,dtype=int)).tolist()
    groupC = np.array(np.loadtxt(fileC,dtype=int)).tolist()
    print (groupA)
    print (groupB)
    print (groupC)
    #groupA, groupB = get_contact_atoms('crystal_structure-openmmawsem.pdb', chainA=cA, chainB=cB)
    #pull_d = CustomCentroidBondForce(2, 'distance(g1,g2)-R0') # 为什么这里给了R0，实际distance变成了2倍？
    #pull_d.addGlobalParameter("R0", 0.0*angstroms)
    pull_d = CustomCentroidBondForce(3, f"r1*cos(theta);\
                                r1=distance(p1,p2);\
                                theta=angle(p1, p2, p3);")

    #pull_d = CustomCompoundBondForce(3, f"r1*cos(theta);\
     #                           r1=distance(p1,p2);\
      #                          theta=angle(p1, p2, p3);")

    g1=pull_d.addGroup(groupA)
    g2=pull_d.addGroup(groupB) # addGroup(groupB)
    #pull_d.addBond([0,1,2])
    g3=pull_d.addGroup(groupC)
    #g4=pull_d.addGroup(groupA)
    #g5=pull_d.addGroup(groupB)


    #print (g1,g2,g3,g4)
    #pull_d.addBond([2800,3435,3780])
    pull_d.addBond([0,1,2],[])
    #pull_d.addBond([g0])
    #pull_d.addBond([1, 2])
    #pull_d.addBond([0, 3])
    pull_d.setForceGroup(forceGroup)
    return pull_d

def create_centroid_system2(oa, fileA="groupA.dat",fileB="groupB.dat",fileC="groupC.dat",k=100,R0=0,forceGroup=26):
    K_pull = k*4.184 * oa.k_awsem
    #R0 = R0*nm
    pull_force = CustomCVForce("0.5*K_pull*(d-R0)^2")
    d = create_dist_vector(oa,fileA=fileA,fileB=fileB,fileC=fileC)#    cA = list(range(68)),
    pull_force.addCollectiveVariable("d", d)
    pull_force.addGlobalParameter("K_pull", K_pull)
    pull_force.addGlobalParameter("R0", R0)
    pull_force.setForceGroup(forceGroup)
    return pull_force

def rg_term(oa, convertToAngstrom=True):
    rg_square = CustomBondForce("1/normalization*r^2")
    # rg = CustomBondForce("1")
    rg_square.addGlobalParameter("normalization", oa.nres*oa.nres)
    for i in range(oa.nres):
        for j in range(i+1, oa.nres):
            rg_square.addBond(oa.ca[i], oa.ca[j], [])
    if convertToAngstrom:
        unit = 10
    else:
        unit = 1
    rg = CustomCVForce(f"{unit}*rg_square^0.5")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(2)
    return rg

def rg_term_res18_to_end(oa, convertToAngstrom=True):
    start = 17  # Python index for residue 18 (0-based indexing)
    end = oa.nres  # last residue index (exclusive in range)

    rg_square = CustomBondForce("1/normalization*r^2")
    n_selected = end - start
    rg_square.addGlobalParameter("normalization", n_selected * n_selected)

    for i in range(start, end):
        for j in range(i + 1, end):
            rg_square.addBond(oa.ca[i], oa.ca[j], [])

    unit = 10 if convertToAngstrom else 1
    rg = CustomCVForce(f"{unit}*rg_square^0.5")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(4)

    return rg


def rg_bias_term(oa, k=1*kilocalorie_per_mole, rg0=0, atomGroup=-1, forceGroup=27):
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_rg = oa.k_awsem * k
    nres, ca = oa.nres, oa.ca
    if atomGroup == -1:
        group = list(range(nres))
    else:
        group = atomGroup     # atomGroup = [0, 1, 10, 12]  means include residue 1, 2, 11, 13.
    n = len(group)
    normalization = n*n
    rg_square = CustomBondForce(f"1.0/{normalization}*r^2")
    # rg = CustomBondForce("1")
    # rg_square.addGlobalParameter("normalization", n*n)
    for i in group:
        for j in group:
            if j <= i:
                continue
            rg_square.addBond(ca[i], ca[j], [])
    rg = CustomCVForce(f"{k_rg}*(rg_square^0.5-{rg0})^2")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(forceGroup)
    return rg

def cylindrical_rg_bias_term(oa, k=1*kilocalorie_per_mole, rg0=0, atomGroup=-1, forceGroup=27):
    k = k.value_in_unit(kilojoule_per_mole)   # convert to kilojoule_per_mole, openMM default uses kilojoule_per_mole as energy.
    k_rg = oa.k_awsem * k
    nres, ca = oa.nres, oa.ca
    if atomGroup == -1:
        group = list(range(nres))
    else:
        group = atomGroup          # atomGroup = [0, 1, 10, 12]  means include residue 1, 2, 11, 13.
    n = len(group)
    normalization = n * n
    rg_square = CustomCompoundBondForce(2, f"1/{normalization}*((x1-x2)^2+(y1-y2)^2)")

    for i in group:
        for j in group:
            if j <= i:
                continue
            rg_square.addBond([ca[i], ca[j]], [])

    rg = CustomCVForce(f"{k_rg}*(rg_square^0.5-{rg0})^2")
    rg.addCollectiveVariable("rg_square", rg_square)
    rg.setForceGroup(forceGroup)
    return rg

def pulling_term(oa, k_pulling=4.184, forceDirect="x", appliedToResidue=1, forceGroup=19):
    # k_m in units of nm^-1, z_m in units of nm.
    # z_m is half of membrane thickness
    # membrane_center is the membrane center plane shifted in z axis.
    # add membrane forces
    # 1 Kcal = 4.184 kJ strength by overall scaling
    k_pulling *= oa.k_awsem
    pulling = CustomExternalForce(f"(-{k_pulling})*({forceDirect})")
    for i in range(oa.natoms):
        if appliedToResidue == "LAST":
            appliedToResidue = oa.nres
        if appliedToResidue == "FIRST":
            appliedToResidue = 1
        if oa.resi[i] == (appliedToResidue-1):
            pulling.addParticle(i)
        # print(oa.resi[i] , oa.seq[oa.resi[i]])
    pulling.setForceGroup(forceGroup)
    return pulling

# Test Pull center of mass on Z 11082024 Rebekah
def pullZ_term(oa,k=100,R0=0,forceGroup=19):
    K_pull = k*4.184 * oa.k_awsem
    pull_Z = CustomCVForce("0.5*K_pull*(d-R0)^2") #harmonic potential
    pull_Z.addGlobalParameter("K_pull", K_pull)
    pull_Z.addGlobalParameter("R0", R0)

    com_z = CustomCentroidBondForce(1, "z1")
    com_z.addGroup([i for i in range(oa.natoms)])  # Use all atoms in the molecule for COM
    com_z.addBond([0], [])  # Adding the only bond between the "group" of atoms
    pull_Z.addCollectiveVariable("d", com_z)
    pull_Z.setForceGroup(forceGroup)
    return pull_Z

# pull center of mass on X 02142025 Rebekah ---Start
def pullX_term(oa,k=20,R0=0,forceGroup=19):
    K_pull = k*4.184 * oa.k_awsem
    pull_X = CustomCVForce("0.5*K_pull*(d-R0)^2") #harmonic potential
    pull_X.addGlobalParameter("K_pull", K_pull)
    pull_X.addGlobalParameter("R0", R0)

    com_x = CustomCentroidBondForce(1, "x1")
    com_x.addGroup([i for i in range(oa.natoms)])  # Use all atoms in the molecule for COM
    com_x.addBond([0], [])  # Adding the only bond between the "group" of atoms
    pull_X.addCollectiveVariable("d", com_x)
    pull_X.setForceGroup(forceGroup)
    return pull_X
# pull center of mass on X 02142025 Rebekah ---End

# Test only pull the COM of n17 on Z 02062025 Rebekah
def pullZn17_term(oa,k=100,R0=0,forceGroup=19):
    K_pull = k * 4.184 * oa.k_awsem
    pull_Z = CustomCVForce("0.5*K_pull*(d-R0)^2")
    pull_Z.addGlobalParameter("K_pull", K_pull)
    pull_Z.addGlobalParameter("R0", R0)

    n17_atoms = []
    for residue in oa.residues[:17]:  # Select first 17 residues
        n17_atoms.extend([atom.index for atom in residue.atoms()])
    # print(n17_atoms)
    com_z_n17 = CustomCentroidBondForce(1, "z1")
    com_z_n17.addGroup(n17_atoms) 
    com_z_n17.addBond([0], [])  
    pull_Z.addCollectiveVariable("d", com_z_n17)
    pull_Z.setForceGroup(forceGroup)

    return pull_Z
# Test only pull the COM of n17 on Z 02062025 Rebekah


def com_z_value(oa,forceGroup=3):
    com_z = CustomCentroidBondForce(1, "z1")
    com_z.addGroup([i for i in range(oa.natoms)])  
    com_z.addBond([0], [])  
    com_z.setForceGroup(forceGroup)
    return com_z
# Test Pull center of mass on Z 11082024 Rebekah


def com_x_value(oa,forceGroup=3):
# Compute center of mass on X 02142025 Rebekah
    com_x = CustomCentroidBondForce(1, "x1")
    com_x.addGroup([i for i in range(oa.natoms)])  
    com_x.addBond([0], [])  
    com_x.setForceGroup(forceGroup)
    return com_x

# Test COM of n17 on Z 01132025 Rebekah
def com_z_value_N17(oa,forceGroup=3):
    # atom_indices=list(range(0,101)) #the atom range of HTT N17 is 0 to 100 # should be 0 to 99 05132025
    # atom_indices=list(range(0,101))+list(range(220,320)) #this one is for NQ20 dimer
    atom_indices=list(range(0,100))+list(range(269,369)) #this one is for NQ20P10 dimer
    com_z_N17 = CustomCentroidBondForce(1, "z1")
    com_z_N17.addGroup(atom_indices) 
    com_z_N17.addBond([0], [])  
    com_z_N17.setForceGroup(forceGroup)
    return com_z_N17
# Test COM of n17 on Z 01132025 Rebekah

# Test COM of n17 on Z 02062025 Rebekah
def com_z_value_N17_2(oa,forceGroup=3):
    n17_atoms = []
    for residue in oa.residues[:17]:  # Select first 17 residues
        n17_atoms.extend([atom.index for atom in residue.atoms()])
    # atom_indices=list(range(0,101))
    com_z_N17 = CustomCentroidBondForce(1, "z1")
    com_z_N17.addGroup(n17_atoms) 
    com_z_N17.addBond([0], [])  
    com_z_N17.setForceGroup(forceGroup)
    return com_z_N17
# Test COM of n17 on Z 02062025 Rebekah

# Test pull two molecules apart 11222024 Rebekah
def pulling_two(oa,d0=0, k=20, groupA='groupA.dat', groupB='groupB.dat',forceGroup=19):
    k_pull = k*4.184 * oa.k_awsem
    print('pulling_two term bias is :',k_pull)
    groupA = np.loadtxt(groupA, dtype=int)
    groupB = np.loadtxt(groupB, dtype=int)
    # print("Max atom index in groupA:", max(groupA))
    # print("Max atom index in groupB:", max(groupB))
    # print("Total atoms in system:", oa.natoms)

    pulling_force = CustomCentroidBondForce(2, f"0.5*{k_pull}*(distance(g1,g2)-d0)^2")
    pulling_force.addGlobalParameter('d0', d0)
    pulling_force.addGroup(groupA)
    pulling_force.addGroup(groupB)
    pulling_force.addBond([0, 1])
    pulling_force.setForceGroup(forceGroup)

    return pulling_force

def distance_two(oa, groupA='groupA.dat', groupB='groupB.dat',forceGroup=4):
    groupA = np.loadtxt(groupA, dtype=int)
    groupB = np.loadtxt(groupB, dtype=int)

    distance_force = CustomCentroidBondForce(2, "distance(g1,g2)")
    distance_force.addGroup(groupA)
    distance_force.addGroup(groupB)
    distance_force.addBond([0, 1], [])

    # distance_force.addCollectiveVariable("d", distance_force)
    distance_force.setForceGroup(forceGroup)

    return distance_force

# 03282025 pin protein on top of vesicle
def pinXY_term(oa, k=100, x0=0, y0=0, forceGroup=19):
    K_pull = k * 4.184 * oa.k_awsem  # kcal/mol/nm^2 -> OpenMM units

    pull_XY = CustomCVForce("0.5*K_pull_x*(x_com - x0)^2 + 0.5*K_pull_y*(y_com - y0)^2")
    pull_XY.addGlobalParameter("K_pull_x", K_pull)
    pull_XY.addGlobalParameter("K_pull_y", K_pull)
    pull_XY.addGlobalParameter("x0", x0)
    pull_XY.addGlobalParameter("y0", y0)

    # COM in x-direction
    com_x = CustomCentroidBondForce(1, "x1")
    com_x.addGroup([i for i in range(oa.natoms)])
    com_x.addBond([0], [])
    pull_XY.addCollectiveVariable("x_com", com_x)

    # COM in y-direction
    com_y = CustomCentroidBondForce(1, "y1")
    com_y.addGroup([i for i in range(oa.natoms)])
    com_y.addBond([0], [])
    pull_XY.addCollectiveVariable("y_com", com_y)

    pull_XY.setForceGroup(forceGroup)
    return pull_XY

def radial_restraint(oa,k=5.0, R0=2.5):
# def radial_restraint_flat_bottom(atom_indices, k=5.0, R0=2.5):
    """
    Applies a flat-bottom cylindrical restraint along the z-axis.
    Particles are free to move within radius R0, but penalized if they go beyond it.

    Parameters:
        atom_indices (list of int): Atom indices to apply the restraint.
        k (float): Force constant in kcal/mol/nm^2.
        R0 (float): Cutoff radius in nm.

    Returns:
        CustomExternalForce: The restraint force object.
    """
    # Flat-bottom harmonic: no energy if r < R0; harmonic if r > R0
    k_radial = k * 4.184 * oa.k_awsem  # kcal/mol/nm^2 -> OpenMM units

    energy_expression = "step(sqrt(x^2 + y^2) - R0) * 0.5 * k_radial * (sqrt(x^2 + y^2) - R0)^2"
    force = CustomExternalForce(energy_expression)
    force.addGlobalParameter("k_radial", k_radial)
    force.addGlobalParameter("R0", R0)
    for i in range(oa.natoms):
        force.addParticle(i)

    return force