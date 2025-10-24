# Written by Mohd Ibrahim
# Technical University of Munich
# Email: ibrahim.mohd@tum.de

## modified by Jose Cediel-Becerra

import numpy as np
import MDAnalysis as mda
import argparse
import warnings

# Suppress specific warnings from MDAnalysis
warnings.filterwarnings("ignore")#, category=UserWarning, module="MDAnalysis.coordinates.PDB")

def create_universe (n_atoms, name, resname, positions, resids, segid):

    u_new = mda.Universe.empty(n_atoms=n_atoms,
                             n_residues=n_atoms,
                             atom_resindex=np.arange (n_atoms),
                             residue_segindex=np.arange (n_atoms),
                             n_segments=n_atoms,
                             trajectory=True) # necessary for adding coordinate


    u_new.add_TopologyAttr('name',   name)
    u_new.add_TopologyAttr('resid', resids)
    u_new.add_TopologyAttr('resname', resname)
    u_new.atoms.positions = positions
    u_new.add_TopologyAttr('segid', n_atoms*[segid])
    u_new.add_TopologyAttr('chainID', n_atoms * [segid])


    return u_new
    
def get_nme_pos (end_residue):

    if "OXT" in end_residue.names:
        index = np.where (end_residue.names == "OXT")[0][0]
        N_position = end_residue.positions [index]
        index_c = np.where (end_residue.names == "C")[0][0]
        carbon_position = end_residue.positions [index_c]
        vector = N_position - carbon_position
        vector /= np.sqrt (sum(vector**2))
        
        C_position = N_position + vector*1.36

        return N_position, C_position

    else:
        # find midpoint of O and CA
        index_o = np.where (end_residue.names == "O")[0][0]
        index_ca = np.where (end_residue.names == "CA")[0][0]

        mid_point = (end_residue.positions [index_o] + end_residue.positions [index_ca] )/2

        # find vector connecting mid_point and C
        index_c = np.where (end_residue.names == "C")[0][0]
        vector  = end_residue.positions [index_c] - mid_point
        vector /= np.sqrt (sum(vector**2))
        N_position = end_residue.positions [index_c] + 1.36* vector
        ##
        C_position = N_position + 1.36*vector
    
    return N_position, C_position 

def get_ace_pos (end_residue):
    
    index_ca = np.where (end_residue.names == "CA")[0][0]
    index_n  = np.where (end_residue.names == "N")[0][0]
    vector   = end_residue.positions [index_n] - end_residue.positions [index_ca] 
    vector  /= np.sqrt (sum(vector**2))

    C1_position = end_residue.positions [index_n] + 1.36*vector

    xa, ya, za =  end_residue.positions [index_ca] 
    xg, yg, zg = C1_position

    # arbritray unit vector
    # create an arbritray orientaiton for the ACE residue
    # does not really matter
    orientation  = np.array([2*np.random.rand () -1, 2*np.random.rand () -1,2*np.random.rand () -1])
    nx, ny, nz =  orientation/np.sqrt (sum(orientation**2))

    ## The carbon and oxygen are placed on the vertices of an equilatrel triangle
    # with another vertex as the Nitrogen atom and the C as the centroid
    # The plane of the triangle is placed in an arbritrary orientation as defined before
    # The orientation does not matter
    ######################################
    x1 = xg - (xa-xg)/2 + np.sqrt (3)*(ny*(za-zg) - nz*(ya-yg))/2
    y1 = yg - (ya-yg)/2 + np.sqrt (3)*(nz*(xa-xg) - nx*(za-zg))/2
    z1 = zg - (za-zg)/2 + np.sqrt (3)*(nx*(ya-yg) - ny*(xa-xg))/2
    
    ## second coordinate
    x2 = xg - (xa-xg)/2 - np.sqrt (3)*(ny*(za-zg) - nz*(ya-yg))/2
    y2 = yg - (ya-yg)/2 - np.sqrt (3)*(nz*(xa-xg) - nx*(za-zg))/2
    z2 = zg - (za-zg)/2 - np.sqrt (3)*(nx*(ya-yg) - ny*(xa-xg))/2

    C2_position = np.array ([x1,y1,z1])
    O_position = np.array ([x2,y2,z2])

    ### rescale distances, the above points may be a bit far apart like 2.1 angstrom but usual bonds are 1.4 or so
    ## Therefore we shrink it
    #  C positinos
    
    vector = C2_position - C1_position
    vector /= np.sqrt (sum (vector**2))
    
    C2_position = C1_position + 1.36*vector

    # O positions
    vector = O_position - C1_position
    vector /= np.sqrt (sum (vector**2))
    
    O_position = C1_position + 1.36*vector
    
    return C1_position, C2_position, O_position


def add_TER_between_chains(pdb_in, pdb_out):
    """Insert 'TER' between NME and ACE residues in multi-chain capped PDBs."""
    with open(pdb_in, "r") as f:
        lines = f.readlines()

    new_lines = []
    nme_found = False

    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            resname = line[17:20].strip()
            if resname == "NME":
                nme_found = True
            elif resname == "ACE" and nme_found:
                new_lines.append("TER\n")
                nme_found = False
        new_lines.append(line)

    with open(pdb_out, "w") as f:
        f.writelines(new_lines)

    return pdb_out

def cap_protein(in_file, out_file):
    u = mda.Universe(in_file)

    # Access each fragment separately
    res_start = 0
    segment_universes = []

    for seg in u.segments:
        chain = u.select_atoms(f"segid {seg.segid}")

        # --- Add ACE ---
        resid_c = chain.residues.resids[0]
        end_residue = u.select_atoms(f"segid {seg.segid} and resid {resid_c}")
        ace_positions = get_ace_pos(end_residue)
        ace_names = ["C", "CH3", "O"]
        resid = chain.residues.resids[0]
        kwargs = dict(
            n_atoms=len(ace_positions),
            name=ace_names,
            resname=len(ace_names)*["ACE"],
            positions=ace_positions,
            resids=resid*np.ones(len(ace_names)),
            segid=chain.segids[0]
        )

        ace_universe = create_universe(**kwargs)

        # --- Add NME ---
        resid_c = chain.residues.resids[-1]
        end_residue = u.select_atoms(f"segid {seg.segid} and resid {resid_c}")
        nme_positions = get_nme_pos(end_residue)
        nme_names = ["N", "CH3"]
        resid = chain.residues.resids[-1] + 2

        kwargs = dict(
            n_atoms=len(nme_names),
            name=nme_names,
            resname=len(nme_names)*["NME"],
            positions=nme_positions,
            resids=resid*np.ones(len(nme_names)),
            segid=chain.segids[0]
        )

        nme_universe = create_universe(**kwargs)

        # --- Handle OXT atom if present ---
        if "OXT" in end_residue.names:
            index = np.where(end_residue.names == "OXT")[0][0]
            OXT = end_residue[index]
            Chain = u.select_atoms(f"segid {seg.segid} and not index {OXT.index}")
        else:
            Chain = u.select_atoms(f"segid {seg.segid}")

        # --- Merge ACE, Protein, and NME ---
        u_all = mda.Merge(ace_universe.atoms, Chain, nme_universe.atoms)

        # --- Renumber residues ---
        resids_ace = [res_start+1, res_start+1, res_start+1]
        resids_pro = np.arange(resids_ace[0]+1, Chain.residues.n_residues+resids_ace[0]+1)
        resids_nme = [resids_pro[-1]+1, resids_pro[-1]+1]

        u_all.atoms.residues.resids = np.concatenate([resids_ace, resids_pro, resids_nme])
        res_start = u_all.atoms.residues.resids[-1]

        segment_universes.append(u_all)

    # --- Join all the universes ---
    all_uni = mda.Merge(*(seg.atoms for seg in segment_universes))
    all_uni.atoms.write(out_file)
    ter_output = out_file.replace(".pdb", "_TER.pdb")
    add_TER_between_chains(out_file, ter_output)


# Only parse CLI if run as script
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Add ACE/NME caps to a PDB")
    parser.add_argument("-i", dest="in_file", required=True)
    parser.add_argument("-o", dest="out_file", required=True)
    args = parser.parse_args()

    cap_protein(args.in_file, args.out_file)
