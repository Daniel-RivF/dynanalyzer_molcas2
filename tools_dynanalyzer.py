import numpy as np
import os
# Read lists of output names in all subdirectories using os or glob, functions work with one file.

def readfile(filename):
    with open(filename,'r') as f:
     #content = f.readlines()
     content = f.read().splitlines()
    for i in range(len(content)-1, -1, -1):
        if '       New Coordinates' in content[i]:
            a = i
            if isinstance(a,int):
                break
    workfile = content[:a]
    return  workfile

def chunk_between_tokens(content_list,token1,token2):
    token1_indexes = [i for i,k in enumerate(content_list) if token1 in k]
    token2_indexes = [j for j,l in enumerate(content_list) if token2 in l]
    if len(token1_indexes) == len(token2_indexes):
        u = zip(token1_indexes,token2_indexes[:len(token1_indexes)])
    else:
        print 'ERROR SHIT CODE'
#    chunks = []
#    for ind1,ind2 in u:
#        chunks.append(content_list[ind1:ind2])
    chunks = [content_list[ind1:ind2] for ind1,ind2 in u]
    return chunks

########################################################
#################### ENERGIES ##########################
########################################################

def listof_2RAS_Energies(workfile):
    # First to energies to read from RASSCF
    # The remaining from Tully (OOLGnuplot)
    token1ras,token2ras = 'Wave function printout:','Pseudonatural active orbitals and approximate occupation numbers'
    RASdata = chunk_between_tokens(workfile,token1ras,token2ras)[:2]
    token1E = 'Final state energy(ies):'
    token2E = 'Molecular orbitals:'
    ras_E = []
    for RAS in RASdata:
        #ras_E.append(chunk_between_tokens(RAS,token1E,token2E))
        a = chunk_between_tokens(RAS,token1E,token2E)
        for i in a:
            listaras = [line.split() for line in i]
            n = len(listaras)-1
            enerRAS = [ p[-1] for p in listaras[3:n]]
            ras_E.append(enerRAS)
            # ras_E as nested list of strings.
    return ras_E


def listof_tully_energies(workfile):
    # Is NStates really necessary?
    # Only works if energies are printed in the same line that OOLGnplt token.
    a = [ i.split() for i in workfile if 'OOLgnuplt:' in i]
    eners = []
    for j in a:
        n = len(j)
        eners.append(j[3:n])
    # eners as list of strings.
    hops_index_E = []
    OOL_indexes = [i for i,k in enumerate(workfile) if 'OOLgnuplt:' in k]
    HOP_indexes = [i for i, x in enumerate(workfile) if 'ALLOWED' in x]
    tohopindex = [ i-1 for i in HOP_indexes]
    afterhop = [OOL_indexes.index(i) for i in tohopindex]
    hop = [OOL_indexes[i-1] for i in afterhop]
    hop_E = [workfile[i] for i in hop]
    all_E = [workfile[j] for j in OOL_indexes]
    N_of_points = len(all_E) + 2 
    return (all_E, hop_E)


def just_energies(workfile):
    first2RAS = listof_2RAS_Energies(workfile)
    first2RAS_L = [i[-1] for i in first2RAS]
    # Includes the first two RAS energies ( before Tulyy starts).
    OOL = [ i.split()[1:] for i in listof_tully_energies(workfile)[0] ]
    OOL_hop = [ j.split()[1:] for j in listof_tully_energies(workfile)[1] ]
    NStates = ( len(OOL[0]) -1 ) / 2
    E_states = [p[NStates:NStates*2] for p in OOL]
    E_states_hop = [ p[NStates:NStates*2] for p in OOL_hop ] 
    Ecurrent = [ p[NStates*2] for p in OOL ]
    return ( first2RAS + E_states,E_states_hop,first2RAS_L + Ecurrent)

def times(workfile ,time_step_fs):
    tully_Es = listof_tully_energies(workfile)[0]
    tully_points = len(tully_Es)
    n_of_steps = tully_points + 2
    times = [ (i+1)*time_step_fs for i in  range(n_of_steps) ]
    ### Times HOP (to compare with Alessio)
    tully_hop_Es = listof_tully_energies(workfile)[1]
    times_hop =  [ ( tully_Es.index(j) + 3)*time_step_fs for j in tully_hop_Es]
    return (times, times_hop)
    #return (times, [times_hop[0]])


def populations(workfile):
    # Depends on NStates
    a = [ i.split() for i in workfile if 'OOLgnuplt:' in i]
    b = [ map(float,i[1:]) for i in a]
    NStates = ( len(b[0]) - 1 ) / 2
    # Filter pops (should be between 0 and 1):
    pop = [ j[0:NStates] for j in b]
    return pop



###################################################################
########################### STRUCTURAL STUFF ######################
###################################################################

def trajectory_extractor(filename):
    #trajectory as nested list (or np.array)
    workfile = readfile(filename)
    token1 = 'Cartesian coordinates in Angstrom:'
    token2 = 'Nuclear repulsion energy'
    geoms_info = chunk_between_tokens(workfile,token1,token2)
    geoms_a = [i[4:-1] for i in geoms_info]
    Natom = len(geoms_a[0])
    geoms = []
    for geom in geoms_a:
        geom_i = [g.split()[1:] for g in geom]
        geoms.append(geom_i)
    ## Removing index from Z:
    for ge in geoms:
        # Make more general, only takes fist character of atom symbol, wont work for two char names like Na, Fe, Sr, etc.
        for k in ge:
            k[0] = k[0][0]
            ##if isinstance(k[0][1],int):
            ##    k[0] = k[0][0]
            ##    Zs.append(k[0])
            ##elif isinstance(k[0][1],str):
            ##    k[0] = k[0][::1]
            ##    Zs.append(k[0])

    return geoms


#### Writer for trajectories ######
def writer_xyz(filename):
    geoms = trajectory_extractor(filename)
    Natom = len(geoms[0])
    fileout = os.path.splitext(filename)[0] + '.MD.xyz'
    with open(fileout,'w') as fo:
     for ge_w in geoms:
         fo.write('%s \n' % Natom)
         fo.write('chiripitiflautico \n')
         for u in ge_w:
             fo.write('%s \n' % ' '.join(u) ) 
    
    return 











    














   #def rindex(lis, item):
   #for i in range(len(lis)-1, -1, -1):
   #if item == lis[i]:
   #return i
   #else:
   #raise ValueError("rindex(lis, item): item not in lis")









     #def parse_molcas_output(filename,tokens):
     #    return result
     #
     #
     #
     #def parse_input_options(user_input):
     #    return options
     #
     #
     #def extract_xyz(filename):
     #    # Extract the xyz trajectory
     #    # To read: Natom, atomic numbers (Z), symbols of the atoms , xyz_coordinates
     #    return xyz_traj_data
     #
     #def writer_xyz(file_xyz_name):
     #    # Takes output of the previous function
     #    return 
     #
     #
     #
     #def enepop_data(filename,fileout):
     #    # format columns: 
     #    # filename step time energy1 energy2 energyN relaxed_stateEnergy population
     #    # igual es mejor crear funciones pequenas que parseen cada dato, mirar sobre la marcha
     #    return something
     #
     #def internal_coordfile(filename):
     #    # Nombre del fileout igual que el filename pero sin extension .out.
     #    return data_internal
     #
     #def extract_FC_data(filename):
     #    # Extracts, energies, excitation energies, and internals (can use the previous function for internal coordinates.
     #    return FC_data
     #




