import tools_dynanalyzer as tools
import numpy as np
import math, os
from numpy.linalg import norm
# Input options

class AngleGeometryError(Exception): pass 
class DihedralGeometryError(Exception): pass


# Internals as nested list, check if distance, angle or dihedral.
# Make folders: Trajectories, EnePop, DATA, DATA_swarm
# write the EnePoP file, and one data (and dataHOP) file for each user-defined coordinate. Write one all_data file for each defined coordinate (columns: point, current state, populations, energies, current energy and coordinate value).
# Write a file with info about all hops.
# join all the data files for the swarm of trajs.

def transform2array(geoms):
    o = []
    for geom in geoms:
        n = [map(float,j[1:]) for j in geom]
        o.append(n)
    geom_matrix = [ np.array(i) for i in o ]
    return geom_matrix


def unit_vect(v):
    return v / norm(v)



def mybool(v1_u,v2_u):
   acc = True
   for x,y in zip(v1_u, v2_u):
       acc & (x ==y)
   return acc

def getgeoms(filename):
    g = tools.trajectory_extractor(filename)
    return g

#geoms = tools.trajectory_extractor(filename)

###########################################################################
################# Distances #############################################
##########################################################################


def get_Distance(geoms, index1, index2):
    geom_matrix = transform2array(geoms)
    dists = []
    for g in geom_matrix:
        dist = norm(g[index2-1] - g[index1-1])
        dists.append(dist)
    return dists

#############################################################################
################# Angles ####################################################
############################################################################

def angle(v1,v2):
    v1u = unit_vect(v1)
    v2u = unit_vect(v2)
    #angle = np.arccos(np.dot(v1u, v2u))
    angle = math.acos(np.dot(v1u, v2u))
    if np.isnan(angle):
       if (v1u == v2u).all():
           return 0.0
       else:
           return np.pi
    return angle



def get_angles(geoms, index1,index2,index3):
    coords = transform2array(geoms)
    angles = []
    for ge in coords:
        v12 =  ge[index1 -1] - ge[index2-1] 
        v32 = ge[index3-1] -ge[index2 - 1]
        ang = angle(v12,v32)
        angles.append(ang)
    angles_deg = [ i*(180/np.pi) for i in angles]
    return angles_deg


################################################################################
######################### Dihedral angles ######################################
################################################################################


#def dihedral(v21,v32,v43):
#    normal1u = unit_vect(np.cross(v21,v32))
#    normal2u = unit_vect(np.cross(v32,v43))
#    torsion = angle(normal1u,normal2u) * 180 / np.pi
#    return torsion
      

def dihedral(v21,v32,v43):
    normal1u = unit_vect(np.cross(v21,v32))
    normal2u = unit_vect(np.cross(v32,v43))
    v32n = unit_vect(v32)
    m1 = np.cross(normal1u,v32n)
    x = np.dot(normal1u,normal2u)
    y = np.dot(m1,normal2u)
    torsion = math.atan2(y,x) * 180 / math.pi
    return torsion

#def correction_dihedrals(uncorr_dih_list):
#    ref0 = uncorr_dih_list[0]
#    ref = uncorr_dih_list[0]
#    dlist = []
#    for d in uncorr_dih_list[1:]:
#        if d > 0:
#            p = d - 360
#        elif d < 0:
#            p = d + 360
#        wa = d - ref
#        wb = p - ref
#        if math.fabs(wa) <= math.fabs(wb):
#            dihed = d
#        else: 
#            dihed = p
#        dlist.append(dihed)
#        ref = dihed
#    return [ref0] + dlist

def correction_dihedrals(uncorr_dih_list):
    ref0 = uncorr_dih_list[0]
    ref = uncorr_dih_list[0]
    if ref > 0:
        ref = ref - 360
        ref0 = ref0 - 360 
    dlist = []
    for d in uncorr_dih_list[1:]:
        if d > 0:
            p = d - 360
        elif d < 0:
            p = d + 360
        wa = d - ref
        wb = p - ref
        if math.fabs(wa) <= math.fabs(wb):
            dihed = d
        else: 
            dihed = p
        dlist.append(dihed)
        ref = dihed
    return [ref0] + dlist

def get_dihedrals(geoms, i1, i2, i3, i4):
    # Returns dihedral angles in degrees, corrected 
    coords = transform2array(geoms)
    dihedrals_deg = []
    for ge in coords:
        v21 = ge[i2 -1] - ge[i1-1]
        v32 = ge[i3-1] - ge[i2-1]
        v43 = ge[i4-1] - ge[i3-1]
        dih = dihedral(v21,v32,v43)
        dihedrals_deg.append(dih)
    corr_dihs = correction_dihedrals(dihedrals_deg)
#    if corr_dihs[0] < 0:
#        a = [i+360 for i in corr_dihs]
#        return a
#    else:
    return corr_dihs
    

###### WRITERS ###########################
### Write files to plot using Gnuplot ####
##########################################


########## Energy, populations vs t ####################

def Enepop_data(filename,time_step):
    workfile = tools.readfile(filename)
    ener_states = tools.just_energies(workfile)[0]
    ener_hops = tools.just_energies(workfile)[1]
    current_eners = tools.just_energies(workfile)[2]
    populations = tools.populations(workfile)
    first2pops = [[0,1.0],[0,1.0]]
    allpops = [ map(str,i) for i in first2pops + populations]
    times = tools.times(workfile ,time_step)[0]
    timeshop = tools.times(workfile ,time_step)[1]
    # Same name as the molcas .out but different extension.
    fileout = os.path.splitext(filename)[0] + '.EnePop'
    with open(fileout,'w') as fo:
        for i ,j, k, l in zip(map(str,times),allpops,ener_states,current_eners):
            fo.write("{tim} {pop} {eners} {currE} \n".format(tim=i,pop=' '.join(j)  ,eners=' '.join(k), currE= l))
    return




##############  Coordinates vs t (and hop points) ###############################################

#def dataFile(filename,internal_args):
    # Internalcoordinates_args 

def check_internal(geoms,argument):
    class ErrorArgument(Exception): pass
    #Takes a list as argument and checks its length ( 2 for distance, 3 for angles, 4 for dihedrals)
    if all(isinstance(item,int) for item in argument):
        if len(argument) == 2:
            i,j = argument[0],argument[1]
            d = get_Distance(geoms,i,j)
            return d
        elif len(argument) == 3:
            i,j,k = argument[0],argument[1],argument[2]
            angles = get_angles(geoms, i,j,k)
            return angles
        elif len(argument) == 4:
            i,j,k,l = argument[0],argument[1],argument[2],argument[3]
            dihs = get_dihedrals(geoms,i,j,k,l)
            return dihs
        elif   len(argument) >= 5 or len(argument) == 0:
            raise ErrorArgument('Error, list must have less than five elements')
    raise ErrorArgument('List must contain the indexes of atoms i.e. [1,2] for the distance between atoms 1 and 2')




####################################################################################
############################## DATA WRITER #########################################
###################################################################################


def main_user_coordinates(filename,delta_t,arguments):
    #arguments will be a nested list with the desired coordinates.
    workfile = tools.readfile(filename)
    geoms = tools.trajectory_extractor(filename)
    internals = [ check_internal(geoms,argument) for argument in arguments ] 
    coordinates = zip(*internals)
    times_all = tools.times(workfile ,delta_t)[0]
    times_hop = tools.times(workfile ,delta_t)[1]
    time_internals = zip(times_all,coordinates)
    inds = [times_all.index(i) for i in times_hop]
    time_internals_hop = [time_internals[j] for j in inds]
    return (time_internals , time_internals_hop)

          
   
def internals_writer(filename_inp,delta_t,arguments):
    data = main_user_coordinates(filename_inp,delta_t,arguments)
    fdata = os.path.splitext(filename_inp)[0] + '.data'
    with open(fdata,'w') as fd:
        for i in data[0]:
            fd.write('%.3f %s \n' % ( i[0] , ' '.join(map(str,i[1]))))
    if data[1] != []:
        fd_hop = os.path.splitext(filename_inp)[0] + '.hopdata'
        with open(fd_hop,'w') as fdh:
            for j in data[1]:
                fdh.write('%.3f %s \n' % ( j[0] , ' '.join(map(str,j[1]))))
        return
    else:
        return


#########################################################################################
####################### COMPLETE SWARM OF TRAJECTORIES ######################
#########################################################################################

def joiner_data(list_data,fileout_swarm):
    with open(fileout_swarm , 'w') as fsw:
        for filedata in list_data:
            with open(filedata,'r') as f:
                a = f.readlines()
                fa = [i.split() for i in a]
                for j in fa:
                    fsw.write('%s \n' % ' '.join(j) )
                fsw.write('  \n')
    return

            


#        # parse individual data files (readlines or some stud like that).
#        for filename in list_outs:
#            file_data = os.path.splitext(filename_inp)[0]
#            data_i = main_user_coordinates(filename,delta_t,arguments)
#

    






    
















