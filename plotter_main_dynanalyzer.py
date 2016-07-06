#!/usr/bin/python

import os, glob
import tools_dynanalyzer as tools_dyn
import tools_dynanalyzer_user as tools_usr

arguments = [ [5,8],[1,7,13,8],[4,8,7,14],[4,5,6,11],[5,6,12,1],[6,7] ]
time_step = 0.499741501
fileouts = sorted(glob.glob('*.out'))

### ORGANIZING STUFF INTO FOLDERS #######

curr_dir = os.getcwd()
resultsdir = curr_dir + '/RESULTS'
coords = resultsdir + '/coords/'
enepopf = resultsdir+'/EnePop/'
alltrajs = resultsdir+'/ALL_TRAJS/'


if not os.path.exists(resultsdir):
    os.makedirs(resultsdir)
    os.makedirs(coords)
    os.makedirs(enepopf)
    os.makedirs(alltrajs)

############ GENERATES THE MAIN ANALYSIS FILES ################

#for filename in fileouts:
#    print 'Extracting data from file %s' % filename
#    tools_dyn.writer_xyz(filename)
#    tools_usr.Enepop_data(filename,time_step)
#    tools_usr.internals_writer(filename,time_step,arguments)

###################JOIN DATAFILES TO PLOT THE SWARM OF TRAJECTORIES #############

list_data = sorted(glob.glob('*.data'))
list_datahop = sorted(glob.glob('*.hopdata')) 

fileout_swarm = 'alltrajs.data'
fileout_swarmhop = 'alltrajs.hopdata'

tools_usr.joiner_data(list_data,fileout_swarm)
tools_usr.joiner_data(list_datahop,fileout_swarmhop)

###################  GNUPLOT SCRIPTS ##########################


ncols = len(arguments) + 1
cols = [ i+1 for i in range(ncols)[1:]]

for col in map(str,cols):
    basenametrajs = os.path.splitext(fileout_swarm)[0]
    filegp = basenametrajs + 'coord' + col + '.gp'
    with open(filegp, 'w') as fgp:
        fgp.write('set title "Trajectories" \n')
        fgp.write('set output "%s" \n' % (basenametrajs + 'coord' + col + '.png'))
        fgp.write('set term pngcairo size 1440,900 enhanced font ", 15 " \n')
        fgp.write('set ylabel "coordinate %s (set label manually in the script)" \n' % col)
        fgp.write('set xlabel "t (fs)" \n')
        fgp.write('set xrange [0:1000] \n')
        fgp.write('plot "%s" u 1:%s w l lc rgb "black", "%s" u 1:%s pt 7 ps 0.7 lc rgb "#0000FF" \n' % (fileout_swarm,col,fileout_swarmhop,col))
    os.system('gnuplot < %s' % filegp)
    namepng = basenametrajs + 'coord' + col + '.png'
    os.rename(namepng , alltrajs + '/' + namepng)

for col in map(str,cols):
    a = [os.path.splitext(i)[0] for i in list_data if 'geom' in i]
    b = [os.path.splitext(i)[0] for i in list_datahop if 'geom' in i]
    for fdata in a:
        fgpil = fdata + '.gp'
        with open(fgpil,'w') as fgpi:
            fgpi.write('set title "Trajectory "  \n')
            fgpi.write('set output "%s" \n' % (fdata + 'coord' + col + '.png'))
            fgpi.write('set term pngcairo size 1440,900 enhanced font ", 15 " \n')
            fgpi.write('set ylabel "coordinate %s (set label manually in the script)" \n' % col)
            fgpi.write('set xlabel "t (fs)" \n')
            fgpi.write('set xrange [0:1000] \n')
            if fdata in b:
                fgpi.write('plot "%s" u 1:%s w l lc rgb "black", "%s" u 1:%s pt 7 ps 0.7 lc rgb "#0000FF" \n' % (fdata + '.data' , col, fdata + '.hopdata',col))
            else:
                fgpi.write('plot "%s" u 1:%s w l lc rgb "black" \n' % (fdata + '.data' , col))
        os.system('gnuplot < %s' % fgpil)
        namepng = fdata + 'coord' + col + '.png'
        os.rename(namepng,coords + '/' + namepng)



########################## GNUPLOT ENERGIES AND PPULATIONS ############################


files_enepop = sorted(glob.glob('*.EnePop'))

for enepopfile in files_enepop:
    basename_enepop = os.path.splitext(enepopfile)[0]
    filegp_enepop = basename_enepop + '.Enepop.gp'
    with open(filegp_enepop,'w') as fgp_e:
        fgp_e.write('set title "Trajectory %s" \n' % enepopfile)
        fgp_e.write('set output "%s" \n' % ( basename_enepop + 'EnePop.png') )
        fgp_e.write('set term pngcairo size 1440,900 enhanced font ", 15 " \n')
        fgp_e.write('set ylabel "E (a.u.)"\n')
        fgp_e.write('set xrange [0:1200] \n')
        fgp_e.write('set y2range [0:1] \n')
        # Only for 2 states:
        fgp_e.write('plot "%s" u 1:2 axes x1y2 w filledcurves x1 lt 1 lc rgb "#E5E5E5" t "S0 Population", "" u  1:4 w lines linecolor rgb "green" t "S0" , "" u 1:5  w lines linecolor rgb "red" t "S1" , "" u 1:6 w lines linecolor rgb "black"  lt "dashed" lw 3 t "RlxRoot" ' % enepopfile)
    os.system('gnuplot < %s' % filegp_enepop)
    namepng = basename_enepop + 'EnePop.png'
    os.rename(namepng, enepopf + '/' + namepng)


