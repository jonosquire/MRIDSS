#!/usr/bin/python

import os,shutil,math,sys
from WriteMRIDSS_Aux import *

MRIdir='../MRIDSS/'

Rm=1000.
Pm=2.
noise = 4.
filename='FullQl_Noise' + n2s(noise) + 'Rm' + n2s(Rm)+ 'Pm' + n2s(Pm)

class RD: # Passed to function to write input file
    equation_type='MHD_BQlin'
    L=[1.,math.pi,1.]
    N=[24,32,32]
    q=1.5
    nu=Pm / Rm
    eta=1 / Rm
    f_noise = noise
    dt=0.1
    time_interval = [0,5.]
    tv_save= 5*dt
    full_save=time_interval[1] + 1.
    save_enQ=1
    save_amQ=1
    save_dissQ=1
    save_ReyQ=1
    save_mfQ=1
    init_By = -1e-15
    remapQ=1
    QuasiLinearQ=1
    StartFromSavedQ=0


DELETE_OLD_FILES = 1 # Wether or not to delete old .DSSinput files

if DELETE_OLD_FILES: # If desired, move .DSSinput files to .OldInputFiles
    file_list=os.listdir(MRIdir)
    for fname in file_list:
        if fname.find('.DSSinput') != -1:
            shutil.move(MRIdir+fname, MRIdir+'.OldInputFiles/'+fname)
            print('Moved file ' + fname+' to .OldInputFiles/')

write_input_file(MRIdir + filename + '.DSSinput', RD)

print '<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>'
print 'Done writing '+ filename + '.DSSinput'
print '<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>'

# Run mpiexec
print ''
#os.system('mpiexec -np 8 ../mridss_prog '+ filename)
print 'Done case '+filename



