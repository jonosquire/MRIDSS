# scons build file for building MRIDSS outside of xcode

import os

# Compiler flags, libraries etc.   
env = Environment(CXX = 'mpic++', CCFLAGS = '-Wall -Ofast -ffast-math -march=native')
print "CXX is:", env['CXX']
env.Append(CPPPATH = ['/usr/local/include/'])
env.Append(LIBPATH = ['/usr/local/lib/'])
env.Append(LIBS = ['fftw3']) #'mpi_cxx','mpi',

# Source files

base = 'MRIDSS/'
models = base+'Models/'
aux = base+'Auxiliary/'
integs = base+'Integrators/'
s = '   '

basesrc = Split(base + 'main.cc')
modsrc = Split(models+'MHD_FullUBQlin.cc'     +s+    models+'MHD_fullBQlin.cc'     +s+    models+'MHD_BQlin.cc'     +s+     models+'HD_fullU.cc'     +s+     models+'MHD_fullBQlin_noMaxwell.cc'       +s+    models+'Constant_Damping.cc')
intsrc = Split(integs+'intEuler.cc'   +s+     integs+'intEulerCN.cc'  +s+   integs+'intRK2CN.cc'     +s+      integs+'intRK3CN.cc')
auxsrc = Split(aux+'Input_parameters.cc'  +s+   aux+'General_Auxiliary.cc'   +s+      aux+'Initialization_routines.cc'     +s+     aux+'MPIdata.cc'    +s+    aux+'fftwPlans.cc'    +s+    aux+'TimeVariables.cc'    +s+    aux+'InterfaceOutput.cc')

# Compile
env.Program('mridss_prog',basesrc+modsrc+intsrc+auxsrc)



