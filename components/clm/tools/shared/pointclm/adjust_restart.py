#!/usr/bin/env/python
import os, sys, csv, time, math
from optparse import OptionParser
import netcdf_functions as nffun

parser = OptionParser()

parser.add_option('--rundir', dest='rundir', default="", \
                    help = 'location of run directory')
parser.add_option('--casename', dest='casename', default="", \
                    help='Name of case')
parser.add_option('--restart_year', dest='restart_year', default="", \
                    help='Year of restart file to modify')
parser.add_option('--BGC', dest='bgc', default=False, action="store_true", \
                    help='Flag to set for BGC compsets')

(options,args)=parser.parse_args()

#The purpose of this code is to replace the restart value with the average value of the last met cycle
#This will bring the accelerated pools more into line with the expected results.

casename = options.casename
year     = options.restart_year

fname_restart = options.rundir+'/'+casename+'.clm2.r.0'+str(year)+ \
    '-01-01-00000.nc'
fname_hist    = options.rundir+'/'+casename+'.clm2.h1.0'+str(year)+ \
    '-01-01-00000.nc'

#save original restart file
#os.system('cp '+fname_restart+' '+fname_restart+'.orig')

#Accelerated pools:  Coarse woody debris, dead stem, coarse root, soil3 and soil4 C and N
var_names   = ['DEADSTEMC', 'DEADSTEMN', 'DEADSTEMP', 'DEADCROOTC', 'DEADCROOTN', 'DEADCROOTP']
if (options.bgc):
  var_names2d = ['CWDC_vr', 'CWDN_vr', 'CWDP_vr', 'SOIL3C_vr', 'SOIL3N_vr', 'SOIL3P_vr', 'SOIL2C_vr', \
                   'SOIL2N_vr', 'SOIL2P_vr']
else:
  var_names2d = ['CWDC_vr', 'CWDN_vr', 'CWDP_vr', 'SOIL4C_vr', 'SOIL4N_vr', 'SOIL4P_vr', 'SOIL3C_vr', \
                   'SOIL3N_vr', 'SOIL3P_vr']

for v in range(0,len(var_names)):
  hist_vals = nffun.getvar(fname_hist, var_names[v])
  rest_vals = nffun.getvar(fname_restart, var_names[v].lower())
  #Loop through all valid values in the restart file.
  n_rest = len(rest_vals)
  for i in range(0,n_rest):
    if (float(rest_vals[i]) > 0.0 and float(hist_vals[0][i]) < 1.0e10 and float(hist_vals[0][i]) > 0.001):
      print v, i, rest_vals[i], hist_vals[0][i]
      rest_vals[i] = hist_vals[0][i]
  ierr = nffun.putvar(fname_restart, var_names[v].lower(), rest_vals)
 
for v in range(0,len(var_names2d)):
  hist_vals = nffun.getvar(fname_hist, var_names2d[v])
  rest_vals = nffun.getvar(fname_restart, var_names2d[v].lower())
  #Loop through all valid values in the restart file.
  for i in range(0,1020642):
    print v, i
    for j in range(0,10):
      if (float(rest_vals[i][j]) > 0.0 and float(hist_vals[0][j][i]) < 1.0e10 and float(hist_vals[0][j][i]) > 1e-10):
        #print i,j, rest_vals[i][j], hist_vals[0][j][i]
        rest_vals[i][j] = hist_vals[0][j][i]
  ierr = nffun.putvar(fname_restart, var_names2d[v].lower(), rest_vals)
