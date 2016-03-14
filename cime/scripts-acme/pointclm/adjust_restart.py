#!/usr/bin/env/python
import os, sys, csv, time, math
from optparse import OptionParser

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

def getvar(fname, varname):
    usescipy = False
    try:
    	import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"r")
        var = nffile.variables[varname]
        varvals = var.copy()    #works for vector only?
        nffile.close()
    else:
        varvals=[]
    	nffile = netcdf.NetCDFFile(fname,"r")
    	var = nffile.variables[varname]
    	varvals = var.getValue()
    	nffile.close()
    return varvals

def putvar(fname, varname, varvals):
    usescipy = False
    try:
        import Scientific.IO.NetCDF as netcdf
    except ImportError:
        import scipy
        from scipy.io import netcdf
        usescipy = True
    if (usescipy):
        nffile = netcdf.netcdf_file(fname,"a")
        var = nffile.variables[varname]
        var[:] = varvals[:]
        nffile.close()
    else:
        nffile = netcdf.NetCDFFile(fname,"a")
        var = nffile.variables[varname]
        var.assignValue(varvals)
        nffile.close()
    ierr = 0
    return ierr

casename = options.casename
year     = options.restart_year

fname_restart = options.rundir+'/'+casename+'.clm2.r.0'+str(year)+ \
    '-01-01-00000.nc'
fname_hist    = options.rundir+'/'+casename+'.clm2.h1.0'+str(year)+ \
    '-01-01-00000.nc'

#save original restart file
os.system('cp '+fname_restart+' '+fname_restart+'.orig')

#Accelerated pools:  Coarse woody debris, dead stem, coarse root, soil3 and soil4 C and N
var_names   = ['DEADSTEMC', 'DEADSTEMN', 'DEADSTEMP', 'DEADCROOTC', 'DEADCROOTN', 'DEADCROOTP']
if (options.bgc):
  var_names2d = ['CWDC_vr', 'CWDN_vr', 'CWDP_vr', 'SOIL3C_vr', 'SOIL3N_vr', 'SOIL3P_vr', 'SOIL2C_vr', \
                   'SOIL2N_vr', 'SOIL2P_vr']
else:
  var_names2d = ['CWDC_vr', 'CWDN_vr', 'CWDP_vr', 'SOIL4C_vr', 'SOIL4N_vr', 'SOIL4P_vr', 'SOIL3C_vr', \
                   'SOIL3N_vr', 'SOIL3P_vr']

for v in range(0,len(var_names)):
  hist_vals = getvar(fname_hist, var_names[v])
  rest_vals = getvar(fname_restart, var_names[v].lower())
  #Loop through all valid values in the restart file.
  for i in range(0,len(rest_vals)):
    if (float(rest_vals[i]) > 0.0 and float(hist_vals[0][i]) < 1.0e10 and float(hist_vals[0][i]) > 0.001):
      print rest_vals[1], hist_vals[0][i]
      rest_vals[i] = hist_vals[0][i]
      ierr = putvar(fname_restart, var_names[v].lower(), rest_vals)
 
for v in range(0,len(var_names2d)):
  hist_vals = getvar(fname_hist, var_names2d[v])
  rest_vals = getvar(fname_restart, var_names2d[v].lower())
  #Loop through all valid values in the restart file.
  for i in range(0,2):
    for j in range(0,10):
      if (float(rest_vals[i][j]) > 0.0 and float(hist_vals[0][j][i]) < 1.0e10 and float(hist_vals[0][j][i]) > 1e-10):
        print i,j, rest_vals[i][j], hist_vals[0][j][i]
        rest_vals[i][j] = hist_vals[0][j][i]
        ierr = putvar(fname_restart, var_names2d[v].lower(), rest_vals)
       

  
