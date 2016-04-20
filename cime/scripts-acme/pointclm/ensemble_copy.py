#!/usr/bin/env python

import netcdf_functions as nffun
import os, sys, csv, time, math, numpy
from optparse import OptionParser

#Create, run and process a CLM/ALM model ensemble member
#  given specified case and parameters (see parm_list and parm_data files)
#  Parent case must be pre-built and all namelists in run directory.
#  Post-processing calcuates normalized sum of squared errors (SSE) given
#  data constraints specified in "constraints" directory"
#  DMRicciuto 12/1/2015
#
#  Note:  This will only work for single-point CLM/ALM compiled with MPI_SERIAL

#-------------------Parse options-----------------------------------------------

parser = OptionParser()

parser.add_option("--runroot", dest="runroot", default="../../run", \
                  help="Directory where the run would be created")
parser.add_option("--ens_num", dest="ensnum", default=1, \
                  help="Ensemble member number")
parser.add_option("--case", dest="casename", default="", \
                  help="Name of case")

(options, args) = parser.parse_args()


parm_names=[]
parm_indices=[]
parm_values=[]
myinput = open('./parm_list', 'r')
casename = options.casename

#get parameter names and PFT information
for s in myinput:
   pdata = s.split()
   parm_names.append(pdata[0])
   parm_indices.append(int(pdata[1]))
myinput.close()

#get parameter values
myinput = open('./parm_data', 'r')
for s in myinput:    
    parm_values.append(float(s))
myinput.close()
os.system('rm ./parm_data')

n_parameters = len(parm_names)
gst=str(100000+int(options.ensnum))

#create ensemble directory from original case 
orig_dir = str(os.path.abspath(options.runroot)+'/'+casename+'/run')
ens_dir  = os.path.abspath(options.runroot)+'/UQ/'+casename+'/g'+gst[1:]
est = str(100000+int(options.ensnum))
		
os.system('mkdir -p '+options.runroot+'/UQ/'+casename+'/g'+gst[1:]+'/timing/checkpoints')
os.system('cp  '+orig_dir+'/*_in* '+ens_dir)
os.system('cp  '+orig_dir+'/*nml '+ens_dir)
if (not ('CB' in casename)):
    os.system('cp  '+orig_dir+'/*stream* '+ens_dir)
os.system('cp  '+orig_dir+'/*.rc '+ens_dir)
os.system('cp  '+orig_dir+'/surf*.nc '+ens_dir)
os.system('cp  '+orig_dir+'/domain*.nc '+ens_dir)
os.system('cp  '+orig_dir+'/*para*.nc '+ens_dir)

#loop through all filenames, change directories in namelists, change parameter values
for f in os.listdir(ens_dir):
    if (os.path.isfile(ens_dir+'/'+f) and (f[-2:] == 'in' or f[-3:] == 'nml' or 'streams' in f)):
        myinput=open(ens_dir+'/'+f)
        myoutput=open(ens_dir+'/'+f+'.tmp','w')
        for s in myinput:
            if ('paramfile' in s):
                paramfile_orig = ((s.split()[2]).strip("'"))
                if (paramfile_orig[0:2] == './'):
                   paramfile_orig = orig_dir+'/'+paramfile_orig[2:]
                paramfile_new  = './clm_params_'+est[1:]+'.nc'
                os.system('cp '+paramfile_orig+' '+ens_dir+'/'+paramfile_new)
                myoutput.write(" paramfile = '"+paramfile_new+"'\n")
                pftfile = ens_dir+'/clm_params_'+est[1:]+'.nc'
                pnum = 0
                for p in parm_names:
                   if ('INI' not in p):
                      param = nffun.getvar(pftfile, p)
                      if (parm_indices[pnum] > 0):
                         param[parm_indices[pnum]-1] = parm_values[pnum]
                      elif (parm_indices[pnum] == 0):
                         param = parm_values[pnum]
                      else:
                         param[:] = parm_values[pnum]
                      ierr = nffun.putvar(pftfile, p, param)
                      pnum = pnum+1
            elif ('finidat = ' in s):
                finidat_file_orig = ((s.split()[2]).strip("'"))
                finidat_file_new  = ens_dir+'/'+(finidat_file_orig.split('/')[-1:])[0]
                if (finidat_file_orig[0:2] == './'):
                    finidat_file_orig = orig_dir+'/'+finidat_file_orig[2:]
                os.system('cp '+finidat_file_orig+' '+finidat_file_new)
                pnum = 0
                #Apply scaling factor for soil organic carbon (slowest  pool only)
                if ('BGC' in casename):
                    scalevars = ['soil3c_vr','soil3n_vr','soil3p_vr']
                else:
                    scalevars = ['soil4c_vr','soil4n_vr','soil4p_vr']
                sumvars = ['totsomc','totsomp','totcolc','totcoln','totcolp']
                for p in parm_names:
                    if ('.nc' in finidat_file_new and 'INI_somfac' in p):
                        for v in scalevars:
                            myvar = nffun.getvar(finidat_file_new, v)
                            myvar = parm_values[pnum] * myvar
                            ierr = nffun.putvar(finidat_file_new, v, myvar)
                        #TEMPORARY - add 3 gN to npool
                        myvar = nffun.getvar(finidat_file_new,'npool')
                        myvar = myvar+3.
                        ierr = nffun.putvar(finidat_file_new,'npool',myvar)
                    pnum=pnum+1
                myoutput.write(" finidat = '"+finidat_file_new+"'\n")
            elif ('logfile =' in s):
                os.system('date +%y%m%d-%H%M%S > mytime')
                mytinput=open('./mytime','r')
                for st in mytinput:
                    timestr = st.strip()
                mytinput.close()
                os.system('rm mytime')
                myoutput.write(s.replace('`date +%y%m%d-%H%M%S`',timestr))
            else:
                myoutput.write(s.replace(orig_dir,ens_dir))
        myoutput.close()
        myinput.close()
        os.system(' mv '+ens_dir+'/'+f+'.tmp '+ens_dir+'/'+f)
