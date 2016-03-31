.SUFFIXES : .F .f .c .o .a .f90 .f95 

CASEDIR := /global/cscratch1/sd/yxl232/acme_scratch/gpm_capt.v0.1.0.ne30_ne30

# CRTM compiled without -g
CRTMFLAGS       := -I$(HOME)/local/CRTM/include -I$(CASEDIR)/bld/atm/obj/cosp -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share 
CRTMLINKS       := -L$(CASEDIR)/bld/atm/obj/cosp -lcosp -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share -lcsm_share  -L$(HOME)/local/CRTM/lib -lCRTM   

# CRTM compiled with -g
#CRTMFLAGS       := -I$(HOME)/local/CRTM_debug/include  -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/atm/obj/cosp -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share
#CRTMLINKS       := -L$(HOME)/local/CRTM_debug/lib -lCRTM -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/atm/obj/cosp -lcosp -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share -lcsm_share

NCFLAG := -I/opt/cray/netcdf/4.3.3.1/INTEL/14.0/include

NCLINK := -L/opt/cray/netcdf/4.3.3.1/INTEL/14.0/lib -lnetcdff -L/opt/cray/hdf5/1.8.14/INTEL/14.0/lib -lnetcdf


CAMHISTLINKS :=  -I/global/cscratch1/sd/yxl232/acme_scratch/SMS_D.ne16_ne16.FC5COSP.corip1_intel.cospfix.latest/bld/atm/obj -L/global/cscratch1/sd/yxl232/acme_scratch/SMS_D.ne16_ne16.FC5COSP.corip1_intel.cospfix.latest/bld/lib -latm  -L/global/cscratch1/sd/yxl232/acme_scratch/SMS_D.ne16_ne16.FC5COSP.corip1_intel.cospfix.latest/bld/intel/mpt/debug/nothreads/lib -lpio

DEBUGFLAGS := -free  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source  -O0 -g -check uninit -check bounds -check pointers -fpe1 -check noarg_temp_created  -warn all  -gen-interfaces 


OBJS = gmi_offline.o 

gmi.exe : $(OBJS)
	ftn $(OBJS)    -o gmi.exe $(DEBUGFLAGS) $(CRTMLINKS) 

%.o: %.F90
	ftn $(CRTMFLAGS)   -c $<  $(DEBUGFLAGS) $(CRTMFLAGS) 

clean:
	rm *.o *.mod gmi.exe

test: ./gpm.exe
	./gpm.exe

