CRTMFLAGS       := -I$(HOME)/local/CRTM/include  -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/atm/obj/cosp -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share
CRTMLINKS       := -L$(HOME)/local/CRTM/lib -lCRTM -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/atm/obj/cosp -lcosp -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share -lcsm_share


gpm.exe : gpm.o gpm_gmi.o
	ftn gpm_gmi.o  gpm.o -o gpm.exe $(CRTMFLAGS) $(CRTMLINKS)  

gpm.o : gpm.F90 gpm_gmi.o
	ftn $(CRTMFLAGS) -c gpm.F90

gpm_gmi.o: gpm_gmi.F90
	ftn $(CRTMFLAGS) -c gpm_gmi.F90



