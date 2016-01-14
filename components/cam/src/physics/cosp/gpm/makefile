CRTMFLAGS       := -I$(HOME)/local/CRTM/include  -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/atm/obj/cosp -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share
CRTMLINKS       := -L$(HOME)/local/CRTM/lib -lCRTM -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/atm/obj/cosp -lcosp -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share -lcsm_share

CAMHISTLINKS :=  -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/atm/obj -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/lib -latm -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/lib -lpio 


gpm.exe : gpm.o gpm_gmi.o gpm_gmi_sensor.o
	ftn gpm_gmi.o gpm_gmi_sensor.o  gpm.o   $(CAMHISTLINKS)  -o gpm.exe $(CRTMFLAGS) $(CRTMLINKS)  

gpm.o : gpm.F90 gpm_gmi.o gpm_gmi_sensor.o
	ftn $(CRTMFLAGS)  -c gpm.F90

gpm_gmi.o: gpm_gmi.F90 gpm_gmi_sensor.o
	ftn $(CRTMFLAGS) $(CAMHISTLINKS)  -c gpm_gmi.F90

cmparray_mod.o: cmparray_mod.F90
	ftn $(CRTMFLAGS) $(CAMHISTLINKS)  -c cmparray_mod.F90

gpm_gmi_sensor.o: gpm_gmi_sensor.F90
	ftn $(CRTMFLAGS) $(CAMHISTLINKS)   -c gpm_gmi_sensor.F90

