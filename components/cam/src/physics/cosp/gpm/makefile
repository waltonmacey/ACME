.SUFFIXES : .F .f .c .o .a .f90 .f95 
CRTMFLAGS       := -I$(HOME)/local/CRTM/include  -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/atm/obj/cosp -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share
CRTMLINKS       := -L$(HOME)/local/CRTM/lib -lCRTM -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/atm/obj/cosp -lcosp -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/MCT/noesmf/a1l1r1i1o1g1w1/csm_share -lcsm_share

CAMHISTLINKS :=  -I/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/atm/obj -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/lib -latm -L/global/cscratch1/sd/yxl232/acme_scratch/cosp_gpm.v0.1.2.test/bld/intel/mpt/debug/nothreads/lib -lpio 


OBJS = gpm.o  gpm_crtm_sensor.o   \
   			gpm_crtm.o gpm_crtm_simulator.o \
				gpm_crtm_result.o  gpm_crtm_constants.o \
				gpm_utility.o  gpmsimulator_intr.o

gpm.exe : $(OBJS)
	ftn $(OBJS) $(CAMHISTLINKS)  -o gpm.exe $(CRTMFLAGS) $(CRTMLINKS)  

%.o: %.F90
	ftn $(CRTMFLAGS) $(CAMHISTLINKS) -c $<

gpm.o :  gpm_crtm_sensor.o gpm_crtm.o gpm_crtm_simulator.o gpm_utility.o gpmsimulator_intr.o
#gpm_gmi.o:  gpm_gmi_sensor.o
gpm_crtm_sensor.o: gpm_crtm_constants.o
gpm_crtm_simulator.o: gpm_crtm_result.o gpm_crtm.o
gpmsimulator.o: gpm_utility.o gpm_crtm_sensor.o
clean:
	rm *.o *.mod gpm.exe

test: ./gpm.exe
	./gpm.exe

