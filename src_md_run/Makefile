all : install

MODOBJ = compar.o combas.o comerg.o commis.o compsc.o comcmm.o \
         comcmmc.o commol.o phycns.o calc_time.o commpi.o

OBJS = 	bstft.o calcm.o  caldih.o calrsf.o caltmp.o chklin.o chkpha.o \
	combin.o conj.o  covang.o cutoff_PB.o cutoff_CB.o \
	dnergy.o dnecap.o dneres_PB.o dneres_CB.o dnergv_PB.o dnergv_CB.o \
	eigjac.o factor.o psr.o  infjob.o inpbou.o \
	inpcrb.o inpinp.o inpmd.o  inpmin.o \
	inpout.o inppdb.o inprst.o inptpb.o inptpl.o allocate_arrays.o \
	input.o  inpvar.o md.o mini.o mkabsr.o mksrdg.o monmd.o  \
	monmin.o outcpu.o outcpu2.o outcrb.o outpdb.o output.o \
	outrst.o outtpb.o outtpl.o ranno.o ranun.o  dhr.o dsr.o remotv.o \
	rescrd.o scalet.o setshk.o setvel.o shake_PB.o shake_CB.o \
	simulu.o sryelm.o staver.o steep.o verlet.o verlethe.o \
	ellpbc.o ellprp.o inpcls.o rdlstatm.o setmnt.o inmcmd.o mnstmd.o inprep.o \
	repprp.o celprp.o incel.o celpar_PB.o celpar_CB.o celoth.o celpol.o neartb.o \
	cmmdo.o enear_CMM_redu.o enear_CMM_high.o enear_CMM_doub.o \
	enear_ACT_redu_PB.o enear_ACT_high_PB.o enear_ACT_doub_PB.o \
	enear_ACT_redu_CB.o enear_ACT_high_CB.o enear_ACT_doub_CB.o \
	enear_ZD_redu_PB.o enear_ZD_high_PB.o enear_ZD_doub_PB.o \
	enear_ZD_redu_CB.o enear_ZD_high_CB.o enear_ZD_doub_CB.o \
        enear_extCMM.o setextCMM.o PBcord.o SIMD_check.o atom_specifier.o

COBJS =	flopen.o flclos.o efcol.o  rdfree.o

# Linux Intel Fortran
FC = ifort

ifeq ($(FC), ifort)
    LFLAGS = -O2 -fpp -ipo -no-prec-div -no-prec-sqrt -xHost -qopenmp
#    LFLAGS = -CB -fpp -traceback -pg -O2 -ipo -no-prec-div -xHost -qopenmp
#    LFLAGS = -check all -fpp -traceback -pg -O2 -ipo -no-prec-div -xHost -qmkl=sequential -qopenmp
else ifeq ($(FC), gfortran)
    LFLAGS = -O2 -cpp -march=native -fopenmp
else
    $(error Unsupported compiler: $(FC))
endif

CFLAGS = -c $(LFLAGS)

#       -lm -> Option to link math carnel library
#LIB = -lm -lpthread
LIB = -lm

install : md_run mkshkl

clean :
	rm *.o *.mod md_run mkshkl enear_*_*.f90 dneres_*.f90 dnergv_*.f90 shake_*.f90 cutoff_*.f90 celpar_*.f90

md_run : $(MODOBJ) main.o $(OBJS) $(COBJS) 
	$(FC) $(LFLAGS) $(MODOBJ) main.o $(OBJS) $(COBJS) -o md_run $(LIB)

main.o : main.f90 COMPAR.f90 COMBAS.f90 COMCMM.f90 COMCMMC.f90 COMMIS.f90 COMMPI.f90
	$(FC) $(CFLAGS) main.f90 $(LIB)

mkshkl : COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 PHYCNS.f90 \
         mkshkl.f90 inptpl.o allocate_arrays.o $(COBJS)
	$(FC) $(LFLAGS) COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 \
              PHYCNS.f90 mkshkl.f90 inptpl.o allocate_arrays.o $(COBJS) -o mkshkl $(LIB)

compar.o : COMPAR.f90
	$(FC) $(CFLAGS) COMPAR.f90 -o compar.o $(LIB)

combas.o : COMPAR.f90 COMBAS.f90
	$(FC) $(CFLAGS) COMBAS.f90 -o combas.o $(LIB)

comerg.o : COMPAR.f90 COMERG.f90
	$(FC) $(CFLAGS) COMERG.f90 -o comerg.o $(LIB)

commis.o : COMPAR.f90 COMMIS.f90
	$(FC) $(CFLAGS) COMMIS.f90 -o commis.o $(LIB)

compsc.o : COMPAR.f90 COMPSC.f90
	$(FC) $(CFLAGS) COMPSC.f90 -o compsc.o $(LIB)

comcmm.o : COMPAR.f90 COMCMM.f90
	$(FC) $(CFLAGS) COMCMM.f90 -o comcmm.o $(LIB)

comcmmc.o : COMCMMC.f90
	$(FC) $(CFLAGS) COMCMMC.f90 -o comcmmc.o $(LIB)
 
commol.o : COMMOL.f90
	$(FC) $(CFLAGS) COMMOL.f90 -o commol.o $(LIB)

phycns.o : PHYCNS.f90
	$(FC) $(CFLAGS) PHYCNS.f90 -o phycns.o $(LIB)

calc_time.o : CALC_TIME.f90
	$(FC) $(CFLAGS) CALC_TIME.f90 -o calc_time.o $(LIB)

commpi.o : COMMPI.f90
	$(FC) $(CFLAGS) COMMPI.f90 -o commpi.o $(LIB)

bstft.o : bstft.f90
	$(FC) $(CFLAGS) bstft.f90 $(LIB)

calcm.o : calcm.f90 COMCMMC.f90
	$(FC) $(CFLAGS) calcm.f90 $(LIB)

caldih.o : caldih.f90 
	$(FC) $(CFLAGS) caldih.f90 $(LIB)

calrsf.o : calrsf.f90
	$(FC) $(CFLAGS) calrsf.f90 $(LIB)

caltmp.o : caltmp.f90 COMMOL.f90
	$(FC) $(CFLAGS) caltmp.f90 $(LIB)

chklin.o : chklin.f90 COMCMMC.f90
	$(FC) $(CFLAGS) chklin.f90 $(LIB)

chkpha.o : chkpha.f90 PHYCNS.f90
	$(FC) $(CFLAGS) chkpha.f90 $(LIB)

combin.o : combin.f90
	$(FC) $(CFLAGS) combin.f90 $(LIB)

conj.o : conj.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) conj.f90 $(LIB)

covang.o : covang.f90 PHYCNS.f90
	$(FC) $(CFLAGS) covang.f90 $(LIB)

cutoff_PB.o : cutoff.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat cutoff.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > cutoff_PB.f90
	$(FC) $(CFLAGS) cutoff_PB.f90 $(LIB)

cutoff_CB.o : cutoff.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat cutoff.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > cutoff_CB.f90
	$(FC) $(CFLAGS) cutoff_CB.f90 $(LIB)

dneres_PB.o : dneres.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 \
           PHYCNS.f90 COMCMM.f90 COMCMMC.f90
	cat dneres.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > dneres_PB.f90
	$(FC) $(CFLAGS) dneres_PB.f90 $(LIB)

dneres_CB.o : dneres.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 \
           PHYCNS.f90 COMCMM.f90 COMCMMC.f90
	cat dneres.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > dneres_CB.f90
	$(FC) $(CFLAGS) dneres_CB.f90 $(LIB)

dnecap.o : dnecap.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 \
           PHYCNS.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) dnecap.f90 $(LIB)

dnergv_PB.o : dnergv.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 PHYCNS.f90 \
           COMCMM.f90 COMCMMC.f90
	cat dnergv.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > dnergv_PB.f90
	$(FC) $(CFLAGS) dnergv_PB.f90 $(LIB)

dnergv_CB.o : dnergv.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 PHYCNS.f90 \
           COMCMM.f90 COMCMMC.f90
	cat dnergv.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > dnergv_CB.f90
	$(FC) $(CFLAGS) dnergv_CB.f90 $(LIB)

dnergy.o : dnergy.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 \
           COMCMMC.f90
	$(FC) $(CFLAGS) dnergy.f90 $(LIB)

efcol.o : efcol.f90
	$(FC) $(CFLAGS) efcol.f90 $(LIB)

eigjac.o : eigjac.f90
	$(FC) $(CFLAGS) eigjac.f90 $(LIB)

factor.o : factor.f90
	$(FC) $(CFLAGS) factor.f90 $(LIB)

flclos.o : flclos.f90 
	$(FC) $(CFLAGS) flclos.f90 $(LIB)

flopen.o : flopen.f90
	$(FC) $(CFLAGS) flopen.f90 $(LIB)

psr.o : psr.f90 COMPAR.f90 COMBAS.f90 COMMIS.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) psr.f90 $(LIB)

infjob.o : infjob.f90
	$(FC) $(CFLAGS) infjob.f90 $(LIB)

inmcmd.o : inmcmd.f90 COMPAR.f90 COMBAS.f90 COMPSC.f90 COMCMM.f90
	$(FC) $(CFLAGS) inmcmd.f90 $(LIB)

inpbou.o : inpbou.f90 COMPAR.f90 COMBAS.f90 COMMIS.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) inpbou.f90 $(LIB)

inpcrb.o : inpcrb.f90
	$(FC) $(CFLAGS) inpcrb.f90 $(LIB)

inpinp.o : inpinp.f90 COMCMM.f90 COMERG.f90
	$(FC) $(CFLAGS) inpinp.f90 $(LIB)

inpmd.o : inpmd.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 \
          COMCMM.f90 COMPSC.f90 COMMOL.f90 PHYCNS.f90
	$(FC) $(CFLAGS) inpmd.f90 $(LIB)

inpmin.o : inpmin.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 \
           PHYCNS.f90
	$(FC) $(CFLAGS) inpmin.f90 $(LIB)

inpout.o : inpout.f90
	$(FC) $(CFLAGS) inpout.f90 $(LIB)

inppdb.o : inppdb.f90
	$(FC) $(CFLAGS) inppdb.f90 $(LIB)

inprst.o : inprst.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) inprst.f90 $(LIB)

inptpb.o : inptpb.f90 COMPAR.f90 COMBAS.f90 COMERG.f90
	$(FC) $(CFLAGS) inptpb.f90 $(LIB)

inptpl.o : inptpl.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 PHYCNS.f90
	$(FC) $(CFLAGS) inptpl.f90 $(LIB)

allocate_arrays.o : allocate_arrays.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90
	$(FC) $(CFLAGS) allocate_arrays.f90 $(LIB)

input.o : input.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 \
          PHYCNS.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) input.f90 $(LIB)

inpvar.o : inpvar.f COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMMC.f90
	$(FC) $(CFLAGS) inpvar.f $(LIB)

md.o : md.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 COMCMM.f90 \
       COMPSC.f90 PHYCNS.f90 COMMOL.f90 COMCMMC.f90 CALC_TIME.f90
	$(FC) $(CFLAGS) md.f90 $(LIB)

mini.o : mini.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 COMCMM.f90\
         COMCMMC.f90
	$(FC) $(CFLAGS) mini.f90 $(LIB)

mkabsr.o : mkabsr.f90 COMBAS.f90 COMCMM.f90
	$(FC) $(CFLAGS) mkabsr.f90 $(LIB)

mksrdg.o : mksrdg.f90 COMBAS.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) mksrdg.f90 $(LIB)

mnstmd.o : mnstmd.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMPSC.f90 \
           PHYCNS.f90
	$(FC) $(CFLAGS) mnstmd.f90 $(LIB)

monmd.o : monmd.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 \
          PHYCNS.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) monmd.f90 $(LIB)

monmin.o : monmin.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 \
           COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) monmin.f90 $(LIB)

outcpu.o : outcpu.f90
	$(FC) $(CFLAGS) outcpu.f90 $(LIB)

outcpu2.o : outcpu2.f90
	$(FC) $(CFLAGS) outcpu2.f90 $(LIB)

outcrb.o : outcrb.f90 COMCMMC.f90
	$(FC) $(CFLAGS) outcrb.f90 $(LIB)

outpdb.o : outpdb.f90
	$(FC) $(CFLAGS) outpdb.f90 $(LIB)

output.o : output.f90 COMPAR.f90 COMBAS.f90 COMCMM.f90
	$(FC) $(CFLAGS) output.f90 $(LIB)

outrst.o : outrst.f90
	$(FC) $(CFLAGS) outrst.f90 $(LIB)

outtpb.o : outtpb.f90 COMPAR.f90 COMBAS.f90 COMERG.f90
	$(FC) $(CFLAGS) outtpb.f90 $(LIB)

outtpl.o : outtpl.f COMPAR.f90 COMBAS.f90 COMERG.f90 PHYCNS.f90
	$(FC) $(CFLAGS) outtpl.f $(LIB)

ranno.o : ranno.f90
	$(FC) $(CFLAGS) ranno.f90 $(LIB)

ranun.o : ranun.f90
	$(FC) $(CFLAGS) ranun.f90 $(LIB)

dhr.o : dhr.f COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 \
           PHYCNS.f90
	$(FC) $(CFLAGS) dhr.f $(LIB)

dsr.o : dsr.f90 COMPAR.f90 COMBAS.f90 COMMIS.f90
	$(FC) $(CFLAGS) dsr.f90 $(LIB)

rdfree.o : rdfree.f90
	$(FC) $(CFLAGS) rdfree.f90 $(LIB)

remotv.o : remotv.f90 COMCMMC.f90
	$(FC) $(CFLAGS) remotv.f90 $(LIB)

rescrd.o : rescrd.f90 COMMIS.f90
	$(FC) $(CFLAGS) rescrd.f90 $(LIB)

scalet.o : scalet.f90 COMCMMC.f90
	$(FC) $(CFLAGS) scalet.f90 $(LIB)

inpcls.o : inpcls.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 \
           PHYCNS.f90
	$(FC) $(CFLAGS) inpcls.f90 $(LIB)

rdlstatm.o : rdlstatm.f90 COMBAS.f90 COMCMM.f90
	$(FC) $(CFLAGS) rdlstatm.f90 $(LIB)

setmnt.o : setmnt.f90 COMPAR.f90 COMBAS.f90 COMPSC.f90 COMERG.f90
	$(FC) $(CFLAGS) setmnt.f90 $(LIB)

inprep.o : inprep.f90 COMPAR.f90 COMBAS.f90 COMERG.f90
	$(FC) $(CFLAGS) inprep.f90 $(LIB)

repprp.o : repprp.f90 COMPAR.f90 COMBAS.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) repprp.f90 $(LIB)

setshk.o : setshk.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90
	$(FC) $(CFLAGS) setshk.f90 $(LIB)

setvel.o : setvel.f90 COMCMMC.f90
	$(FC) $(CFLAGS) setvel.f90 $(LIB)

shake_PB.o : shake.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 PHYCNS.f90
	cat shake.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > shake_PB.f90
	$(FC) $(CFLAGS) shake_PB.f90 $(LIB)

shake_CB.o : shake.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMMIS.f90 PHYCNS.f90
	cat shake.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > shake_CB.f90
	$(FC) $(CFLAGS) shake_CB.f90 $(LIB)

simulu.o : simulu.f90
	$(FC) $(CFLAGS) simulu.f90 $(LIB)

sryelm.o : sryelm.f90
	$(FC) $(CFLAGS) sryelm.f90 $(LIB)

staver.o : staver.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 PHYCNS.f90 COMCMM.f90\
           COMCMMC.f90
	$(FC) $(CFLAGS) staver.f90 $(LIB)

steep.o : steep.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 PHYCNS.f90 COMCMM.f90\
          COMCMMC.f90
	$(FC) $(CFLAGS) steep.f90 $(LIB)

verlet.o : verlet.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 PHYCNS.f90 COMCMM.f90 \
           COMCMMC.f90 CALC_TIME.f90
	$(FC) $(CFLAGS) verlet.f90 $(LIB)

verlethe.o : verlethe.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMPSC.f90 \
             COMMIS.f90 COMMOL.f90 PHYCNS.f90 COMCMM.f90 COMCMMC.f90 \
             CALC_TIME.f90
	$(FC) $(CFLAGS) verlethe.f90 $(LIB)

ellpbc.o : ellpbc.f90 COMBAS.f90 COMCMMC.f90
	$(FC) $(CFLAGS) ellpbc.f90 $(LIB)

ellprp.o : ellprp.f90 COMCMMC.f90
	$(FC) $(CFLAGS) ellprp.f90 $(LIB)

celprp.o : celprp.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) celprp.f90 $(LIB)

incel.o : incel.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) incel.f90 $(LIB)

celpar_PB.o : celpar.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat celpar.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > celpar_PB.f90
	$(FC) $(CFLAGS) celpar_PB.f90 $(LIB)

celpar_CB.o : celpar.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat celpar.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > celpar_CB.f90
	$(FC) $(CFLAGS) celpar_CB.f90 $(LIB)

celoth.o : celoth.f90 COMPAR.f90 COMBAS.f90 COMCMM.f90
	$(FC) $(CFLAGS) celoth.f90 $(LIB)

celpol.o : celpol.f90 COMPAR.f90 COMBAS.f90 COMCMM.f90
	$(FC) $(CFLAGS) celpol.f90 $(LIB)

neartb.o : neartb.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90
	$(FC) $(CFLAGS) neartb.f90 $(LIB)

cmmdo.o : cmmdo.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) cmmdo.f90 $(LIB)

enear_CMM_redu.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#CMM/ && $$0 ~ /\#REDU/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#CMM/ || $$0 ~ /\#REDU/) {print $$0} } else {print $$0}}' > enear_CMM_redu.f90
	$(FC) $(CFLAGS) enear_CMM_redu.f90 $(LIB)

enear_CMM_high.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#CMM/ && $$0 ~ /\#HIGH/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#CMM/ || $$0 ~ /\#HIGH/) {print $$0} } else {print $$0}}' > enear_CMM_high.f90
	$(FC) $(CFLAGS) enear_CMM_high.f90 $(LIB)

enear_CMM_doub.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#CMM/ && $$0 ~ /\#DOUB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#CMM/ || $$0 ~ /\#DOUB/) {print $$0} } else {print $$0}}' > enear_CMM_doub.f90
	$(FC) $(CFLAGS) enear_CMM_doub.f90 $(LIB)

enear_ACT_redu_PB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ACT/ && $$0 ~ /\#REDU/ && $$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ACT/ || $$0 ~ /\#REDU/ || $$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > enear_ACT_redu_PB.f90
	$(FC) $(CFLAGS) enear_ACT_redu_PB.f90 $(LIB)

enear_ACT_high_PB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ACT/ && $$0 ~ /\#HIGH/ && $$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ACT/ || $$0 ~ /\#HIGH/ || $$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > enear_ACT_high_PB.f90
	$(FC) $(CFLAGS) enear_ACT_high_PB.f90 $(LIB)

enear_ACT_doub_PB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ACT/ && $$0 ~ /\#DOUB/ && $$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ACT/ || $$0 ~ /\#DOUB/ || $$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > enear_ACT_doub_PB.f90
	$(FC) $(CFLAGS) enear_ACT_doub_PB.f90 $(LIB)

enear_ACT_redu_CB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ACT/ && $$0 ~ /\#REDU/ && $$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ACT/ || $$0 ~ /\#REDU/ || $$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > enear_ACT_redu_CB.f90
	$(FC) $(CFLAGS) enear_ACT_redu_CB.f90 $(LIB)

enear_ACT_high_CB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ACT/ && $$0 ~ /\#HIGH/ && $$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ACT/ || $$0 ~ /\#HIGH/ || $$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > enear_ACT_high_CB.f90
	$(FC) $(CFLAGS) enear_ACT_high_CB.f90 $(LIB)

enear_ACT_doub_CB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ACT/ && $$0 ~ /\#DOUB/ && $$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ACT/ || $$0 ~ /\#DOUB/ || $$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > enear_ACT_doub_CB.f90
	$(FC) $(CFLAGS) enear_ACT_doub_CB.f90 $(LIB)

enear_ZD_redu_PB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ZD/ && $$0 ~ /\#REDU/ && $$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ZD/ || $$0 ~ /\#REDU/ || $$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > enear_ZD_redu_PB.f90
	$(FC) $(CFLAGS) enear_ZD_redu_PB.f90 $(LIB)

enear_ZD_high_PB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ZD/ && $$0 ~ /\#HIGH/ && $$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ZD/ || $$0 ~ /\#HIGH/ || $$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > enear_ZD_high_PB.f90
	$(FC) $(CFLAGS) enear_ZD_high_PB.f90 $(LIB)

enear_ZD_doub_PB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ZD/ && $$0 ~ /\#DOUB/ && $$0 ~ /\#PB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ZD/ || $$0 ~ /\#DOUB/ || $$0 ~ /\#PB/) {print $$0} } else {print $$0}}' > enear_ZD_doub_PB.f90
	$(FC) $(CFLAGS) enear_ZD_doub_PB.f90 $(LIB)

enear_ZD_redu_CB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ZD/ && $$0 ~ /\#REDU/ && $$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ZD/ || $$0 ~ /\#REDU/ || $$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > enear_ZD_redu_CB.f90
	$(FC) $(CFLAGS) enear_ZD_redu_CB.f90 $(LIB)

enear_ZD_high_CB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ZD/ && $$0 ~ /\#HIGH/ && $$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ZD/ || $$0 ~ /\#HIGH/ || $$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > enear_ZD_high_CB.f90
	$(FC) $(CFLAGS) enear_ZD_high_CB.f90 $(LIB)

enear_ZD_doub_CB.o : enear.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	cat enear.f90 | gawk '{ if ($$0 ~ /\#SLCT/) { if ($$0 ~ /\#ZD/ && $$0 ~ /\#DOUB/ && $$0 ~ /\#CB/) {print $$0} } else if ($$0 ~ /\#SLC2/) { if ($$0 ~ /\#ZD/ || $$0 ~ /\#DOUB/ || $$0 ~ /\#CB/) {print $$0} } else {print $$0}}' > enear_ZD_doub_CB.f90
	$(FC) $(CFLAGS) enear_ZD_doub_CB.f90 $(LIB)

enear_extCMM.o : enear_extCMM.f90 COMPAR.f90 COMBAS.f90 COMERG.f90 COMCMM.f90 COMCMMC.f90
	$(FC) $(CFLAGS) enear_extCMM.f90 $(LIB)

setextCMM.o : setextCMM.f90 COMCMM.f90
	$(FC) $(CFLAGS) setextCMM.f90 $(LIB)

PBcord.o : PBcord.f90 COMBAS.f90 COMCMMC.f90
	$(FC) $(CFLAGS) PBcord.f90 $(LIB)

SIMD_check.o : SIMD_check.f90 COMCMM.f90
	$(FC) $(CFLAGS) SIMD_check.f90 $(LIB)

atom_specifier.o : atom_specifier.f90
	$(FC) $(CFLAGS) atom_specifier.f90

