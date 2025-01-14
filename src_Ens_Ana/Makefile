all : install

install : Ens_Ana

clean : 
	rm *.o *.mod Ens_Ana

OBJS = main.o error.o output_time.o infjob.o inpinp.o intmp.o store_element.o atom_specifier.o \
       inref.o pdb_chk.o inweight.o inrefpdbL.o inrefbinL.o ranf.o Rg.o RMSDprep.o RMSDcalc.o \
       mkNOElst.o mkJCClst.o Xplor.o DISCOVER.o \
       PCA.o inptpl.o rdfree.o allocate_arrays.o efcol.o mk_vdWrad_list.o \
       Qvalue.o distribution.o contact.o measure.o pocket_search.o SOAP.o \
       inpdb.o inbin.o inrest.o \
       infil.o noe_calc.o geminal.o RMSF.o distance_calc.o angle_calc.o dihedral_calc.o outpdb.o DSSP.o
#       inotherfmt.o vmdfio.o hoqrvd_new.o \

FC = ifort
ifeq ($(FC), ifort)
    LFLAGS = -O2 -fpp -ipo -no-prec-div -no-prec-sqrt -xHost -qmkl=sequential -qopenmp
#    LFLAGS = -CB -fpp -traceback -pg -O2 -ipo -no-prec-div -xHost -qmkl=sequential -qopenmp
#    LFLAGS = -check all -fpp -traceback -pg -O2 -ipo -no-prec-div -xHost -qmkl=sequential -qopenmp
else ifeq ($(FC), gfortran)
    LFLAGS = -O2 -cpp -march=native -fopenmp
else
    $(error Unsupported compiler: $(FC))
endif

#CC = gcc


#place_lib = /usr/local/lib
#place_libmorph = libmorph_er.so
#place_lib = /home/usr3/m70293a/local/bin/lib
#place_libmorph = /home/usr3/m70293a/library/morph_area_fortran/libmorph_er.so

CFLAGS = -c $(LFLAGS)
#CFLAGS2 = -c

LIB = -ldl
#LIB = ${place_libmorph} -L${place_lib} -lvoro++ -ldl -lstdc++ -qmkl=sequential
#ERMOD_PLUGINS = vmdplugins/libexec

#Ens_Ana : COMINP.f90 COMIFN.f90 COMVAL.f90 morph_f.f90 trajectory.f90 $(OBJS) 
#	$(FC) $(LFLAGS) -o Analy.exe COMINP.f90 COMIFN.f90 COMVAL.f90 \
#	COMPAR.f90 morph_f.f90 trajectory.f90 $(OBJS) $(LIB)
Ens_Ana : COMINP.f90 COMIFN.f90 COMVAL.f90 COMTPL.f90 $(OBJS) 
	$(FC) $(LFLAGS) -o Ens_Ana COMINP.f90 COMIFN.f90 COMVAL.f90 COMTPL.f90 \
	$(OBJS) $(LIB)
main.o : main.f90 COMINP.f90 COMIFN.f90 COMVAL.f90 COMTPL.f90
	$(FC) $(CFLAGS) COMINP.f90 COMIFN.f90 COMVAL.f90 COMTPL.f90 main.f90

error.o : error.f90
	$(FC) $(CFLAGS) error.f90

output_time.o : output_time.f90
	$(FC) $(CFLAGS) output_time.f90

infjob.o : infjob.f90
	$(FC) $(CFLAGS) infjob.f90

inpinp.o : inpinp.f90 COMINP.f90
	$(FC) $(CFLAGS) inpinp.f90

intmp.o : intmp.f90 COMVAL.f90
	$(FC) $(CFLAGS) intmp.f90

store_element.o : store_element.f90 COMINP.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) store_element.f90

atom_specifier.o : atom_specifier.f90
	$(FC) $(CFLAGS) atom_specifier.f90

inref.o : inref.f90 COMVAL.f90
	$(FC) $(CFLAGS) inref.f90

pdb_chk.o : pdb_chk.f90 COMVAL.f90
	$(FC) $(CFLAGS) pdb_chk.f90

inweight.o : inweight.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) inweight.f90

inrefpdbL.o : inrefpdbL.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) inrefpdbL.f90

inrefbinL.o : inrefbinL.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) inrefbinL.f90

ranf.o : ranf.f90
	$(FC) $(CFLAGS) ranf.f90

Rg.o : Rg.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) Rg.f90

RMSDprep.o : RMSDprep.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) RMSDprep.f90

RMSDcalc.o : RMSDcalc.f90
	$(FC) $(CFLAGS) RMSDcalc.f90

mkNOElst.o : mkNOElst.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) mkNOElst.f90

mkJCClst.o : mkJCClst.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) mkJCClst.f90

Xplor.o : Xplor.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) Xplor.f90

DISCOVER.o : DISCOVER.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) DISCOVER.f90

PCA.o : PCA.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) PCA.f90

inptpl.o : inptpl.f90 COMVAL.f90 COMTPL.f90
	$(FC) $(CFLAGS) inptpl.f90

rdfree.o : rdfree.f90
	$(FC) $(CFLAGS) rdfree.f90

allocate_arrays.o : allocate_arrays.f90 COMTPL.f90
	$(FC) $(CFLAGS) allocate_arrays.f90

efcol.o : efcol.f90
	$(FC) $(CFLAGS) efcol.f90

mk_vdWrad_list.o : COMTPL.f90 mk_vdWrad_list.f90
	$(FC) $(CFLAGS) mk_vdWrad_list.f90

Qvalue.o : Qvalue.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) Qvalue.f90

distribution.o : distribution.f90 COMIFN.f90 COMVAL.f90 COMTPL.f90
	$(FC) $(CFLAGS) distribution.f90

ifdef ASA
contact.o : contact.f90 COMIFN.f90 COMVAL.f90 morph_f.f90
	$(FC) $(CFLAGS) morph_f.f90 contact.f90

else
contact.o : contact.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) contact.f90

endif
measure.o : measure.f90 COMVAL.f90
	$(FC) $(CFLAGS) measure.f90

pocket_search.o : pocket_search.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) pocket_search.f90

SOAP.o : SOAP.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) SOAP.f90

inbin.o : inbin.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) inbin.f90

inpdb.o : inpdb.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) inpdb.f90

inrest.o : inrest.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) inrest.f90

infil.o : infil.f90
	$(FC) $(CFLAGS) infil.f90

noe_calc.o : noe_calc.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) noe_calc.f90

geminal.o : geminal.f90
	$(FC) $(CFLAGS) geminal.f90

RMSF.o : RMSF.f90 COMIFN.f90 COMVAL.f90
	${FC} ${CFLAGS} RMSF.f90

distance_calc.o : distance_calc.f90 COMVAL.f90 COMIFN.f90
	$(FC) $(CFLAGS) distance_calc.f90

angle_calc.o : angle_calc.f90 COMVAL.f90 COMIFN.f90
	$(FC) $(CFLAGS) angle_calc.f90

dihedral_calc.o : dihedral_calc.f90 COMVAL.f90 COMIFN.f90
	$(FC) $(CFLAGS) dihedral_calc.f90

outpdb.o : outpdb.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) outpdb.f90

DSSP.o : DSSP.f90 COMIFN.f90 COMVAL.f90
	$(FC) $(CFLAGS) DSSP.f90

#inotherfmt.o : inotherfmt.f90 COMIFN.f90 COMVAL.f90 trajectory.f90
#	$(FC) $(CFLAGS) trajectory.f90 inotherfmt.f90
#
#vmdfio.o : vmdfio.c config.h
#	$(CC) $(CFLAGS2) vmdfio.c
#
#config.h: 
#	touch $@
#
#hoqrvd_new.o : hoqrvd_new.f90
#	$(FC) $(CFLAGS) hoqrvd_new.f90
#
#hoqrvd.o : hoqrvd.f90
#	$(FC) $(CFLAGS) hoqrvd.f90
#
