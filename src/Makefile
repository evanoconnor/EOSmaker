include ../make.inc

MODINC=$(HDF5INCS)

SOURCES= eos_table_module.F90 derivatives_production.F90 electron_eos.F90 eos_table.F90 fix_compositions.F90 fix_table_units.F90 helpers.F90 interpolate_initial_table_cubic_hermite.F90 loweos.F90 read_table.F90 setup_final_table.F90 trilin_interp.F90 findNANs.F90 write_table.F90

FSOURCES= modified_timmes.f

OBJECTS=$(SOURCES:.F90=.o)
FOBJECTS=$(FSOURCES:.f=.o)

all: eos_table 

eos_table: $(OBJECTS) $(FOBJECTS)
	$(F95) $(LDFLAGS) -o ../eos_table $(OBJECTS) $(FOBJECTS) $(HDF5LIBS)

$(OBJECTS): %.o: %.F90
	$(F95) $(F95FLAGS) $(MODINC) -c $< -o $@

$(FOBJECTS): %.o: %.f
	$(F95) $(F95FLAGS) $(MODINC) -c $< -o $@

clean:	
	rm ../eos_table *.o *.mod
