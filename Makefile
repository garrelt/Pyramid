#F90 = gfortran

#GFORTFLAGS = -O

#F90FLAGS1 = $(GFORTFLAGS) 

#F90FLAGS = $(F90FLAGS1)

#OPTIONS = $(F90FLAGS)

#LDFLAGS = $(OPTIONS)
#LIBS = 

F90 = ifort

#GFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -DIFORT -shared-intel -check all -traceback

GFORTFLAGS = -g -u -fpe0 -DIFORT -check all -traceback

F90FLAGS1 = $(GFORTFLAGS) 

F90FLAGS = $(F90FLAGS1)

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS)
LIBS = 

main:	precision.o input.o cgsconstants.o cgsastroconstants.o type.o array.o romberg.o output.o field.o radiation.o pyramid_analytical_column_density.o pyramid_analytical_photo_rate.o pyramid_global_source_transformation.o cartesian_analytical_column_density.o cartesian_analytical_photo_rate.o integration.o solid_angle.o pyramid_initialization.o pyramid_column_density.o pyramid_photo_rate.o pyramid_source_domain_transformation.o case_division.o geometry_transformation.o short_coordinate_transformation.o short.o short_column_density.o short_photo_rate.o long_coordinate_transformation.o long.o long_column_density.o long_photo_rate.o trilinear.o main.o

	$(F90) $(OPTIONS) -o $@ precision.o input.o cgsconstants.o cgsastroconstants.o type.o array.o romberg.o output.o field.o radiation.o pyramid_analytical_column_density.o pyramid_analytical_photo_rate.o pyramid_global_source_transformation.o cartesian_analytical_column_density.o cartesian_analytical_photo_rate.o integration.o solid_angle.o pyramid_initialization.o pyramid_column_density.o pyramid_photo_rate.o pyramid_source_domain_transformation.o case_division.o geometry_transformation.o short_coordinate_transformation.o short.o short_column_density.o short_photo_rate.o long_coordinate_transformation.o long.o long_column_density.o long_photo_rate.o trilinear.o main.o

clean: 
	rm -f *.o *.mod *.l *.il

.f.o:
	$(F90) -c $(OPTIONS) $<

.f90.o:
	$(F90) -c $(OPTIONS) $<

.F90.o:
	$(F90) -c $(OPTIONS) $<

.f90.mod:
	$(F90) -c $(OPTIONS) $<

.SUFFIXES: .f90 .F90 .mod .o


