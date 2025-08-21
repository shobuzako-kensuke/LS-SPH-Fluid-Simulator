# Intel Fortran ===============================================================#
COMPILER = ifx
FFLAGS   = -O2 -xHost -ipo -fiopenmp -convert big_endian -static-intel \
           -fpp -qmkl=sequential -static-intel -I"${MKLROOT}/include" \
		   -D$(SPH_model) -D$(kernel_function) -D$(variable_density) \
		   -D$(virtual_marker) -D$(wall_accuracy) -D$(wall_bottom) \
		   -D$(wall_top) -D$(wall_left) -D$(wall_right)

# Link for Intel MKL ==========================================================#
FLIBS = -Wl,--start-group \
	    ${MKLROOT}/lib/libmkl_intel_lp64.a \
		${MKLROOT}/lib/libmkl_intel_thread.a \
		${MKLROOT}/lib/libmkl_core.a \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl

# Compilation Options =========================================================#
SPH_model        = LS_SPH_2ND
kernel_function  = WENDLAND_C2
variable_density = CONSTANT_DENSITY
virtual_marker   = VM_ON
wall_accuracy    = WL_1ST
wall_bottom      = NEUMANN_BOTTOM
wall_top         = NEUMANN_TOP
wall_left        = DIRICHLET_LEFT
wall_right       = DIRICHLET_RIGHT


#==============================================================================#
#=========================== NOT CHANGE BELOW =================================#
#==============================================================================#

#=========================#
#  file name              #
#=========================#
TARGET = start_calculation

#=========================#
# write all files as .f90 #
#=========================#
COMMON_MOD = \
	input.f90 \
	global_variables.f90 \
	lib_file_operations.f90 \
	lib_kernel_functions.f90 \
	check_options.f90 \
	set_VM.f90 \
	set_system.f90 \
	set_background_cell.f90 \
	cal_background_cell.f90 \
	cal_CSPH_0th.f90 \
	cal_classical_Laplacian.f90 \
	cal_LSSPH_typeA.f90 \
	cal_LSSPH_typeB.f90 \
	cal_density.f90 \
	cal_VM_to_WL.f90 \
	cal_diffusion.f90 \
	out_input.f90 \
	out_system_info.f90 \
	out_data.f90 \
	out_progress.f90 \
	main.f90
	
#=========================#
#  write all files as .o  #
#=========================#
OBJECTS = \
	input.o \
	global_variables.o \
	lib_file_operations.o \
	lib_kernel_functions.o \
	check_options.o \
	set_VM.o \
	set_system.o \
	set_background_cell.o \
	cal_background_cell.o \
	cal_CSPH_0th.o \
	cal_classical_Laplacian.o \
	cal_LSSPH_typeA.o \
	cal_LSSPH_typeB.o \
	cal_density.o \
	cal_VM_to_WL.o \
	cal_diffusion.o \
	out_input.o \
	out_system_info.o \
	out_data.o \
	out_progress.o \
	main.o
	
#=========================#
#  define rules: f90 >> o #
#=========================#
.SUFFIXES:
.SUFFIXES: .o .f90

#=========================#
#  compile .f90 >> .o     #
#=========================#
.f90.o:
	$(COMPILER) -c $(FFLAGS) $<

$(TARGET): $(OBJECTS)
	$(COMPILER) $(FFLAGS) $(FLIBS) -o $(TARGET) $(OBJECTS)
	@echo + -------------------------------------------------------- +
	@echo [message] All programs have been successfully compiled.
	@echo [message] Please input: ./start_calculation
	@echo + -------------------------------------------------------- +

#=========================#
#  if module changed      #
#=========================#
$(OBJECTS): $(COMMON_MOD)

#=========================#
#  make clean             #
#=========================#
# Linux=@rm, Command Prompt=@del
.PHONY: clean
clean:
	@rm -rf $(TARGET)     $(OBJECTS) *.mod
#	@del    $(TARGET).exe $(OBJECTS) *.mod
	@echo + -------------------------------------------------------- +
	@echo [message] make clean
	@echo + -------------------------------------------------------- +

# END #
