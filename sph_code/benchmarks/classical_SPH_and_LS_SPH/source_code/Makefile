#=========================#
#  intel fortran          # 
#=========================#
COMPILER = ifx

# For debug, add -check bounds -traceback
FFLAGS   = -O2 -xHost -fiopenmp -ipo -convert big_endian -fpp -heap-arrays \
		   -qmkl=sequential	-static-intel -I"${MKLROOT}/include" \
		   -D$(problem_name) -D$(kernel_function) -D$(SPH_model) -D$(particle_shift) \
		   -D$(wall_model) \
		   -D$(wall_uvw_top) -D$(wall_uvw_bottom) -D$(wall_uvw_left) -D$(wall_uvw_right)

# link for Intel MKL
FLIBS = -Wl,--start-group \
	    ${MKLROOT}/lib/libmkl_intel_lp64.a \
		${MKLROOT}/lib/libmkl_intel_thread.a \
		${MKLROOT}/lib/libmkl_core.a \
		-Wl,--end-group -liomp5 -lpthread -lm -ldl

#=========================#
#  compile options        #
#=========================#
problem_name    = TAYLOR_GREEN
kernel_function = WENDLAND_C2
SPH_model       = LSSPH_2ND
particle_shift  = PS_2ND

wall_model      = WL_2ND
wall_uvw_top    = FREE_SLIP_TOP
wall_uvw_bottom = FREE_SLIP_BOTTOM
wall_uvw_left   = FREE_SLIP_LEFT
wall_uvw_right  = FREE_SLIP_RIGHT


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
	lib_kernel_function.f90 \
	lib_solve_equation.f90 \
	check_compile_options.f90 \
	set_system_info.f90 \
	set_xyz_kind.f90 \
	set_VM_xyz_kind.f90 \
	set_physical_quantity.f90 \
	set_recalculation.f90 \
	set_background_cell.f90 \
	cal_EOS.f90 \
	cal_WL_0th.f90 \
	cal_WL_1st.f90 \
	cal_WL_LSSPH.f90 \
	cal_classical.f90 \
	cal_CSPH.f90 \
	cal_LSSPH.f90 \
	cal_cell_info.f90 \
	cal_VM_to_WL.f90 \
	cal_EOC_EOM_EOE.f90 \
	cal_Runge_Kutta.f90 \
	cal_outside.f90 \
	cal_PS_vec.f90 \
	cal_PS.f90 \
	out_input.f90 \
	out_system_info.f90 \
	out_progress.f90 \
	out_data.f90 \
	main.f90
	
#=========================#
#  write all files as .o  # 
#=========================#
OBJECTS = \
	input.o \
	global_variables.o \
	lib_file_operations.o \
	lib_kernel_function.o \
	lib_solve_equation.o \
	check_compile_options.o \
	set_system_info.o \
	set_xyz_kind.o \
	set_VM_xyz_kind.o \
	set_physical_quantity.o \
	set_recalculation.o \
	set_background_cell.o \
	cal_EOS.o \
	cal_WL_0th.o \
	cal_WL_1st.o \
	cal_WL_LSSPH.o \
	cal_classical.o \
	cal_CSPH.o \
	cal_LSSPH.o \
	cal_cell_info.o \
	cal_VM_to_WL.o \
	cal_EOC_EOM_EOE.o \
	cal_Runge_Kutta.o \
	cal_outside.o \
	cal_PS_vec.o \
	cal_PS.o \
	out_input.o \
	out_system_info.o \
	out_progress.o \
	out_data.o \
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
	@echo [message] compilation has finished
	@echo [message] input ./start_calculation
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
