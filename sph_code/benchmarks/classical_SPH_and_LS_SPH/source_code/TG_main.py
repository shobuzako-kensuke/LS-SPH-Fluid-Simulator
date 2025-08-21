#=========================#
#  file name              # 
#=========================#
save_name = 'TG_test'
#=========================#
#  analysis file          # 
#=========================#
analysis_step = 'n'  # if specifying analysis step, input file ID
#=========================#
#  program switch         # 
#=========================#
exe_fig_ini = 'y'    # initial settings
exe_uv_plot = 'y'    # uv profile
exe_p_time  = 'y'    # time vs pressure
exe_energy  = 'y'    # time vs kinematic energy
exe_fig_u   = 'n'    # snapshot of u
exe_fig_v   = 'n'    # snapshot of v
exe_fig_UV  = 'y'    # snapshot of speed
exe_fig_rho = 'n'    # snapshot of density
exe_fig_pre = 'y'    # snapshot of pressure
add_wall    = 'n'    # add wall to snapshots
#=========================#
#  SPH model              # 
#=========================#
kernel_switch = 3    # 1=cubic, 2=quintic, 3=wenC2, 4=wenC4, 5=wenC6
SPH_model     = 2    # 0=classical, 1=CSPH, 2=LSSPH


#==============================================================================#
#=========================== NOT CHANGE BELOW =================================#
#==============================================================================#

#=========================#
#  module                 # 
#=========================#
import os
import numpy as np
import time
import math
import ana_read_fortran_files as read_f90
import TG_initial_settings as ini
import TG_energy_decay as energy_decay
import TG_uv_plot as uv_plot
import TG_p_vs_time as p_vs_time
import TG_snapshot as snap

# cpu time
start_time = time.perf_counter()
print('+ -------------------------------------------------------- +')
print('[message] main.py has started')

# mkdir
save_path = '../fig'
os.makedirs(save_path, exist_ok=True)

#=========================#
#  read variables         # 
#=========================#
parent_path   = '../output/{}'.format(save_name)                    # path of the saved data
variable_path = parent_path + '/system_info/system_info_VALUES.dat' # file path

start_step, end_step, write_step, total_step, \
N_sys, N_inn, N_out, N_VM, \
SP_mass, Delta_x, h, h_eff, WL_thick, W_ave, \
zeta, xi, N_x, N_y, N_h, width, height, \
Delta_t_CFL, Delta_t_vis, Delta_t_th, Delta_t_CFL_relax, Delta_t_vis_relax, Delta_t_eff, \
rho_ref, vis_ref, K_ref, k_th, c_p, alpha, \
U_top, T_top, T_bot, gravity, \
Delta_temp, temp_ave, sound_speed, kinematic_vis, thermal_dif, \
Re, Ra, Pr, \
= read_f90.read_variables(variable_path) # read

#=========================#
#  count                  # 
#=========================#
N_file = 0
if (analysis_step == 'n'):
    last = end_step + 1 # including initial conditions
else:
    last = analysis_step + 1

for i in range(0, last, write_step):
    file_path = parent_path + '/data/{}/SP_xyz.dat'.format(i) # ex) SP_xyz
    if os.path.exists(file_path): # if file exists
        N_file += 1               # the number of the saved files
        last_file = i             # save last file name  
    else:                         # if file does not exist
        break
print('[message] number of read files: {}'.format(N_file))
print('+ -------------------------------------------------------- +')
#==============================================================================#
#                            main program below                                #
#==============================================================================#

#=========================#
#  initial conditions     # 
#=========================#
if (exe_fig_ini == 'y'):
    print('[message] exe_fig_ini has started')
    save_path = '../fig/{}/initial_settings'.format(save_name)
    os.makedirs(save_path, exist_ok=True) # mkdir
    ini.xyz_kind(N_sys, N_inn, WL_thick, parent_path, save_path)
    print('          >> xyz_kind')
    ini.cell(N_sys, N_inn, WL_thick, h_eff, parent_path, save_path)
    print('          >> cell')
    ini.xyz_kind_VM(N_VM, N_inn, WL_thick, parent_path, save_path)
    print('          >> xyz_kind_VM')
    ini.cell_VM(N_VM, N_inn, WL_thick, h_eff, parent_path, save_path)
    print('          >> cell_VM')
    print('[message] exe_fig_ini has finished')
    print('+ -------------------------------------------------------- +')

#=========================#
#  exe_uv_plot            # 
#=========================#
if (exe_uv_plot=='y'):
    print('[message] exe_uv_plot has started')
    save_path = '../fig/{}'.format(save_name)
    os.makedirs(save_path, exist_ok=True) # mkdir
    uv_plot.main(kernel_switch, SPH_model, N_sys, Delta_t_eff, WL_thick, last_file, SP_mass, h, Re, parent_path, save_path)
    print('[message] exe_uv_plot has finished')
    print('+ -------------------------------------------------------- +')

#=========================#
#  exe_p_time             # 
#=========================#
if (exe_p_time=='y'):
    print('[message] exe_p_time has started')
    save_path = '../fig/{}'.format(save_name)
    os.makedirs(save_path, exist_ok=True) # mkdir
    p_vs_time.main(Re, WL_thick, SP_mass, kernel_switch, h, N_sys, N_file, last_file, write_step, Delta_t_eff, parent_path, save_path)
    print('')
    print('[message] exe_p_time has finished')
    print('+ -------------------------------------------------------- +')

#=========================#
#  exe_energy             # 
#=========================#
if (exe_energy=='y'):
    print('[message] exe_energy has started')
    save_path = '../fig/{}'.format(save_name)
    os.makedirs(save_path, exist_ok=True) # mkdir
    energy_decay.main(Re, N_sys, N_inn, N_file, last_file, write_step, Delta_t_eff, parent_path, save_path)
    print('[message] exe_energy has finished')
    print('+ -------------------------------------------------------- +')

#=========================#
#  exe_snapshots          # 
#=========================#
if ((exe_fig_u=='y') or (exe_fig_v=='y') or (exe_fig_UV=='y') or (exe_fig_rho=='y') or (exe_fig_pre=='y')):
    print('[message] exe_snapshots has started')
    snap.main(N_sys, N_inn, WL_thick, N_file, last_file, write_step, \
              Delta_t_eff, U_top, rho_ref, Re, parent_path, save_name, \
              exe_fig_u, exe_fig_v, exe_fig_UV, exe_fig_rho, exe_fig_pre, add_wall)

print('+ -------------------------------------------------------- +')
end_time = time.perf_counter()  # cpu time
print('[message] program has finished: {:.2f} [s]'.format(end_time - start_time))
print('')

# END #