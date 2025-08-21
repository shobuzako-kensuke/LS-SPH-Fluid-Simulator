#=========================#
#  file name              # 
#=========================#
save_name = 'BC_test1'
#=========================#
#  analysis file          # 
#=========================#
analysis_step = 'n'  # if specifying analysis step, input file ID
#=========================#
#  program switch         # 
#=========================#
exe_fig_ini = 'n'    # initial settings
exe_V_rms   = 'y'    # time vs V_rms
exe_Nu      = 'y'    # Nusselt number
exe_Mach    = 'y'    # time vs Mach number
exe_fig_u   = 'n'    # snapshot of u
exe_fig_v   = 'n'    # snapshot of v
exe_fig_UV  = 'n'    # snapshot of speed
exe_fig_rho = 'n'    # snapshot of density
exe_fig_pre = 'n'    # snapshot of pressure
exe_fig_tem = 'y'    # snapshot of temperature
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
import BC_initial_settings as ini
import BC_V_rms as V_rms
import BC_Nu as Nu
import BC_Mach_number as mach
import BC_snapshot as snap

start_time = time.perf_counter() # cpu time
print('+ -------------------------------------------------------- +')
print('[message] main.py has started')

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
    ini.xyz_kind(N_sys, N_inn, WL_thick, width, parent_path, save_path)
    print('          >> xyz_kind')
    ini.cell(N_sys, N_inn, WL_thick, width, h_eff, parent_path, save_path)
    print('          >> cell')
    ini.xyz_kind_VM(N_VM, N_inn, WL_thick, width, parent_path, save_path)
    print('          >> xyz_kind_VM')
    ini.cell_VM(N_VM, N_inn, WL_thick, width, h_eff, parent_path, save_path)
    print('          >> cell_VM')
    print('[message] exe_fig_ini has finished')
    print('+ -------------------------------------------------------- +')

#=========================#
#  exe_V_rms              # 
#=========================#
if (exe_V_rms=='y'):
    print('[message] exe_V_rms has started')
    save_path = '../fig/{}'.format(save_name)
    os.makedirs(save_path, exist_ok=True) # mkdir
    V_rms.main(N_sys, N_inn, N_file, last_file, write_step, Delta_t_eff, \
               Pr, gravity, alpha, height, Delta_temp, thermal_dif, parent_path, save_path)
    print('[message] exe_V_rms has finished')
    print('+ -------------------------------------------------------- +')

#=========================#
#  exe_Nu                 # 
#=========================#
if (exe_Nu=='y'):
    print('[message] exe_Nu has started')
    save_path = '../fig/{}'.format(save_name)
    os.makedirs(save_path, exist_ok=True) # mkdir
    Nu.main(kernel_switch, SPH_model, N_sys, WL_thick, last_file, SP_mass, \
            h, temp_ave, height, Delta_temp, parent_path, save_path)
    print('[message] exe_Nu has finished')
    print('+ -------------------------------------------------------- +')

#=========================#
#  exe_Mach               # 
#=========================#
if (exe_Mach=='y'):
    print('[message] exe_Mach has started')
    save_path = '../fig/{}'.format(save_name)
    os.makedirs(save_path, exist_ok=True) # mkdir
    mach.main(N_sys, N_inn, N_file, last_file, write_step, Delta_t_eff, \
                sound_speed, zeta, xi, parent_path, save_path)
    print('[message] exe_Mach has finished')
    print('+ -------------------------------------------------------- +')

#=========================#
#  exe_snapshots          # 
#=========================#
if ((exe_fig_u=='y') or (exe_fig_v=='y') or (exe_fig_UV=='y') or (exe_fig_rho=='y') or (exe_fig_pre=='y') or (exe_fig_tem=='y')):
    print('[message] exe_snapshots has started')
    snap.main(N_sys, N_inn, WL_thick, N_file, last_file, write_step, \
              Delta_t_eff, rho_ref, alpha, Delta_temp, temp_ave, gravity, thermal_dif, height, \
              Ra, Pr, parent_path, save_name, \
              exe_fig_u, exe_fig_v, exe_fig_UV, exe_fig_rho, exe_fig_pre, exe_fig_tem, add_wall)

print('+ -------------------------------------------------------- +')
end_time = time.perf_counter()  # cpu time
print('[message] program has finished: {:.2f} [s]'.format(end_time - start_time))
print('')

# END #