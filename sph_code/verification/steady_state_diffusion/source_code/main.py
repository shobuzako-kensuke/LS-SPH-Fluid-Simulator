#=========================#
#  file name              #
#=========================#
save_name = 'test'
#=========================#
#  program switch         #
#=========================#
fig_ini = False   # figures for initial settings
fig_f   = True   # snapshots of function value f

# メモ：解析解の図，最終的な解析解とのL2ノルム計測（内部と境界で分ける）

#=========================#
#  scatter size           #
#=========================#
mark_size = 40
err_min   = 1e-3
err_max   = 1

#==============================================================================#
#=========================== NOT CHANGE BELOW =================================#
#==============================================================================#

#=========================#
#  module                 #
#=========================#
import os
import glob
import time
import ana_read_fortran as read_f90
import ana_initial_setting as ini
import ana_snapshot_movie as snap_movie

start_time = time.perf_counter()
print('+ -------------------------------------------------------- +')

#=========================#
#  read variables         #
#=========================#
parent_path   = '../output/{}'.format(save_name)                     # path to the saved data
variable_path = parent_path + '/system_info/system_info_VALUES.dat'  # file path
 
write_step, Nx, Nh, rho_ref, x_rand, threshold, \
D, m, h, h_eff, WL_thick, dt, N_WL, N_inn, N_out, N_sys \
= read_f90.read_variables(variable_path)

#=========================#
#  file count             #
#=========================#
file_list = glob.glob(parent_path + '/data/*/SP_f.dat')  # get the number of files
N_file = len(file_list)  # total number of files

file_number = []
for i in range(N_file):
    parts = file_list[i].split('/')    # split by '/'
    file_number.append(int(parts[4]))  # append each file number

file_number.sort()                     # sort
last_file = max(file_number)           # last file number

print('[message] Number of read files: {}'.format(N_file))
print('[message] Last file number    : {}'.format(last_file))
print('+ -------------------------------------------------------- +')
#==============================================================================#
#                            main program below                                #
#==============================================================================#

#=========================#
#  initial settings       #
#=========================#
if fig_ini:
    print('[message] Figures for initial settings are being made.')
    save_path = '../fig/{}'.format(save_name)
    os.makedirs(save_path, exist_ok=True)
    ini.main(N_inn, N_out, N_sys, WL_thick, D, mark_size, parent_path, save_path)
    print('[message] Figures for initial settings have been made.')
    print('+ -------------------------------------------------------- +')

#=========================#
#  exe_snapshots          #
#=========================#
if fig_f:
    print('[message] Snapshots of function value f are being made.')
    tmp_name = 'f'
    snap_movie.fig_cmap(tmp_name, file_number, last_file, N_file, dt, N_sys, N_inn, WL_thick, err_min, err_max, mark_size, parent_path, save_name)
    
print('+ -------------------------------------------------------- +')
end_time = time.perf_counter()
print('[message] Program has finished: {:.2f} [s]'.format(end_time - start_time))
print('')

# END #