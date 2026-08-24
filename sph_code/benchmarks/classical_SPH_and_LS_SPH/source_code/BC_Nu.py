#=========================#
#  module                 # 
#=========================#
import math
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import ana_read_fortran_files as read_f90
import ana_kernel_function_s as KFC
plt.rcParams['mathtext.fontset'] = 'cm' # mathfont for figure

#=========================#
#  main                   # 
#=========================#
def main(kernel_switch, SPH_model, N_sys, WL_thick, last_file, SP_mass, \
         h, temp_ave, height, Delta_temp, parent_path, save_path):
    
    if (SPH_model != 2):
        print('[error] input: SPH_model = 2')
        sys.exit()
    
    #=========================#
    #  read                   # 
    #=========================#
    file_name = parent_path + '/data/{}/SP_xyz.dat'.format(last_file) # SP_xyz
    SP_xyz    = read_f90.binary_files(file_name, N_sys, 2)
    # SP_xyz    = SP_xyz[:,:] - WL_thick                                # translation
    # SP_xyz   /= height                                                # normalize

    file_name = parent_path + '/data/{}/SP_tem.dat'.format(last_file) # SP_tem
    SP_tem    = read_f90.binary_files(file_name, N_sys, 1)
    # SP_tem[:,0] = (SP_tem[:,0] - temp_ave) / Delta_temp               # normalize
   
    file_name = parent_path + '/data/{}/SP_rho.dat'.format(last_file) # SP_rho
    SP_rho    = read_f90.binary_files(file_name, N_sys, 1)

    #=========================#
    #  calculation point      # 
    #=========================#
    N_cal = 50
    MP_xyz = np.zeros((2*N_cal, 2)) # (name, (x,y)), :N_cal are at bottom, N_cal: are at top
    MP_Nu  = np.zeros((2*N_cal   )) # (name)
    MP_xyz[:N_cal, 1] = WL_thick
    MP_xyz[N_cal:, 1] = WL_thick + height

    d_MP = height / (float(N_cal) - 1.0) # delta_x for MP
    MP_xyz[0,    0] = WL_thick
    MP_xyz[N_cal,0] = WL_thick
    for i in range(1, N_cal):
        MP_xyz[i,0] = MP_xyz[i-1, 0] + d_MP
    for i in range(N_cal+1, 2*N_cal):
        MP_xyz[i,0] = MP_xyz[i-1, 0] + d_MP

    #=========================#
    #  interpolation          # 
    #=========================#
    for me in range(2*N_cal):
        a_vec = np.zeros((6))
        M_mat = np.zeros((6,6))
        T_vec = np.zeros((6))

        for you in range(N_sys):
            x_ji = SP_xyz[you,0] - MP_xyz[me,0]
            y_ji = SP_xyz[you,1] - MP_xyz[me,1]
            r_ij = math.sqrt(x_ji**2.0 + y_ji**2.0)
            V_j  = SP_mass/SP_rho[you,0]
            if (kernel_switch == 1):
                W_ij = KFC.cubic_spline_d0W(r_ij, h)
            elif (kernel_switch == 2):
                W_ij = KFC.quintic_spline_d0W(r_ij, h)
            elif (kernel_switch == 3):
                W_ij = KFC.wendland_C2_d0W(r_ij, h)
            elif (kernel_switch == 4):
                W_ij = KFC.wendland_C4_d0W(r_ij, h)
            elif (kernel_switch == 5):
                W_ij = KFC.wendland_C6_d0W(r_ij, h)
            # LSSPH
            a_vec[0] = 1.0
            a_vec[1] = x_ji           / h
            a_vec[2] = y_ji           / h
            a_vec[3] = x_ji**2.0 /2.0 / h**2.0
            a_vec[4] = x_ji*y_ji      / h**2.0
            a_vec[5] = y_ji**2.0 /2.0 / h**2.0


            for i in range(6):
                for j in range(6):
                    M_mat[i,j] += V_j * W_ij * a_vec[i] * a_vec[j]
                T_vec[i] += V_j * W_ij * a_vec[i] * SP_tem[you,0]
        
        ans = np.linalg.solve(M_mat, T_vec)
        MP_Nu[me] = ans[2] / h

    MP_Nu /= (-Delta_temp/height)

    Nu_ave_bot = np.mean(MP_Nu[:N_cal])
    Nu_ave_top = np.mean(MP_Nu[N_cal:])
    #=========================#
    #  figure                 # 
    #=========================#
    fig = plt.figure(figsize=(14,5.5), facecolor='white')
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.18, top=0.93, wspace=0.4, hspace=0.4)
    
    ax = fig.add_subplot(121, facecolor='white')
    # plot
    ax.plot((MP_xyz[:N_cal, 0]-WL_thick)/height, MP_Nu[:N_cal], c='r', marker='s')
    # axis label
    ax.set_xlabel('$x$', fontsize=30, labelpad=10)
    ax.set_ylabel('$Nu$', fontsize=30, labelpad=16)
    # lim
    ax.set_xlim(-0.1, 1.1)
    # ticks
    ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=22)
    ax.minorticks_off() # minor ticks
    # grid
    ax.grid(which='major', color='gray', linewidth=0.5)  # None
    ax.set_axisbelow(True) # grid back
    # outer frame
    ax.spines["bottom"].set_linewidth(1.2)
    ax.spines["top"].set_linewidth(1.2)
    ax.spines["right"].set_linewidth(1.2)
    ax.spines["left"].set_linewidth(1.2)
    #title
    ax.set_title('Nu (bot)={:2.5e}'.format(Nu_ave_bot), c='k', y=1.02, fontsize=18)

    ax = fig.add_subplot(122, facecolor='white')
    # plot
    ax.plot((MP_xyz[N_cal:, 0]-WL_thick)/height, MP_Nu[N_cal:], c='r', marker='s')
    # axis label
    ax.set_xlabel('$x$', fontsize=30, labelpad=10)
    ax.set_ylabel('$Nu$', fontsize=30, labelpad=16)
     # lim
    ax.set_xlim(-0.1, 1.1)
    # ticks
    ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=22)
    ax.minorticks_off() # minor ticks
    # grid
    ax.grid(which='major', color='gray', linewidth=0.5)  # None
    ax.set_axisbelow(True) # grid back
    # outer frame
    ax.spines["bottom"].set_linewidth(1.2)
    ax.spines["top"].set_linewidth(1.2)
    ax.spines["right"].set_linewidth(1.2)
    ax.spines["left"].set_linewidth(1.2)
    #title
    ax.set_title('Nu (top)={:2.5e}'.format(Nu_ave_top), c='k', y=1.02, fontsize=18)
    # save
    fig.savefig('{}/Nu.png'.format(save_path), format='png', dpi=300, transparent=False)
    plt.close()

    # save
    np.savez(parent_path + '/Nu.npz', xyz=MP_xyz, Nu=MP_Nu)         # save

# END #