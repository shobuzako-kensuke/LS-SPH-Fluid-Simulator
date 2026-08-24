#=========================#
#  module                 # 
#=========================#
import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import ana_read_fortran_files as read_f90
import ana_kernel_function_s as KFC
plt.rcParams['mathtext.fontset'] = 'cm' # mathfont for figure

#=========================#
#  main                   # 
#=========================#
def main(Re, WL_thick, SP_mass, kernel_switch, h, N_sys, N_file, last_file, write_step, Delta_t_eff, parent_path, save_path):
    count  = 0
    time   = np.zeros((N_file)) # for time record
    p_time = np.zeros((N_file)) # for p_time record

    #=========================#
    #  read & cal             # 
    #=========================#
    for t in range(0, last_file+1, write_step):
        time[count] = Delta_t_eff * float(t)
        
        file_name = parent_path + '/data/{}/SP_xyz.dat'.format(t) # SP_xyz
        SP_xyz    = read_f90.binary_files(file_name, N_sys, 2)
        SP_xyz    = SP_xyz[:,:] - WL_thick                        # translation
        file_name = parent_path + '/data/{}/SP_rho.dat'.format(t) # SP_rho
        SP_rho    = read_f90.binary_files(file_name, N_sys, 1)    # already normalized
        file_name = parent_path + '/data/{}/SP_pre.dat'.format(t) # SP_pre
        SP_pre    = read_f90.binary_files(file_name, N_sys, 1)    # already normalized
        
        #=========================#
        #  interpolation          # 
        #=========================#
        a_vec = np.zeros((3))
        M_mat = np.zeros((3,3))
        p_vec = np.zeros((3))

        for you in range(N_sys):
            x_ji = SP_xyz[you,0] - 0.5
            y_ji = SP_xyz[you,1] - 0.5
            r_ij = math.sqrt(x_ji**2.0 + y_ji**2.0)
            V_j  = SP_mass/SP_rho[you]
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
            a_vec[1] = x_ji / h
            a_vec[2] = y_ji / h
            for i in range(3):
                for j in range(3):
                    M_mat[i,j] += V_j * W_ij * a_vec[i] * a_vec[j]
                p_vec[i] += V_j * W_ij * a_vec[i] * SP_pre[you,0]

        ans = np.linalg.solve(M_mat, p_vec)
        p_time[count] = ans[0]

        count += 1

        #=========================#
        #  progress               # 
        #=========================# 
        if (i==0):
            print('          {}/{} step is completed'.format(count, N_file), end='')
        else:
            print('\r          {}/{} step is completed'.format(count, N_file), end='')

    #=========================#
    #  figure                 # 
    #=========================#
    ana_p = np.zeros((N_file))
    ana_p[:] = -0.5*np.exp(-4*np.pi**2/Re*time[:])

    fig = plt.figure(figsize=(6.7,6), facecolor='white')
    plt.subplots_adjust(left=0.25, right=0.95, bottom=0.18, top=0.95, wspace=0.35, hspace=0.4)
    ax = fig.add_subplot(111, facecolor='white')
    # plot
    ax.plot(time, p_time, c='r', label='result', linewidth=1.5)
    ax.plot(time, ana_p, c='k', label='Analytical', linewidth=1.5)
    # axis label
    ax.set_xlabel('$t$', fontsize=30, labelpad=10)
    ax.set_ylabel('$p$', fontsize=30, labelpad=16)
    # ticks
    ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=22)
    ax.minorticks_off() # minor ticks
    # grid
    ax.grid(which='major', color='none', linewidth=0.5) # None
    ax.set_axisbelow(True) # grid back
    # outer frame
    ax.spines["bottom"].set_linewidth(1.2)
    ax.spines["top"].set_linewidth(1.2)
    ax.spines["right"].set_linewidth(1.2)
    ax.spines["left"].set_linewidth(1.2)
    # lagend
    ax.legend(fontsize=16, frameon=True, fancybox=True, edgecolor='silver')
    # save
    fig.savefig('{}/p_vs_time.png'.format(save_path), format='png', dpi=300, transparent=False)
    plt.close()

    # save
    np.savez(parent_path + '/ana_p_vs_time.npz', t=time, p=p_time)  # save

# END #