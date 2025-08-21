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
def main(kernel_switch, SPH_model, N_sys, Delta_t_eff, WL_thick, last_file, SP_mass, h, Re, parent_path, save_path):
    
    #=========================#
    #  read                   # 
    #=========================#
    file_name = parent_path + '/data/{}/SP_xyz.dat'.format(last_file) # SP_xyz
    SP_xyz    = read_f90.binary_files(file_name, N_sys, 2)
    SP_xyz    = SP_xyz[:,:] - WL_thick                                # translation
    file_name = parent_path + '/data/{}/SP_uvw.dat'.format(last_file) # SP_uvw
    SP_uvw    = read_f90.binary_files(file_name, N_sys, 2)            # already normalized
    file_name = parent_path + '/data/{}/SP_rho.dat'.format(last_file) # SP_rho
    SP_rho    = read_f90.binary_files(file_name, N_sys, 1)            # already normalized
    file_name = parent_path + '/data/{}/SP_pre.dat'.format(last_file) # SP_pre
    SP_pre    = read_f90.binary_files(file_name, N_sys, 1)            # already normalized

    #=========================#
    #  calculation point      # 
    #=========================#
    N_cal = 50
    MP_xyz = np.zeros((2*N_cal, 2)) # (name, (x,y)), :N_cal are x=0.5, N_cal: are y =0.5
    MP_uvw = np.zeros((2*N_cal, 2)) # (name, (x,y))
    MP_xyz[:N_cal, 0] = 0.5
    MP_xyz[N_cal:, 1] = 0.5

    d_MP = 1.0 / (float(N_cal) - 1.0) # delta_x for MP
    for i in range(1, N_cal):
        MP_xyz[i,1] = MP_xyz[i-1, 1] + d_MP
    for i in range(N_cal+1, 2*N_cal):
        MP_xyz[i,0] = MP_xyz[i-1, 0] + d_MP

    #=========================#
    #  interpolation          # 
    #=========================#
    MP_KC = np.zeros((2*N_cal)) # Kernel Correction for CSPH
    for me in range(2*N_cal):
        if (SPH_model == 2):
            a_vec = np.zeros((3))
            M_mat = np.zeros((3,3))
            u_vec = np.zeros((3))
            v_vec = np.zeros((3))

        for you in range(N_sys):
            x_ji = SP_xyz[you,0] - MP_xyz[me,0]
            y_ji = SP_xyz[you,1] - MP_xyz[me,1]
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
            # SPH
            if (SPH_model==0):
                MP_uvw[me,:] += V_j * W_ij * SP_uvw[you,:]
            # CSPH
            elif (SPH_model==1):
                MP_uvw[me,:] += V_j * W_ij * SP_uvw[you,:]
                MP_KC[me]    += V_j * W_ij
            # LSSPH
            elif (SPH_model==2):
                a_vec[0] = 1.0
                a_vec[1] = x_ji / h
                a_vec[2] = y_ji / h
                for i in range(3):
                    for j in range(3):
                        M_mat[i,j] += V_j * W_ij * a_vec[i] * a_vec[j]
                    u_vec[i] += V_j * W_ij * a_vec[i] * SP_uvw[you,0]
                    v_vec[i] += V_j * W_ij * a_vec[i] * SP_uvw[you,1]
        
        if (SPH_model==2):
            ans = np.linalg.solve(M_mat, u_vec)
            MP_uvw[me,0] = ans[0]
            ans = np.linalg.solve(M_mat, v_vec)
            MP_uvw[me,1] = ans[0]

    if (SPH_model==1):
        MP_uvw[:,0] /= MP_KC[:]
        MP_uvw[:,1] /= MP_KC[:]

    #=========================#
    #  figure                 # 
    #=========================#
    t_ana = Delta_t_eff* last_file
    ana_X = np.linspace(0, 1, 1000)
    ana_U = np.zeros((1000, 2))  # (points, (u,v))
    ana_U[:,0] =  np.exp(-2* math.pi**2* t_ana/ Re)* np.sin(math.pi* 0.5)* np.cos(math.pi* ana_X[:])  # at x=0.5
    ana_U[:,1] = -np.exp(-2* math.pi**2* t_ana/ Re)* np.cos(math.pi* ana_X[:])* np.sin(math.pi* 0.5)  # at y=0.5

    fig = plt.figure(figsize=(12,5.5), facecolor='white')
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.18, top=0.93, wspace=0.4, hspace=0.4)
    
    ax = fig.add_subplot(121, facecolor='white')
    # plot
    ax.scatter(MP_uvw[:N_cal, 0], MP_xyz[:N_cal, 1], label='result', c='r', marker='s', s=20)
    ax.plot(ana_U[:,0], ana_X, label='Analytical', c='k', linewidth=1.5)
    # axis label
    ax.set_xlabel('$u$', fontsize=30, labelpad=10)
    ax.set_ylabel('$y$', fontsize=30, labelpad=16)
    # ticks
    # ax.set_xticks([0.0, 0.5, 1.0])
    # ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # ax.set_xticklabels(['$0.0$', '$0.5$', '$1.0$'])
    # ax.set_yticklabels(['$0.0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'])
    ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=22)
    ax.minorticks_off() # minor ticks
    # lim
    # ax.set_xlim(0, 1.0)
    ax.set_ylim(0.0, 1.0)
    # grid
    ax.grid(which='major', color='none', linewidth=0.5)  # None
    ax.set_axisbelow(True) # grid back
    # outer frame
    ax.spines["bottom"].set_linewidth(1.2)
    ax.spines["top"].set_linewidth(1.2)
    ax.spines["right"].set_linewidth(1.2)
    ax.spines["left"].set_linewidth(1.2)
    # lagend
    ax.legend(fontsize=16, frameon=True, fancybox=True, edgecolor='silver')

    ax = fig.add_subplot(122, facecolor='white')
    # plot
    ax.scatter(MP_xyz[N_cal:, 0], MP_uvw[N_cal:, 1], label='result', c='r', marker='s', s=20)
    ax.plot(ana_X, ana_U[:,1], label='Analytical', c='k', linewidth=1.5)

    # axis label
    ax.set_xlabel('$x$', fontsize=30, labelpad=10)
    ax.set_ylabel('$v$', fontsize=30, labelpad=16)
    ax.set_xlim(0, 1)
    # ticks & lim
    # ax.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # ax.set_xticklabels(['$0.0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'])
    # ax.set_xlim(0.0, 1.0)
    # if (Re==100):
    #     ax.set_yticks([-0.4, -0.2, 0.0, 0.2])
    #     ax.set_yticklabels(['$-0.4$', '$-0.2$', '$0.0$', '$0.2$'])
    #     ax.set_ylim(-0.4, 0.2)
    # elif (Re==1000):
    #     ax.set_yticks([-0.6, -0.4, -0.2, 0.0, 0.2, 0.4])
    #     ax.set_yticklabels(['$-0.6$', '$-0.4$', '$-0.2$', '$0.0$', '$0.2$', '$0.4$'])
    #     ax.set_ylim(-0.6, 0.4)
    # elif (Re==10000):
    #     ax.set_yticks([-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6])
    #     ax.set_yticklabels(['$-0.6$', '$-0.4$', '$-0.2$', '$0.0$', '$0.2$', '$0.4$' '$0.6$'])
    #     ax.set_ylim(-0.6, 0.6)
    ax.tick_params(axis='both', which='major', direction='out', length=6, width=1, labelsize=22)
    ax.minorticks_off() # minor ticks
    # grid
    ax.grid(which='major', color='none', linewidth=0.5)  # None
    ax.set_axisbelow(True) # grid back
    # outer frame
    ax.spines["bottom"].set_linewidth(1.2)
    ax.spines["top"].set_linewidth(1.2)
    ax.spines["right"].set_linewidth(1.2)
    ax.spines["left"].set_linewidth(1.2)
    # lagend
    ax.legend(fontsize=16, frameon=True, fancybox=True, edgecolor='silver')
    # save
    fig.savefig('{}/uv_plot.png'.format(save_path), format='png', dpi=300, transparent=False)
    plt.close()

    # save
    np.savez(parent_path + '/ana_uv_plot.npz', xyz=MP_xyz, uvw=MP_uvw)  # save

# END #