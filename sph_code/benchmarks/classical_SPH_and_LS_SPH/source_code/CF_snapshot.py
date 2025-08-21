#=========================#
#  module                 # 
#=========================#
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import ana_read_fortran_files as read_f90
import ana_mk_movie as mk_movie
plt.rcParams['mathtext.fontset'] = 'cm' # mathfont for figure

#=========================#
#  fig_cbar_yticks        # 
#=========================#
def fig_cbar_yticks(Re, val, cbar, ax):
    if (Re==100):
        if (val=='u'):
            cbar.ax.set_yticks([-1, -0.5, 0, 0.5, 1])
            cbar.ax.set_yticklabels(['$-1.0$', '$-0.5$', '$0.0$', '$0.5$', '$1.0$'])
        if (val=='v'):
            cbar.ax.set_yticks([-1, -0.5, 0, 0.5, 1])
            cbar.ax.set_yticklabels(['$-1.0$', '$-0.5$', '$0.0$', '$0.5$', '$1.0$'])
        if (val=='UV'):
            cbar.ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
            cbar.ax.set_yticklabels(['$0.0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'])
        if (val=='rho'):
            cbar.ax.set_yticks([0.9, 1, 1.1])
            cbar.ax.set_yticklabels(['$0.9$', '$1.0$', '$1.1$'])
        if (val=='pre'):
            cbar.ax.set_yticks([-0.3, -0.15, 0, 0.15, 0.3])
            cbar.ax.set_yticklabels(['$-0.3$', '$-0.15$', '$0.0$', '$0.15$', '$0.3$'])

    if (Re==1000):
        if (val=='u'):
            cbar.ax.set_yticks([-1, -0.5, 0, 0.5, 1])
            cbar.ax.set_yticklabels(['$-1.0$', '$-0.5$', '$0.0$', '$0.5$', '$1.0$'])
        if (val=='v'):
            cbar.ax.set_yticks([-1, -0.5, 0, 0.5, 1])
            cbar.ax.set_yticklabels(['$-1.0$', '$-0.5$', '$0.0$', '$0.5$', '$1.0$'])
        if (val=='UV'):
            cbar.ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
            cbar.ax.set_yticklabels(['$0.0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1.0$'])
        if (val=='rho'):
            cbar.ax.set_yticks([0.9, 1, 1.1])
            cbar.ax.set_yticklabels(['$0.9$', '$1.0$', '$1.1$'])
        if (val=='pre'):
            cbar.ax.set_yticks([-0.5, -0.25, 0, 0.25, 0.5])
            cbar.ax.set_yticklabels(['$-0.5$', '$-0.25$', '$0.0$', '$0.25$', '$0.5$'])

#=========================#
#  fig_detail             # 
#=========================#
def fig_detail(Re, val, fig, lim_switch, boudary_switch, time, step, SP_xyz, SP_val, val_min, val_max, sca_size):
    ax = fig.add_subplot(111, facecolor='white')
    cfig = ax.scatter(SP_xyz[:,0], SP_xyz[:,1], c=SP_val, marker='.', cmap='jet', \
                      vmin=val_min, vmax=val_max, s=sca_size)
    # color bar
    cbar = plt.colorbar(cfig, aspect=30, shrink=0.9, ax=ax, orientation='vertical', pad=0.05, location='right')
    cbar.ax.tick_params(direction='out', length=4, width=1, labelsize=20)
    fig_cbar_yticks(Re, val, cbar, ax)
    cbar.ax.set_ylim(val_min, val_max)
    # axis label
    ax.set_xlabel('$x$', fontsize=30, labelpad=10)
    ax.set_ylabel('$y$', fontsize=30, labelpad=16)
    # lim
    if (lim_switch == 1):
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        #ax.axis('equal')
    # boundary
    if (boudary_switch == 1):
        ax.vlines(x=[0, 1], ymin=0, ymax=1, color='k', linestyle='--')
        ax.hlines(y=[0, 1], xmin=0, xmax=1, color='k', linestyle='--')
    # ticks
    ax.set_xticks([0.0, 0.5, 1.0])
    ax.set_yticks([0.0, 0.5, 1.0])
    ax.set_xticklabels(['$0.0$', '$0.5$', '$1.0$'])
    ax.set_yticklabels(['$0.0$', '$0.5$', '$1.0$'])
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
    # title
    ax.set_title('{:3.2e} s  ({} step)'.format(time, step), c='k', y=1.02, fontsize=18)

#=========================#
#  main                   # 
#=========================#
def main(N_sys, N_inn, WL_thick, N_file, last_file, write_step, \
         Delta_t_eff, U_top, rho_ref, Re, parent_path, save_name, \
         exe_fig_u, exe_fig_v, exe_fig_UV, exe_fig_rho, exe_fig_pre, add_wall):
    
    #=========================#
    #  make snapshots         # 
    #=========================#
    count = 0
    for i in range(0, last_file+1, write_step):
        count += 1
        #=========================#
        #  read                   # 
        #=========================#   
        time  = Delta_t_eff * float(i) # time in simulation

        file_name = parent_path + '/data/{}/SP_xyz.dat'.format(i)      # SP_xyz
        SP_xyz    = read_f90.binary_files(file_name, N_sys, 2)
        SP_xyz    = SP_xyz[:,:] - WL_thick                             # translation

        if ((exe_fig_u=='y') or (exe_fig_v=='y') or (exe_fig_UV=='y')):
            file_name = parent_path + '/data/{}/SP_uvw.dat'.format(i)  # read SP_uvw
            SP_uvw    = read_f90.binary_files(file_name, N_sys, 2)
            SP_uvw   /= U_top                                          # normalized by U_top
            if (count==1):
                if (exe_fig_u=='y'):
                    save_path = '../fig/{}/snapshot_u'.format(save_name)
                    os.makedirs(save_path, exist_ok=True)              # mkdir for u
                if (exe_fig_v=='y'):
                    save_path = '../fig/{}/snapshot_v'.format(save_name)
                    os.makedirs(save_path, exist_ok=True)              # mkdir for v
                if (exe_fig_UV=='y'):
                    save_path = '../fig/{}/snapshot_UV'.format(save_name)
                    os.makedirs(save_path, exist_ok=True)              # mkdir for UV

        if (exe_fig_rho=='y'):
            file_name = parent_path + '/data/{}/SP_rho.dat'.format(i)  # SP_rho
            SP_rho    = read_f90.binary_files(file_name, N_sys, 1)
            SP_rho   /= rho_ref                                        # normalized by rho_ref
            if (count==1):
                save_path = '../fig/{}/snapshot_rho'.format(save_name)
                os.makedirs(save_path, exist_ok=True)                  # mkdir for rho

        if (exe_fig_pre=='y'):
            file_name = parent_path + '/data/{}/SP_pre.dat'.format(i)  # SP_pre
            SP_pre    = read_f90.binary_files(file_name, N_sys, 1)
            SP_pre   /= (rho_ref*U_top**2.0)                          # normalized by (rho_ref*U_top**2)
            if (count==1):
                save_path = '../fig/{}/snapshot_pre'.format(save_name)
                os.makedirs(save_path, exist_ok=True)                  # mkdir for pre

        #=========================#
        #  scatter size           # 
        #=========================#   
        if (add_wall=='y'):
            if (N_inn==50**2):
                sca_size = 80
            elif (N_inn==100**2):
                sca_size = 15
            elif (N_inn==200**2):
                sca_size = 1
        else:
            if (N_inn==50**2):
                sca_size = 110
            elif (N_inn==100**2):
                sca_size = 25
            elif (N_inn==200**2):
                sca_size = 3

        #=========================#
        #  min/max values         # 
        #=========================#
        if (Re==100):
            u_min   = - 1.0
            u_max   = 1.0
            v_min   = - 1.0
            v_max   = 1.0
            UV_min  = 0.0
            UV_max  = 1.0
            rho_min = 0.9
            rho_max = 1.1
            pre_min = - 0.3
            pre_max = 0.3

        if (Re==1000):
            u_min   = - 1.0
            u_max   = 1.0
            v_min   = - 1.0
            v_max   = 1.0
            UV_min  = 0.0
            UV_max  = 1.0
            rho_min = 0.9
            rho_max = 1.1
            pre_min = - 0.3
            pre_max = 0.3

        #=========================#
        #  figure                 # 
        #=========================#   
        if (exe_fig_u=='y'):
            fig = plt.figure(figsize=(7.5,6.5), facecolor='white')
            plt.subplots_adjust(left=0.18, right=0.95, bottom=0.18, top=0.90, wspace=0.35, hspace=0.4)
            # details
            if (add_wall=='y'):
                fig_detail(Re, 'u', fig, 0, 1, time, i, SP_xyz, SP_uvw[:,0], u_min, u_max, sca_size)
            else:
                fig_detail(Re, 'u', fig, 1, 0, time, i, SP_xyz, SP_uvw[:,0], u_min, u_max, sca_size)
            # save
            fig.savefig('../fig/{}/snapshot_u/{}.png'.format(save_name, i), format='png', dpi=300, transparent=False)
            plt.close()

        if (exe_fig_v=='y'):
            fig = plt.figure(figsize=(7.5,6.5), facecolor='white')
            plt.subplots_adjust(left=0.18, right=0.95, bottom=0.18, top=0.90, wspace=0.35, hspace=0.4)
            # details
            if (add_wall=='y'):
                fig_detail(Re, 'v', fig, 0, 1, time, i, SP_xyz, SP_uvw[:,1], v_min, v_max, sca_size)
            else:
                fig_detail(Re, 'v', fig, 1, 0, time, i, SP_xyz, SP_uvw[:,1], v_min, v_max, sca_size)
            # save
            fig.savefig('../fig/{}/snapshot_v/{}.png'.format(save_name, i), format='png', dpi=300, transparent=False)
            plt.close()

        if (exe_fig_UV=='y'):
            fig = plt.figure(figsize=(7.5,6.5), facecolor='white')
            plt.subplots_adjust(left=0.18, right=0.95, bottom=0.18, top=0.90, wspace=0.35, hspace=0.4)
            SP_UV = np.sqrt(SP_uvw[:,0]**2.0 + SP_uvw[:,1]**2.0)
            # details
            if (add_wall=='y'):
                fig_detail(Re, 'UV', fig, 0, 1, time, i, SP_xyz, SP_UV, UV_min, UV_max, sca_size)
            else:
                fig_detail(Re, 'UV', fig, 1, 0, time, i, SP_xyz, SP_UV, UV_min, UV_max, sca_size)
            # save
            fig.savefig('../fig/{}/snapshot_UV/{}.png'.format(save_name, i), format='png', dpi=300, transparent=False)
            plt.close()

        if (exe_fig_rho=='y'):
            fig = plt.figure(figsize=(7.5,6.5), facecolor='white')
            plt.subplots_adjust(left=0.18, right=0.95, bottom=0.18, top=0.90, wspace=0.35, hspace=0.4)
            # details
            if (add_wall=='y'):
                fig_detail(Re, 'rho', fig, 0, 1, time, i, SP_xyz, SP_rho, rho_min, rho_max, sca_size)
            else:
                fig_detail(Re, 'rho', fig, 1, 0, time, i, SP_xyz, SP_rho, rho_min, rho_max, sca_size)
            # save
            fig.savefig('../fig/{}/snapshot_rho/{}.png'.format(save_name, i), format='png', dpi=300, transparent=False)
            plt.close()

        if (exe_fig_pre=='y'):
            fig = plt.figure(figsize=(7.5,6.5), facecolor='white')
            plt.subplots_adjust(left=0.18, right=0.95, bottom=0.18, top=0.90, wspace=0.35, hspace=0.4)
            # details
            if (add_wall=='y'):
                fig_detail(Re, 'pre', fig, 0, 1, time, i, SP_xyz, SP_pre, pre_min, pre_max, sca_size)
            else:
                fig_detail(Re, 'pre', fig, 1, 0, time, i, SP_xyz, SP_pre, pre_min, pre_max, sca_size)
            # save
            fig.savefig('../fig/{}/snapshot_pre/{}.png'.format(save_name, i), format='png', dpi=300, transparent=False)
            plt.close()
        
        #=========================#
        #  progress               # 
        #=========================# 
        if (i==0):
            print('          >> {:.2f} %'.format(float(count)/float(N_file)*100), end='')
        else:
            print('\r          >> {:.2f} %'.format(float(count)/float(N_file)*100), end='')

    #=========================#
    #  make movie             #
    #=========================#
    print('')
    if (exe_fig_u=='y'):
        print('          making a movie for u ... ', end='', flush=True)
        mk_movie.func_animation(last_file, write_step, save_name, 'u')
        print('finish')
    if (exe_fig_v=='y'):
        print('          making a movie for v ... ', end='', flush=True)
        mk_movie.func_animation(last_file, write_step, save_name, 'v')
        print('finish')
    if (exe_fig_UV=='y'):
        print('          making a movie for UV ... ', end='', flush=True)
        mk_movie.func_animation(last_file, write_step, save_name, 'UV')
        print('finish')
    if (exe_fig_rho=='y'):
        print('          making a movie for rho ... ', end='', flush=True)
        mk_movie.func_animation(last_file, write_step, save_name, 'rho')
        print('finish')
    if (exe_fig_pre=='y'):
        print('          making a movie for pre ... ', end='', flush=True)
        mk_movie.func_animation(last_file, write_step, save_name, 'pre')
        print('finish')
        
# END #