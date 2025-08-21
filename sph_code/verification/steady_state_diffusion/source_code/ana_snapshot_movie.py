#=========================#
#  module                 #
#=========================#
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ana_read_fortran as read_f90
import ana_mk_movie as mk_movie
plt.rcParams['mathtext.fontset'] = 'cm' # mathfont for figure

#=========================#
#  fig_cmap               #
#=========================#
def fig_cmap(tmp_name, file_number, last_file, N_file, dt, N_sys, N_inn, WL_thick, err_min, err_max, mark_size, parent_path, save_name):

    #=========================#
    #  make directory         #
    #=========================#
    save_path = '../fig/{}/snapshot_{}'.format(save_name, tmp_name)
    os.makedirs(save_path, exist_ok=True)

    #=========================#
    #  make snapshots         #
    #=========================#
    count = 0

    file_name = parent_path + '/data/0/SP_x.dat'
    SP_x = read_f90.binary_files(file_name, N_sys, 2)
    SP_x = SP_x - WL_thick

    for i in file_number:
        count += 1
        #=========================#
        #  read                   #
        #=========================#   
        time  = dt * float(i)  # simulation time

        file_name = parent_path + '/data/{}/SP_f.dat'.format(i)
        SP_f = read_f90.binary_files(file_name, N_sys, 1)

        SP_err = np.abs(SP_f[:N_inn,0] - np.sin(math.pi* SP_x[:N_inn,0]) * np.cos(math.pi* SP_x[:N_inn,1]))

        #=========================#
        #  figure                 #
        #=========================#
        fig = plt.figure(figsize=(12,6), facecolor='white', constrained_layout=True)        
        ax1 = fig.add_subplot(121, facecolor='white')
        ax2 = fig.add_subplot(122, facecolor='white')

        #=========================#
        #  ax1                    #
        #=========================#
        cfig = ax1.scatter(SP_x[:N_inn, 0], SP_x[:N_inn, 1], c=SP_f[:N_inn], ec='k', marker='.', \
                           cmap='jet', s=mark_size, linewidth=0.1, vmin=-1, vmax=1)
        cbar = plt.colorbar(cfig, aspect=30, shrink=0.7, ax=ax1, orientation='vertical', pad=0.05, location='right', \
                            ticks=[-1.0, -0.5, 0, 0.5, 1.0])
    
        # color bar
        cbar.ax.tick_params(direction='out', length=4, width=1, labelsize=14)
        cbar.ax.set_ylim(-1, 1)

        # text
        ax1.text(0.04, -0.3, '(a) Calculated function value $f(x,y)$', fontsize=16)

        #=========================#
        #  ax2                    #
        #=========================#
        cfig = ax2.scatter(SP_x[:N_inn, 0], SP_x[:N_inn, 1], c=SP_err, ec='k', marker='.', \
                           cmap='binary', s=mark_size, linewidth=0.1, \
                           norm=colors.LogNorm(vmin=err_min, vmax=err_max))
        cbar = plt.colorbar(cfig, aspect=30, shrink=0.7, ax=ax2, orientation='vertical', pad=0.05, location='right')

        # color bar
        cbar.ax.tick_params(which='major', direction='out', length=4, width=1, labelsize=14)
        cbar.ax.tick_params(which='minor', direction='out')
        cbar.ax.set_ylim(err_min, err_max)

        # text
        ax2.text(0.25, -0.3, '(b) Absolute error', fontsize=16)

        ax_list = [ax1, ax2]
        for ax in ax_list:
            # axis label
            ax.set_xlabel(r'$x$', fontsize=20, labelpad=10)
            ax.set_ylabel(r'$y$', fontsize=20, labelpad=16)
            
            # lim
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_aspect('equal', 'box')

            # ticks
            ax.tick_params(axis='both', which='major', direction='out', length=4, width=1, labelsize=14)
            ax.minorticks_off()

            # grid
            ax.grid(which='major', color='none', linewidth=0.5)
            ax.set_axisbelow(True) # grid back

            # title
            ax.set_title('{:3.2e} [s]  ({} step)'.format(time, i), c='k', y=1.02, fontsize=14)

        # save
        fig.savefig('../fig/{}/snapshot_{}/{}.png'.format(save_name, tmp_name, i), format='png', dpi=300, transparent=False)
        plt.close()
        #=========================#
        #  progress               #
        #=========================#
        if (i==0):
            print('          >> {:.2f} %'  .format(float(count)/float(N_file)*100), end='')
        else:
            print('\r          >> {:.2f} %'.format(float(count)/float(N_file)*100), end='')

    #=========================#
    #  make movie             #
    #=========================#
    print('')
    print('          Movie is being made ... ', end='', flush=True)
    mk_movie.func_animation(file_number, save_name, tmp_name)
    print('finish')

# END #