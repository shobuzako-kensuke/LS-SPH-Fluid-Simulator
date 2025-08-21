#=========================#
#  module                 #
#=========================#
import numpy as np
import matplotlib.pyplot as plt
import ana_read_fortran as read_f90
plt.rcParams['mathtext.fontset'] = 'cm' # mathfont for figure

#=========================#
#  main                   #
#=========================#
def main(N_inn, N_out, N_sys, WL_thick, D, mark_size, parent_path, save_path):

    #=========================#
    #  read                   #
    #=========================#
    file_name = parent_path + '/data/0/SP_x.dat'
    SP_x = read_f90.binary_files(file_name, N_sys, 2)
    SP_x = SP_x - WL_thick

    file_name = parent_path + '/data/0/SP_kind.dat'
    SP_kind = read_f90.binary_files(file_name, N_sys, 1)
    SP_kind = SP_kind.astype(np.int64)  # float >> int

    file_name = parent_path + '/data/0/VM_x.dat'
    VM_x = read_f90.binary_files(file_name, N_out, 2)
    VM_x = VM_x - WL_thick

    file_name = parent_path + '/data/0/VM_kind.dat'
    VM_kind = read_f90.binary_files(file_name, N_out, 1)
    VM_kind = VM_kind.astype(np.int64)  # float >> int
    
    #==============================================================================#
    #  figure 1                                                                    #
    #==============================================================================#
    fig = plt.figure(figsize=(12,11), facecolor='white', constrained_layout=True)
    ax1 = fig.add_subplot(221, facecolor='white')
    ax2 = fig.add_subplot(222, facecolor='white')
    ax3 = fig.add_subplot(223, facecolor='white')
    ax4 = fig.add_subplot(224, facecolor='white')

    #=========================#
    #  ax1                    #
    #=========================#
    ax1.scatter(SP_x[:N_inn, 0], SP_x[:N_inn, 1], c='w', ec='dimgray', marker='.', s=mark_size)

    for i in range(N_inn, N_sys):
        if (SP_kind[i] == 1):
            ax1.scatter(SP_x[i,0], SP_x[i,1], c='w', ec='r', marker='.', s=mark_size, label='bottom')
        elif (SP_kind[i] == 2):
            ax1.scatter(SP_x[i,0], SP_x[i,1], c='w', ec='b', marker='.', s=mark_size, label='top')
        elif (SP_kind[i] == 3):
            ax1.scatter(SP_x[i,0], SP_x[i,1], c='w', ec='g', marker='.', s=mark_size, label='left')
        elif (SP_kind[i] == 4):
            ax1.scatter(SP_x[i,0], SP_x[i,1], c='w', ec='orange', marker='.', s=mark_size, label='right')
        elif (SP_kind[i] == 5):
            ax1.scatter(SP_x[i,0], SP_x[i,1], c='w', ec='magenta', marker='.', s=mark_size, label='bottom left')
        elif (SP_kind[i] == 6):
            ax1.scatter(SP_x[i,0], SP_x[i,1], c='w', ec='purple', marker='.', s=mark_size, label='bottom right')
        elif (SP_kind[i] == 7):
            ax1.scatter(SP_x[i,0], SP_x[i,1], c='w', ec='cyan', marker='.', s=mark_size, label='top left')
        elif (SP_kind[i] == 8):
            ax1.scatter(SP_x[i,0], SP_x[i,1], c='w', ec='lime', marker='.', s=mark_size, label='top left')

    #=========================#
    #  ax2, ax3, and ax4      #
    #=========================#
    for i in range(N_out):
        if (VM_kind[i] == 1):
            ax2.scatter(VM_x[i,0], VM_x[i,1], c='r', ec='none', marker='.',s=mark_size)
        elif (VM_kind[i] == 2):
            ax2.scatter(VM_x[i,0], VM_x[i,1], c='b', ec='none', marker='.',s=mark_size)
        elif (VM_kind[i] == 3):
            ax3.scatter(VM_x[i,0], VM_x[i,1], c='g', ec='none', marker='.',s=mark_size)
        elif (VM_kind[i] == 4):
            ax3.scatter(VM_x[i,0], VM_x[i,1], c='orange', ec='none', marker='.',s=mark_size)
        elif (VM_kind[i] == 5):
            ax4.scatter(VM_x[i,0], VM_x[i,1], c='magenta', ec='none', marker='.',s=mark_size)
        elif (VM_kind[i] == 6):
            ax4.scatter(VM_x[i,0], VM_x[i,1], c='purple', ec='none', marker='.',s=mark_size)
        elif (VM_kind[i] == 7):
            ax4.scatter(VM_x[i,0], VM_x[i,1], c='cyan', ec='none', marker='.',s=mark_size)
        elif (VM_kind[i] == 8):
            ax4.scatter(VM_x[i,0], VM_x[i,1], c='lime', ec='none', marker='.',s=mark_size)

    ax_list = [ax1, ax2, ax3, ax4]
    for ax in ax_list:
        # lines
        ax.plot([0,0], [0,1], color='k', linewidth=1.2, linestyle='-')
        ax.plot([0,1], [0,0], color='k', linewidth=1.2, linestyle='-')
        ax.plot([0,1], [1,1], color='k', linewidth=1.2, linestyle='-')
        ax.plot([1,1], [0,1], color='k', linewidth=1.2, linestyle='-')

        # axis label
        ax.set_xlabel(r'$x$', fontsize=24, labelpad=10)
        ax.set_ylabel(r'$y$', fontsize=24, labelpad=16)

        # lim
        ax.set_xlim(0-WL_thick-D/2, 1+WL_thick+D/2)
        ax.set_ylim(0-WL_thick-D/2, 1+WL_thick+D/2)
        
        # ticks
        ax.tick_params(axis='both', which='major', direction='out', length=4, width=1, labelsize=14)
        ax.minorticks_off()
    
        # grid
        ax.grid(which='major', color='none', linewidth=0.5)
        ax.set_axisbelow(True) # grid back

    # save
    fig.savefig('{}/initial_setting_1.png'.format(save_path), format='png', dpi=800, transparent=False)
    plt.close()


    #==============================================================================#
    #  figure 2                                                                    #
    #==============================================================================#
    fig = plt.figure(figsize=(12,6), facecolor='white', constrained_layout=True)
    ax1 = fig.add_subplot(121, facecolor='white')
    ax2 = fig.add_subplot(122, facecolor='white')

    #=========================#
    #  ax1                    #
    #=========================#
    def ana_f(x,y):
        return np.sin(np.pi * x) * np.cos(np.pi * y)

    x = np.linspace(0, 1, 200)
    y = np.linspace(0, 1, 200)
    X, Y = np.meshgrid(x, y)  # make mesh grids
    Z = ana_f(X, Y)           # cal f (steady solution)

    cfig = ax1.contourf(X, Y, Z, levels=200, cmap='jet', vmin=-1, vmax=1)
    cbar = plt.colorbar(cfig, aspect=30, shrink=0.8, ax=ax1, orientation='vertical', pad=0.05, location='right', \
                        ticks=[-1.0, -0.5, 0, 0.5, 1.0])

    # color bar
    cbar.ax.tick_params(direction='out', length=4, width=1, labelsize=14)
    cbar.ax.set_ylim(-1, 1)

    # text
    ax1.text(0.08, -0.25, '(a) Analytical steady-state solution', fontsize=16)

    #=========================#
    #  ax2                    #
    #=========================#
    ax2.scatter(SP_x[:N_inn, 0], SP_x[:N_inn, 1], c='silver', ec='k', marker='.', s=mark_size)

    # text
    ax2.text(0.15, -0.25, '(b) Initial particle distribution', fontsize=16)

    ax_list = [ax1, ax2]
    for ax in ax_list:
        # lines
        ax.plot([0,0], [0,1], color='k', linewidth=1.2, linestyle='-')
        ax.plot([0,1], [0,0], color='k', linewidth=1.2, linestyle='-')
        ax.plot([0,1], [1,1], color='k', linewidth=1.2, linestyle='-')
        ax.plot([1,1], [0,1], color='k', linewidth=1.2, linestyle='-')

        # axis label
        ax.set_xlabel(r'$x$', fontsize=24, labelpad=10)
        ax.set_ylabel(r'$y$', fontsize=24, labelpad=16)

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

    # save
    fig.savefig('{}/initial_setting_2.png'.format(save_path), format='png', dpi=800, transparent=False)
    plt.close()
    
# END #