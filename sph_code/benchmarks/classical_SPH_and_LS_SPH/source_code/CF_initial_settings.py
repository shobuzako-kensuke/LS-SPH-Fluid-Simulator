#=========================#
#  module                 # 
#=========================#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import ana_read_fortran_files as read_f90
plt.rcParams['mathtext.fontset'] = 'cm' # mathfont for figure

#=========================#
#  draw circle            # 
#=========================#
def draw_circle(ax, radius, center=(0,0), n=300, mycolor='blue'):
    theta = np.linspace(0, 2*np.pi, n)
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)
    ax.plot(x, y, color=mycolor, zorder=5)

#=========================#
#  fig_detail             # 
#=========================#
def fig_detail(ax, lim_switch, boudary_switch):
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
        ax.vlines(x=[0, 1], ymin=0, ymax=1, color='k', linestyle='-')
        ax.hlines(y=[0, 1], xmin=0, xmax=1, color='k', linestyle='-')
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

#=========================#
#  xyz_kind               # 
#=========================#
def xyz_kind(N_sys, N_inn, WL_thick, parent_path, save_path):
    #=========================#
    #  read                   # 
    #=========================#
    file_name = parent_path + '/data/0/SP_xyz.dat'          # SP_xyz
    SP_xyz    = read_f90.binary_files(file_name, N_sys, 2)
    SP_xyz    = SP_xyz[:,:] - WL_thick                      # translation
    file_name = parent_path + '/data/0/SP_kind.dat'         # SP_kind
    SP_kind   = read_f90.binary_files(file_name, N_sys, 1)
    SP_kind   = SP_kind.astype(np.int64)                    # float >> int
    #=========================#
    #  #1. ONLY inners        # 
    #=========================#
    fig = plt.figure(figsize=(6,6), facecolor='white')
    plt.subplots_adjust(left=0.2, right=0.95, bottom=0.18, top=0.95, wspace=0.35, hspace=0.4)
    ax = fig.add_subplot(111, facecolor='white')
    # scatter
    if (N_inn==50**2):
        sca_size = 160
    elif (N_inn==100**2):
        sca_size = 40
    elif (N_inn==200**2):
        sca_size = 15
    ax.scatter(SP_xyz[:N_inn,0], SP_xyz[:N_inn,1], c='dimgray', ec='none', marker='.', s=sca_size)
    # axis label
    ax.set_xlabel('$x$', fontsize=30, labelpad=10)
    ax.set_ylabel('$y$', fontsize=30, labelpad=16)
    # detail
    fig_detail(ax, 1, 0)
    # save
    fig.savefig('{}/xyz_kind.png'.format(save_path), format='png', dpi=300, transparent=False)
    fig.savefig('{}/xyz_kind.pdf'.format(save_path), format='pdf', transparent=True)
    plt.close()
    #=========================#
    #  #2. inners + outers    # 
    #=========================#
    fig = plt.figure(figsize=(6,6), facecolor='white')
    plt.subplots_adjust(left=0.2, right=0.95, bottom=0.18, top=0.95, wspace=0.35, hspace=0.4)
    ax = fig.add_subplot(111, facecolor='white')
    # scatter
    if (N_inn==50**2):
        sca_size = 110
    elif (N_inn==100**2):
        sca_size = 30
    elif (N_inn==200**2):
        sca_size = 10
    ax.scatter(SP_xyz[:N_inn,0], SP_xyz[:N_inn,1], c='dimgray', ec='none', marker='.', s=sca_size, label='inner')
    for i in range(N_inn, N_sys):
        if (SP_kind[i] == 1):
            ax.scatter(SP_xyz[i,0], SP_xyz[i,1], c='r', ec='none', marker='.', s=sca_size, label='bottom')
        elif (SP_kind[i] == 2):
            ax.scatter(SP_xyz[i,0], SP_xyz[i,1], c='b', ec='none', marker='.', s=sca_size, label='top')
        elif (SP_kind[i] == 3):
            ax.scatter(SP_xyz[i,0], SP_xyz[i,1], c='g', ec='none', marker='.', s=sca_size, label='left')
        elif (SP_kind[i] == 4):
            ax.scatter(SP_xyz[i,0], SP_xyz[i,1], c='orange', ec='none', marker='.', s=sca_size, label='right')
        elif (SP_kind[i] == 5):
            ax.scatter(SP_xyz[i,0], SP_xyz[i,1], c='magenta', ec='none', marker='.', s=sca_size, label='bottom left')
        elif (SP_kind[i] == 6):
            ax.scatter(SP_xyz[i,0], SP_xyz[i,1], c='purple', ec='none', marker='.', s=sca_size, label='bottom right')
        elif (SP_kind[i] == 7):
            ax.scatter(SP_xyz[i,0], SP_xyz[i,1], c='cyan', ec='none', marker='.', s=sca_size, label='top left')
        elif (SP_kind[i] == 8):
            ax.scatter(SP_xyz[i,0], SP_xyz[i,1], c='silver', ec='none', marker='.', s=sca_size, label='top left')
    # detail
    fig_detail(ax, 0, 1)
    # save
    fig.savefig('{}/xyz_kind_All.png'.format(save_path), format='png', dpi=300, transparent=False)
    plt.close()
    
#=========================#
#  cell number            # 
#=========================#
def cell(N_sys, N_inn, WL_thick, h_eff, parent_path, save_path):
    #=========================#
    #  read                   # 
    #=========================#
    file_name = parent_path + '/data/0/SP_xyz.dat'          # SP_xyz
    SP_xyz    = read_f90.binary_files(file_name, N_sys, 2)
    SP_xyz    = SP_xyz[:,:] - WL_thick                      # translation
    file_name = parent_path + '/cell/initial_cell.dat'      # initial cell
    SP_cell   = read_f90.binary_files(file_name, N_sys, 1)
    SP_cell   = SP_cell.astype(np.int64)                    # float >> int
    #=========================#
    #  cell                   # 
    #=========================#
    fig = plt.figure(figsize=(7.5,6), facecolor='white')
    plt.subplots_adjust(left=0.18, right=0.95, bottom=0.18, top=0.95, wspace=0.35, hspace=0.4)
    ax = fig.add_subplot(111, facecolor='white')
    # scatter
    if (N_inn==50**2):
        sca_size = 80
    elif (N_inn==100**2):
        sca_size = 10
    elif (N_inn==200**2):
        sca_size = 1
    cfig = ax.scatter(SP_xyz[:,0], SP_xyz[:,1], c=SP_cell, marker='.', cmap='prism', \
                      vmin=SP_cell.min(), vmax=SP_cell.max(), s=sca_size)
    draw_circle(ax, h_eff, center=(SP_xyz[0,0], SP_xyz[0,1])) # first inner particles
    # color bar
    cbar = plt.colorbar(cfig, aspect=30, shrink=0.9, ax=ax, orientation='vertical', pad=0.05, location='right')
    cbar.ax.tick_params(direction='out', length=4, width=1, labelsize=20)
    # detail
    fig_detail(ax, 0, 1)
    # save
    fig.savefig('{}/cell.png'.format(save_path), format='png', dpi=300, transparent=False)
    plt.close()

#=========================#
#  VM_kind                # 
#=========================#
def xyz_kind_VM(N_VM, N_inn, WL_thick, parent_path, save_path):
    #=========================#
    #  read                   # 
    #=========================#
    file_name = parent_path + '/data/0/VM_xyz.dat'          # VM_xyz
    VM_xyz    = read_f90.binary_files(file_name, N_VM, 2)
    VM_xyz    = VM_xyz[:,:] - WL_thick                      # translation
    file_name = parent_path + '/data/0/VM_kind.dat'         # SP_kind
    VM_kind   = read_f90.binary_files(file_name, N_VM, 1)
    VM_kind   = VM_kind.astype(np.int64)                    # float >> int
    #=========================#
    #  VM kind                # 
    #=========================#
    fig = plt.figure(figsize=(12,12), facecolor='white')
    plt.subplots_adjust(left=0.12, right=0.95, bottom=0.1, top=0.93, wspace=0.35, hspace=0.4)
    ax = fig.add_subplot(221, facecolor='white')
    # scatter
    if (N_inn==50**2):
        sca_size = 100
    elif (N_inn==100**2):
        sca_size = 30
    elif (N_inn==200**2):
        sca_size = 10

    for i in range(N_VM):
        if (VM_kind[i] == 1):
            ax.scatter(VM_xyz[i,0], VM_xyz[i,1], c='r', ec='none', marker='.',s=sca_size)
        elif (VM_kind[i] == 2):
            ax.scatter(VM_xyz[i,0], VM_xyz[i,1], c='b', ec='none', marker='.',s=sca_size)
    # detail
    fig_detail(ax, 1, 0)

    ax = fig.add_subplot(222, facecolor='white')
    for i in range(N_VM):
        if (VM_kind[i] == 3):
            ax.scatter(VM_xyz[i,0], VM_xyz[i,1], c='g', ec='none', marker='.',s=sca_size)
        elif (VM_kind[i] == 4):
            ax.scatter(VM_xyz[i,0], VM_xyz[i,1], c='orange', ec='none', marker='.',s=sca_size)
    # detail
    fig_detail(ax, 1, 0)

    ax = fig.add_subplot(223, facecolor='white')
    for i in range(N_VM):
        if (VM_kind[i] == 5):
            ax.scatter(VM_xyz[i,0], VM_xyz[i,1], c='magenta', ec='none', marker='.',s=sca_size)
        elif (VM_kind[i] == 6):
            ax.scatter(VM_xyz[i,0], VM_xyz[i,1], c='purple', ec='none', marker='.',s=sca_size)
    # detail
    fig_detail(ax, 1, 0)

    ax = fig.add_subplot(224, facecolor='white')
    for i in range(N_VM):
        if (VM_kind[i] == 7):
            ax.scatter(VM_xyz[i,0], VM_xyz[i,1], c='cyan', ec='none', marker='.',s=sca_size)
        elif (VM_kind[i] == 8):
            ax.scatter(VM_xyz[i,0], VM_xyz[i,1], c='silver', ec='none', marker='.',s=sca_size)
    # detail
    fig_detail(ax, 1, 0)
    # save
    fig.savefig('{}/xy_kind_VM.png'.format(save_path), format='png', dpi=300, transparent=False)
    plt.close()

#=========================#
#  cell number for VM     # 
#=========================#
def cell_VM(N_VM, N_inn, WL_thick, h_eff, parent_path, save_path):
    #=========================#
    #  read                   # 
    #=========================#
    file_name = parent_path + '/data/0/VM_xyz.dat'          # VM_xyz
    VM_xyz    = read_f90.binary_files(file_name, N_VM, 2)
    VM_xyz    = VM_xyz[:,:] - WL_thick                      # translation
    file_name = parent_path + '/cell/initial_VM_cell.dat'   # initial cell for VM
    VM_cell   = read_f90.binary_files(file_name, N_VM, 1)
    VM_cell   = VM_cell.astype(np.int64)                    # float >> int
    #=========================#
    #  VM kind                # 
    #=========================#
    fig = plt.figure(figsize=(7.5,6), facecolor='white')
    plt.subplots_adjust(left=0.18, right=0.95, bottom=0.18, top=0.95, wspace=0.35, hspace=0.4)
    ax = fig.add_subplot(111, facecolor='white')
    # scatter
    if (N_inn==50**2):
        sca_size = 100
    elif (N_inn==100**2):
        sca_size = 20
    elif (N_inn==200**2):
        sca_size = 0.5
    cfig = ax.scatter(VM_xyz[:,0], VM_xyz[:,1], c=VM_cell, marker='.', cmap='prism', vmin=VM_cell.min(), vmax=VM_cell.max(), s=sca_size)
    draw_circle(ax, h_eff, center=(VM_xyz[0,0], VM_xyz[0,1])) # first inner particles
    # color bar
    cbar = plt.colorbar(cfig, aspect=30, shrink=0.9, ax=ax, orientation='vertical', pad=0.05, location='right')
    cbar.ax.tick_params(direction='out', length=4, width=1, labelsize=20)
    # detail
    fig_detail(ax, 0, 1)
    # save
    fig.savefig('{}/cell_VM.png'.format(save_path), format='png', dpi=300, transparent=False)
    plt.close()

# END #