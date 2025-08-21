#=========================#
#  module                 # 
#=========================#
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import ana_read_fortran_files as read_f90
plt.rcParams['mathtext.fontset'] = 'cm' # mathfont for figure

#=========================#
#  main                   # 
#=========================#
def main(N_sys, N_inn, N_file, last_file, write_step, Delta_t_eff, \
        Pr, gravity, alpha, height, Delta_temp, thermal_dif, parent_path, save_path):
    
    count = 0
    time  = np.zeros((N_file)) # for time record
    V_rms = np.zeros((N_file)) # for V_rms record

    #=========================#
    #  read & cal             # 
    #=========================#
    for i in range(0, last_file+1, write_step):
        time[count] = Delta_t_eff * float(i)
        
        file_name = parent_path + '/data/{}/SP_uvw.dat'.format(i)  # read SP_uvw
        SP_uvw    = read_f90.binary_files(file_name, N_sys, 2)
        if (Pr==0.71):
            SP_uvw   /= math.sqrt(gravity* alpha* height* Delta_temp) # normalize
        elif (Pr>1e10):
            SP_uvw   /= (thermal_dif/ height)                         # normalize
            
        V_rms_tmp = 0.0
        for me in range(0, N_inn):
            V_rms_tmp += SP_uvw[me,0]**2.0 + SP_uvw[me,1]**2.0
        V_rms[count] = math.sqrt(V_rms_tmp / float(N_inn))

        count += 1

    #=========================#
    #  figure                 # 
    #=========================#
    fig = plt.figure(figsize=(7,6), facecolor='white')
    plt.subplots_adjust(left=0.25, right=0.95, bottom=0.18, top=0.9, wspace=0.35, hspace=0.4)
    ax = fig.add_subplot(111, facecolor='white')
    # plot
    ax.plot(time, V_rms, c='k', linewidth=1.5)
    # axis label
    ax.set_xlabel('$t$', fontsize=30, labelpad=10)
    ax.set_ylabel('$V_{\mathrm{rms}}$', fontsize=30, labelpad=16)
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
    # title
    ax.set_title(r'$V_{}={:.5e}$'.format('\mathrm{rms}', V_rms[-1]), y=1.03, fontsize=24)
    # save
    fig.savefig('{}/V_rms.png'.format(save_path), format='png', dpi=300, transparent=False)
    plt.close()

# END #