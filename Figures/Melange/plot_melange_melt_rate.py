


import os
import numpy as np
import matplotlib.pyplot as plt


def read_melange_melt_rate(config_dir):
    file_names = []
    for file_name in os.listdir(config_dir + '/run_melange/diags/BRGflx'):
        if file_name[0]!='.' and file_name.endswith('.data'):
            file_names.append(config_dir + '/run_melange/diags/BRGflx/'+file_name)
    file_names = sorted(file_names)

    for file_name in file_names:
        grid = np.fromfile(file_name, '>f4').reshape((2*61, 375, 450))

        print(np.min(grid), np.max(grid))

        # plt.pcolormesh(np.sum(grid,axis=0))
        # plt.title(file_name)
        # plt.colorbar()
        # plt.show()

    return grid





config_dir = '/Users/mike/Documents/Research/Projects/Greenland Model Analysis/Fjord/Upernavik'

read_melange_melt_rate(config_dir)


