
import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4

file_name = '/Users/mike/Documents/Research/Projects/Iceberg Modeling/Github/iceberg/glacier_visualization/' \
            'src/data/calving_schedules/1992/calving_schedule_180'


schedule = np.fromfile(file_name, '>f8').reshape((4,100000)).T

# sort by the first column
# schedule = schedule[schedule[:, 0].argsort()]
print(schedule)

plt.subplot(3,1,1)
plt.plot(schedule[:,0], schedule[:,1])
plt.title('Width')

plt.subplot(3,1,2)
plt.plot(schedule[:,0], schedule[:,2])
plt.title('Length')

plt.subplot(3,1,3)
plt.plot(schedule[:,0], schedule[:,3])
plt.title('Thickness')

plt.show()





