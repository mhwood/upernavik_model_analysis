
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4


total_volume_flux = 8e9
size_categories = 50
beta = -8/5
max_size = 9
min_size = 4

interval_size = (max_size-min_size)/size_categories
v = np.linspace(min_size,max_size,size_categories+1)
v = v[:-1]+interval_size/2
v_edges = np.linspace(min_size,max_size,size_categories+1)

v = 10**v
v_edges = 10**v_edges

# n(v) = c v^-beta
# this is the density of fragments
# integrate over size interval to get total number of bergs
n = (v**beta)

# compute total volume to solve for c
iceberg_volumes_test = np.zeros_like(n)
for vi in range(len(v)):
    iceberg_volumes_test[vi] = n[vi]*v[vi]*(v_edges[vi+1]-v_edges[vi])
c = total_volume_flux/np.sum(iceberg_volumes_test)

number_of_bergs = np.zeros_like(n)
for vi in range(len(v)):
    number_of_bergs[vi] = c*n[vi]*(v_edges[vi+1]-v_edges[vi])
number_of_bergs = np.round(number_of_bergs)

iceberg_volumes = np.zeros_like(n)
for vi in range(len(v)):
    iceberg_volumes[vi] = number_of_bergs[vi]*v[vi]

print('Total number of bergs: '+str(np.sum(number_of_bergs)))
print('Total iceberg volume: '+str(np.sum(iceberg_volumes))+' (check: '+str(total_volume_flux)+')')

fig = plt.figure(figsize=(12,4))

plt.subplot(1,3,1)
plt.loglog(v,c*(v**beta),'b-', label = '$n(v) = $'+str(round(c))+'$\cdot v^{'+'{:.2f}'.format(beta)+'}$')
plt.ylabel('$n(v)$ [1/m$^3$]')
plt.legend()
# plt.xlabel('$v [m$^3$]')

plt.subplot(1,3,2)
plt.loglog(v,number_of_bergs,'b.',label='Total Bergs: '+str(int(np.sum(number_of_bergs))))
plt.ylabel('Number of Icebergs')
plt.xlabel('Volumetric Size Categories [m$^3$]')
plt.legend()
print(number_of_bergs)

plt.subplot(1,3,3)
plt.loglog(v,iceberg_volumes,'b.',label='Total Volume: '+'{:.2e}'.format(np.sum(iceberg_volumes)))
plt.ylabel('Iceberg Volumes [m$^3$]')
plt.legend()

output_file = '/Users/mike/Documents/Research/Projects/Iceberg Modeling/Figures/Glacier/Calving Function Distribution.png'
plt.savefig(output_file)
plt.close(fig)


