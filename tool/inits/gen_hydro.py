from pyFC import LogNormalFractalCube, plot_field_stats, write_cube
#import pyFC
import matplotlib.pyplot as pl
import matplotlib.cm as cm

N = 3

fc = LogNormalFractalCube(ni=N, nj=N, nk=N, kmin=1, kmax=None, mean=1, sigma=2.23606797749979, beta=-1.66666666666)

fc.gen_cube()

plot_field_stats(fc, scaling='log', vmin=-2.1, vmax=2.1, cmap=cm.jet)
write_cube(fc=fc, fname='data.dbl', app=True, prec='double')
