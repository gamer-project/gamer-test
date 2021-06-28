import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()



idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

colormap    = 'arbre'
dpi         = 150
field       = 'particle_mass'


yt.enable_parallelism()

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
p = []
for idx in range(idx_start,idx_end+1,1):
   ds = yt.load('Data_%06d'%idx)
   p = yt.ParticlePlot(ds, 'particle_position_x', 'particle_position_y','particle_mass', center='c')
   field       = 'particle_mass'         
   p.set_background_color( field )
   p.set_unit('particle_mass', 'Msun')
   p.set_axes_unit( 'code_length' )
   p.set_unit( field, 'code_mass' )
   p.annotate_timestamp( time_unit='code_time', corner='upper_right', text_args={'color':'k'} )
   p.save( mpl_kwargs={"dpi":dpi} )
   p.save()
