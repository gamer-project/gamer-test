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

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

field       = 'density'
#field       = [ 'density', 'Cloud0', 'Cloud1' ]
dpi         = 150


yt.enable_parallelism()

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   pz_dens = yt.ProjectionPlot( ds, 'z', field )
   pz_dens.set_zlim( field, 1.0e-5, 5.0e-2 )
   pz_dens.set_font( {'size':16} )
   pz_dens.set_axes_unit( 'code_length' )
   pz_dens.set_unit( field, 'code_mass/code_length**2' )
   pz_dens.annotate_timestamp( time_unit='code_time', corner='upper_right' )
   pz_dens.save( mpl_kwargs={"dpi":dpi} )

