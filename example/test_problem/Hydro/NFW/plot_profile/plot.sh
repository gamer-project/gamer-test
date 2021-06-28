# $1 : Input    file__Initial   Density
# $2 : Input    file__Final     Density
# $3 : Output   Profile name
cp ../$1 $1
cp ../$2 $2

#export LD_LIBRARY_PATH= YOUR_HDF5_LIBRARY_PATH :$LD_LIBRARY_PATH #export hdf5 library path if it is needed
./GAMER_ExtractProfile  -D -i $1 -S -r 1 -n 128 -x 1.5 -y 1.5 -z 1.5 -p # Where we assume the coordinate of the center of the box is (1.5,1.5,1.5). You may change 1.5 to other number x if the coordinate of the center of the box is (x,x,x).
mv AveParDens AveParDens_init
./GAMER_ExtractProfile  -D -i $2 -S -r 1 -n 128 -x 1.5 -y 1.5 -z 1.5 -p
mv AveParDens AveParDens_final
gnuplot plot_profile.gpt
mv Density_Profile.png $3