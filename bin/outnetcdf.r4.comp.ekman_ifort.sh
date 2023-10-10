module purge
module load intel

cmd="cd ../ftn"
echo $cmd
eval $cmd
EXTINC="-I/srv/lib/netcdf-last/include"
EXTLIB="-L/srv/lib/netcdf-last/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lcurl"

cmd="mpiifort -O2 -assume byterecl -align all -automatic -c outnetcdf.f90 $EXTINC $EXTLIB "
echo $cmd
eval $cmd

cmd="mv outnetcdf.o ../obj/ ; mv outnetcdf.mod ../mod/"
echo $cmd
eval $cmd
