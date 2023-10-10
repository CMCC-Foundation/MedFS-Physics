cmd="cd ../ftn"
echo $cmd
eval $cmd

#EXTINC="-I${NETCDF}/include $CFLAGS"
#EXTLIB="-L${NETCDF}/lib -lnetcdff -lnetcdf $LDFLAGS"

EXTINC="$CFLAGS"
EXTLIB="$LDFLAGS"

cmd="mpiifort -g -r8 -O2 -fp-model source -assume byterecl -c outnetcdf.f90 $EXTINC $EXTLIB "
echo $cmd
eval $cmd

cmd="mv outnetcdf.o ../obj/ ; mv outnetcdf.mod ../mod/"
echo $cmd
eval $cmd
