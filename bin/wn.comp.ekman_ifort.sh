module purge
module load intel

cmd="cd ../ftn"
echo $cmd
eval $cmd

cmd="mpiifort -r8 -O2 -assume byterecl -align all -automatic -c nemowwcoupling.f90 " #-i4 -double-size 64
echo $cmd
eval $cmd

cmd="mv nemowwcoupling.o ../obj/ ; mv nemowwcoupling.mod ../mod/"
echo $cmd
eval $cmd
