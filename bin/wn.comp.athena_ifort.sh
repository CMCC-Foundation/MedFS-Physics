cmd="cd ../ftn"
echo $cmd
eval $cmd

cmd="mpiifort -g -r8 -O2 -fp-model source -assume byterecl -c nemowwcoupling.f90 "
echo $cmd
eval $cmd

cmd="mv nemowwcoupling.o ../obj/ ; mv nemowwcoupling.mod ../mod/"
echo $cmd
eval $cmd
