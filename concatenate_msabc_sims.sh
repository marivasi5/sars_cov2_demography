#bash command to merge all the independent msABC simulations that were performed in separate 'mini' files into one representation
#in this case, we are merging the asian simulations 
#i.e. the msAS_* files

head -1 msAS_1.out > run_AS_3R_ms.out
#prosekse posa arxeia exeis gia to seq
for i in `seq 1 20` ; do grep -vi Tajima msAS_$i.out; done >> run_AS_3R_ms.out

