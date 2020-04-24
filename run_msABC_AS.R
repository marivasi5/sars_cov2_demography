library(foreach)
library(doMC)

N <- 20
offset <- 20
registerDoMC(N)
setwd('/home/pavlos/synology/pavlos/cov/analysis/ABC')

#400xil sims



rand1 <- sample(1:1000000, N)
rand2 <- sample(1:1000000, N)
rand3 <- sample(1:1000000, N)

#europe_gisaid_cov2020_sequences.20200330.ns.ms.SIM
cmds <- foreach(i=1:N) %dopar% {
    cmd <- paste("msABC 385 15000 -t -R 100 3000 -G -R 1 5000  -seeds ", rand1[i], " ", rand2[i], " ", rand3[i], " > simulations/maya_msAS_",offset+i,".out", sep="")
    system(cmd)
    return(cmd)
}

