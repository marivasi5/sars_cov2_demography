library(foreach)
library(doMC)
registerDoMC(20)
setwd('/home/mariav/synology/pavlos/cov/analysis/ABC')

#400xil sims
N <- 16


rand1 <- sample(1:1000000, N)
rand2 <- sample(1:1000000, N)
rand3 <- sample(1:1000000, N)

#europe_gisaid_cov2020_sequences.20200330.ns.ms.SIM
cmds <- foreach(i=1:N) %dopar% {
    cmd <- paste("msABC 811 25000 -G -R 1 2000 -t -R 100 2000  -seeds ", rand1[i], " ", rand2[i], " ", rand3[i], " > simulations/msEU_",i,".out", sep="")
    system(cmd)
    return(cmd)
}

