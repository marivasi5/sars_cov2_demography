# ABC parameter inference for the north_american population
#input: 
#sim_file: file containing the simulations performed with msABC
#both parameters were drawn from log-uniform priors (using the -R argument)
#obs_file: the vector of the observed summary statistics


#in this case the independent simulations were created in two different files that require merging

#INPUT
#ta dika mou simulations exoun orisei anapota tis parametrous apo tou pavlou
#opote gia na ta enwsw:
pavlos_file= '/home/mariav/synology/pavlos/cov/analysis/ABC/20200402/no_gisaid_cov2020_sequences_hc_20200402_oneline_fullseq_ns_20200402.fasta.MSABC'
my_file= '/home/mariav/synology/pavlos/cov/analysis/ABC/20200402/sims/run_NO_3R_ms.out'

obs_file='/home/mariav/synology/pavlos/cov/analysis/ABC/20200402/no_gisaid_cov2020_sequences_hc_20200402_oneline_fullseq_ns_20200402.fasta.stats'
run='no_FINAL'
tolerance=0.005
path='/home/mariav/synology/pavlos/cov/analysis/ABC/plots'

#===========================================================================================================================================
library(abc)
setwd(path)

#_______________________________________SIMULATED DATA________________________________________________________________________

#merging: ennonw ta dika mou data me tou pavlou:
my_simulations  <- read.table(my_file, header=T)
pavlos_simulations  <- read.table(pavlos_file, header=T)
pavlos_simulations =pavlos_simulations[1:100000,]
pavlos_reordered=pavlos_simulations[,c(2,1,3:ncol(pavlos_simulations))]

simulations=rbind(my_simulations, pavlos_reordered)

params=simulations[,c(1:2)]
stats=simulations[,c(3:ncol(simulations))]
stats=stats[,-c(1,6)] #<-----DIWXNW SEGS KAI FAY---------

#vlepw  sta ss an kapoio column exei mideniko variance      
filtered_exp_stats <- Filter(function(x) var(x) != 0, stats)
ncol(filtered_exp_stats)==ncol(stats)

any(is.na(params))
any(is.na(stats))

#_____________________________________________OBSERVED DATA________________________________________________________________________________________________
obs = read.table(obs_file, h=T)
obs=obs[,-c(1,6)]#DIWXNW SEGS KAI FAY
names(obs)==colnames(stats)    

#for ploting
names(obs)=c('Theta pi','Theta W',"Tajima's D", "ZnS", "DVK","DVH","Thomson estimation", 'Thomson variance')
colnames(stats)=c('Theta pi','Theta W',"Tajima's D", "ZnS", "DVK","DVH","Thomson estimation", 'Thomson variance')

#========================================================================================================================================
#________________________________________________ABC_____________________________________________________________________________________
#========================================================================================================================================

#===============================            LOCLINEAR logit               =====================================================================================================
#logit boundaries
mins <- apply(params, 2, min)
maxs <- apply(params, 2, max)
oria <- cbind(mins, maxs)

abc.res.lin<-abc(target=obs, param=params, sumstat=stats, tol=tolerance, method='loclinear', hcorr=TRUE, transf='logit', logit.bounds=oria)
params_linear_bounds=summary(abc.res.lin)                       
mallon= data.frame(abc.res.lin$ss)  
#write.table(round(params_linear_bounds,2), file='final3/table_NO_lin_logit.txt', quote = F)

#===========  PARAMS CURVE ===============
axonasx=c(expression(log[10](alpha)), expression(log[10](theta)))

onoma_pdf=paste('final3/run_',run,'_params_CURVE_linear_logit_t_',tolerance, '.pdf', sep='')   
onoma_png= gsub('.pdf','.png', onoma_pdf )

pdf(onoma_pdf, height = 4, width = 8)
a<-dev.cur()  #kolpa
png(onoma_png, height = 4, width = 8, units = 'in', res=300)
dev.control("enable")
#plot
par(mfrow=c(1,2))
for(i in 1:ncol(params)){
  par(mar=c(5,4.1,2,2))
  plot(density(log10(params[,i])), 
       xlab=axonasx[i],
       ylab='Density', main='', col='azure2', cex.lab=1.2, cex.axis=0.8, lwd=1.4, frame.plot=F ,
       ylim=c(0,max(density(log10(abc.res.lin$unadj.values[,i]))$y,density(log10(abc.res.lin$adj.values[,i]))$y))
  )  
  polygon(density(log10(params[,i])), col=alpha('azure2', 0.7), border = 'azure2')
  points(density(log10(abc.res.lin$unadj.values[,i])), type='l', col='#4393C3', lwd=1.4)
  polygon(density(log10(abc.res.lin$unadj.values[,i])), col=alpha('#4393C3', 0.5), border = '#4393C3')
  points(density(log10(abc.res.lin$adj.values[,i])), type='l', col='#D6604D', lwd=1.4)
  polygon(density(log10(abc.res.lin$adj.values[,i])), col=alpha('#D6604D', 0.5), border = '#D6604D')
}
#legend("topleft", inset=.02, title="Probability distributions of thleta",
#       c("Prior distribution","Density distribution of the ABC accepted simulations","Posterior distribution after regression adjustment"), fill=c('azure2','#4393C3','#D6604D' ), horiz=F, cex=0.8, border = 'white')

#mtext('Prior and posterior distributions of expansion parameters', side = 3, line = -2, outer = TRUE, cex=1.8)

#kolpa trela kopla
dev.copy(which=a)
dev.off()
dev.off()

#__________________________________________________________________________________________________
#==================== NEW CONTOUR with LOG10(stats)  ==============================================
set.seed(2020)
k=7                                                                 
library(RColorBrewer); my.cols <- rev(brewer.pal(k, "RdYlBu"))        #colors
axonasx=c(expression(log[10](alpha)), expression(log[10](theta)))
#xrisimopoiw ligotera sims gia na fanoun ws koukides sto plot
ns <- tolerance * nrow(params)
#ss <- sample(1:nrow(abc.res.neur$ss),  2000) 
ss <- sample(1:nrow(abc.res.lin$ss),  2500) #gonna use all sims as dots afterall

onoma_pdf=paste('final3/run_',run, '_contour_','linear_logit_t', tolerance , '.pdf', sep='')              
onoma_png= gsub('.pdf','.png', onoma_pdf )

pdf(onoma_pdf, height=12, width=12)
a<-dev.cur()  #kolpa
png(onoma_png, height = 12, width = 12, units = 'in', res=300)
dev.control("enable")
#plot
layout(matrix(1:16, nrow=4, byrow=T))
for(i in 1:2){
  par(mar=c(4,4.1,2,2))
  print(i)
  for(j in 1:8){
    offset <- 3
    yy <- log10(abc.res.lin$ss[ss,j] + offset)
    yall <- log10(abc.res.lin$ss[,j] + offset)
    #line: observed value
    obsyy <- log10(obs[j]+ offset)
    
    plot(log10(abc.res.lin$unadj.values[ss,i]), yy,
         ylim=c(min(yy), max(obsyy, yy)),
         xlab=axonasx[i], ylab=paste('log10(',  colnames(abc.res.lin$ss)[j], ')', sep=''), 
         cex.lab=1.2, cex.axis=0.8, cex=0.5) 
    #contour
    kati= kde2d(log10(abc.res.lin$unadj.values[,i]), yall)
    contour(kati, drawlabels=FALSE, nlevels=k, add=TRUE ,col=my.cols, lwd=2)
    
    abline(h=obsyy, col="red")
  }  
}  
#kolpa trela kopla
dev.copy(which=a)
dev.off()
dev.off()

#===================  simulation stats =======================================
onoma_pdf=paste('final3/run_',run,'_simulationsZOOM_linear_logit.pdf', sep='')
onoma_png= gsub('.pdf','.png', onoma_pdf )

pdf(onoma_pdf, height=6, width=15)
a<-dev.cur()  #kolpa
png(onoma_png, height = 6, width = 15, units = 'in', res=300)
dev.control("enable")
#plot
layout(matrix(1:8, nrow=2, byrow=TRUE))
for (col in 1:ncol(stats)) {
  par(mar = c(5, 4.1, 2, 2))
  
  #xrisimopoiw to 0.90 qualtile 
  if(col==1 | col==2 | col==7 | col==8 ){
    plot(density(stats[,col]), xlab=colnames(stats)[col], ylab="Density", main= ' ',col='antiquewhite4', frame.plot=F,
         cex.lab=1.2, cex.axis=0.8,
         ylim=c(0,max(density(mallon[,col])$y, density(stats[,col])$y)),
         xlim=c(quantile(stats[,col], 0.01), quantile(stats[,col], 0.90))
    )
    points(density(mallon[,col]), type='l', col='deepskyblue4')
    abline(v=obs[col], col="darkorange1")
  }
  # gia to dvh pou einai left skewed tha valw to katw lim sto 0.05
  else if(col==6 ){
    print(colnames(stats)[col]);print(quantile(stats[,col], 0.95))
    plot(density(stats[,col]), xlab=colnames(stats)[col], ylab="Density", main= ' ',col='antiquewhite4', frame.plot=F,
         cex.lab=1.2, cex.axis=0.8,
         ylim=c(0,max(density(mallon[,col])$y, density(stats[,col])$y)),
         xlim=c(quantile(stats[,col], 0.05), quantile(stats[,col], 1))
    )
    points(density(mallon[,col]), type='l', col='deepskyblue4')
    abline(v=obs[col], col="darkorange1")
  }
  
  # ta 0.01 kai 0.99 quantile
  else {
    print(colnames(stats)[col]);print(quantile(stats[,col], 0.99))
    plot(density(stats[,col]), xlab=colnames(stats)[col], ylab="Density", main= ' ',col='antiquewhite4', frame.plot=F,
         cex.lab=1.2, cex.axis=0.8,
         ylim=c(0,max(density(mallon[,col])$y, density(stats[,col])$y)),
         xlim=c(quantile(stats[,col], 0.01), quantile(stats[,col], 0.99))
    )
    points(density(mallon[,col]), type='l', col='deepskyblue4')
    abline(v=obs[col], col="darkorange1")
  }
}
#kolpa trela kopla
dev.copy(which=a)
dev.off()
dev.off()

#parameters correlation ___

pdf('final3/no_correlation.pdf')
plot(x=abc.res.lin$unadj.values[,1], y=abc.res.lin$unadj.values[,2])
dev.off()

cor.test(abc.res.lin$unadj.values[,1], abc.res.lin$unadj.values[,2])

