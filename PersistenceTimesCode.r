## R code for "Biodiversity maintenance may be lower under partial niche differentiation than under neutrality", by Rafael D'Andrea* and Annette Ostling
## *contacting author: rdandrea@umich.edu


library(plyr)

ind=as.numeric(commandArgs(TRUE))[1]

## ---------------- List of scenarios --------------------
## baseline --> pi = const, ri = const, m = 0.01
## A 		--> pi ~ logseries
## B 		--> pi ~ logseries, m = 0.1
## C 		--> ri = xi(1-xi)
## D 		--> ri ~ U(0,1)
## E 		--> pi ~ logseries, w various
## F 		--> pi ~ logseries, w various, circular axis
## G 		--> pi ~ logseries, w various, m = 0.1
## -------------------------------------------------------

d1=data.frame(scen=c('baseline','A','B','C','D'))
d2=data.frame(mode=c('H0','H1'),w=c(Inf,0.063))
d3=data.frame(run=1:10)
scenarios=merge(merge(d1,d2),d3)
d4=data.frame(w=c(.2,.0266,.015,.01,.0075,.005,.0035,.001))
scenarios=rbind(scenarios,merge(merge(data.frame(scen=c('E','F','G'),mode='H1'),d3),d4))
scenarios=scenarios[order(scenarios$scen,scenarios$mode),]
scenarios=rbind(scenarios,merge(data.frame(scen=c('F','G'),mode='H1',w=.063),d3))
mode=scenarios$mode[ind]; scenario=scenarios$scen[ind]; w=scenarios$w[ind]; run=scenarios$run[ind]
axis=ifelse(scenario=='F','circular','linear')
r=ifelse(scenario=='C','parabolic',ifelse(scenario=='D','random','flat'))


## --------------------------------------------------- Metacommunity ----------------------------------------------------------
S=400; mtrait=0:(S-1)/S
p=.99966; k=seq(150e3) #the chosen value of p gives approx. 150,000 individuals in the regional pool
met=if(scenario%in%c('A','B','E','F','G')){ set.seed(run); sample(k,size=S,replace=TRUE,prob=-1/log(1-p)*p^k/k)} else rep(1,S)
##-----------------------------------------------------------------------------------------------------------------------------
	
## --------------------------------------- Parameters ------------------------------------------------
J=21e3
numsteps=100e6
d=as.matrix(dist(mtrait)); if(scenario=='F') d=pmin(d,max(d)-d)
Amet=exp(-d^4/w^4)
r0met=if(scenario=='C') mtrait*(1-mtrait) else if(scenario=='D'){ set.seed(run); runif(S)} else rep(1,S)
m=if(scenario%in%c('B','G')) .1 else .01
## ---------------------------------------------------------------------------------------------------

## -------------- Persistence times list ---------------
outlist=physoutlist=NULL
threshhold=50e6
## -----------------------------------------------------
	
## ----------- Initial Community -------------------
com=sample(mtrait,size=J,prob=met,replace=TRUE)
trait=plyr::count(com)$x
N=plyr::count(com)$freq
A=Amet[mtrait%in%trait,mtrait%in%trait]
r0=r0met[mtrait%in%trait]
B=r0met[mtrait%in%trait]
if(mode=='H1') D=A%*%N else D=sum(N)
##--------------------------------------------------

## ----------------------------- Dynamics --------------------------------------------------------------------------------------------
steps=time=0; set.seed(run)
while(steps<numsteps){
	if(steps==threshhold){ 
		inlist=cbind(steps,which(mtrait%in%trait))
		physinlist=cbind(time,which(mtrait%in%trait))
	}
	if(flux==0 & steps%%1e5==0) plot(trait,N,t='h',main=paste0(mode,' ',scenario,' ',steps/1e4))
	steps=steps+1; time=time+rexp(1,rate=sum(D*N))
	death=sample(seq_along(trait),size=1,prob=D*N)
	N[death]=N[death]-1
	if(mode=='H1') D=D-as.numeric(A[,death])
	if(N[death]==0){
		if(steps>threshhold){ 
			outlist=rbind(outlist,c(steps,which(mtrait==trait[death])))
			physoutlist=rbind(physoutlist,c(time,which(mtrait==trait[death])))
		}
		trait=trait[-death]
		B=B[-death]
		if(mode=='H1'){ 
			D=D[-death]
			A=as.matrix(A[-death,-death])
		}
		N=N[-death]
	}
	if(runif(1)>m){
		birth=sample(seq_along(trait),size=1,prob=B*N)
		N[birth]=N[birth]+1
		if(mode=='H1') D=D+as.numeric(A[,birth])
	}else{
		migrant=sample(mtrait,size=1,prob=met)
		if(migrant%in%trait){ 
			N[trait==migrant]=N[trait==migrant]+1
			if(mode=='H1') D=D+as.numeric(A[,trait==migrant])
		}else{
			if(steps>threshhold) inlist=rbind(inlist,c(steps,which(mtrait==migrant)))
			trait=c(trait,migrant)
			N=c(N,1)
			B=c(B,r0met[mtrait==migrant])
			B=B[order(trait)]
			N=N[order(trait)]
			trait=sort(trait)
			if(mode=='H1'){ 
				A=Amet[mtrait%in%trait,mtrait%in%trait]
				D=A%*%N 
			}
		} 
	}
}
## -----------------------------------------------------------------------------------------------------------------------------------------

## ----------------------------------------- Create Data Frame ---------------------------------------------------------------------
outlist=rbind(outlist,cbind(steps,which(mtrait%in%trait)))
lifetimes=cbind(inlist,outlist)
colnames(lifetimes)=c('in.steps','in.species','out.steps','out.species')
lifetimes=data.frame(lifetimes)
lt=sapply(seq(S),function(s){ i=subset(lifetimes,in.species==s)$in.steps; o=subset(lifetimes,out.species==s)$out.steps; mean(o-i)})

df=ddply(subset(lifetimes,out.steps<1e8),.(out.species),function(v) nrow(v))
ext=rep(0,S); ext[seq(S)%in%df$out.species]=df$V1

abun=rep(0,S); abun[mtrait%in%trait]=N

dat=data.frame(mode=mode,scenario=scenario,w=w,axis=axis,m=m,r=r,run=run,trait=mtrait,N=abun,meta=met,lifetime=lt,extrate=ext)
## ---------------------------------------------------------------------------------------------------------------------------------

