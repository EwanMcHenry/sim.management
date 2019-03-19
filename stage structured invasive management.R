# script file name "stage structured invasive management.R"
# Run this code

set.seed(121212121)
##########################################################################################################
##                                  load libraries ----
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#

library(partitions)
library(gtools)
library(compiler)
#library(scatterplot3d)
#library(plot3D)
#library(rgl)
library(abind)
library(fitdistrplus)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(viridis)
library(partykit) # for ctree

# funciton to find midpoint of range
midpoint = function(x){min(x) + (range(x)[2]- range(x)[1])/2}

#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
##                                  define working directories  -----
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
this.wd = "..."
perfect.reaction.simulations.directory = paste(this.wd, "\\reactive and constant simulations.R", sep = "")
imperfect.monitoring.source.code = paste(this.wd, "\\example run of imperfect monitoring v3.0.R", sep = "")
monitoirng.source.code = paste(this.wd, "\\Imperfect monitoring simulation.R", sep = "")
# ----
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#




#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
##                                      SIMULAITON SPECIFICATION     ------
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
            #####  parameters relating to time ----
            burn.in.years = 1000 ########################### 19.10.18
            sim.years = 1000
            years = burn.in.years + sim.years            # number of years to simulate over
            year.length = 4         # number of timesteps in each year
            max.time = years * year.length +1      # maximum timestep
            step.length = 1
        
            season.prop = 4                     #number of seasons per a year
            season.length = (year.length / season.prop) # number of timesteps in
            
            start.s1 = seq(from = step.length, to = max.time , by = year.length)
            start.s2 = seq(
                from = round(year.length / season.prop) + step.length,
                to = max.time ,
                by = year.length
            )
            start.s3 = seq(from = round(2 * (year.length / season.prop)) + step.length,
                           to = max.time ,
                           by = year.length)
            start.s4 = seq(from = round(3 * (year.length / season.prop)) + step.length,
                           to = max.time ,
                           by = year.length)
            
            is.it.winter = rep(0, times = max.time)
            is.it.winter[start.s4] = 1   # 1 when winter to kill(1- is it winter is multiplied by juvinile population)
            is.it.birth = rep(0, times = max.time)
            is.it.birth[start.s3] = 1 # 1 when it is season 3
            
            which.year = c(0, rep(rep(1:(
                years + burn.in.years
            ) , each = year.length), length = max.time - 1)) # note of year numbers
            which.season.is.it = c(4, rep(rep(1:4, each = year.length / 4), length = max.time -
                                              1)) # note of season numbers (within years)
            which.step = 1:max.time                                                             # note of timestep number
            analysable.years = ((burn.in.years + 1):(years))
            analysable.steps = which.step[which.year %in% analysable.years]
            analysible.s2 = which.step[which.season.is.it == 2 &
                                           which.step %in% analysable.steps] #
            
            # life history timing ----
            single.pulse.maturation.times = start.s3
            first.maturation.start = single.pulse.maturation.times[1]
            maturation.length = 2
            first.maturation.end = first.maturation.start + (maturation.length-1)
            first.maturation.season =maturation.seasons = first.maturation.start:first.maturation.end
            while (max.time > max(maturation.seasons)+ year.length){
                maturation.seasons = c(first.maturation.season,maturation.seasons+year.length)}
            
            maturation.times = maturation.seasons
            maturing.time = rep(0, times = max.time)
            maturing.time[maturation.seasons] = 1
            
            is.it.winter = rep(0,times = max.time) 
            is.it.winter[start.s4] =1   # 1 when winter to kill(1- is it winter is multiplied by juvinile population)
            
            birth.times = start.s3
            is.it.birth = rep(0, times = max.time)
            is.it.birth[birth.times] = 1
        
            # monitoring times
            season.of.monitoring = matrix(NA, nrow = max.time, ncol = year.length)
            for (y in 1:dim(season.of.monitoring)[2]) {
                for (x in 1:dim(season.of.monitoring)[1]) {
                    season.of.monitoring[x, y] = max(which.step[c(2:4, 1)[which.season.is.it] == y &
                                                                    which.step <= x])
                }
            }
            season.of.monitoring[season.of.monitoring == "-Inf"] = NA
            
            ##### parameters relating to the invasive populaiton and its trapping -----
            # BIRTH RATES - DEPENDING ON LEVEL OF ENVIRONEMENTAL STOCHASTICITY ----
            mean.b = 2 # mean rate of dispersing juvenile production by the adult populaiton
            
            # stochastic birth rates----
            dummy.b.params = seq(0, 10, by = 0.1) # range of potential beta vals
            dummy.less.than.1 = pbeta(0.25, dummy.b.params, dummy.b.params)# find proportion of the aproximatley normal beata distribution parameterised using b(dummy.b.params,dummy.b.params) that is less than 0.25 (1/2 of the mean)
            beta.vals = c(dummy.b.params[which.min(abs(dummy.less.than.1 - c(0.05)))] ,
                          dummy.b.params[which.min(abs(dummy.less.than.1 - c(0.1)))] ,
                          dummy.b.params[which.min(abs(dummy.less.than.1 - c(0.2)))]) # find which is the closest of the beta parameters to produce a distribution where 5%, 10% and 20% of the distribution is less than 0.25 (1/2 of the mean)
            # draw a random vector of birth rates from the distributions definied by the beta parameters selected above
            stochb.05less1 = 4 * rbeta(max.time, beta.vals[1], beta.vals[1])
            stochb.10less1 = 4 * rbeta(max.time, beta.vals[2], beta.vals[2])
            stochb.20less1 = 4 * rbeta(max.time, beta.vals[3], beta.vals[3])
            
            # NATURAL MORTALITY ----
            yearly.z.A =  0.4      #natural adult predator death rate
        #    yearly.z =  0.2      #natural predator death rate
            death.chances.for.juv = 5   # number of mortality events experienced by juveniles - to simulate higher mortality for dispersing individuals
            mean.z.A = 1 - (1 - yearly.z.A) ^ (1 / year.length)
            z.A = rep (mean.z.A, times = max.time)
            # z.A = z.A.stoch ("constant")
            z.J = 1 - (1 - z.A) ^ death.chances.for.juv  # here juvinile mortality is modeled as multiple adult mortality events
            yearly.z.J =  0.2     #natural predator death rate
            yearly.survival = 1-yearly.z.A
            death.chances.for.juv = 5
        
            # SETTLEMENT ----
            settleable =  0.75 # proportion of floting juveniles able to find territory
            pk = pk1= 20 # total amount of territory availiblility (carrying capacity)
            m=  0.9       #predator maturation rate: assume that all non-maturing juviniles die
            yearly.fill.prop = rep(0,max.time)
            yearly.fill.prop[maturation.seasons] = fill.prop = 0.1
            proportion.to.fill = yearly.fill.prop
            ratio.for.saturation = 1
            fill.prop = 0.7
            proportion.to.fill = rep(0,max.time)
            proportion.to.fill[maturation.seasons] = fill.prop 
            
            #  TRAPPING ----
            prop.caught.w.1.trap = 0.1 #  proportion of the adult population removed by a trapping effort of 1
            juv.mods = c(0, 1, 2, 3, 4) # variations of jvenile trappability -- can be thoght of as the number of trapping events encountered by juveniles for every one encountered by the adult population
            yearly.trapping.effort = 1 # mean trappign effort to be spent each year
            
            # for constnat prpororiton #####
            h = 0.2
            #for DE effort+ efficiency ########
            trap.b0 = 0.2 # must be >= 0
            trap.b1 = 1
            # trap.effort = constant.effort = rep(1, times = max.time)
            trap.b2 = 0
            trap.effic = rep(1, times = max.time)
            effort.multiplier = 1
            
            # have a look at the trapping function
            #for effort
            dummy.effort = seq(from = 0, to = 10, length = 100)
            dummy.v = (trap.b0 * (dummy.effort^trap.b1))/(1+(trap.b0 * (dummy.effort^trap.b1)))
            #plot(dummy.v~dummy.effort, xlab = "Effort(N trap nights)", ylab = "propotion of populaiton harvested"  , ylim = c(0,1))
            #for "raft" catching function ######
            largest.effort = 20
            # effort budget ----
            horizon.effort.budget = yearly.trapping.effort*length(analysable.years)
        #-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
#  INITIAL SYSTEM STATE ############
    prey.init = 350
    adult.init = 20
    juv.init = 2
    float.init = 0
    A.harvest.init = 0
    J.harvest.init = 0
    prey.eaten.init = 0
    proporiton.harvest = 0
    n.settle = 0
    n.born = 0
    winter.dieoff = 0
    trapping.effort = 0
    init.state = c(
        prey.init ,
        adult.init ,
        juv.init  ,
        float.init,
        A.harvest.init,
        J.harvest.init,
        prey.eaten.init ,
        proporiton.harvest ,
        n.settle,
        n.born ,
        winter.dieoff ,
        trapping.effort
    )
    state = matrix(c(init.state , rep(rep(
        NA, times = length(init.state)
    ), times = max.time - 1)),
    ncol = length(init.state),
    byrow = T)
    state.col.names = c(
        "prey",
        "A",
        "J",
        "F",
        "A.h",
        "J.h",
        "eat",
        "h.prop",
        "set",
        "born",
        "winter",
        "trapping.effort"
    )
    real.state.col.names = c(
        "Prey",
        "Adult",
        "Juv",
        "Floaters",
        "Adult harv",
        "Juv harv",
        "Eaten",
        "Harvest proportion",
        "Settling",
        "Born",
        "Winter",
        "Trapping.effort"
    )
    colnames(state) = state.col.names
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
##                                      SIMULATION FUNCTIONS                                           ## 
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
# function to add juveniles to the population, produced by the adult population
J.pred.di.birth= function (pops= state[t+1,]){
    #  Juv.Pred density independent birth
    # b = number of females born per adult female inthe population
    #is.it.birth = binary vector, 1 = birth timestep, 0 = not
    birth = pops[2]*b[t]*is.it.birth[t] # number born into populaiton. = 0 is not birth season
    z1 <- pops[3] + birth
    return (c(pops[1], pops[2] , z1  , pops[4]  , pops[5] , pops[6] , pops[7] , pops[8] , pops[9] , birth , pops[11] , pops[12] ))  
}
comp.J.pred.di.birth = cmpfun (J.pred.di.birth)

# function to simulate natural adult mortality
A.pred.di.mort = function (pops= state[t+1,]){
    # computes pops after density independent adult predator mortality
    # z.A - adult mortality rate per timestep
    y1 <- pops[2] - pops[2]*z.A[t] # can easily alter here to add age specific death
    return (c(pops[1], y1 , pops[3] , pops[4] , pops[5] , pops[6] , pops[7] , pops[8] , pops[9] , pops[10] , pops[11] , pops[12] ))  
}
comp.A.pred.di.mort = cmpfun (A.pred.di.mort)

J.pred.di.mort = function (pops= state[t+1,]){
    # computes pops after density independent Juvinile predator mortality
    #parameters
    # z.J - juvinile mortality rate per time step
    
    z1 <- pops[3] - pops[3]*z.J[t]
    return (c(pops[1], pops[2] , z1  , pops[4] , pops[5] , pops[6] , pops[7] , pops[8] , pops[9] , pops[10] , pops[11] , pops[12] ))  
}

comp.J.pred.di.mort = cmpfun (J.pred.di.mort)


# function to move a proportion of the juvenile population to the adult population -- simulating settlement into limited adult territory 
simplest.dd.settlement = function(pops= state[t+1,], able.to.settle = settleable){
    empty.n = max(c(pk-pops[2],0)) # computes amount of space in adult territory. 
    #pk is the total availible adult territory
    # if adult popualiton > territory avialible empty space = 0
    settling = min(c(pops[3]*able.to.settle * empty.n/pk, empty.n))    
    # only a proportion (able.to.settle) of the juvenile populaiton is able to settle, simulating imperfect territory searching
    # of that abel to settle a proporiton equal to the proportion of empty habiat will settle, up to a maximum of hte amount of empty territory
    y1 = pops[2] + settling
    z1 = pops[3] - settling
    return (c(pops[1], y1 , z1 , pops[4]  , pops[5] , pops[6] , pops[7] , pops[8] , settling , pops[10] , pops[11] , pops[12] ))
}
comp.simplest.dd.settlement = cmpfun (simplest.dd.settlement)

# funciton to remove all juveniles from the population at the end of the winter season
kill.all.J= function (pops= state[t+1,]){
    # kills all remaning juviniles if IN WINTER   
    # is.it.winter = vecotr of winter if its winter ==1
    
    z1 = pops[3] * (1-is.it.winter[t])
    return (c(pops[1], pops[2] , z1 , pops[4]  , pops[5] , pops[6] , pops[7] , pops[8] , pops[9] , pops[10] , pops[3] *(is.it.winter[t]) , pops[12] )) 
}

comp.kill.all.J = cmpfun (kill.all.J)


# funciton to make all populations = 0
        reset.trap.eat = function (pops = state[t+1,]){
            return (c( pops[1:4] , rep(0, times = 8)) )
        }
        comp.reset.trap.eat = cmpfun(reset.trap.eat) # compile function
        
# function to remove adults and juveniles by trapping
        better.trapping= function(pops= state[t+1,], juv.trap.chances = juvinile.traping.modifier, current.effort = trap.effort[t], prob.catch1 = 0.1){
            adult.trap.prop = v1 =  1-((1-prob.catch1)^current.effort)
            juv.trap.prop = 1- ((1- adult.trap.prop )^ juv.trap.chances )
            a.h = pops[2] * adult.trap.prop
            j.h = pops[3] * juv.trap.prop
            y1 = pops[2] - a.h
            z1 = pops[3] - j.h
            return (c( pops[1], y1 , z1 , pops[4] , a.h , j.h , pops[7] , v1 , pops[9], pops[10] ,pops[11] , current.effort )) 
            #juvinile.traping.modifier is the trappability of juveniles relative to adults
            # prob.catch1 is the proportion of hte adult population removed with a trapping effort of 1
        }        
        comp.better.trapping = cmpfun(better.trapping)
        
        
# function to simulate monitoing of the adult populaiton from a number of detection events and return the proportion of detection events that return positive detections 
        raft.monitor= function(pops= state[t+1,], n.rafts = 1, prob.1adult = 0.1){
            detectable.pop = sum(pops[2] )
            detect.prob = v1 =  1-((1-prob.1adult)^detectable.pop) # probability of a positive detection increases with populaiton density
            n.detections = rbinom(1, n.rafts, detect.prob)
            prop.positive = n.detections/n.rafts
            return (prop.positive) 
        }
        comp.raft.monitor = cmpfun(raft.monitor)
        


#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
##              STRATEDGIES -- PARAMETARISING SEASONAL SPLIT OF TRAPPING EFFORT ----
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
# split effort by season in 1/4s, w' all possible combinations of each proportion
eff.partition.a<- parts(4)/4 # first enumerate partitions of your effort
eff.partition.b<- perms(4) # get permutations of effort 

eff.partition.ab<- lapply(data.frame(as.matrix(eff.partition.b)), function(x)eff.partition.a[x,]) # get all permutations of each partition
factorial.season.effort.partition1<- unique(do.call(rbind, lapply(eff.partition.ab, t))) # put all together and remove duplicates

factorial.season.effort.partition.reorder = factorial.season.effort.partition1[c(1 ,12 , 16, 19, 2, 7, 9, 13, 24, 26, 17, 28, 30, 20, 29, 32, 3, 8, 10, 25, 27, 31, 4, 6, 11, 14, 15, 33, 18, 22, 34, 21, 35, 23, 5),]
factorial.season.effort.partition =factorial.season.effort.partition.reorder  
rownames(factorial.season.effort.partition) = 1:35
colnames(factorial.season.effort.partition) = c("S1", "S2", "S3", "S4")

# repeat each seasonal effort the number of steps in that season
seasonal.step.effort.mod = factorial.season.effort.partition[,(rep(1:4, rep(season.length,4)))]/season.length
n.seasonal.eff.combos = dim(seasonal.step.effort.mod)[1]

#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
##              STRATEDGIES -- REACTIVE TRAPPING WITH CONCENTRATED EFFORT ----
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#

eff.prop.if.over = seq(from = 2, to = 10, by = 1)
eff.prop.if.over = c(2,3,4,5,6)

eff.prop.if.less = rep(0, length(eff.prop.if.over))

#monitoring.season = 1:year.length
monitoring.season = 1

#juv.misID.rate = c(0,0.25,0.5,1)
juv.misID.rate = 1
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
# explore stoat stochasticity usign paramteres form King et al ----

stoat.brates = c(rnorm(10000, 4.4, 1.35), rnorm(10000,1.64 , 0.62) , rnorm(10000, 0.75 , 0.58 ), rnorm(10000, 3.95 , 1.28 ))
hist(stoat.brates)
minkier = stoat.brates/ (mean(stoat.brates)/2)
hist(minkier) # can be really variable

single.phase.stoat.brates = rnorm(10000, 4.4, 1.35)/(4.4/mean.b)
par(mfrow = c(1,1)) ;  hist(rnorm(10000, 4.4, 1.35)/(4.4/mean.b))
pnorm((4.4/2), 4.4, 1.35)
pnorm((1.64/2),1.64 , 0.62)
pnorm((0.75/2), 0.75 , 0.58)
pnorm((3.95/2), 3.95 , 1.28)
# ----
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
# specification of example traces ----
which.trigs = c(1,2,4)
examp.which.ss = c(35, 2,3) # seasonal stratedgies for monitoring sim
the.ylim = c(0,10)
which.ss = 35 # effort eveningly spread across years
examp.juv.mod = 1 # for equally trappable adults and juveniles 
which.c = which(    juv.mods %in% examp.juv.mod)
trace.length =50        # length of example traces (years)
trace.start = 0
show.time = trace.start:(trace.start + trace.length) # years included in trace of fig 2

which.interestin.eff.overs = 4
interesting.nrafts = 10 # number of rafts for monitoring trace
# ----
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#






#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
##                          SELECT LEVEL OF ENVIRONEMENTAL STOCHASTICITY AND RUN ---- 

#constant
b = rep (mean.b , time = max.time)
source (perfect.reaction.simulations.directory)
determ.sim = sims.for.analysis
save(determ.sim, file = paste(this.wd, "\\determ.sim.RData", sep = ""))
determ.analysible.s2.q90 = apply(sims.for.analysis[analysible.s2,2,,,,,], c(2:4), FUN = quantile, probs = 0.9)
save(determ.analysible.s2.q90, file = paste(this.wd, "\\determ.analysible.s2.q90.RData", sep = ""))
determ.trigs = sim5.triggered.aprox.triggers
save(determ.trigs, file = paste(this.wd, "\\determ.trigs.RData", sep = ""))

# environementally stochastic b
# low stochasticity
b = stochb.05less1 
source (perfect.reaction.simulations.directory)
s05.sim = sims.for.analysis
save(s05.sim, file = paste(this.wd, "\\s05.sim", sep = ""))
s05.analysible.s2.q90 = apply(sims.for.analysis[analysible.s2,2,,,,,], c(2:4), FUN = quantile, probs = 0.9)
save(s05.analysible.s2.q90, file = paste(this.wd, "\\s05.analysible.s2.q90.RData", sep = ""))
s05.trigs = sim5.triggered.aprox.triggers
save(s05.trigs, file = paste(this.wd, "\\s05.trigs", sep = ""))

# high stochasticity
b = stochb.20less1 
source (perfect.reaction.simulations.directory)
s20.sim = sims.for.analysis
save(s20.sim, file = paste(this.wd, "\\s20.sim", sep = ""))
s20.analysible.s2.q90 = apply(sims.for.analysis[analysible.s2,2,,,,,], c(2:4), FUN = quantile, probs = 0.9)
save(s20.analysible.s2.q90, file = paste(this.wd, "\\s20.analysible.s2.q90.RData", sep = ""))
s20.trigs = sim5.triggered.aprox.triggers
save(s20.trigs, file = paste(this.wd, "\\s20.trigs", sep = ""))
##########################################################################################################
##                          EXAMPLE IMPERFECT MONITORING ----
# specifying example  traces
# 3 lines, seasonality of trapping: even, pre-juv and post juv recruitment
# deterministic & medium stochasticity KD of 3 as example, add line for threshold

b = rep (mean.b , time = max.time)
examp.trigs = determ.trigs
source(imperfect.monitoring.source.code) # makes e.g.trace.save
determ.monitor.trace = e.g.trace.save[, which.ss ,  , which.interestin.eff.overs , ]
save(determ.monitor.trace, file = paste(this.wd, "\\e.g.determ.monitor.trace", sep = ""))

b = stochb.05less1 
examp.trigs = s05.trigs
source(imperfect.monitoring.source.code) # makes e.g.trace.save
stoch05.monitor.trace = e.g.trace.save[, which.ss ,  , which.interestin.eff.overs , ]
save(stoch05.monitor.trace, file = paste(this.wd, "\\e.g.s05.monitor.trace", sep = ""))

b = stochb.20less1 
examp.trigs = s20.trigs
source(imperfect.monitoring.source.code) # makes e.g.trace.save
stoch20.monitor.trace = e.g.trace.save[, which.ss ,  , which.interestin.eff.overs , ]
save(stoch20.monitor.trace, file = paste(this.wd, "\\e.g.s20.monitor.trace", sep = ""))
##########################################################################################################
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#






#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
##                      LOAD SIM OBJECTS ------
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
load(paste(this.wd, "\\determ.sim.RData", sep = ""))
load(paste(this.wd, "\\determ.analysible.s2.q90.RData", sep = ""))
load(paste(this.wd, "\\determ.trigs.RData", sep = ""))
load(paste(this.wd, "\\e.g.determ.monitor.trace", sep = ""))

load(paste(this.wd, "\\s05.sim", sep = ""))
load(paste(this.wd, "\\s05.analysible.s2.q90.RData", sep = ""))
load(paste(this.wd, "\\s05.trigs", sep = ""))
load(paste(this.wd, "\\e.g.s05.monitor.trace", sep = ""))

load(paste(this.wd, "\\s20.sim", sep = ""))
load(paste(this.wd, "\\s20.analysible.s2.q90.RData", sep = ""))
load(paste(this.wd, "\\s20.trigs", sep = ""))
load(paste(this.wd, "\\e.g.s20.monitor.trace", sep = ""))
##  ------

#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#


##########################################################################################################
##                          MAKE FIG 2 - example traces                              ####

the.trig = which.interestin.eff.overs

determ.trace.const = determ.sim[analysible.s2,2,,which.ss, which.c ,1,]
determ.trace.perf = determ.sim[analysible.s2,2,,which.ss, which.c ,the.trig,]
determ.trace.imperf = determ.monitor.trace[analysible.s2]
determ.trace.all= matrix(NA,nrow = length(determ.trace.const)*3,ncol = 3)
determ.trace.all = data.frame(s2.dens = c(determ.trace.perf,  determ.trace.imperf , determ.trace.const) , 
                              kd.strat = c(rep ("perf", times = length (determ.trace.perf)), rep ("imperf", times = length (determ.trace.imperf)), rep("const", times = length (determ.trace.const))   ),
                              year = (analysible.s2-3)/4 - burn.in.years + 1)
determ.eg.traces = determ.trace.all[determ.trace.all$year %in% show.time,]

s05.trace.const = s05.sim[analysible.s2,2,,which.ss, which.c ,1,]
s05.trace.perf = s05.sim[analysible.s2,2,,which.ss, which.c ,the.trig,]
s05.trace.imperf = stoch05.monitor.trace[analysible.s2]
s05.trace.all= matrix(NA,nrow = length(s05.trace.const)*3,ncol = 3)
s05.trace.all = data.frame(s2.dens = c( s05.trace.perf,  s05.trace.imperf , s05.trace.const) , 
                           kd.strat = c( rep ("perf", times = length (s05.trace.perf)), rep ("imperf", times = length (s05.trace.imperf)) , rep("const", times = length (s05.trace.const))  ),
                           year = (analysible.s2-3)/4 - burn.in.years + 1 )
s05.eg.traces = s05.trace.all[s05.trace.all$year %in% show.time,]

s20.trace.const = s20.sim[analysible.s2,2,,which.ss, which.c ,1,]
s20.trace.perf = s20.sim[analysible.s2,2,,which.ss, which.c ,the.trig,]
s20.trace.imperf = stoch20.monitor.trace[analysible.s2]
s20.trace.all= matrix(NA,nrow = length(s20.trace.const)*3,ncol = 3)
s20.trace.all = data.frame(s2.dens = c(s20.trace.perf,  s20.trace.imperf , s20.trace.const) , 
                           kd.strat = c( rep ("perf", times = length (s20.trace.perf)), rep ("imperf", times = length (s20.trace.imperf)) , rep("const", times = length (s20.trace.const))  ),
                           year = (analysible.s2-3)/4 - burn.in.years + 1 )
s20.eg.traces = s20.trace.all[s20.trace.all$year %in% show.time,]


trace.range = c(0,10)

determ.marg.const = ggdensity(determ.trace.all[determ.trace.all$kd.strat == "imperf" | determ.trace.all$ kd.strat == "perf",], "s2.dens" , fill = "kd.strat" , palette = "jco" )+ 
    clean_theme() +
    xlim (trace.range)+ 
    coord_flip()#+
#    theme(plot.margin = unit(c(1,0,1,0), "lines"))


s05.marg.const = ggdensity(s05.trace.all, "s2.dens", fill = "kd.strat", palette = "jco")+
    clean_theme() +
    xlim (trace.range) +
    coord_flip()#+
#    theme(plot.margin = unit(c(1,0,1,0), "lines"))


s20.marg.const = ggdensity(s20.trace.all, "s2.dens", fill = "kd.strat" , palette = "jco"
)+
    clean_theme() +
    xlim (trace.range) +
    coord_flip()#+
#    theme(plot.margin = unit(c(1,0,1,0), "lines"))



determ.trace.plot = ggplot() + 
    geom_line(aes( y = s2.dens, x = year , colour = kd.strat ),
                                         data = determ.eg.traces , stat = "identity"
                                         ) + 
    scale_y_continuous (lim = trace.range, name = "Adult breeding density
")+
    scale_x_continuous(breaks = c(0,25,50))#+
#    theme(plot.margin = unit(c(0.5,0,1,1), "lines"))


s05.trace.plot = ggplot() + geom_line(aes( y = s2.dens, x = year , colour = kd.strat),
                                      data = s05.eg.traces# , stat = "identity"
) + scale_y_continuous (lim = trace.range, name = "")+
    scale_x_continuous(breaks = c(0,25,50))#+
#    theme(plot.margin = unit(c(0.5,0,1,1), "lines"))


s20.trace.plot = ggplot() + geom_line(aes( y = s2.dens, x = year , colour = kd.strat),
                                      data = s20.eg.traces# , stat = "identity"
) + scale_y_continuous (lim = trace.range, name = "")+
    scale_x_continuous(breaks = c(0,25,50))#+
    #theme(plot.margin = unit(c(0.5,0,1,1), "lines"))



ggarrange( determ.trace.plot , determ.marg.const  , s05.trace.plot , s05.marg.const  , s20.trace.plot , s20.marg.const , 
           ncol = 6, nrow = 1, # align = "hv", 
           widths = rep(c(4, 2), 3), heights = c(1),
           common.legend = TRUE)





##########################################################################################################
##                          MAKE FIG 3 - binary partitioning
# 1- create explanitory vairables and plot config----
s1.eff.expl = array(NA, dim = dim(determ.sim)[3: length(dim(determ.sim)) ])
for (i in 1:dim(factorial.season.effort.partition)[1]){
    s1.eff.expl[,i,,,] = factorial.season.effort.partition[i,1]
}
s2.eff.expl = array(NA, dim = dim(s1.eff.expl))
for (i in 1:dim(factorial.season.effort.partition)[1]){
    s2.eff.expl[,i,,,] = factorial.season.effort.partition[i,2]
}
s3.eff.expl = array(NA, dim = dim(s1.eff.expl))
for (i in 1:dim(factorial.season.effort.partition)[1]){
    s3.eff.expl[,i,,,] = factorial.season.effort.partition[i,3]
}
s4.eff.expl = array(NA, dim = dim(s1.eff.expl))
for (i in 1:dim(factorial.season.effort.partition)[1]){
    s4.eff.expl[,i,,,] = factorial.season.effort.partition[i,4]
}
over.trig.expl = array(NA, dim = dim(s1.eff.expl))
over.trig4analysis = c(1,eff.prop.if.over)
for (i in 1:length(over.trig4analysis)){
    over.trig.expl[,,,i,] = over.trig4analysis[i]
}

#### ### ###
ctrees1 = list(NA)
juv0.minlim = 4
juv1.minlim = 4
juv2.minlim = 2
juv3.minlim = 0
juv4.minlim = 0
juv5.minlim = 0

juv0.maxlim = 11
juv1.maxlim = 8
juv2.maxlim = 8
juv3.maxlim = 12
juv4.maxlim = 12
juv5.maxlim = 6.3
# deterministic ctree analysis ----
    # juv mod = 0 dum.s3.eff. <=0.25, dumm.trig <= 3, dum.s3.eff. <= 0, dumm.trig <= 2 , dum.s4.eff. <=0.5,  ----
i = 1
dtt = data.frame(
    dum.q90. = c(determ.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
    
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. <=0.25 ,] # n = 150
ndtt2 = dtt1[!dtt1$dum.s3.eff. <=0.25 ,] # n = 60
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig , data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dumm.trig <= 3 ,] #  ,] # n = 75
ndtt3 = dtt2[!dtt2$dumm.trig <= 3 ,] # n = 75
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dum.s3.eff. <= 0,] # n = 45
ndtt4 = dtt3[!dtt3$dum.s3.eff. <= 0 ,] # n = 30
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dumm.trig <= 2  ,] # n = 30
ndtt5 = dtt4[!dtt4$dumm.trig <= 2  ,] # n = 15
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt6 = dtt5[dtt5$dum.s4.eff. <=0.5 ,] # n = 24
ndtt6 = dtt5[!dtt5$dum.s4.eff. <=0.5 ,] # n = 6
just.dum.trig = dtt6$dumm.trig.
dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt7 = dtt6[dtt6$dumm.trig. <= 1 ,] # n = 12
ndtt7 = dtt6[!dtt6$dumm.trig. <= 1 ,] # n = 12
just.dum.trig = dtt7$dumm.trig.
dtt7$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt7, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt8 = dtt7[dtt7$dum.s4.eff. <=0 ,] # n = 5
ndtt8 = dtt7[!dtt7$dum.s4.eff. <=0,] # n = 4
just.dum.trig = dtt8$dumm.trig.
dtt8$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt8, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plot ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(ndtt6[,1]) ),
    cbind("FF", as.numeric(dtt6[,1])  )
)

determ.juv0.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(determ.juv0.parted.data) = c( "split", "q")
determ.juv0.parted.data$split = as.factor(determ.juv0.parted.data$split)
determ.juv0.parted.data$q = as.numeric(determ.juv0.parted.data$q)

cut.determ.juv0.parted.data = determ.juv0.parted.data
cut.determ.juv0.parted.data.w.tot = rbind(cut.determ.juv0.parted.data , cut.determ.juv0.parted.data) 
cut.determ.juv0.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.determ.juv0.parted.data)[1]), cut.determ.juv0.parted.data$split))

cut.determ.juv0.part = ggplot(cut.determ.juv0.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv0.minlim,juv0.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

    # juv mod = 1 dumm.trig <=3, dumm.trig <=2, dum.s3.eff. <=0.25, dumm.trig <= 1, dum.s3.eff. <= 0  ----

i = 2
dtt = data.frame(
    dum.q90. = c(determ.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
    
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dumm.trig <=3 ,] # n = 105
ndtt2 = dtt1[!dtt1$dumm.trig <=3 ,]# not when!  # n = 105
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dumm.trig <=2 ,] # n = 70
ndtt3 = dtt2[!dtt2$dumm.trig <=2 ,] # n = 35
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dum.s3.eff. <=0.25 ,] # n = 50
ndtt4 = dtt3[!dtt3$dum.s3.eff. <=0.25 ,] # n = 20
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dumm.trig <= 1 ,] # n = 25
ndtt5 = dtt4[!dtt4$dumm.trig <= 1,]# not when! # n = 25
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt6 = dtt5[dtt5$dum.s3.eff. <= 0 ,] # n = 15
ndtt6 = dtt5[!dtt5$dum.s3.eff. <= 0 ,] # n = 10
just.dum.trig = dtt6$dumm.trig.
dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt7 = dtt6[dtt6$dum.s4.eff.  <= 0.25 ,] # n = 9
ndtt7 = dtt6[!dtt6$dum.s4.eff. <= 0.25 ,] # n = 6
just.dum.trig = dtt7$dumm.trig.
dtt7$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt7, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt8 = dtt7[dtt7$dum.s4.eff. <=0 ,] # n = 5
ndtt8 = dtt7[!dtt7$dum.s4.eff. <=0 ,] # n = 4
just.dum.trig = dtt8$dumm.trig.
dtt8$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt8, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))



        # violin plot ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(ndtt6[,1]) ),
    cbind("FF", as.numeric(dtt6[,1])  )
)

determ.juv1.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(determ.juv1.parted.data) = c( "split", "q")
determ.juv1.parted.data$split = as.factor(determ.juv1.parted.data$split)
determ.juv1.parted.data$q = as.numeric(determ.juv1.parted.data$q)

cut.determ.juv1.parted.data = determ.juv1.parted.data
cut.determ.juv1.parted.data.w.tot = rbind(cut.determ.juv1.parted.data , cut.determ.juv1.parted.data) 
cut.determ.juv1.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.determ.juv1.parted.data)[1]), cut.determ.juv1.parted.data$split))

cut.determ.juv1.part = ggplot(cut.determ.juv1.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv1.minlim,juv1.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

    # juv mod = 2 dum.s3.eff. > 0.25, dumm.trig <=3, dum.s3.eff. > 0.5 , dum.s3.eff. > 0.75, NA----
i = 3
dtt = data.frame(
    dum.q90. = c(determ.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. > 0.25 ,] # n = 60
ndtt2 = dtt1[!dtt1$dum.s3.eff. > 0.25 ,]# 150
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dumm.trig. <= 3 ,] # n = 30
ndtt3 = dtt2[!dtt2$dumm.trig. <= 3 ,] # n = 30
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dum.s3.eff. > 0.5  ,] # n = 12
ndtt4 = dtt3[!dtt3$dum.s3.eff. > 0.5 ,] # n = 18 
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dum.s3.eff. > 0.75 ,] # n = 3
ndtt5 = dtt4[!dtt4$dum.s3.eff. > 0.75 ,] # n = 9
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))
# 
# dtt6 = dtt5[dtt5$dum.s3.eff. > 0.5 ,] # n = 4
# ndtt6 = dtt5[!dtt5$dum.s3.eff. > 0.5 ,] # n = 6
# just.dum.trig = dtt6$dumm.trig.
# dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
# ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
# plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plots ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(dtt5[,1]) ),
    cbind("FF", as.numeric(NA)  )
)

determ.juv2.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(determ.juv2.parted.data) = c( "split", "q")
determ.juv2.parted.data$split = as.factor(determ.juv2.parted.data$split)
determ.juv2.parted.data$q = as.numeric(determ.juv2.parted.data$q)

cut.determ.juv2.parted.data = determ.juv2.parted.data
cut.determ.juv2.parted.data.w.tot = rbind(cut.determ.juv2.parted.data , cut.determ.juv2.parted.data) 
cut.determ.juv2.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.determ.juv2.parted.data)[1]), cut.determ.juv2.parted.data$split))

cut.determ.juv2.part = ggplot(cut.determ.juv2.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv2.minlim,juv2.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

    # juv mod = 3 dum.s3.eff. > 0.25 , dum.s3.eff. > 0.5 , dumm.trig <=3, dum.s3.eff. > 0.75, NA ----

i = 4
dtt = data.frame(
    dum.q90. = c(determ.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. > 0.25 ,] # n = 60
ndtt2 = dtt1[!dtt1$dum.s3.eff. > 0.25 ,]# 150
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dum.s3.eff. > 0.5 ,] # n = 24
ndtt3 = dtt2[!dtt2$dum.s3.eff. > 0.5 ,] # n = 36
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dumm.trig. <= 3 ,] # n = 12
ndtt4 = dtt3[!dtt3$dumm.trig. <= 3 ,] # n = 12
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dum.s3.eff. > 0.75 ,] # n = 3
ndtt5 = dtt4[!dtt4$dum.s3.eff. > 0.75 ,] # n = 9
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))


        # violin plots ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(dtt5[,1]) ),
    cbind("FF", as.numeric(NA ))
)

determ.juv3.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(determ.juv3.parted.data) = c( "split", "q")
determ.juv3.parted.data$split = as.factor(determ.juv3.parted.data$split)
determ.juv3.parted.data$q = as.numeric(determ.juv3.parted.data$q)

cut.determ.juv3.parted.data = determ.juv3.parted.data
cut.determ.juv3.parted.data.w.tot = rbind(cut.determ.juv3.parted.data , cut.determ.juv3.parted.data) 
cut.determ.juv3.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.determ.juv3.parted.data)[1]), cut.determ.juv3.parted.data$split))

cut.determ.juv3.part = ggplot(cut.determ.juv3.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv3.minlim,juv3.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

    # juv mod = 4 dum.s3.eff. > 0.25, dum.s3.eff. > 0.5, dumm.trig <=4, dumm.trig. <= 3 , dumm.trig. <= 1 ----

i = 5
dtt = data.frame(
    dum.q90. = c(determ.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. > 0.25 ,] # n = 60
ndtt2 = dtt1[!dtt1$dum.s3.eff. > 0.25 ,] # 150
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dum.s3.eff. > 0.5 ,] # n = 24
ndtt3 = dtt2[!dtt2$dum.s3.eff. > 0.5 ,] # n = 36
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dumm.trig. <= 4  ,] # n = 16
ndtt4 = dtt3[!dtt3$dumm.trig. <= 4 ,] # n = 8 
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dumm.trig. <= 3 ,] # n = 12
ndtt5 = dtt4[!dtt4$dumm.trig. <= 3 ,] # n = 4
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt6 = dtt5[dtt5$dumm.trig. <= 1 ,] # n = 4
ndtt6 = dtt5[!dtt5$dumm.trig. <= 1 ,] # n = 8
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plots ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(ndtt6[,1]) ),
    cbind("FF", as.numeric(dtt6[,1])  )
)

determ.juv4.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(determ.juv4.parted.data) = c( "split", "q")
determ.juv4.parted.data$split = as.factor(determ.juv4.parted.data$split)
determ.juv4.parted.data$q = as.numeric(determ.juv4.parted.data$q)

cut.determ.juv4.parted.data = determ.juv4.parted.data
cut.determ.juv4.parted.data.w.tot = rbind(cut.determ.juv4.parted.data , cut.determ.juv4.parted.data) 
cut.determ.juv4.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.determ.juv4.parted.data)[1]), cut.determ.juv4.parted.data$split))

cut.determ.juv4.part = ggplot(cut.determ.juv4.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv4.minlim,juv4.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

# low stoch ctree analysis ----
    # juv mod = 0 dum.s3.eff. <=0.25, mid.trig <= 2.25, dum.s3.eff. <=0, dum.s4.eff. <= 0.5 , dumm.trig. <= 3 ----
i = 1
dtt = data.frame(
    dum.q90. = c(s05.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
    
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. <=0.25 ,] # n = 150
ndtt2 = dtt1[!dtt1$dum.s3.eff. <=0.25 ,] # n = 60
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig , data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$mid.trig <= 2.25 ,] # n = 100
ndtt3 = dtt2[!dtt2$mid.trig <= 2.25 ,] # n = 50
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dum.s3.eff. <=0 ,] # n = 60
ndtt4 = dtt3[!dtt3$dum.s3.eff. <=0 ,] # n = 40
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dum.s4.eff. <= 0.5 ,] # n = 48
ndtt5 = dtt4[!dtt4$dum.s4.eff. <= 0.5 ,] # n = 12
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt6 = dtt5[dtt5$dumm.trig. <= 3 ,] # n = 24
ndtt6 = dtt5[!dtt5$dumm.trig. <= 3 ,] # n = 24
just.dum.trig = dtt6$dumm.trig.
dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt7 = dtt6[dtt6$dum.s4.eff. <= 0.25 ,] # n = 18
ndtt7 = dtt6[!dtt6$dum.s4.eff. <= 0.25 ,] # n = 6
just.dum.trig = dtt7$dumm.trig.
dtt7$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt7, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt8 = dtt7[dtt7$dumm.trig. <= 2 ,] # n = 9
ndtt8 = dtt7[!dtt7$dumm.trig. <= 2 ,] # n = 9
just.dum.trig = dtt8$dumm.trig.
dtt8$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt8, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt9 = dtt8[dtt8$dum.s4.eff. <= 0 ,] # n = 5
ndtt9 = dtt8[!dtt8$dum.s4.eff. <= 0 ,] # n = 4
just.dum.trig = dtt9$dumm.trig.
dtt9$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt9, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plot ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(ndtt6[,1]) ),
    cbind("FF", as.numeric(dtt6[,1])  )
)

stoch05.juv0.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(stoch05.juv0.parted.data) = c( "split", "q")
stoch05.juv0.parted.data$split = as.factor(stoch05.juv0.parted.data$split)
stoch05.juv0.parted.data$q = as.numeric(stoch05.juv0.parted.data$q)

cut.stoch05.juv0.parted.data = stoch05.juv0.parted.data
cut.stoch05.juv0.parted.data.w.tot = rbind(cut.stoch05.juv0.parted.data , cut.stoch05.juv0.parted.data) 
cut.stoch05.juv0.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.stoch05.juv0.parted.data)[1]), cut.stoch05.juv0.parted.data$split))

cut.stoch05.juv0.part = ggplot(cut.stoch05.juv0.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv0.minlim,juv0.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

    # juv mod = 1 mid.trig <= 2.25, dum.s3.eff. <=0.5, dum.s3.eff. <= 0.25, dum.s4.eff. <= 0.5, dumm.trig. <= 3 ----

i = 2
dtt = data.frame(
    dum.q90. = c(s05.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
    
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$mid.trig <= 2.25 ,] # n = 140
ndtt2 = dtt1[!dtt1$mid.trig <= 2.25 ,]# not when!  # n = 70
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dum.s3.eff. <= 0.5 ,] # n = 124
ndtt3 = dtt2[!dtt2$dum.s3.eff. <= 0.5 ,] # n = 16
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dum.s3.eff. <= 0.25 ,] # n = 100
ndtt4 = dtt3[!dtt3$dum.s3.eff. <= 0.25 ,] # n = 24
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dum.s4.eff. <= 0.5 ,] # n = 84
ndtt5 = dtt4[!dtt4$dum.s4.eff. <= 0.5 ,]# n = 16
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt6 = dtt5[dtt5$dumm.trig. <= 3 ,] # n = 42
ndtt6 = dtt5[!dtt5$dumm.trig. <= 3 ,] # n = 42
just.dum.trig = dtt6$dumm.trig.
dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt7 = dtt6[dtt6$dum.s3.eff. <= 0 ,] # n = 24
ndtt7 = dtt6[!dtt6$dum.s3.eff. <= 0 ,] # n = 18
just.dum.trig = dtt7$dumm.trig.
dtt7$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt7, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt8 = dtt7[dtt7$dum.s4.eff. <= 0.25 ,] # n = 18
ndtt8 = dtt7[!dtt7$dum.s4.eff. <= 0.25 ,] # n = 6
just.dum.trig = dtt8$dumm.trig.
dtt8$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt8, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))


        # violin plot ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(ndtt6[,1]) ),
    cbind("FF", as.numeric(dtt6[,1])  )
)

stoch05.juv1.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(stoch05.juv1.parted.data) = c( "split", "q")
stoch05.juv1.parted.data$split = as.factor(stoch05.juv1.parted.data$split)
stoch05.juv1.parted.data$q = as.numeric(stoch05.juv1.parted.data$q)

cut.stoch05.juv1.parted.data = stoch05.juv1.parted.data
cut.stoch05.juv1.parted.data.w.tot = rbind(cut.stoch05.juv1.parted.data , cut.stoch05.juv1.parted.data) 
cut.stoch05.juv1.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.stoch05.juv1.parted.data)[1]), cut.stoch05.juv1.parted.data$split))

cut.stoch05.juv1.part = ggplot(cut.stoch05.juv1.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv1.minlim,juv1.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()


    # juv mod = 2 dum.s3.eff. > 0.25, mid.trig <= 2.25, dum.s3.eff. > 0.5 , dum.s3.eff. > 0.75 , NA  ----

i = 3
dtt = data.frame(
    dum.q90. = c(s05.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. > 0.25 ,] # n = 60
ndtt2 = dtt1[!dtt1$dum.s3.eff. > 0.25 ,]# n = 150
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$mid.trig <= 2.25 ,] # n = 40
ndtt3 = dtt2[!dtt2$mid.trig <= 2.25 ,] # n = 20
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dum.s3.eff. > 0.5 ,] # n = 16
ndtt4 = dtt3[!dtt3$dum.s3.eff. > 0.5 ,] # n = 24 
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dum.s3.eff. > 0.75 ,] # n = 12
ndtt5 = dtt4[!dtt4$dum.s3.eff. > 0.75 ,] # n = 4
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

# dtt6 = dtt5[dtt5$dum.s4.eff. <= 0.25 ,] # n = 50
# ndtt6 = dtt5[!dtt5$dum.s4.eff. <= 0.25 ,] # n = 12
# just.dum.trig = dtt6$dumm.trig.
# dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
# ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
# plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plots ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(dtt5[,1]) ),
    cbind("FF", as.numeric(NA)  )
)

stoch05.juv2.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(stoch05.juv2.parted.data) = c( "split", "q")
stoch05.juv2.parted.data$split = as.factor(stoch05.juv2.parted.data$split)
stoch05.juv2.parted.data$q = as.numeric(stoch05.juv2.parted.data$q)

cut.stoch05.juv2.parted.data = stoch05.juv2.parted.data
cut.stoch05.juv2.parted.data.w.tot = rbind(cut.stoch05.juv2.parted.data , cut.stoch05.juv2.parted.data) 
cut.stoch05.juv2.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.stoch05.juv2.parted.data)[1]), cut.stoch05.juv2.parted.data$split))

cut.stoch05.juv2.part = ggplot(cut.stoch05.juv2.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv2.minlim,juv2.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

    # juv mod = 3 dum.s3.eff. > 0.25, dum.s3.eff. > 0.5, dum.s3.eff. > 0.75, dumm.trig. <= 4 , NA  ----

i = 4
dtt = data.frame(
    dum.q90. = c(s05.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
) 
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. > 0.25 ,] # n = 60
ndtt2 = dtt1[!dtt1$dum.s3.eff. > 0.25 ,]# n = 150
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dum.s3.eff. > 0.5 ,] # n = 24
ndtt3 = dtt2[!dtt2$dum.s3.eff. > 0.5 ,] # n = 36
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dum.s3.eff. > 0.75 ,] # n = 6
ndtt4 = dtt3[!dtt3$dum.s3.eff. > 0.75 ,] # n = 18 
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dumm.trig. <= 4 ,] # n = 4
ndtt5 = dtt4[!dtt4$dumm.trig. <= 4 ,] # n = 2
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))
# 
# dtt6 = dtt5[dtt5$mid.trig <= 0.25 ,] # n = 2
# ndtt6 = dtt5[!dtt5$mid.trig <= 0.25 ,] # n = 2
# just.dum.trig = dtt6$dumm.trig.
# dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
# ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
# plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

            # violin plots ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(dtt5[,1]) ),
    cbind("FF", as.numeric(NA)  )
)

stoch05.juv3.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(stoch05.juv3.parted.data) = c( "split", "q")
stoch05.juv3.parted.data$split = as.factor(stoch05.juv3.parted.data$split)
stoch05.juv3.parted.data$q = as.numeric(stoch05.juv3.parted.data$q)

cut.stoch05.juv3.parted.data = stoch05.juv3.parted.data
cut.stoch05.juv3.parted.data.w.tot = rbind(cut.stoch05.juv3.parted.data , cut.stoch05.juv3.parted.data) 
cut.stoch05.juv3.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.stoch05.juv3.parted.data)[1]), cut.stoch05.juv3.parted.data$split))

cut.stoch05.juv3.part = ggplot(cut.stoch05.juv3.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv3.minlim,juv3.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

    ## juv mod = 4 dum.s3.eff. > 0.25, dum.s3.eff. > 0.5, dumm.trig. <= 5 , dumm.trig. <= 4 , dumm.trig. <= 1 ----

i = 5
dtt = data.frame(
    dum.q90. = c(s05.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
) 
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. > 0.25 ,] # n = 60
ndtt2 = dtt1[!dtt1$dum.s3.eff. > 0.25 ,]# n = 150
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dum.s3.eff. > 0.5 ,] # n = 24
ndtt3 = dtt2[!dtt2$dum.s3.eff. > 0.5 ,] # n = 36
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dumm.trig. <= 5 ,] # n = 20
ndtt4 = dtt3[!dtt3$dumm.trig. <= 5 ,] # n = 4 
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dumm.trig. <= 4 ,] # n = 16
ndtt5 = dtt4[!dtt4$dumm.trig. <= 4 ,] # n = 4
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt6 = dtt5[dtt5$dumm.trig. <= 1 ,] # n = 4
ndtt6 = dtt5[!dtt5$dumm.trig. <= 1 ,] # n = 12
just.dum.trig = dtt6$dumm.trig.
dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plots ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(dtt5[,1]) ),
    cbind("FF", as.numeric(NA)  )
)

stoch05.juv4.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(stoch05.juv4.parted.data) = c( "split", "q")
stoch05.juv4.parted.data$split = as.factor(stoch05.juv4.parted.data$split)
stoch05.juv4.parted.data$q = as.numeric(stoch05.juv4.parted.data$q)

cut.stoch05.juv4.parted.data = stoch05.juv4.parted.data
cut.stoch05.juv4.parted.data.w.tot = rbind(cut.stoch05.juv4.parted.data , cut.stoch05.juv4.parted.data) 
cut.stoch05.juv4.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.stoch05.juv4.parted.data)[1]), cut.stoch05.juv4.parted.data$split))

cut.stoch05.juv4.part = ggplot(cut.stoch05.juv4.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv4.minlim,juv4.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

# high stoch ctree analysis ----
    # juv mod = 0 dum.s3.eff. <=0.25, mid.trig <= 2.25 , dum.s3.eff. <= 0 , dum.s4.eff. <= 0.5 , dum.s4.eff. <= 0.25 ----
i = 1
dtt = data.frame(
    dum.q90. = c(s20.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
    
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. <=0.25 ,] # n = 150
ndtt2 = dtt1[!dtt1$dum.s3.eff. <=0.25 ,] # n = 60
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig , data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$mid.trig <= 2.25 ,] # n = 100
ndtt3 = dtt2[!dtt2$mid.trig <= 2.25 ,] # n = 50
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dum.s3.eff. <= 0 ,] # n = 60
ndtt4 = dtt3[!dtt3$dum.s3.eff. <= 0 ,] # n = 40
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dum.s4.eff. <= 0.5 ,] # n = 48
ndtt5 = dtt4[!dtt4$dum.s4.eff. <= 0.5 ,] # n = 12
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt6 = dtt5[dtt5$dum.s4.eff. <= 0.25 ,] # n = 36
ndtt6 = dtt5[!dtt5$dum.s4.eff. <= 0.25 ,] # n = 12
just.dum.trig = dtt6$dumm.trig.
dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt7 = dtt6[dtt6$mid.trig <= 0.25 ,] # n = 18
ndtt7 = dtt6[!dtt6$mid.trig <= 0.25 ,] # n = 18
just.dum.trig = dtt7$dumm.trig.
dtt7$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt7, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt8 = dtt7[dtt7$dumm.trig. <= 3 ,] # n = 9
ndtt8 = dtt7[!dtt7$dum.s4.eff. <= 3 ,] # n = 9
just.dum.trig = dtt8$dumm.trig.
dtt8$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt8, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt9 = dtt8[dtt8$dum.s4.eff. <= 0 ,] # n = 5
ndtt9 = dtt8[!dtt8$dum.s4.eff. <= 0 ,] # n = 4
just.dum.trig = dtt9$dumm.trig.
dtt9$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt9, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plot ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(ndtt6[,1]) ),
    cbind("FF", as.numeric(dtt6[,1])  )
)

stoch20.juv0.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(stoch20.juv0.parted.data) = c( "split", "q")
stoch20.juv0.parted.data$split = as.factor(stoch20.juv0.parted.data$split)
stoch20.juv0.parted.data$q = as.numeric(stoch20.juv0.parted.data$q)

cut.stoch20.juv0.parted.data = stoch20.juv0.parted.data
cut.stoch20.juv0.parted.data.w.tot = rbind(cut.stoch20.juv0.parted.data , cut.stoch20.juv0.parted.data) 
cut.stoch20.juv0.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.stoch20.juv0.parted.data)[1]), cut.stoch20.juv0.parted.data$split))

cut.stoch20.juv0.part = ggplot(cut.stoch20.juv0.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv0.minlim,juv0.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()



    # juv mod = 1 mid.trig <= 2.25, dum.s3.eff. <=0.25, dum.s4.eff. <= 0.5 , dum.s3.eff. <= 0 , mid.trig <= 0.25 ----

i = 2
dtt = data.frame(
    dum.q90. = c(s20.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
    
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$mid.trig <= 2.25 ,] # n = 140
ndtt2 = dtt1[!dtt1$mid.trig <= 2.25 ,]# not when!  # n = 70
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dum.s3.eff. <= 0.25 ,] # n = 100
ndtt3 = dtt2[!dtt2$dum.s3.eff. <= 0.25 ,] # n = 40
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dum.s4.eff. <= 0.5 ,] # n = 84
ndtt4 = dtt3[!dtt3$dum.s4.eff. <= 0.5 ,] # n = 16
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dum.s3.eff. <= 0 ,] # n = 48
ndtt5 = dtt4[!dtt4$dum.s3.eff. <= 0 ,]# n = 36
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt6 = dtt5[dtt5$mid.trig <= 0.25 ,] # n = 24
ndtt6 = dtt5[!dtt5$mid.trig <= 0.25 ,] # n = 24
just.dum.trig = dtt6$dumm.trig.
dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt7 = dtt6[dtt6$dumm.trig.  <= 3 ,] # n = 12
ndtt7 = dtt6[!dtt6$dumm.trig. <= 3 ,] # n = 12
just.dum.trig = dtt7$dumm.trig.
dtt7$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt7, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt8 = dtt7[dtt7$dum.s4.eff. <= 0.25 ,] # n = 9
ndtt8 = dtt7[!dtt7$dum.s4.eff. <= 0.25 ,] # n = 3
just.dum.trig = dtt8$dumm.trig.
dtt8$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt8, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plot ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(ndtt6[,1]) ),
    cbind("FF", as.numeric(dtt6[,1])  )
)

stoch20.juv1.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(stoch20.juv1.parted.data) = c( "split", "q")
stoch20.juv1.parted.data$split = as.factor(stoch20.juv1.parted.data$split)
stoch20.juv1.parted.data$q = as.numeric(stoch20.juv1.parted.data$q)

cut.stoch20.juv1.parted.data = stoch20.juv1.parted.data
cut.stoch20.juv1.parted.data.w.tot = rbind(cut.stoch20.juv1.parted.data , cut.stoch20.juv1.parted.data) 
cut.stoch20.juv1.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.stoch20.juv1.parted.data)[1]), cut.stoch20.juv1.parted.data$split))

cut.stoch20.juv1.part = ggplot(cut.stoch20.juv1.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv1.minlim,juv1.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()



    # juv mod = 2  dum.s3.eff. > 0.25, dum.s3.eff. > 0.5 , mid.trig <= 2.25, dum.s3.eff. > 0.75 , NA ----

i = 3
dtt = data.frame(
    dum.q90. = c(s20.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. > 0.25 ,] # n = 60
ndtt2 = dtt1[!dtt1$dum.s3.eff. > 0.25 ,]# n = 150
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dum.s3.eff. > 0.5 ,] # n = 24
ndtt3 = dtt2[!dtt2$dum.s3.eff. > 0.5 ,] # n = 36
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$mid.trig <= 2.25 ,] # n = 16
ndtt4 = dtt3[!dtt3$mid.trig <= 2.25 ,] # n = 8 
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dum.s3.eff. > 0.75 ,] # n = 4
ndtt5 = dtt4[!dtt4$dum.s3.eff. > 0.75 ,] # n = 12
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

# dtt6 = dtt5[dtt5$dum.s4.eff. <= 0 ,] # n = 15
# ndtt6 = dtt5[!dtt5$dum.s4.eff. <= 0 ,] # n = 10
# just.dum.trig = dtt6$dumm.trig.
# dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
# ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
# plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plots ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(dtt5[,1]) ),
    cbind("FF", as.numeric(NA)  )
)


stoch20.juv2.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(stoch20.juv2.parted.data) = c( "split", "q")
stoch20.juv2.parted.data$split = as.factor(stoch20.juv2.parted.data$split)
stoch20.juv2.parted.data$q = as.numeric(stoch20.juv2.parted.data$q)

cut.stoch20.juv2.parted.data = stoch20.juv2.parted.data
cut.stoch20.juv2.parted.data.w.tot = rbind(cut.stoch20.juv2.parted.data , cut.stoch20.juv2.parted.data) 
cut.stoch20.juv2.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.stoch20.juv2.parted.data)[1]), cut.stoch20.juv2.parted.data$split))

cut.stoch20.juv2.part = ggplot(cut.stoch20.juv2.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv2.minlim,juv2.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

    # juv mod = 3 dum.s3.eff. > 0.25, dum.s3.eff. > 0.5, dumm.trig. <= 4 , dumm.trig. <= 3 , dumm.trig. <= 2 ----

i = 4
dtt = data.frame(
    dum.q90. = c(s20.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. >0.25 ,] # n = 60
ndtt2 = dtt1[!dtt1$dum.s3.eff. >0.25 ,]# n = 150
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dum.s3.eff. > 0.5 ,] # n = 24
ndtt3 = dtt2[!dtt2$dum.s3.eff. > 0.5 ,] # n = 36
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dumm.trig. <= 4 ,] # n = 16
ndtt4 = dtt3[!dtt3$dumm.trig. <= 4 ,] # n = 8 
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dumm.trig. <= 3 ,] # n = 12
ndtt5 = dtt4[!dtt4$dumm.trig. <= 3 ,] # n = 4
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt6 = dtt5[dtt5$dumm.trig. <= 2 ,] # n = 8
ndtt6 = dtt5[!dtt5$dumm.trig. <= 2 ,] # n = 4
just.dum.trig = dtt6$dumm.trig.
dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plots ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(ndtt6[,1]) ),
    cbind("FF", as.numeric(dtt6[,1])  )
)

stoch20.juv3.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(stoch20.juv3.parted.data) = c( "split", "q")
stoch20.juv3.parted.data$split = as.factor(stoch20.juv3.parted.data$split)
stoch20.juv3.parted.data$q = as.numeric(stoch20.juv3.parted.data$q)

cut.stoch20.juv3.parted.data = stoch20.juv3.parted.data
cut.stoch20.juv3.parted.data.w.tot = rbind(cut.stoch20.juv3.parted.data , cut.stoch20.juv3.parted.data) 
cut.stoch20.juv3.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.stoch20.juv3.parted.data)[1]), cut.stoch20.juv3.parted.data$split))

cut.stoch20.juv3.part = ggplot(cut.stoch20.juv3.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv3.minlim,juv3.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()

    # juv mod = 4 dum.s3.eff. >0.25 , dumm.trig. <= 4 , dumm.trig. <= 3 , dumm.trig. <= 2 , dumm.trig. <= 1----

i = 5
dtt = data.frame(
    dum.q90. = c(s20.analysible.s2.q90[,i,])
    ,dum.s1.eff. = c(s1.eff.expl[,,i,,])
    ,dum.s2.eff. = c(s2.eff.expl[,,i,,])
    ,dum.s3.eff. = c(s3.eff.expl[,,i,,])
    ,dum.s4.eff. = c(s4.eff.expl[,,i,,])
    ,dumm.trig. = c(over.trig.expl[,,i,,])
)
dtt1 = dtt
just.dum.trig = dtt1$dumm.trig.
dtt1$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt1, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt2 = dtt1[dtt1$dum.s3.eff. >0.25 ,] # n = 60
ndtt2 = dtt1[!dtt1$dum.s3.eff. >0.25 ,]# n = 150
just.dum.trig = dtt2$dumm.trig.
dtt2$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt2, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt3 = dtt2[dtt2$dumm.trig. <= 4 ,] # n = 40
ndtt3 = dtt2[!dtt2$dumm.trig. <= 4 ,] # n = 20
just.dum.trig = dtt3$dumm.trig.
dtt3$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt3, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt4 = dtt3[dtt3$dumm.trig. <= 3 ,] # n = 30
ndtt4 = dtt3[!dtt3$dumm.trig. <= 3 ,] # n = 10 
just.dum.trig = dtt4$dumm.trig.
dtt4$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt4, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt5 = dtt4[dtt4$dumm.trig. <= 2 ,] # n = 20
ndtt5 = dtt4[!dtt4$dumm.trig. <= 2 ,] # n = 10
just.dum.trig = dtt5$dumm.trig.
dtt5$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt5, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

dtt6 = dtt5[dtt5$dumm.trig. <= 1 ,] # n = 10
ndtt6 = dtt5[!dtt5$dumm.trig. <= 1 ,] # n = 10
just.dum.trig = dtt6$dumm.trig.
dtt6$mid.trig = (just.dum.trig-midpoint(just.dum.trig))^2
ctrees1[[i]] =  ctree(dum.q90. ~ dum.s1.eff. + dum.s2.eff. + dum.s3.eff. + dum.s4.eff. + dumm.trig. + mid.trig ,data = dtt6, control = ctree_control (mincriterion = 0, maxdepth = 1,minsplit = 1, minbucket =1) )
plot(ctrees1[[i]] , main = paste("C = ", juv.mods[i], "; B = const"))

        # violin plots ----
aaa= rbind(
    cbind("AA", as.numeric(ndtt2[,1]) ),
    cbind("BB", as.numeric(ndtt3[,1]) ),
    cbind("CC", as.numeric(ndtt4[,1]) ),
    cbind("DD", as.numeric(ndtt5[,1]) ),
    cbind("EE", as.numeric(ndtt6[,1]) ),
    cbind("FF", as.numeric(dtt6[,1])  )
)

stoch20.juv4.parted.data = data.frame (aaa[,1], as.numeric(aaa[,2]) )
names(stoch20.juv4.parted.data) = c( "split", "q")
stoch20.juv4.parted.data$split = as.factor(stoch20.juv4.parted.data$split)
stoch20.juv4.parted.data$q = as.numeric(stoch20.juv4.parted.data$q)

cut.stoch20.juv4.parted.data = stoch20.juv4.parted.data
cut.stoch20.juv4.parted.data.w.tot = rbind(cut.stoch20.juv4.parted.data , cut.stoch20.juv4.parted.data) 
cut.stoch20.juv4.parted.data.w.tot$split = factor(c(rep("0.1", dim(cut.stoch20.juv4.parted.data)[1]), cut.stoch20.juv4.parted.data$split))

cut.stoch20.juv4.part = ggplot(cut.stoch20.juv4.parted.data.w.tot, aes(x=split, y=q)) + 
    ylim(juv4.minlim,juv4.maxlim) +
    geom_violin(trim=T, scale = "width", fill= "grey") + geom_boxplot(width=0.3) +   theme_classic()


# arrange plot ----

#cut.determ.juv0.part cut.stoch05.juv0.part, cut.stoch20.juv0.part
#cut.determ.juv1.part, cut.stoch05.juv1.part, cut.stoch20.juv1.part
#cut.determ.juv2.part, cut.stoch05.juv2.part, cut.stoch20.juv2.part
#cut.determ.juv3.part, cut.stoch05.juv3.part, cut.stoch20.juv3.part
#cut.determ.juv4.part, cut.stoch05.juv4.part, cut.stoch20.juv4.part


grid.arrange(cut.determ.juv0.part, cut.stoch05.juv0.part, cut.stoch20.juv0.part,
             cut.determ.juv1.part, cut.stoch05.juv1.part, cut.stoch20.juv1.part,
             cut.determ.juv2.part, cut.stoch05.juv2.part, cut.stoch20.juv2.part)#,ncol=3)
#####
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#

# #


####################################################################################################
##                          MAKE FIG 4 - imperfect monitoring
# boxplots of density for different monitoring effort for best juv stratedgy -----
n.iterations = 100
setwd(this.wd)
interesting.nrafts = seq(from = 2, to = 40, by = 2) #1:10
nu.percentil90s2.Apop = array (NA, dim = c( n.seasonal.eff.combos, length(juv.mods), length(eff.prop.if.over) , length(interesting.nrafts), n.iterations))
set.seed(121212121)
# run monitoring sims
# determ ---- 
    # juvmod 0   ----
interested.juv.mods = 0
which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
which(determ.analysible.s2.q90[,which.interesting.juv.mods,] == sort(determ.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
which.interesting.ss = c(17,35)
interested.eff.overs = 2:3
which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)

b = rep (mean.b , time = max.time)
unresponsive.q90 = determ.analysible.s2.q90[ which.interesting.ss , which.interesting.juv.mods , 1]
sim.kd.tiggers = determ.trigs

source(monitoirng.source.code)
save(nu.percentil90s2.Apop, file = paste(this.wd, "\\determ.juv0.monitor", sep = ""))

    # juvmod 1 ----
interested.juv.mods = 1
which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
which(determ.analysible.s2.q90[,which.interesting.juv.mods,] == sort(determ.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
which.interesting.ss = c(17,35)
interested.eff.overs = 2:3
which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)

b = rep (mean.b , time = max.time)
unresponsive.q90 = determ.analysible.s2.q90[ which.interesting.ss , which.interesting.juv.mods , 1]
sim.kd.tiggers = determ.trigs

source(monitoirng.source.code)
save(nu.percentil90s2.Apop, file = paste(this.wd, "\\determ.juv1.monitor", sep = ""))

    # juvmod 2 ----
interested.juv.mods = 2
which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
which(determ.analysible.s2.q90[,which.interesting.juv.mods,] == sort(determ.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
which.interesting.ss = c(3,35)
interested.eff.overs = 2:3
which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)

b = rep (mean.b , time = max.time)
unresponsive.q90 = determ.analysible.s2.q90[ which.interesting.ss , which.interesting.juv.mods , 1]
sim.kd.tiggers = determ.trigs

source(monitoirng.source.code)
save(nu.percentil90s2.Apop, file = paste(this.wd, "\\determ.juv2.monitor", sep = ""))

# s05 ----
    # juvmod 0   ----
interested.juv.mods = 0
which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
which(s05.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s05.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
which.interesting.ss = c(17,35)
interested.eff.overs = 2:3
which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)

b = stochb.05less1
unresponsive.q90 = s05.analysible.s2.q90[ which.interesting.ss , which.interesting.juv.mods , 1]
sim.kd.tiggers = s05.trigs

source(monitoirng.source.code)
save(nu.percentil90s2.Apop, file = paste(this.wd, "\\s05.juv0.monitor", sep = ""))

    # juvmod 1 ----
interested.juv.mods = 1
which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
which(s05.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s05.analysible.s2.q90[,which.interesting.juv.mods,])[3], arr.ind = T)
which.interesting.ss = c(17,35)
interested.eff.overs = 2:3
which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)

b = stochb.05less1
unresponsive.q90 = s05.analysible.s2.q90[ which.interesting.ss , which.interesting.juv.mods , 1]
sim.kd.tiggers = s05.trigs

source(monitoirng.source.code)
save(nu.percentil90s2.Apop, file = paste(this.wd, "\\s05.juv1.monitor", sep = ""))

    # juvmod 2 ----
interested.juv.mods = 2
which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
which(s05.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s05.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
which.interesting.ss = c(3,35)
interested.eff.overs = 2:5
which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)

b = stochb.05less1
unresponsive.q90 = s05.analysible.s2.q90[ which.interesting.ss , which.interesting.juv.mods , 1]
sim.kd.tiggers = s05.trigs

source(monitoirng.source.code)
save(nu.percentil90s2.Apop, file = paste(this.wd, "\\s05.juv2.monitor", sep = ""))

# s20 ----
    # juvmod 0   ----
interested.juv.mods = 0
which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
which(s20.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s20.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
which.interesting.ss = c(17,35)
interested.eff.overs = 2:3
which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)

b = stochb.05less1
unresponsive.q90 = s20.analysible.s2.q90[ which.interesting.ss , which.interesting.juv.mods , 1]
sim.kd.tiggers = s20.trigs

source(monitoirng.source.code)
save(nu.percentil90s2.Apop, file = paste(this.wd, "\\s20.juv0.monitor", sep = ""))

    # juvmod 1  ----
interested.juv.mods = 1
which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
which(s20.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s20.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
which.interesting.ss = c(17,35)
interested.eff.overs = 2:3
which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)

b = stochb.05less1
unresponsive.q90 = s20.analysible.s2.q90[ which.interesting.ss , which.interesting.juv.mods , 1]
sim.kd.tiggers = s20.trigs

source(monitoirng.source.code)
save(nu.percentil90s2.Apop, file = paste(this.wd, "\\s20.juv1.monitor", sep = ""))

    # juvmod 2 ----
interested.juv.mods = 2
which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
which(s20.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s20.analysible.s2.q90[,which.interesting.juv.mods,])[2], arr.ind = T)
which.interesting.ss = c(3,35)
interested.eff.overs = 2:5
which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)

b = stochb.05less1
unresponsive.q90 = s20.analysible.s2.q90[ which.interesting.ss , which.interesting.juv.mods , 1]
sim.kd.tiggers = s20.trigs

source(monitoirng.source.code)
save(nu.percentil90s2.Apop, file = paste(this.wd, "\\s20.juv2.monitor", sep = ""))

#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-# ----
#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#--#-#-#-#-#-#-#-#-#-#--#-#-#-#-#--#-#-#-#-#-#-#--#-#-#-#-#
# load and plot monitoring sims ----
which.interesting.ss.to.plot = 2 # 1 = juv mod specialised ss, 2 = even spread of effort
x.axis.labels = c(0, rep("",4),10, rep("",4),20, rep("",4),30, rep("",4),40)
x.axis.breaks = seq(0,40, by = 2)
pd <- position_dodge(0.5) # move them .05 to the left and right

# determ ----
    # juvmod 0   ----
interested.juv.mods = 0
which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
which(determ.analysible.s2.q90[,which.interesting.juv.mods,] == sort(determ.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
which.interesting.ss = c(17,35)
interested.eff.overs = 2:3
which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)

b = rep (mean.b , time = max.time)
unresponsive.q90 = determ.analysible.s2.q90[ which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , 1]
sim.kd.tiggers = determ.trigs

load(paste(this.wd, "\\determ.juv0.monitor", sep = ""))

monitoring.df = data.frame (q90 = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                            eff =  rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                            n.rafts = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))))  

monitoring.df$eff = as.factor(rep(interested.eff.overs, times = length(interesting.nrafts)*n.iterations))
monitoring.df$n.rafts = rep(rep(interesting.nrafts, each = length(interested.eff.overs)), times =n.iterations )
monitoring.df$q90 = c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,])

monitoring.df$eff = as.factor(monitoring.df$eff)
monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
monitoring.df$q90 = as.numeric(monitoring.df$q90)

determ.juv0.monitor.plot = 
    ggplot(monitoring.df, aes( x = n.rafts , y = q90, fill = eff )) + 
    #  ylim(4.7,8.5) +
    geom_boxplot(outlier.size = 0.07) +
    geom_hline(yintercept = unresponsive.q90, lwd = 1, lty = 1) +      
    theme_classic()


    # juvmod 1 ----
interested.juv.mods = 1
    which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
    which(determ.analysible.s2.q90[,which.interesting.juv.mods,] == sort(determ.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
    which.interesting.ss = c(17,35)
    interested.eff.overs = 2:3
    which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)
    
    b = rep (mean.b , time = max.time)
    unresponsive.q90 = determ.analysible.s2.q90[ which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , 1]
    sim.kd.tiggers = determ.trigs
    
    load(paste(this.wd, "\\determ.juv1.monitor", sep = ""))
    
    monitoring.df = data.frame (q90 = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                eff =  rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                n.rafts = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))))  
    
    monitoring.df$eff = as.factor(rep(interested.eff.overs, times = length(interesting.nrafts)*n.iterations))
    monitoring.df$n.rafts = rep(rep(interesting.nrafts, each = length(interested.eff.overs)), times =n.iterations )
    monitoring.df$q90 = c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,])
    
    monitoring.df$eff = as.factor(monitoring.df$eff)
    monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
    monitoring.df$q90 = as.numeric(monitoring.df$q90)
    
    plot.df = data.frame(n.rafts = rep (unique(monitoring.df$n.rafts), length(unique(monitoring.df$eff))),
                         eff = rep(unique(monitoring.df$eff), each = length(unique(monitoring.df$n.rafts))),
                         q.med =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), median )),
                         q.97.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.975 )),
                         q.02.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.025 ))
                         
    )
    

    determ.juv1.monitor.plot =ggplot(plot.df, aes(x=n.rafts, y=q.med, colour=eff)) + 
        geom_errorbar(aes(ymin=q.02.5, ymax=q.97.5), width=.1, position=pd) +
        geom_hline(yintercept = unresponsive.q90, lwd = 1, lty = 1) +      
        geom_point(position=pd) 
        
 
    # juvmod 2 ----
    interested.juv.mods = 2
    which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
    which(determ.analysible.s2.q90[,which.interesting.juv.mods,] == sort(determ.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
    which.interesting.ss = c(3,35)
    interested.eff.overs = 2:3
    which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)
    
    b = rep (mean.b , time = max.time)
    unresponsive.q90 = determ.analysible.s2.q90[ which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , 1]
    sim.kd.tiggers = determ.trigs
    
    load(paste(this.wd, "\\determ.juv2.monitor", sep = ""))
    
    monitoring.df = data.frame (q90 = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                eff =  rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                n.rafts = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))))  
    
    monitoring.df$eff = as.factor(rep(interested.eff.overs, times = length(interesting.nrafts)*n.iterations))
    monitoring.df$n.rafts = rep(rep(interesting.nrafts, each = length(interested.eff.overs)), times =n.iterations )
    monitoring.df$q90 = c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,])
    
    monitoring.df$eff = as.factor(monitoring.df$eff)
    monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
    monitoring.df$q90 = as.numeric(monitoring.df$q90)
    
    determ.juv2.monitor.plot = 
        ggplot(monitoring.df, aes( x = n.rafts , y = q90, fill = eff )) + 
        #  ylim(4.7,8.5) +
        geom_boxplot(outlier.size = 0.07) +
        geom_hline(yintercept = unresponsive.q90, lwd = 1, lty = 1) +      
        theme_classic()
    
# s05 ----
    # juvmod 0   ----
    interested.juv.mods = 0
    which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
    which(s05.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s05.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
    which.interesting.ss = c(17,35)
    interested.eff.overs = 2:3
    which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)
    
    b = stochb.05less1
    unresponsive.q90 = s05.analysible.s2.q90[ which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , 1]
    sim.kd.tiggers = s05.trigs
    
    load(paste(this.wd, "\\s05.juv0.monitor", sep = ""))
    
    monitoring.df = data.frame (q90 = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                eff =  rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                n.rafts = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))))  
    
    monitoring.df$eff = as.factor(rep(interested.eff.overs, times = length(interesting.nrafts)*n.iterations))
    monitoring.df$n.rafts = rep(rep(interesting.nrafts, each = length(interested.eff.overs)), times =n.iterations )
    monitoring.df$q90 = c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,])
    
    monitoring.df$eff = as.factor(monitoring.df$eff)
    monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
    monitoring.df$q90 = as.numeric(monitoring.df$q90)
    
    plot.df = data.frame(n.rafts = rep (unique(monitoring.df$n.rafts), length(unique(monitoring.df$eff))),
                         eff = rep(unique(monitoring.df$eff), each = length(unique(monitoring.df$n.rafts))),
                         q.med =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), median )),
                         q.97.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.975 )),
                         q.02.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.025 ))
    )

    s05.juv0.monitor.plot = ggplot(plot.df, aes(x=n.rafts, y=q.med, colour=eff)) + 
        geom_errorbar(aes(ymin=q.02.5, ymax=q.97.5), width=.1, position=pd) +
        geom_hline(yintercept = unresponsive.q90, lwd = 0.1, lty = 2) +      
        scale_x_discrete( labels = x.axis.labels, breaks = x.axis.breaks, name = "Number of observations")+
                geom_point(position=pd, size = 1) +
        scale_color_viridis(discrete = T) +
        theme_classic()

    # juvmod 1 ----
    interested.juv.mods = 1
    which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
    which(s05.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s05.analysible.s2.q90[,which.interesting.juv.mods,])[3], arr.ind = T)
    which.interesting.ss = c(17,35)
    interested.eff.overs = 2:3
    which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)
    
    b = stochb.05less1
    unresponsive.q90 = s05.analysible.s2.q90[ which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , 1]
    sim.kd.tiggers = s05.trigs
    
    load(paste(this.wd, "\\s05.juv1.monitor", sep = ""))
    
    monitoring.df = data.frame (q90 = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                eff =  rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                n.rafts = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))))  
    
    monitoring.df$eff = as.factor(rep(interested.eff.overs, times = length(interesting.nrafts)*n.iterations))
    monitoring.df$n.rafts = rep(rep(interesting.nrafts, each = length(interested.eff.overs)), times =n.iterations )
    monitoring.df$q90 = c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,])
    
    monitoring.df$eff = as.factor(monitoring.df$eff)
    monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
    monitoring.df$q90 = as.numeric(monitoring.df$q90)
  
    plot.df = data.frame(n.rafts = rep (unique(monitoring.df$n.rafts), length(unique(monitoring.df$eff))),
                         eff = rep(unique(monitoring.df$eff), each = length(unique(monitoring.df$n.rafts))),
                         q.med =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), median )),
                         q.97.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.975 )),
                         q.02.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.025 ))
                         
    )
    
    s05.juv1.monitor.plot =ggplot(plot.df, aes(x=n.rafts, y=q.med, colour=eff)) + 
        geom_errorbar(aes(ymin=q.02.5, ymax=q.97.5), width=.1, position=pd) +
        geom_hline(yintercept = unresponsive.q90, lwd = 0.1, lty = 2) +      
        scale_x_discrete( labels = x.axis.labels, breaks = x.axis.breaks, name = "Number of observations")+
        scale_color_viridis(discrete = T) +
        geom_point(position=pd, size = 1) +
        theme_classic()
      
    
    # juvmod 2 ----
    interested.juv.mods = 2
    which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
    which(s05.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s05.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
    which.interesting.ss = c(3,35)
    interested.eff.overs = 2:3
    which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)
    
    b = stochb.05less1
    unresponsive.q90 = s05.analysible.s2.q90[ which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , 1]
    sim.kd.tiggers = s05.trigs
    
    load(paste(this.wd, "\\s05.juv2.monitor", sep = ""))
    
    monitoring.df = data.frame (q90 = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                eff =  rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                n.rafts = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))))  
    
    monitoring.df$eff = as.factor(rep(interested.eff.overs, times = length(interesting.nrafts)*n.iterations))
    monitoring.df$n.rafts = rep(rep(interesting.nrafts, each = length(interested.eff.overs)), times =n.iterations )
    monitoring.df$q90 = c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,])
    
    monitoring.df$eff = as.factor(monitoring.df$eff)
    monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
    monitoring.df$q90 = as.numeric(monitoring.df$q90)

    plot.df = data.frame(n.rafts = rep (unique(monitoring.df$n.rafts), length(unique(monitoring.df$eff))),
                         eff = rep(unique(monitoring.df$eff), each = length(unique(monitoring.df$n.rafts))),
                         q.med =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), median )),
                         q.97.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.975 )),
                         q.02.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.025 ))
                         
    )
    
    s05.juv2.monitor.plot =ggplot(plot.df, aes(x=n.rafts, y=q.med, colour=eff)) + 
        geom_errorbar(aes(ymin=q.02.5, ymax=q.97.5), width=.1, position=pd) +
        geom_hline(yintercept = unresponsive.q90, lwd = 0.1, lty = 2) +      
        scale_color_viridis(discrete = T) +
        scale_x_discrete( labels = x.axis.labels, breaks = x.axis.breaks, name = "Number of observations")+
                geom_point(position=pd, size = 1) +
        theme_classic()
    
# s20 ----
    # juvmod 0   ----
    interested.juv.mods = 0
    which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
    which(s20.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s20.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
    which.interesting.ss = c(17,35)
    interested.eff.overs = 2:3
    which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)
    
    b = stochb.05less1
    unresponsive.q90 = s20.analysible.s2.q90[ which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , 1]
    sim.kd.tiggers = s20.trigs
    
    load(paste(this.wd, "\\s20.juv0.monitor", sep = ""))
    
    monitoring.df = data.frame (q90 = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                eff =  rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                n.rafts = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))))  
    
    monitoring.df$eff = as.factor(rep(interested.eff.overs, times = length(interesting.nrafts)*n.iterations))
    monitoring.df$n.rafts = rep(rep(interesting.nrafts, each = length(interested.eff.overs)), times =n.iterations )
    monitoring.df$q90 = c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,])
    
    monitoring.df$eff = as.factor(monitoring.df$eff)
    monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
    monitoring.df$q90 = as.numeric(monitoring.df$q90)
    
     plot.df = data.frame(n.rafts = rep (unique(monitoring.df$n.rafts), length(unique(monitoring.df$eff))),
                         eff = rep(unique(monitoring.df$eff), each = length(unique(monitoring.df$n.rafts))),
                         q.med =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), median )),
                         q.97.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.975 )),
                         q.02.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.025 ))
                         
    )
    
    s20.juv0.monitor.plot =ggplot(plot.df, aes(x=n.rafts, y=q.med, colour=eff)) + 
        geom_errorbar(aes(ymin=q.02.5, ymax=q.97.5), width=.1, position=pd) +
        geom_hline(yintercept = unresponsive.q90, lwd = 0.1, lty = 2) +      
        scale_color_viridis(discrete = T) +
        scale_x_discrete( labels = x.axis.labels, breaks = x.axis.breaks, name = "Number of observations")+
                geom_point(position=pd, size = 1) +
        theme_classic()
    
    # juvmod 1  ----
    interested.juv.mods = 1
    which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
    which(s20.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s20.analysible.s2.q90[,which.interesting.juv.mods,])[1], arr.ind = T)
    which.interesting.ss = c(17,35)
    interested.eff.overs = 2:3
    which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)
    
    b = stochb.05less1
    unresponsive.q90 = s20.analysible.s2.q90[ which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , 1]
    sim.kd.tiggers = s20.trigs
    
    load(paste(this.wd, "\\s20.juv1.monitor", sep = ""))
    
    monitoring.df = data.frame (q90 = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                eff =  rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                n.rafts = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))))  
    
    monitoring.df$eff = as.factor(rep(interested.eff.overs, times = length(interesting.nrafts)*n.iterations))
    monitoring.df$n.rafts = rep(rep(interesting.nrafts, each = length(interested.eff.overs)), times =n.iterations )
    monitoring.df$q90 = c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,])
    
    monitoring.df$eff = as.factor(monitoring.df$eff)
    monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
    monitoring.df$q90 = as.numeric(monitoring.df$q90)

    plot.df = data.frame(n.rafts = rep (unique(monitoring.df$n.rafts), length(unique(monitoring.df$eff))),
                         eff = rep(unique(monitoring.df$eff), each = length(unique(monitoring.df$n.rafts))),
                         q.med =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), median )),
                         q.97.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.975 )),
                         q.02.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.025 ))
                         
    )
    
    
    s20.juv1.monitor.plot =ggplot(plot.df, aes(x=n.rafts, y=q.med, colour=eff)) + 
        geom_errorbar(aes(ymin=q.02.5, ymax=q.97.5), width=.1, position=pd) +
        geom_hline(yintercept = unresponsive.q90, lwd = 0.1, lty = 2) +      
        scale_color_viridis(discrete = T) +
        scale_x_discrete( labels = x.axis.labels, breaks = x.axis.breaks, name = "Number of observations")+
                geom_point(position=pd, size = 1) +
        theme_classic()
    
    
    # juvmod 2 ----
    interested.juv.mods = 2
    which.interesting.juv.mods = which(juv.mods %in% interested.juv.mods)
    which(s20.analysible.s2.q90[,which.interesting.juv.mods,] == sort(s20.analysible.s2.q90[,which.interesting.juv.mods,])[2], arr.ind = T)
    which.interesting.ss = c(3,35)
    interested.eff.overs = 2:3
    which.interesting.eff.overs = which(eff.prop.if.over %in% interested.eff.overs)
    
    b = stochb.05less1
    unresponsive.q90 = s20.analysible.s2.q90[ which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , 1]
    sim.kd.tiggers = s20.trigs
    
    load(paste(this.wd, "\\s20.juv2.monitor", sep = ""))
    
    monitoring.df = data.frame (q90 = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                eff =  rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))), 
                                n.rafts = rep (NA, times = length (c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,]))))  
    
    monitoring.df$eff = as.factor(rep(interested.eff.overs, times = length(interesting.nrafts)*n.iterations))
    monitoring.df$n.rafts = rep(rep(interesting.nrafts, each = length(interested.eff.overs)), times =n.iterations )
    monitoring.df$q90 = c(nu.percentil90s2.Apop[which.interesting.ss[which.interesting.ss.to.plot] , which.interesting.juv.mods , which.interesting.eff.overs,,])
    
    monitoring.df$eff = as.factor(monitoring.df$eff)
    monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
    monitoring.df$q90 = as.numeric(monitoring.df$q90)

    plot.df = data.frame(n.rafts = rep (unique(monitoring.df$n.rafts), length(unique(monitoring.df$eff))),
                         eff = rep(unique(monitoring.df$eff), each = length(unique(monitoring.df$n.rafts))),
                         q.med =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), median , na.rm =T )),
                         q.97.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.975 , na.rm =T )),
                         q.02.5 =c(tapply(monitoring.df$q90, list(monitoring.df$n.rafts,monitoring.df$eff), quantile, probs = 0.025 , na.rm =T ))
                         )
    
    s20.juv2.monitor.plot =ggplot(plot.df, aes(x=n.rafts, y=q.med, colour=eff)) + 
        geom_errorbar(aes(ymin=q.02.5, ymax=q.97.5), width=.1, position=pd) +
        geom_hline(yintercept = unresponsive.q90, lwd = 0.1, lty = 2) +      
        scale_color_viridis(discrete = T) +
        scale_x_discrete( labels = x.axis.labels, breaks = x.axis.breaks, name = "Number of observations")+
                geom_point(position=pd, size = 1) +
        theme_classic()
    
    
    
    
    
    
#######################  

# plot the lot ----
# determ.juv0.monitor.plot
# determ.juv1.monitor.plot
# determ.juv2.monitor.plot
# 
# s05.juv0.monitor.plot
# s05.juv1.monitor.plot
# s05.juv2.monitor.plot
# 
# 
# s20.juv0.monitor.plot
# s20.juv1.monitor.plot
# s20.juv2.monitor.plot

    grid.arrange(s05.juv0.monitor.plot, s20.juv0.monitor.plot,
                 s05.juv1.monitor.plot, s20.juv1.monitor.plot,
                 s05.juv2.monitor.plot, s20.juv2.monitor.plot, ncol=2 , common.legend = TRUE)
    grid.arrange(s05.juv1.monitor.plot, s20.juv1.monitor.plot,
                  ncol=2 , common.legend = TRUE)
    
    ggarrange( s05.juv1.monitor.plot , s20.juv1.monitor.plot, 
               ncol = 2, nrow = 1, heights = c(1),
               common.legend = TRUE)
    
# lookign at detection probabilities of interest

prob.detect.1indiv = prop.caught.w.1.trap
s05.detections.trigger = 1-((1-prob.detect.1indiv)^s05.trigs) 
s20.detections.trigger = 1-((1-prob.detect.1indiv)^s20.trigs) 

s05.detections.trigger[,c(35), 1:3,1:2,]
s20.detections.trigger[,c(35), 1:3,1:2,]


# ----
# some statistics for paper -----


juv5.determ.tot.median = median(cut.determ.juv5.parted.data.w.tot$q[cut.determ.juv5.parted.data.w.tot$split==0.1])
juv5.determ.halfsummer.median = median(cut.determ.juv5.parted.data.w.tot$q[as.numeric(cut.determ.juv5.parted.data.w.tot$split)>2])
juv5.determ.halfsummer.median/juv5.determ.tot.median

juv5.stoch05.tot.median = median(cut.stoch05.juv5.parted.data.w.tot$q[cut.stoch05.juv5.parted.data.w.tot$split==0.1])
juv5.stoch05.halfsummer.median = median(cut.stoch05.juv5.parted.data.w.tot$q[as.numeric(cut.stoch05.juv5.parted.data.w.tot$split)>3])
juv5.stoch05.halfsummer.median/juv5.stoch05.tot.median


juv5.stoch20.tot.median = median(cut.stoch20.juv5.parted.data.w.tot$q[cut.stoch20.juv5.parted.data.w.tot$split==0.1])
juv5.stoch20.halfsummer.median = median(cut.stoch20.juv5.parted.data.w.tot$q[as.numeric(cut.stoch20.juv5.parted.data.w.tot$split)>3])
juv5.stoch20.halfsummer.median/juv5.stoch20.tot.median

#standard deviation of determ scenarios by juv mod
sd(cut.determ.juv0.parted.data$q)
sd(cut.determ.juv1.parted.data$q)
sd(cut.determ.juv2.parted.data$q, na.rm = T)

range(cut.determ.juv0.parted.data$q)
range(cut.determ.juv1.parted.data$q)
range(cut.determ.juv2.parted.data$q, na.rm = T)

max(cut.determ.juv0.parted.data$q) - min(cut.determ.juv0.parted.data$q)
max(cut.determ.juv1.parted.data$q) - min(cut.determ.juv1.parted.data$q)
max(cut.determ.juv2.parted.data$q, na.rm = T) - min(cut.determ.juv2.parted.data$q, na.rm = T)


########## makign figure 4
 



setwd(wd.s05)

load("analysable.s2.q90.data.rda") 
# juv mod = 1
load("1000percentil90s2.Apop.juvmod1.rda") # percentil90s2.Apop.juvmod1

unresponsive.q90 = analysible.s2.q90[2,2,1]


workinq90 = percentil90s2.Apop.juvmod1[,2,2,2:3,,,]

monitor.dat.eff = array(NA, dim = dim(workinq90)) # ss, juvmod, eff, nrafts, iterations
monitor.dat.eff[1,,] = 2
monitor.dat.eff[2,,] = 3

monitor.dat.rafts = array(NA, dim = dim(workinq90)) # ss, juvmod, eff, nrafts, iterations
for ( i in 1:(dim(monitor.dat.rafts)[2])){
    monitor.dat.rafts[,i,] = interesting.nrafts[i]
}

aaa= cbind(workinq90 , monitor.dat.eff, monitor.dat.rafts )
monitoring.df = data.frame (aaa )
names(monitoring.df) = c( "q90", "eff", "n.rafts")

monitoring.df$eff = as.factor(monitoring.df$eff)
monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
monitoring.df$q90 = as.numeric(monitoring.df$q90)

monitoring.df.juv1.ss2 = monitoring.df

stoch05.juv1.monitor.plot = 
    ggplot(monitoring.df.juv1.ss2, aes( x = n.rafts , y = q90, fill = eff )) + 
    ylim(4.7,8.5) +
    geom_boxplot(outlier.size = 0.07) +
    geom_hline(yintercept = unresponsive.q90, lwd = 1, lty = 1) +      
    theme_classic()


# STOCHASTIC 20 -----------

setwd(wd.s20)

load("analysable.s2.q90.data.rda") 
# juv mod = 1
load("1000percentil90s2.Apop.juvmod1.rda") # percentil90s2.Apop.juvmod1

unresponsive.q90 = analysible.s2.q90[2,2,1]


workinq90 = percentil90s2.Apop.juvmod1[,2,2,2:3,,,]

monitor.dat.eff = array(NA, dim = dim(workinq90)) # ss, juvmod, eff, nrafts, iterations
monitor.dat.eff[1,,] = 2
monitor.dat.eff[2,,] = 3

monitor.dat.rafts = array(NA, dim = dim(workinq90)) # ss, juvmod, eff, nrafts, iterations
for ( i in 1:(dim(monitor.dat.rafts)[2])){
    monitor.dat.rafts[,i,] = interesting.nrafts[i]
}

aaa= cbind(workinq90 , monitor.dat.eff, monitor.dat.rafts )
monitoring.df = data.frame (aaa )
names(monitoring.df) = c( "q90", "eff", "n.rafts")

monitoring.df$eff = as.factor(monitoring.df$eff)
monitoring.df$n.rafts = as.factor(monitoring.df$n.rafts)
monitoring.df$q90 = as.numeric(monitoring.df$q90)

monitoring.df.juv1.ss2 = monitoring.df

stoch20.juv1.monitor.plot = 
    ggplot(monitoring.df.juv1.ss2, aes( x = n.rafts , y = q90, fill = eff )) + 
    geom_boxplot(outlier.size = 0.07) + 
    ylim(4.7,8.5) +
    geom_hline(yintercept = unresponsive.q90, lwd = 1, lty = 1) +      
    theme_classic()


