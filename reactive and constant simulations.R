# script file name "reactive and constant simulations.R"
# this code is sourced in the script "stage structured invasive management.R"
# save this file as "reactive and constant simulations.R" and run "stage structured invasive management.R"

###### simulation for perfect reaction
##      initial prep simulation to find a starting point for reactive density triggers estimation       ##
##                  and calculate the effort spent over the analysis period for each                    ##
##########################################################################################################
sim5.begin = array(NA, dim = c(max.time, length(init.state), n.seasonal.eff.combos, length(juv.mods) )) #  dime names = (time, state bits , seasonal effort combos, juv.trapping.mod)


# run sims to estiamte starting point triggers  - no reactive trapping
for (which.juvmod in 1:length(juv.mods)){
    #parameterising juvenile trapping chances and trapping coefficients
    juvinile.traping.modifier = juv.mods[which.juvmod]
    for (which.seasonal.strat in 1:n.seasonal.eff.combos){
        trap.effort = (rep(yearly.trapping.effort*seasonal.step.effort.mod[which.seasonal.strat,], times = years+ burn.in.years))   # select seasonal distribution of effort and replicate for years
        state = matrix(c(init.state ,rep(rep(NA, times = length(init.state)), times = max.time-1)),ncol = length(init.state), byrow = T) 
        for (t in 1:(max.time-1)){
            state[t+1,] = state[t,]
            state[t+1,] = comp.reset.trap.eat()
            state[t+1,] = comp.J.pred.di.birth(  )
            state[t+1,] = comp.A.pred.di.mort()
            state[t+1,] = comp.J.pred.di.mort()
            state[t+1,] = comp.better.trapping(juv.trap.chances = juvinile.traping.modifier, current.effort = trap.effort[t], prob.catch1 = prop.caught.w.1.trap)   
            state[t+1,] = comp.simplest.dd.settlement( ) 
            state[t+1,] = comp.kill.all.J()                            # prey consumption
        }
        sim5.begin[,,which.seasonal.strat,which.juvmod] = state                                         # add time serise to output object 
    }
}

# find mean from each to act as startign point for trigger estimation
sim5.tot.pred.init.triggers.byseason = array(NA ,dim = c(length(monitoring.season), n.seasonal.eff.combos, length(juv.mods), length (eff.prop.if.over), length(juv.misID.rate))) # for total predator density trigger
for(i in 1:length(monitoring.season)){
    for (iii in 1:n.seasonal.eff.combos){
        for(iv in 1:length(juv.mods)){
            for(v in 1:length(eff.prop.if.less)){
                for(vi in 1:length(juv.misID.rate)){
                    sim5.tot.pred.init.triggers.byseason[i, iii, iv, v, vi] = mean(apply(
                        cbind(sim5.begin[analysable.steps[which.season.is.it[analysable.steps]==i], c(2), iii,iv], juv.misID.rate[vi] * sim5.begin[analysable.steps[which.season.is.it[analysable.steps]==i], c(3), iii,iv])
                        ,1 , sum))
                }
            }
        }
    }
}





# FIND TOTAL EFFORT RESULTING FROM THE USE OF EACH INITIAL TRIGGER FOR REACTIVE STRATEDGIES
sim5.firsteffort.season.spec.trigger = array(NA, dim = dim(sim5.tot.pred.init.triggers.byseason) )

for( which.misID in 1: length(juv.misID.rate)){                 #juv misID
    print(which.misID)
    for (which.juvmod in 1:length(juv.mods)){                       #juv trapping
        print(which.juvmod/length(juv.mods))
        for (which.seasonal.strat in 1:n.seasonal.eff.combos){          #seasonal stratedgy
            for(which.eff.mod in 1: length(eff.prop.if.over)){              #trigger modifiers
                juvinile.traping.modifier = juv.mods[which.juvmod]
                pre.trigger.effort =   as.numeric (rep(yearly.trapping.effort*seasonal.step.effort.mod[which.seasonal.strat,], times = years + burn.in.years))   # select seasonal distribution of effort and replicate for years
                state = matrix(c(init.state ,rep(rep(NA, times = length(init.state)), times = max.time-1)),ncol = length(init.state), byrow = T)
                for(which.season.trigger in 1: length(monitoring.season)){                     # monitoring season
                    trap.effort = pre.trigger.effort* eff.prop.if.less[which.eff.mod]
                    for (t in 1:(max.time-1)){
                        state[t+1,] = state[t,]
                        # monitoring informaiton form most recent (including current) specified season
                        # monitoring is before trapping/births/deaths/settlement that season
                        ####
                        if (t %in% season.of.monitoring[, monitoring.season[which.season.trigger]] ){ # monitoring at end of s1 
                            if(t >= year.length +1 & !is.na(season.of.monitoring[t, monitoring.season[which.season.trigger]])&
                               sum(state[season.of.monitoring[t, monitoring.season[which.season.trigger]], 2 ], juv.misID.rate[which.misID] * state[season.of.monitoring[t, monitoring.season[which.season.trigger]], 3 ]) >   
                               sim5.tot.pred.init.triggers.byseason[ which.season.trigger, which.seasonal.strat, which.juvmod , which.eff.mod, which.misID]){
                                trap.effort[season.of.monitoring[t, monitoring.season[which.season.trigger]]:(season.of.monitoring[t, monitoring.season[which.season.trigger]] + year.length-1)] = pre.trigger.effort[season.of.monitoring[t, monitoring.season[which.season.trigger]]:(season.of.monitoring[t, monitoring.season[which.season.trigger]]+ year.length-1)]*eff.prop.if.over[which.eff.mod]
                            }
                        }
                        state[t+1,] = comp.reset.trap.eat()
                        state[t+1,] = comp.J.pred.di.birth()
                        state[t+1,] = comp.A.pred.di.mort()
                        state[t+1,] = comp.J.pred.di.mort()
                        state[t+1,] = comp.better.trapping(juv.trap.chances = juvinile.traping.modifier, current.effort = trap.effort[t], prob.catch1 = prop.caught.w.1.trap)    
                        state[t+1,] = comp.simplest.dd.settlement( ) 
                        state[t+1,] = comp.kill.all.J()                           
                    }
                    sim5.firsteffort.season.spec.trigger[which.season.trigger , which.seasonal.strat, which.juvmod , which.eff.mod , which.misID] =   sum(state[analysable.steps, 12])                                       # add time serise to output object 
                }
            }
        }
    }
}

##########################################################################################################
##               SECOND prep simulation to find density triggeres for reactive trapping                 ##
## triggers itterativly changed by decreasing amount untill total effort spent over the analysis period ##
##            approximates the effort budget                                                            ##
##########################################################################################################
sim5.triggered.aprox.triggers = sim5.tot.pred.init.triggers.byseason # start aproximated fair triggers using those calculated in initial simulaiton
sim5.triggered.effort = sim5.firsteffort.season.spec.trigger*10 # effort spent for each stratedgy over analysible period- started usign effort form intial simulaiton


book.keeping.vec = 1
amount.to.change = 0.6 # largest amount by which the dentiy trigger for reactive trapping is changed each iteration of aproximation
# create a vector of exponentially decreasing size, by which the dentity triggeres for reactive trapping are altered, to progressivly approximate a fair (one that uses as close to as possible the total trapping effrt budget) trigger more accuratley
while (amount.to.change[book.keeping.vec]> 0.00001){
    book.keeping.vec = book.keeping.vec+1
    amount.to.change[book.keeping.vec] = amount.to.change[book.keeping.vec-1]/2.9
}


allowable.budget.discrep = 5 # allowable discrepancy between budget and effort spent over analysible period-- speeds up aproximation method by stopping once efffort spent is satisfactorily close (+- 1%) to budget

for (which.changer in 1:(length(amount.to.change))){
    print (paste("fixing", sum((sim5.triggered.effort < horizon.effort.budget -allowable.budget.discrep | sim5.triggered.effort > horizon.effort.budget + allowable.budget.discrep)*1), "of", length (sim5.triggered.effort) ))  # gives an idea of how many triggers to approximate each time the amount to change triggers is reduced 
    for( which.misID in 1: length(juv.misID.rate)){                                 # juv misID
        for (which.juvmod in 1:length(juv.mods)){                                       # juv trapping
            for (which.seasonal.strat in 1:n.seasonal.eff.combos){                          # seasonal trapping
                for(which.eff.mod in 1: length(eff.prop.if.over)){                              # trigger modifiers              
                    for(which.season.trigger in 1: length(monitoring.season)){                                     # monitoring season
                        juvinile.traping.modifier = juv.mods[which.juvmod]
                        pre.trigger.effort =   as.numeric (rep(yearly.trapping.effort*seasonal.step.effort.mod[which.seasonal.strat,], times = years + burn.in.years))   # select seasonal distribution of effort and replicate for years
                        
                        if(sim5.triggered.effort[ which.season.trigger, which.seasonal.strat , which.juvmod, which.eff.mod, which.misID ] < horizon.effort.budget -allowable.budget.discrep | sim5.triggered.effort[which.season.trigger, which.seasonal.strat , which.juvmod, which.eff.mod, which.misID] > horizon.effort.budget + allowable.budget.discrep){
                            #if the effort is out with the acceptable boundaries :
                            
                            #############################
                            # if too little effort is spent reduce the trigger repetitivly and run untill no longer the case
                            #############################
                            if(sim5.triggered.effort[which.season.trigger, which.seasonal.strat , which.juvmod , which.eff.mod, which.misID] < horizon.effort.budget){
                                while (sim5.triggered.effort[which.season.trigger, which.seasonal.strat , which.juvmod , which.eff.mod, which.misID] < horizon.effort.budget){
                                    sim5.triggered.aprox.triggers[ which.season.trigger, which.seasonal.strat, which.juvmod , which.eff.mod, which.misID]  = sim5.triggered.aprox.triggers[ which.season.trigger, which.seasonal.strat, which.juvmod , which.eff.mod, which.misID] - amount.to.change[which.changer] 
                                    state = matrix(c(init.state ,rep(rep(NA, times = length(init.state)), times = max.time-1)),ncol = length(init.state), byrow = T) 
                                    trap.effort = pre.trigger.effort*eff.prop.if.less[which.eff.mod]
                                    for (t in 1:(max.time-1)){
                                        state[t+1,] = state[t,]
                                        # select trapping effort based on previous density.
                                        ####
                                        if (t %in% season.of.monitoring[, monitoring.season[which.season.trigger]] ){
                                            if(t >= year.length +1 & !is.na(season.of.monitoring[t, monitoring.season[which.season.trigger]])&
                                               sum(state[season.of.monitoring[t, monitoring.season[which.season.trigger]], 2 ], juv.misID.rate[which.misID] * state[season.of.monitoring[t, monitoring.season[which.season.trigger]], 3 ]) >   
                                               sim5.triggered.aprox.triggers[ which.season.trigger, which.seasonal.strat, which.juvmod , which.eff.mod, which.misID]){
                                                trap.effort[season.of.monitoring[t, monitoring.season[which.season.trigger]]:(season.of.monitoring[t, monitoring.season[which.season.trigger]]+ year.length-1)] = pre.trigger.effort[season.of.monitoring[t, monitoring.season[which.season.trigger]]:(season.of.monitoring[t, monitoring.season[which.season.trigger]]+ year.length-1)]*eff.prop.if.over[which.eff.mod]
                                            }
                                        }             
                                        state[t+1,] = comp.reset.trap.eat()
                                        state[t+1,] = comp.J.pred.di.birth()
                                        state[t+1,] = comp.A.pred.di.mort()
                                        state[t+1,] = comp.J.pred.di.mort()
                                        state[t+1,] = comp.better.trapping(juv.trap.chances = juvinile.traping.modifier, current.effort = trap.effort[t], prob.catch1 = prop.caught.w.1.trap)    
                                        state[t+1,] = comp.simplest.dd.settlement() 
                                        state[t+1,] = comp.kill.all.J()                            # prey consumption
                                    }
                                    sim5.triggered.effort[ which.season.trigger, which.seasonal.strat , which.juvmod , which.eff.mod, which.misID] = sum(state[analysable.steps,12])
                                }
                            }
                            #############################
                            # if too MUCH effort is spent INCREASE the trigger repetitivly and run untill no longer the case
                            #############################
                            if(sim5.triggered.effort[which.season.trigger, which.seasonal.strat , which.juvmod , which.eff.mod, which.misID] > horizon.effort.budget){
                                while (sim5.triggered.effort[ which.season.trigger, which.seasonal.strat , which.juvmod , which.eff.mod, which.misID] > horizon.effort.budget){
                                    sim5.triggered.aprox.triggers[ which.season.trigger, which.seasonal.strat, which.juvmod , which.eff.mod, which.misID]  = sim5.triggered.aprox.triggers[ which.season.trigger, which.seasonal.strat, which.juvmod , which.eff.mod, which.misID] + amount.to.change[which.changer] 
                                    state = matrix(c(init.state ,rep(rep(NA, times = length(init.state)), times = max.time-1)),ncol = length(init.state), byrow = T) 
                                    trap.effort = pre.trigger.effort*eff.prop.if.less[which.eff.mod]
                                    for (t in 1:(max.time-1)){
                                        state[t+1,] = state[t,]
                                        # select trapping effort based on previous density.
                                        #### 
                                        if (t %in% season.of.monitoring[, monitoring.season[which.season.trigger]] ){
                                            if(t >= year.length +1 & !is.na(season.of.monitoring[t, monitoring.season[which.season.trigger]])&
                                               sum(state[season.of.monitoring[t, monitoring.season[which.season.trigger]], 2 ], juv.misID.rate[which.misID] * state[season.of.monitoring[t, monitoring.season[which.season.trigger]], 3 ]) >   
                                               sim5.triggered.aprox.triggers[ which.season.trigger, which.seasonal.strat, which.juvmod , which.eff.mod, which.misID]){
                                                trap.effort[season.of.monitoring[t, monitoring.season[which.season.trigger]]:(season.of.monitoring[t, monitoring.season[which.season.trigger]]+ year.length-1)] = pre.trigger.effort[season.of.monitoring[t, monitoring.season[which.season.trigger]]:(season.of.monitoring[t, monitoring.season[which.season.trigger]]+ year.length-1)]*eff.prop.if.over[which.eff.mod]
                                            }
                                        }
                                        
                                        state[t+1,] = comp.reset.trap.eat()
                                        state[t+1,] = comp.J.pred.di.birth(  )
                                        state[t+1,] = comp.A.pred.di.mort()
                                        state[t+1,] = comp.J.pred.di.mort()
                                        state[t+1,] = comp.better.trapping(juv.trap.chances = juvinile.traping.modifier, current.effort = trap.effort[t], prob.catch1 = prop.caught.w.1.trap)    
                                        state[t+1,] = comp.simplest.dd.settlement( ) 
                                        state[t+1,] = comp.kill.all.J()                            # prey consumption
                                    }
                                    sim5.triggered.effort[ which.season.trigger , which.seasonal.strat , which.juvmod, which.eff.mod, which.misID] = sum(state[analysable.steps,12])
                                } # while too much effort
                            } # if too much effort
                        } # if the effort is outside +- allowable.budget.discrep of allowed effort
                        
                    } # which season in the previosu year triggers
                } # which.eff.mod
            } # which.seasonal.strat
        } # which.juvmod
    }
} # which.changer


##########################################################################################################    
##          FIRST PROPER SIMULAITON -- CONSTANT ANNUAL EFFORT                                           ##
##########################################################################################################

# first simulation of constatn anual effort stratedgies
ex4.sim.notrig = array(NA, dim = c(max.time, length(init.state), length(monitoring.season) ,n.seasonal.eff.combos, length(juv.mods), 1 , length(juv.misID.rate) )) #  dime names = (time, state bits , seasonal effort combos, juv.trapping.mod)

for( which.misID in 1: length(juv.misID.rate)){                                 # juv misID
    for (which.seasonal.strat in 1:n.seasonal.eff.combos){
        for (which.juvmod in 1:length(juv.mods)){
            juvinile.traping.modifier = juv.mods[which.juvmod]
            state = matrix(c(init.state ,rep(rep(NA, times = length(init.state)), times = max.time-1)),ncol = length(init.state), byrow = T) 
            trap.effort =   as.numeric (rep(yearly.trapping.effort*seasonal.step.effort.mod[which.seasonal.strat,], times = years + burn.in.years))   # select seasonal distribution of effort and replicate for years
            
            for (t in 1:(max.time-1)){
                state[t+1,] = state[t,]
                state[t+1,] = comp.reset.trap.eat()
                state[t+1,] = comp.J.pred.di.birth(  )
                state[t+1,] = comp.A.pred.di.mort()
                state[t+1,] = comp.J.pred.di.mort()
                state[t+1,] = comp.better.trapping(juv.trap.chances = juvinile.traping.modifier, current.effort = trap.effort[t], prob.catch1 = prop.caught.w.1.trap)    
                state[t+1,] = comp.simplest.dd.settlement( ) 
                state[t+1,] = comp.kill.all.J()                            # prey consumption            
            }#  for time
            ex4.sim.notrig [ , ,which.season.trigger,  which.seasonal.strat ,which.juvmod , 1 , which.misID] = state
        }
    }
}

##########################################################################################################    
##          SECOND PROPER SIMULAITON -- REACTIVE ANNUAL EFFORT                                          ##
##########################################################################################################

ex4.sim.state =             array(NA, dim = c(max.time, length(init.state), length(monitoring.season) ,n.seasonal.eff.combos, length(juv.mods), length(eff.prop.if.over), length(juv.misID.rate) )) #  dime names = (time, state bits , seasonal effort combos, juv.trapping.mod)


for( which.misID in 1: length(juv.misID.rate)){                                 # juv misID
    for (which.seasonal.strat in 1:n.seasonal.eff.combos){
        for (which.juvmod in 1:length(juv.mods)){
            for(which.eff.mod in 1: length(eff.prop.if.over)){
                for(which.season.trigger in 1: length(monitoring.season)){                                     # monitoring season
                    juvinile.traping.modifier = juv.mods[which.juvmod]
                    state = matrix(c(init.state ,rep(rep(NA, times = length(init.state)), times = max.time-1)),ncol = length(init.state), byrow = T) 
                    pre.trigger.effort =   as.numeric (rep(yearly.trapping.effort*seasonal.step.effort.mod[which.seasonal.strat,], times = years + burn.in.years))   # select seasonal distribution of effort and replicate for years
                    trap.effort = pre.trigger.effort*eff.prop.if.less[which.eff.mod]
                    for (t in 1:(max.time-1)){
                        state[t+1,] = state[t,]
                        ######################################################
                        # MONITORING
                        ######################################################
                        if (t %in% season.of.monitoring[, monitoring.season[which.season.trigger]] ){
                            if(t >= year.length +1 & !is.na(season.of.monitoring[t, monitoring.season[which.season.trigger]]) & # sum(trap.effort[analysable.steps]) < horizon.effort.budget & # this last bit maybe only needed for imperfect monitoring
                               sum(state[season.of.monitoring[t, monitoring.season[which.season.trigger]], 2 ], juv.misID.rate[which.misID] * state[season.of.monitoring[t, monitoring.season[which.season.trigger]], 3 ]) >   
                               sim5.triggered.aprox.triggers[ which.season.trigger, which.seasonal.strat, which.juvmod , which.eff.mod, which.misID]){
                                trap.effort[season.of.monitoring[t, monitoring.season[which.season.trigger]]:(season.of.monitoring[t, monitoring.season[which.season.trigger]]+ year.length-1)] = pre.trigger.effort[season.of.monitoring[t, monitoring.season[which.season.trigger]]:(season.of.monitoring[t, monitoring.season[which.season.trigger]]+ year.length-1)]*eff.prop.if.over[which.eff.mod]
                            }
                        }
                        #######################################################
                        
                        state[t+1,] = comp.reset.trap.eat()
                        state[t+1,] = comp.J.pred.di.birth(  )
                        state[t+1,] = comp.A.pred.di.mort()
                        state[t+1,] = comp.J.pred.di.mort()
                        state[t+1,] = comp.better.trapping(juv.trap.chances = juvinile.traping.modifier, current.effort = trap.effort[t], prob.catch1 = prop.caught.w.1.trap)    
                        state[t+1,] = comp.simplest.dd.settlement( ) 
                        state[t+1,] = comp.kill.all.J()                            # prey consumption            
                    }#  for time
                    ex4.sim.state [ , ,which.season.trigger,  which.seasonal.strat , which.juvmod , which.eff.mod , which.misID] = state
                }
            }
        }}}




sims.for.analysis = abind(x = ex4.sim.notrig, y = ex4.sim.state, along = 6) # bind together array form unreactive and reactive stratedgies



#########################################################################################################
