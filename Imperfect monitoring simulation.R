# script file name "Imperfect monitoring simulation.R"
# this code is sourced in the script "stage structured invasive management.R"
# save this file as "Imperfect monitoring simulation.R" and run "stage structured invasive management.R"


## #############################################################################################################################
# CALCULATE PROPORTION OF DETECTIONS EXPECTED FOR EACH TRIGGER  set up as array
################################################################################################################################
prob.detect.1indiv = prop.caught.w.1.trap
adult.prop.detections.trigger = 1-((1-prob.detect.1indiv)^sim.kd.tiggers) 

## #############################################################################################################################
# SIMULATION SCENERIO 6: TRIGGERING OF EFFORT BY MONITORING
###############################################################################################################################
interesting.timeframes = c("T", "S1", "S2", "S3", "S4")
interesting.quantities = c("Mean", "Q1", "Median", "Q3")
interesting.timequants = paste (rep(interesting.timeframes, times = length (interesting.quantities)), rep (interesting.quantities, each = length (interesting.timeframes)))

# quick stuff im interested in
total.adult.harvest = array (NA, dim = c( length(monitoring.season), n.seasonal.eff.combos ,length(juv.mods), length(eff.prop.if.over),length(juv.misID.rate) ,length(interesting.nrafts) , n.iterations))
total.harvest =       array (NA, dim = dim(total.adult.harvest))
horizon.effort =                    array (NA, dim = dim(total.adult.harvest))

maxs2.Apop =           array (NA, dim = dim(total.adult.harvest))
medians2.Apop =           array (NA, dim = dim(total.adult.harvest))
maxs3.Apop =           array (NA, dim = dim(total.adult.harvest))
medians3.Apop =           array (NA, dim = dim(total.adult.harvest))
percentil90s3.Apop =           array (NA, dim = dim(total.adult.harvest))
#########################################################################################################################
# prep to save mean and sd for each timestep
#
adult.trace.means = array(NA, dim = c( max.time, length(monitoring.season), n.seasonal.eff.combos ,length(juv.mods), length(eff.prop.if.over),length(juv.misID.rate) ,length(interesting.nrafts)))
adult.trace.sds = array(NA, dim = c( max.time, length(monitoring.season), n.seasonal.eff.combos ,length(juv.mods), length(eff.prop.if.over),length(juv.misID.rate) ,length(interesting.nrafts)))
#
#########################################################################################################################


for(which.rafts in 1:length(interesting.nrafts)){
    for (which.juvmod in 1:length(which.interesting.juv.mods)){
        for (which.seasonal.strat in which.interesting.ss){
            for(which.eff.mod in 1:length(which.interesting.eff.overs)){
                for(which.season.trigger in 1: length(monitoring.season)){
                    for( which.misID in 1: length(juv.misID.rate)){
                        juvinile.traping.modifier = juv.mods[which.interesting.juv.mods[which.juvmod]]
                        state.of.iterations = array(NA, dim = c(dim(state), n.iterations))
                        pre.trigger.effort =   as.numeric (rep(yearly.trapping.effort*seasonal.step.effort.mod[which.seasonal.strat,], times = years))   # select seasonal distribution of effort and replicate for years
                        
                        traces.by.iteration = matrix(NA, nrow = max.time, ncol = n.iterations)
                        
                        for(iteration in 1:n.iterations){
                            state = matrix(c(init.state ,rep(rep(NA, times = length(init.state)), times = max.time-1)),ncol = length(init.state), byrow = T) 
                            trap.effort = pre.trigger.effort*eff.prop.if.less[which.interesting.eff.overs[which.eff.mod]]
                            for (t in 1:(max.time-1)){
                                state[t+1,] = state[t,]
                                # select trapping effort based on monitoring
                                ####
                                if ((t %in% season.of.monitoring[, which.season.trigger]) &
                                    (t >= year.length +1 & !is.na(season.of.monitoring[t, which.season.trigger])& sum(state[analysable.steps,12], na.rm = T) < horizon.effort.budget &
                                     comp.raft.monitor( pops = state[season.of.monitoring[t,monitoring.season[which.season.trigger]],], n.rafts = interesting.nrafts[which.rafts], prob.1adult = prob.detect.1indiv ) >   
                                     adult.prop.detections.trigger[ which.season.trigger, which.seasonal.strat, which.interesting.juv.mods[which.juvmod] , which.interesting.eff.overs[which.eff.mod], which.misID])){
                                    trap.effort[season.of.monitoring[t, which.season.trigger]:(season.of.monitoring[t, which.season.trigger]+ year.length-1)] = 
                                        pre.trigger.effort[season.of.monitoring[t, which.season.trigger]:(season.of.monitoring[t, which.season.trigger]+ year.length-1)]*eff.prop.if.over[which.interesting.eff.overs[which.eff.mod]]
                                }       
                                ####
                                state[t+1,] = comp.reset.trap.eat()
                                #state[t+1,] = comp.prey.bazy()
                                state[t+1,] = comp.J.pred.di.birth(  )
                                state[t+1,] = comp.A.pred.di.mort()
                                state[t+1,] = comp.J.pred.di.mort()
                                state[t+1,] = comp.better.trapping(juv.trap.chances = juvinile.traping.modifier, current.effort = trap.effort[t], prob.catch1 = 0.1)    
                                state[t+1,] = comp.simplest.dd.settlement( ) 
                                state[t+1,] = comp.kill.all.J()                            # prey consumption
                            } #  for time
                            total.adult.harvest [which.season.trigger,which.seasonal.strat, which.interesting.juv.mods[which.juvmod] , which.interesting.eff.overs[which.eff.mod], which.misID ,  which.rafts, iteration] = sum(state[,5])
                            total.harvest [which.season.trigger,which.seasonal.strat, which.interesting.juv.mods[which.juvmod] , which.interesting.eff.overs[which.eff.mod], which.misID ,  which.rafts, iteration] = sum(c(state[,5],state[,6] ))
                            horizon.effort[which.season.trigger,which.seasonal.strat, which.interesting.juv.mods[which.juvmod] , which.interesting.eff.overs[which.eff.mod], which.misID ,  which.rafts, iteration] = sum(state[,12])
                            
                            maxs2.Apop [which.season.trigger,which.seasonal.strat, which.interesting.juv.mods[which.juvmod] , which.interesting.eff.overs[which.eff.mod], which.misID ,  which.rafts, iteration] = max(state[analysible.s2,2])
                            medians2.Apop [which.season.trigger,which.seasonal.strat, which.interesting.juv.mods[which.juvmod] , which.interesting.eff.overs[which.eff.mod], which.misID ,  which.rafts, iteration] = median(state[analysible.s2,2])  
                            nu.percentil90s2.Apop [which.seasonal.strat, which.interesting.juv.mods[which.juvmod] , which.interesting.eff.overs[which.eff.mod], which.rafts, iteration] = quantile(state[analysible.s2,2], probs = 0.9 )

                            traces.by.iteration[,iteration] = state[,2]
                            
                        }
                        
                        adult.trace.means [, which.season.trigger , which.seasonal.strat , which.interesting.juv.mods[which.juvmod]  , which.interesting.eff.overs[which.eff.mod] , which.misID , which.rafts ] = apply(traces.by.iteration , 1, mean )
                        adult.trace.sds [, which.season.trigger , which.seasonal.strat , which.interesting.juv.mods[which.juvmod]  , which.interesting.eff.overs[which.eff.mod] , which.misID , which.rafts ] = apply(traces.by.iteration , 1, sd )
                    }
                }   
            } # which.eff.mod
        } # which.seasonal.strat
    } # which.juvmod
    print(which.rafts/length(interesting.nrafts))
}

