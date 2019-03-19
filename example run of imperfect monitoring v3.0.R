# script file name "example run of imperfect monitoring v3.0.R"
# this code is sourced in the script "stage structured invasive management.R"
# save this file as "example run of imperfect monitoring v3.0.R" and run "stage structured invasive management.R"

## #############################################################################################################################
# CALCULATE PROPORTION OF DETECTIONS EXPECTED FOR EACH TRIGGER  set up as array
################################################################################################################################
prob.detect.1indiv = prop.caught.w.1.trap
adult.prop.detections.trigger = -1*(exp(examp.trigs*(log(1-prob.detect.1indiv)))-1)

## #############################################################################################################################
# SIMULATION SCENERIO 6: TRIGGERING OF EFFORT BY MONITORING
###############################################################################################################################
n.iterations = 1

e.g.trace.save = array(NA, dim = c(max.time, n.seasonal.eff.combos, length(examp.juv.mod), length(eff.prop.if.over), length(interesting.nrafts)))
#########################################################################################################################
for( which.misID in 1: length(juv.misID.rate)){                                 # juv misID
    for(which.rafts in 1:length(interesting.nrafts)){
    for (which.juvmod in 1:length(examp.juv.mod)){
        for (which.seasonal.strat in 1:length(which.ss)){
            for(which.eff.mod in which.interestin.eff.overs){
                juvinile.traping.modifier = examp.juv.mod[which.juvmod]
                state.of.iterations = array(NA, dim = c(dim(state), n.iterations))
                pre.trigger.effort =   as.numeric (rep(yearly.trapping.effort*seasonal.step.effort.mod[which.ss[which.seasonal.strat],], times = years + burn.in.years))   # select seasonal distribution of effort and replicate for years
                
                for(iteration in 1:n.iterations){
                    state = matrix(c(init.state ,rep(rep(NA, times = length(init.state)), times = max.time-1)),ncol = length(init.state), byrow = T) 
                    trap.effort = pre.trigger.effort*eff.prop.if.less[which.eff.mod]
                    for (t in 1:(max.time-1)){
                        state[t+1,] = state[t,]
                        # select trapping effort based on monitoring
                        ####
                        if ((t %in% season.of.monitoring[, monitoring.season[which.season.trigger]] ) & 
                            (t >= year.length +1 & !is.na(season.of.monitoring[t, monitoring.season[which.season.trigger]]) & 
                             sum(state[analysable.steps,12], na.rm = T) < horizon.effort.budget &
                             comp.raft.monitor( pops = state[season.of.monitoring[t,monitoring.season[which.season.trigger]],], n.rafts = interesting.nrafts[which.rafts], prob.1adult = prob.detect.1indiv ) >   
                               adult.prop.detections.trigger[ which.season.trigger, which.ss[which.seasonal.strat], examp.juv.mod[which.juvmod] , which.eff.mod, which.misID])){
                                
                                trap.effort[season.of.monitoring[t, monitoring.season[which.season.trigger]]:(season.of.monitoring[t, monitoring.season[which.season.trigger]]+ year.length-1)] = pre.trigger.effort[season.of.monitoring[t, monitoring.season[which.season.trigger]]:(season.of.monitoring[t, monitoring.season[which.season.trigger]]+ year.length-1)]*eff.prop.if.over[which.eff.mod]
                        
                        }         
                        ####
                        state[t+1,] = comp.reset.trap.eat()
                        #state[t+1,] = comp.prey.bazy()
                        state[t+1,] = comp.J.pred.di.birth(  )
                        state[t+1,] = comp.A.pred.di.mort()
                        state[t+1,] = comp.J.pred.di.mort()
                        state[t+1,] = comp.better.trapping(juv.trap.chances = juvinile.traping.modifier, current.effort = trap.effort[t], prob.catch1 = prop.caught.w.1.trap)    
                        state[t+1,] = comp.simplest.dd.settlement( ) 
                        state[t+1,] = comp.kill.all.J()                            # prey consumption
                    } #  for time
                    e.g.trace.save[, which.ss[which.seasonal.strat] , which.juvmod, which.eff.mod , which.rafts] = state[,2]
                }   
            } # which.eff.mod
        } # which.seasonal.strat
        print(which.juvmod)
    } # which.juvmod
    }
}

