# function for initial treatment randomization
# created 4/8/2025
# updated 12/17/2025

# n_treatments is the # treatments in the initial randomization
# p_treatment is a vector of randomization probabilities associated with each treatment
#### must add to 1
# p_respond is a vector where each element is the probability of responding to initial treatment
##### first element is probability of response to treatment 1, second is prob of response to treatment 2, etc
##### does NOT have to add to 1!

stage1_assign = function(n_treatments_init=2, p_treatment=c(0.5,0.5), p_respond=c(0.5,0.5)) {
  # let A1 be treatment 1, A2 be treatment 2, etc.
  treatment_1 = sample(1:n_treatments_init, size=1, prob = p_treatment)

  # generate true response status
  # responder = 1, nonresponder = 0
  true_response_1 = rbinom(n=1, size=1, prob=p_respond[treatment_1])
  return(list(treatment_1 = treatment_1, true_response_1 = true_response_1))
}


# outputs tailoring variable results (1 = responder, 0 = nonresponder)
# prob_obs_responder
tv_response = function(true_response, sensitivity, specificity, n_tv=3, bp=FALSE) {
  prob_obs_responder = ifelse(true_response==1, sensitivity, 1-specificity)
  if(bp==FALSE) {
    tv_output = rbinom(n=n_tv, size = 1, prob = prob_obs_responder)
  } else {
    tv_output = rep(NA,n_tv)
    tv_output[1] = rbinom(n=1, size=1, prob=prob_obs_responder)
    counter = 1
    while(counter<n_tv & tv_output[counter]==0) {
      counter = counter + 1
      tv_output[counter] = rbinom(n=1, size=1, prob=prob_obs_responder)
    }
  }
  return(tv_output)
}


# Q function for EM alg
Q <- function(params, params_old, dat, n_tv=3) {
  
  sens = params[1]
  spec = params[2]
  prev_arm1 = params[3]
  prev_arm2 = params[4]
  
  sens_old = params_old[1]
  spec_old = params_old[2]
  prev_arm1_old = params_old[3]
  prev_arm2_old = params_old[4]
  
  #prev_old_vec = numeric(nrow(dat)) 
  #for (i in 1:nrow(dat)) {
  #  prev_old_vec[i] = ifelse(dat$treatment_1[i]==1,prev_arm1_old,prev_arm2_old)
  #}
  
  prev_old_vec = data.table::fifelse(dat$treatment_1==1,prev_arm1_old,prev_arm2_old)
  
  y = dat$n_tv_pos
  
  # compute P(Si=1 | data, params_old)
  numerator = prev_old_vec * sens_old^y * (1-sens_old)^(n_tv - y)
  denominator = numerator + (1-prev_old_vec) * (1-spec_old)^y * spec_old^(n_tv - y)
  prob_condit_truepos = numerator / denominator
  
  term1 = log(prev_arm1) * sum(prob_condit_truepos[which(dat$treatment_1==1)])
  
  term2 = log(1-prev_arm1) * sum((1 - prob_condit_truepos[which(dat$treatment_1==1)]))
  
  term3 = log(prev_arm2) * sum(prob_condit_truepos[which(dat$treatment_1==2)])
  
  term4 = log(1-prev_arm2) * sum((1 - prob_condit_truepos[which(dat$treatment_1==2)]))
  
  term5 = log(sens) * sum(prob_condit_truepos * y)
  
  term6 = log(1-sens) * sum(prob_condit_truepos * (n_tv - y))
  
  term7 = log(1-spec) * sum((1-prob_condit_truepos) * y)
  
  term8 = log(spec) * sum((1-prob_condit_truepos) * (n_tv - y))
  
  
  return(term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8)
}

# function that runs EM algorithm, takes in n, parameter values, initial guesses
# params_true = c(sensitivity, specificity, response_prev_arm1, response_prev_arm2)

EM_alg = function(n, params_true, params_init, method=c("optim","manual"), n_tv=3, bp=FALSE) {
  # generate data
  sens = params_true[1]
  spec = params_true[2]
  prev_arm1 = params_true[3]
  prev_arm2 = params_true[4]
  # first stage assignment and true response status
  stage1_results_list = replicate(n,
                                  stage1_assign(
                                    p_respond = c(prev_arm1, prev_arm2)
                                  ),
                                  simplify = FALSE
  )
  
  # Convert list of lists to a data frame
  stage1_results <- do.call(rbind, lapply(stage1_results_list, as.data.frame))

  # now generate tailoring variable values
  # For each patient, generate observed tailoring variables based on sensitivity/specificity and response to first stage treatment
  # Use lapply to generate tailoring variables row-by-row
  tv_list <- lapply(1:nrow(stage1_results), function(i) {
    true_response <- stage1_results$true_response_1[i]
    tv_response(
      true_response = true_response,
      sensitivity = sens,
      specificity = spec,
      n_tv = n_tv,
      bp=bp
    )
  })
  
  # Convert list of vectors to a matrix/data frame
  tv_matrix <- do.call(rbind, tv_list)
  colnames(tv_matrix) <- paste0("tv", 1:ncol(tv_matrix))
  
  # Combine with the original data
  results_with_tv <- cbind(stage1_results, tv_matrix)
  
  # add column for the number of positive tests
  tv_cols = paste0("tv", 1:ncol(tv_matrix))
  results_with_tv$n_tv_pos = rowSums(results_with_tv[,tv_cols], na.rm = TRUE)
  
  
  params <- c("sensitivity" = params_init[1],
              "specificity" = params_init[2],
              "prevalence_arm1" = params_init[3],
              "prevalence_arm2" = params_init[4])
  
  # Intiate loop variables
  tol <- 1e-10
  last_iter <- 5
  diff <- 10
  
  if (method=="optim") {
    while (diff >= tol) {
      params_old <- params
    
      o <- optim(par = params_old, fn = Q, #gr = Qgrad,
                 params_old = params_old, dat = results_with_tv,
                 lower = 0.000001,
                 upper = 0.999999,
                 method = "L-BFGS-B",
                 control = list(fnscale = -1))
    
      params <- o$par
    
      #diff <- abs(o$value - last_iter)
    
      diff <- max(abs(params - params_old))
    
      last_iter <- o$value
    }
  } else {
    if(bp==FALSE) {
      while (diff >= tol) {
        params_old <- params
        sens_old = params_old[1]
        spec_old = params_old[2]
        prev_arm1_old = params_old[3]
        prev_arm2_old = params_old[4]
      
        arm1_ind = which(results_with_tv$treatment_1==1)
        arm2_ind = which(results_with_tv$treatment_1==2)
        n_arm1 = length(arm1_ind)
        n_arm2 = length(arm2_ind)
      
        y = results_with_tv$n_tv_pos
      
        # compute P(Si=1 | data, current est)
        prev_old = ifelse(results_with_tv$treatment_1==1, prev_arm1_old, prev_arm2_old)
        numerator = prev_old * sens_old^y * (1-sens_old)^(n_tv - y)
        denominator = numerator + (1-prev_old) * (1-spec_old)^y * spec_old^(n_tv - y)
        prob_responder_condit_cur = numerator / denominator
      
        prev_arm1_cur = sum(prob_responder_condit_cur[arm1_ind]) / n_arm1
        prev_arm2_cur = sum(prob_responder_condit_cur[arm2_ind]) / n_arm2
        sens_cur = sum(prob_responder_condit_cur * y) / (n_tv * (n_arm1*prev_arm1_cur + n_arm2*prev_arm2_cur))
        spec_cur = 1 - ((sum(y) - sum(prob_responder_condit_cur * y)) / (n_tv * (n - n_arm1*prev_arm1_cur - n_arm2*prev_arm2_cur)))
      
        params <- c("sensitivity" = sens_cur,
                    "specificity" = spec_cur,
                    "prevalence_arm1" = prev_arm1_cur,
                    "prevalence_arm2" = prev_arm2_cur)
      
        diff <- max(abs(params - params_old))
      } 
    } else {
      while (diff >= tol) {
        params_old <- params
        sens_old = params_old[1]
        spec_old = params_old[2]
        prev_arm1_old = params_old[3]
        prev_arm2_old = params_old[4]
        
        arm1_ind = which(results_with_tv$treatment_1==1)
        arm2_ind = which(results_with_tv$treatment_1==2)
        n_arm1 = length(arm1_ind)
        n_arm2 = length(arm2_ind)
        
        y = results_with_tv$n_tv_pos
        R_i = rowSums(!is.na(results_with_tv[,tv_cols]))
        
        # compute P(Si=1 | data, current est)
        prev_old = ifelse(results_with_tv$treatment_1==1, prev_arm1_old, prev_arm2_old)
        numerator = prev_old * sens_old^y * (1-sens_old)^(R_i - y)
        denominator = numerator + (1-prev_old) * (1-spec_old)^y * spec_old^(R_i - y)
        prob_responder_condit_cur = numerator / denominator
        
        prev_arm1_cur = sum(prob_responder_condit_cur[arm1_ind]) / n_arm1
        prev_arm2_cur = sum(prob_responder_condit_cur[arm2_ind]) / n_arm2
        sens_cur = sum(prob_responder_condit_cur * y) / sum(prob_responder_condit_cur * R_i)
        spec_cur = 1 - ((sum(y) - sum(prob_responder_condit_cur * y)) / (sum(R_i) - sum(prob_responder_condit_cur * R_i)))
        
        params <- c("sensitivity" = sens_cur,
                    "specificity" = spec_cur,
                    "prevalence_arm1" = prev_arm1_cur,
                    "prevalence_arm2" = prev_arm2_cur)
        
        diff <- max(abs(params - params_old))
      }
    }
  }
  
  return(params)
}

### create a function that inputs TV assessments for a single person and 
### updated EM alg parameter estimates, then outputs person's prob of response
### the derivation for this is in the Goodnotes file Prob Of Response Derivation
prob_responder = function(tv_vec, prev, sens, spec) {
  
  n_pos = sum(tv_vec==1)
  n_neg = sum(tv_vec==0)
  s1_prod = sens^n_pos * (1-sens)^n_neg
  s0_prod = (1-spec)^n_pos * spec^n_neg

  prob_responder = (prev * s1_prod) / (prev*s1_prod + (1-prev)*s0_prod)
  return(prob_responder)
}


### assign stage 2 treatment based on initial treatment and prob of being a responder
stage2_assign = function(prob_responder) {
  # let A1 be treatment 1, A2 be treatment 2, etc.
  observed_response_1 = rbinom(n=1, size=1, prob=prob_responder)
  if (observed_response_1==1) {
    treatment_2 = 0
  } else {
    treatment_2 = rbinom(n=1, size=1,prob=.5) * 2 - 1
  }

  return(list(observed_response_1 = observed_response_1, treatment_2 = treatment_2))
}





