
#source("~/Box Sync/MelissaMartin_Dissertation/Paper3/code/helper_functions.R")
setwd("/home/martin30/MelissaMartin_Dissertation/Paper3")
source("code/helper_functions.R")

n_sim = 1000
n = 1000 #total number of patients in trial
n_tv = 3
start_EM = 201

prevalence = data.frame(prevalence_arm1 = c(0.7, 0.4), prevalence_arm2 = c(0.6, .3))
sensitivity = c(.8, .5)
specificity = c(.8, .5)
EM_alloc = c(TRUE, FALSE)
setting = merge(
  prevalence,
  expand.grid(sensitivity = sensitivity, specificity = specificity, EM_alloc = EM_alloc)
)

index = as.numeric(commandArgs(trailingOnly = TRUE)[1])
set.seed(index)


#### true parameter values
sensitivity = setting$sensitivity[index] # tv sensitivity
specificity = setting$specificity[index] # tv specificity
response_prev_arm1 = setting$prevalence_arm1[index] # arm 1 prevalence, i.e. prob of being a true responder conditional on arm 1
response_prev_arm2 = setting$prevalence_arm2[index] # arm 2 prevalence, i.e. prob of being a true responder conditional on arm 1
EM_alloc = setting$EM_alloc[index]
#n_tv = setting[index,2] # number of tailoring variable assessments


#### means
outcome_means = list(a1_r_cor = 10,
                     a1_r_mis = 7,
                     a1_nr_cor_a2nr_1 = 9,
                     a1_nr_mis_a2nr_1 = 8,
                     a1_nr_cor_a2nr_2 = 8,
                     a1_nr_mis_a2nr_2 = 7,
                     b1_r_cor = 8,
                     b1_r_mis = 6,
                     b1_nr_cor_b2nr_1 = 8,
                     b1_nr_mis_b2nr_1 = 7,
                     b1_nr_cor_b2nr_2 = 7,
                     b1_nr_mis_b2nr_2 = 6)

sd = 1


for (sim in 1:n_sim) {
  
  # first stage assignment and true response status
  stage1_results_list = replicate(n,
                                  stage1_assign(
                                    p_respond = c(response_prev_arm1, response_prev_arm2)
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
      sensitivity = sensitivity,
      specificity = specificity,
      n_tv = n_tv
    )
  })
  
  # Convert list of vectors to a matrix/data frame
  tv_matrix <- do.call(rbind, tv_list)
  colnames(tv_matrix) <- paste0("tv", 1:ncol(tv_matrix))
  
  # Combine with the original data
  results_with_tv <- cbind(stage1_results, tv_matrix)
  
  # add column for the number of positive tests
  tv_cols = paste0("tv", 1:ncol(tv_matrix))
  results_with_tv$n_tv_pos = rowSums(results_with_tv[,tv_cols])
  
  dat = results_with_tv
  dat$p_i = NA
  dat$observed_response_1 = NA
  dat$treatment_2 = NA
  dat$outcome = NA
  
  
  # go one patient at a time
  param_mat = matrix(nrow=n, ncol=4)
  # create for loop that adds each patient, runs EM algorithm with current patient + previous patients, and generates outcome
  for (i in 1:n) {
    #print(i)
    dat_cur = dat[1:i,]
  
    # hard code initial guess for parameter estimates
    #if (i==1) {
    params = c("sensitivity" = .7,
              "specificity" = .75,
              "prevalence_arm1" = .65,
              "prevalence_arm2" = .65)
    #}
  
    # only run EM if we've accumulated enough data
    if (EM_alloc & i >= start_EM) {
      # Intiate loop variables
      tol <- 1e-5
      #last_iter <- 5
      diff <- 10
    
      iter = 0
      while (diff >= tol) {
        iter=iter+1
        #print(iter)
        params_old <- params
      
        o <- optim(par = params_old, fn = Q, #gr = Qgrad,
                  params_old = params_old, dat = dat_cur, n_tv=n_tv,
                  lower = 0.000001,
                  upper = 0.999999,
                  method = "L-BFGS-B",
                  control = list(fnscale = -1))
      
        params <- o$par
        
        #diff <- abs(o$value - last_iter)
      
        diff <- max(abs(params - params_old))
      
        last_iter <- o$value
        
        #print(params)
        #print(diff)
      }
      
      param_mat[i,] = params
      
      # use EM alg estimates to compute p_i
      prev_i = ifelse(dat[i,1]==1, params[3], params[4])
      p_i = prob_responder(tv_vec = as.numeric(dat[i,tv_cols]), prev = prev_i, sens = params[1], spec=params[2])
      
    } else {
      p_i = sum(as.numeric(dat[i,tv_cols])) / length(tv_cols)
    }

    # truncate if p_i > 0.95 or p_i < 0.05
    p_i = min(c(p_i, 0.95))
    p_i = max(c(p_i, 0.05))
  
    dat$p_i[i] = p_i
  
    stage2 = stage2_assign(p_i)
    dat$observed_response_1[i] = stage2$observed_response_1
    dat$treatment_2[i] = stage2$treatment_2
  
    # generate outcomes
    if (dat$treatment_1[i]==1) {
      if (dat$true_response_1[i]==1 & dat$treatment_2[i]==0) { # correctly classified as a responder
        y = rnorm(1, mean=outcome_means$a1_r_cor, sd=sd)
      } else if (dat$true_response_1[i]==0 & dat$treatment_2[i]==0) { # misclassified as a responder
        y = rnorm(1, mean=outcome_means$a1_r_mis, sd=sd)
      } else if (dat$true_response_1[i]==0 & dat$treatment_2[i]==1) { # correctly classified as nonresponder and receives a2nr = 1
        y = rnorm(1, mean=outcome_means$a1_nr_cor_a2nr_1, sd=sd)
      } else if (dat$true_response_1[i]==1 & dat$treatment_2[i]==1) { # misclassified as nonresponder and received a2nr = 1
        y = rnorm(1, mean=outcome_means$a1_nr_mis_a2nr_1, sd=sd)
      } else if (dat$true_response_1[i]==0 & dat$treatment_2[i]==-1) { # correctly classified as nonresponder and received a2nr = -1
        y = rnorm(1, mean=outcome_means$a1_nr_cor_a2nr_2, sd=sd)
      } else if (dat$true_response_1[i]==1 & dat$treatment_2[i]==-1) { # misclassified as nonresponder and received a2nr = -1
        y = rnorm(1, mean=outcome_means$a1_nr_mis_a2nr_2, sd=sd)
      }
    } else {
      if (dat$true_response_1[i]==1 & dat$treatment_2[i]==0) { # correctly classified as a responder
        y = rnorm(1, mean=outcome_means$b1_r_cor, sd=sd)
      } else if (dat$true_response_1[i]==0 & dat$treatment_2[i]==0) { # misclassified as a responder
        y = rnorm(1, mean=outcome_means$b1_r_mis, sd=sd)
      } else if (dat$true_response_1[i]==0 & dat$treatment_2[i]==1) { # correctly classified as nonresponder and receives a2nr = 1
        y = rnorm(1, mean=outcome_means$b1_nr_cor_b2nr_1, sd=sd)
      } else if (dat$true_response_1[i]==1 & dat$treatment_2[i]==1) { # misclassified as nonresponder and received a2nr = 1
        y = rnorm(1, mean=outcome_means$b1_nr_mis_b2nr_1, sd=sd)
      } else if (dat$true_response_1[i]==0 & dat$treatment_2[i]==-1) { # correctly classified as nonresponder and received a2nr = -1
        y = rnorm(1, mean=outcome_means$b1_nr_cor_b2nr_2, sd=sd)
      } else if (dat$true_response_1[i]==1 & dat$treatment_2[i]==-1) { # misclassified as nonresponder and received a2nr = -1
        y = rnorm(1, mean=outcome_means$b1_nr_mis_b2nr_2, sd=sd)
      }
    }
  
    dat$outcome[i] = y
  }

  sens_label = ifelse(sensitivity>.5, "HIGH", "LOW")
  spec_label = ifelse(specificity>.5, "HIGH", "LOW")
  prev_label = ifelse(response_prev_arm1>.5, "HIGH", "LOW")
  save(dat, file=paste0("smart_simulations/results/EM_alloc_",EM_alloc, "/prev", prev_label, "_sens", sens_label, "_spec", spec_label, "/sim",sim,"/dat_upd.RData"))
  save(param_mat, file=paste0("smart_simulations/results/EM_alloc_",EM_alloc, "/prev", prev_label, "_sens", sens_label, "_spec", spec_label, "/sim",sim,"/param_mat_upd.RData"))
}





