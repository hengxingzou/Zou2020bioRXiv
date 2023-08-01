library(tidyverse)
library(magrittr)
library(reshape2)
library(patchwork)
library(minpack.lm)
library(doParallel)

registerDoParallel(cores = 4)

options(dplyr.summarise.inform = FALSE) # repress dplyr summarise() info

color_scheme = c("#a61b29", "#0eb0c9", "#ECB88A", "#577C8A", "#622954")
shades_xch = c("#d6979d", "#c66d76", "#b6444f", "#a61b29", "#881722", "#6a121b", "#4c0d13")
shades_xch_2 = c("#c66d76", "#6a121b")
shades_kql = c("#91dbe6", "#65ccdc", "#39bed2", "#0eb0c9", "#0c91a5", "#097180", "#07505c")
shades_kql_2 = c("#65ccdc", "#097180")
shades_063 = c("#fdd0b9", "#fcbe9d", "#fbab81", "#fb9966", "#ce7e54", "#a06241", "#73462f")

# Scheme selected from traditional colors of China (zhongguose.com) and Japan (nipponcolors.com)


########## CORE FUNCTIONS ##########


# Interspecific Competition Coefficients

# A function for calculating interspecific competition coefficients that scales sigmoidally with relative stage differences
# Input: parameters for scaling function; all numeric except
#        type: the type of scaling function f with delta p; "v", variable (uses all other input parameters and returns f(delta p)), 
#                                                           "c", constant (uses only B and returns a constant value B/2)
# Output: interspecific competition coefficient; numeric

comp_coeff = function(B, scal, xmid, deltp, type) {
  if (type == "v") {return(B/(1+exp((xmid-deltp)/scal)))} # reciprocal scaling function assuming early-arrival advantage
  else if (type == "c") {return(B/2)} # constant scaling function
}


# Generate Competition Matrix

# A function that generates a matrix with all intra- and interspecific competition coefficients (stage specific)
# Input: s, number of total stages, integer; 
#        alpha, vector of three elements in the order of B (as in comp_coeff()), intraspecific competition of species 1 and 2;
#        note that actual intraspecific competition will be half of the input value in alpha
#        each competition coefficients can be either variable ("v") or constant ("c") regarding delta p; this is manually changed
#        in the function
# Output: a square matrix with dimension 2*s by 2*s
# Note: 1. relative stage difference is calculated as the arrival time of species 1 minus that of species 2,
#          i.e. the stage of species 2 minus that of species 1
#       2. to include adults in all density-dependence, change if ... else if statements

generate_amat = function(s, alpha) {
  
  # Set up the matrix
  alpha_mat = matrix(0, 2*s, 2*s)
  
  # Go through rows
  for (i in 1:(2*s)) {
    
    # Go through columns
    for (j in 1:(2*s)) {
      
      # Calculate intraspecific competition (alpha11)
      if (i<(s-1) & j<(s-1)) {
      # if (i<=(s-1) & j<=(s-1)) { # include adults in density-dependence
        
        # alpha_mat[i, j] = comp_coeff(alpha[2], scal, xmid, (j-i), "v")
        alpha_mat[i, j] = comp_coeff(alpha[2], scal, xmid, (j-i), "c")
      } 
      
      # Calculate interspecific competition (beta12)
      else if (i<(s-1) & j>s & j<(2*s-1)) {
      # else if (i<=(s-1) & j>s & j<=(2*s-1)) { # include adults in density-dependence
        
        alpha_mat[i, j] = comp_coeff(alpha[1], scal, xmid, (j-s-i), "v")
        # alpha_mat[i, j] = comp_coeff(alpha[1], scal, xmid, (j-s-i), "c")
      } 
      
      # Calculate interspecific competition (beta21)
      else if (i>s & i<(2*s-1) & j<(s-1)) {
      # else if (i>s & i<=(2*s-1) & j<=(s-1)) { # include adults in density-dependence
        
        alpha_mat[i, j] = comp_coeff(alpha[1], -scal, xmid, (i-s-j), "v")
        # alpha_mat[i, j] = comp_coeff(alpha[1], -scal, xmid, (i-s-j), "c")
      } 
      
      # Calculate intraspecific competition (alpha22)
      else if (i>s & i<(2*s-1) & j>s & j<(2*s-1)) {
      # else if (i>s & i<=(2*s-1) & j>s & j<=(2*s-1)) { # include adults in density-dependence
        
        # alpha_mat[i, j] = comp_coeff(alpha[3], scal, xmid, (j-i), "v")
        alpha_mat[i, j] = comp_coeff(alpha[3], scal, xmid, (j-i), "c")
      }
    }
  }
  
  return(alpha_mat)
}


# Generate Transition Matrix

# A function that generates the transition matrix of populations between each time step
# Input: spp, species, integer (1 or 2);
#        alpha_mat, competition matrix generated from generate_amat(), matrix;
#        pop, stage-specific population of the two species, matrix, where rows are stages and columns species;
#             stages 1-5 are juveniles, stage 6, adults, stage 7, dormant
#        para, baseline transition rates, in the order of fecundity (R), baseline transition rates of juvenile stages (P),
#              and adult survival (S); matrix, where columns are species
#        mu, mortality rates of stage 1 juveniles at the end of the season by species, vector
#        a, time parameter that indicates species emergence from dormancy, integer, 1 (emergence) or 0
#        b, time parameter that indicates the end of the season, integer, 0 (end of season) or 1
# Output: a square matrix with dimension s by s

generate_tmat = function(spp, alpha_mat, pop, para, mu, a, b) {
  
  # Extract population of species 1 and 2, and number of stages
  N1 = pop[, 1]
  N2 = pop[, 2]
  s = nrow(pop)
  
  # Calculate intraspecific competition 11
  intra11 = alpha_mat[1:s, 1:s] %*% N1
  
  # Calculate interspecific competition 12
  inter12 = alpha_mat[1:s, (s+1):(2*s)] %*% N2
  
  # Calculate interspecific competition 21
  inter21 = alpha_mat[(s+1):(2*s), 1:s] %*% N1
  
  # Calculate intraspecific competition 22
  intra22 = alpha_mat[(s+1):(2*s), (s+1):(2*s)] %*% N2
  
  # Calculate overall density dependence terms
  if (spp == 1) {dd = intra11 + inter12} else {dd = intra22 + inter21}
  
  # Calculate transition probabilities from density dependence
  # - if adult mortality and fecundity are NOT density-independent
  para[2:(s-1), spp] = para[2:(s-1), spp]/(1 + dd[1:(s-2)])
  
  # - if adult mnortality and fecundity are density-dependent
  # para[1:s, spp] = para[1:s, spp]/(1 + dd[1:s])
  
  # Set up transition matrix
  tmat = matrix(0, s, s)
  
  # Set up transition matrix with density dependence
  for (i in 1:(s-2)) {
    tmat[i+1, i] = para[i+1, spp]*b
  }
  
  # Assign within-season reproduction: adults to J1
  tmat[1, s-1] = para[1, spp]*b
  
  # Assign adult survival
  tmat[s-1, s-1] = para[s, spp]*b
  
  # Assign emergence from dormancy
  tmat[1, s] = b
  
  # Assign end-of-season reproduction: adults to dormancy
  tmat[s, s-1] = para[1, spp]*(1-mu)*(1-b)
  
  # Assign remaining in dormancy 
  tmat[s, s] = (1-b)
  
  return(tmat)
}


# Get Emergence Vector

# A function that decides the time of emergence from dormancy for both species
# Input: mean_diff, mean difference in stages, integer, of which absolute value should not exceed s-1;
#        var, variation from mean_diff determining the lower and upper bounds of uniform distribution from which the difference in 
#             arrival time will be drawn; when used, the sum of absolute values of mean_diff and var should not exceed s-1; integer
#        season_len, season length, integer
# Output: a sequence of time parameters a, with which 1 means emergence, of both species, in a matrix with columns being species

get_emergence = function(delta_s, lifecycle_len, season_len) {
  
  # Set up the function for emergence a(t): 1 means emergence; columns are species
  a = matrix(0, nrow = season_len, ncol = 2)

  # Set up the function for dormancy b(t): 1 means no dormancy, 0 means going to dormancy
  b = matrix(rep(ifelse(1:season_len %% season_len == 0, 0, 1), times = 2), 
             ncol = 2, byrow = F)
  
  # If both species arrive at the same time
  if (delta_s == 0) {a[1, ] = 1}
  
  # Species 1 arrives first
  else if (delta_s < 0) {
    
    # Species 1 emerges at the first time step
    a[1, 1] = 1
    # Species 2 emerges late
    a[abs(delta_s)+1, 2] = 1
    
    # Species 1 goes into dormancy earlier
    b[(lifecycle_len+1):season_len, 1] = 0
    # Species 2 emerges late
    b[1:(abs(delta_s)), 2] = 0
  } 
  
  # Species 2 arrives first
  else {
    
    # Species 2 emerges at the first time step
    a[1, 2] = 1
    # Species 1 emerges late
    a[abs(delta_s)+1, 1] = 1
    
    # Species 2 goes into dormancy earlier
    b[(lifecycle_len+1):season_len, 2] = 0
    # Species 1 emerges late
    b[1:(abs(delta_s)), 1] = 0
  }
  
  return(list(a, b))
}


########## SIMULATION ##########


# Seasons

# Do simulation with multiple seasons, allowing either consistent or variable relative stage differences at the beginning of each season
# Input: time, length of the simulation, integer, should be multiples of season_len
#        pop, stage-specific population of the two species, matrix, where rows are stages and columns species;
#             stages 1-5 are juveniles, stage 6, adults, stage 7, dormant
#        para, baseline transition rates, in the order of fecundity (R), baseline transition rates of juvenile stages (P),
#              and adult survival (S); matrix, where columns are species
#        alpha, vector of three elements in the order of B (as in comp_coeff()), intraspecific competition of species 1 and 2;
#        mean_diff (Delta s or mean Delta s), mean difference in stages, integer, of which absolute value should not exceed s-1;
#        var, variation from mean_diff determining the lower and upper bounds of uniform distribution from which the difference in 
#             arrival time will be drawn; when used, the sum of absolute values of mean_diff and var should not exceed s-1; integer
#        season_len, season length, integer
#        return_all, whether to return time series of both population, logical, default to TRUE;
#                    if FALSE, only returns the final population of the two species as a matrix, with columns being species
# Output: if return_all == TRUE (default), returns a list of two matrices with population dynamics: rows, stages; columns, time;
#         if return_all == FALSE, returns the final population of the two species as a matrix, with columns being species

seasons = function(n_seasons, pop, para, alpha, mean_diff, var, n_gens, 
                   return_all = T, periodic = F) {
  
  # Get number of stages from initial population
  s = nrow(pop)
  
  # Set up lists for recording population dynamics
  all_pop_1 = vector("list", length = n_seasons)
  all_pop_2 = vector("list", length = n_seasons)
  
  # Iterate through seasons
  for (i in 1:n_seasons) {
    
    # Get initial stage differences (Delta s) for the current season
    # - Periodic emergence
    if (periodic) {
      delta_s = mean_diff+var
      var = -var
    }
    
    # - Non-periodic emergence
    else {delta_s = round(runif(1, min = mean_diff-var, max = mean_diff+var))}
    
    # Calculate season lengths
    # - Lengths of the actual species life cycle before dormancy
    lifecycle_len = n_gens*(s-1)
    
    # - Extend the season to let the arriver finish development; +1 for the time step of turning dormant
    season_len = lifecycle_len+abs(delta_s)+1
    
    # Get emergence and dormancy times (functions a(t) and b(t))
    a_matrix = get_emergence(delta_s, lifecycle_len, season_len)[[1]]
    b_matrix = get_emergence(delta_s, lifecycle_len, season_len)[[2]]
    
    # Set up matrices for recording population dynamics within season i
    i_pop_1 = matrix(0, ncol = season_len, nrow = s)
    i_pop_2 = matrix(0, ncol = season_len, nrow = s)
    
    # Run the model within season i
    for (t in 1:season_len) {
      
      # Generate transition matrices
      tmat1 = generate_tmat(1, generate_amat(s, alpha), pop, para, mu[1], a_matrix[t, 1], b_matrix[t, 1])
      tmat2 = generate_tmat(2, generate_amat(s, alpha), pop, para, mu[2], a_matrix[t, 2], b_matrix[t, 2])
      
      # Update population
      pop[, 1] = tmat1 %*% pop[, 1]
      pop[, 2] = tmat2 %*% pop[, 2]

      # Record new population
      i_pop_1[, t] = pop[, 1]
      i_pop_2[, t] = pop[, 2]
    }
    
    # Record the population of season i
    all_pop_1[[i]] = i_pop_1
    all_pop_2[[i]] = i_pop_2
    
  }

  # Return all population
  if (return_all) {
    
    all_pop_1 = do.call(cbind, all_pop_1)
    all_pop_2 = do.call(cbind, all_pop_2)
    return(list(all_pop_1, all_pop_2))
    
  }
  
  # Return just the final population
  return(pop)
}


# Convert Output

# Convert output from seasons() into a data frame for further analysis
# Input: all_out, output from seasons(), a list with two matrices containing population dynamics of both species
# Output: a data frame with columns as stages and rows as population of each stage at each time

convert_output = function(all_out) {
  
  pop_df = matrix(ncol = 0, nrow = 0)
  
  for (spp in 1:2) {
    
    output = all_out[[spp]]
    time = ncol(output)
    s = nrow(output)
    all_pop = as.data.frame(t(output))
    colnames(all_pop) = rep(1:s)
    all_pop = cbind(all_pop, time = 1:time, species = factor(spp))
    pop_df = rbind(pop_df, all_pop)
  }
  
  return(pop_df)
}


# Get Stable Distribution Numerically

# Numerically find stable condition of seasonal population by comparing the final population of the previous season to that
# the current season; if the current initial population is close to that of the previous season (< 1e-12), 
# population is considered stable
# Input: time, length of the simulation, integer, should be multiples of season_len
#        max_time: maximum running time, integer; if reached, will return the final population without checking for stability
#        pop, stage-specific population of the two species, matrix, where rows are stages and columns species;
#             stages 1-5 are juveniles, stage 6, adults, stage 7, dormant
#        para, baseline transition rates, in the order of fecundity (R), baseline transition rates of juvenile stages (P),
#              and adult survival (S); matrix, where columns are species
#        alpha, vector of three elements in the order of B (as in comp_coeff()), intraspecific competition of species 1 and 2;
#        mean_diff, mean difference in stages, integer, of which absolute value should not exceed s-1;
#        var, variation from mean_diff determining the lower and upper bounds of uniform distribution from which the difference in 
#             arrival time will be drawn; when used, the sum of absolute values of mean_diff and var should not exceed s-1; integer
#        season_len, season length, integer
# Output: stable stage distributions of both species, in a matrix with columns being species

stabledist = function(n_seasons, max_n_seasons, pop, para, alpha, mean_diff, n_gens, periodic = F) {
  
  # Find the stable stage distribution of both species
  repeat {
    if (n_seasons > max_n_seasons) {
      
      # Warning message if result still not stable after exceeding max running time
      warning(paste(c("exceeds max running time; return value may not be stable at\ ", 
                      mean_diff, "\ with season length\ ", season_len)))
      break}
    
    # Run the simulation
    new_pop = seasons(n_seasons, pop, para, alpha, mean_diff, var = 0, n_gens, return_all = F, periodic = periodic)
    
    # Compare the end population of this round to the end population of the previous round
    if (sqrt(sum((c(new_pop[, 1], new_pop[, 2])-c(pop[, 1], pop[, 2]))^2) < 1e-12)) {break}
    
    # Update population
    pop = new_pop
    
  }
  
  return(new_pop)
}


######### MODEL ANALYSIS ##########


# Deterministic Invasion Analysis Using Frequency/Abundance

# A function that calculates invasion growth rate of a given species by the method from Benaïm and Schreiber (2019) J Math Biol
# Input: time, length of the simulation, integer, should be multiples of season_len
#        pop, stage-specific population of the two species, matrix, where rows are stages and columns species;
#             stages 1-5 are juveniles, stage 6, adults, stage 7, dormant
#        para, baseline transition rates, in the order of fecundity (R), baseline transition rates of juvenile stages (P),
#              and adult survival (S); matrix, where columns are species
#        alpha, vector of three elements in the order of B (as in comp_coeff()), intraspecific competition of species 1 and 2;
#        mean_diff, mean difference in stages, integer, of which absolute value should not exceed s-1;
#        season_len, season length, integer
#        inv_spp, invader species, integer, 1 or 2
# Output: long-term growth rate of invader population

inv_growth = function(n_seasons, pop, para, alpha, mean_diff, var, n_gens, inv_spp, return_all = F, periodic = F) {
  
  # Set up resident and invader
  res_spp = -inv_spp + 3
  
  # Calculate the stable stage distribution of the resident only
  pop[, inv_spp] = 0
  res = stabledist(n_seasons, max_n_seasons = 50*n_seasons, pop, para, alpha, mean_diff, n_gens, periodic)
  
  # Set up population with the invader
  inv = res
  inv[, inv_spp] = rep(1/s, s)
  
  # Set intraspecific competition of the invader to 0
  alpha[inv_spp+1] = 0
  
  lambdas_all = numeric()
  
  # Iterate through seasons
  for (i in 1:n_seasons) {
    
    # Get initial stage differences (Delta s) for the current season
    # - Periodic emergence
    if (periodic) {
      delta_s = mean_diff+var
      # print(paste0("delta s=", delta_s, "var=", var))
      var = -var
    }
    
    # - Non-periodic emergence
    else {delta_s = round(runif(1, min = mean_diff-var, max = mean_diff+var))}
    
    # Calculate season lengths
    # - Lengths of the actual species life cycle before dormancy
    lifecycle_len = n_gens*(s-1)
    
    # - Extend the season to let the arriver finish development; +1 for the time step of turning dormant
    season_len = lifecycle_len+abs(delta_s)+1

    # Get emergence and dormancy times (functions a(t) and b(t))
    a_matrix = get_emergence(delta_s, lifecycle_len, season_len)[[1]]
    b_matrix = get_emergence(delta_s, lifecycle_len, season_len)[[2]]
    
    lambdas = numeric(season_len)
    
    # Run the model within season i
    for (t in 1:season_len) {
      
      # Update resident population (doesn't contain the invader)
      tmat_res = generate_tmat(res_spp, generate_amat(s, alpha), res, para, mu[res_spp], a_matrix[t, res_spp], b_matrix[t, res_spp])
      res[, res_spp] = tmat_res %*% res[, res_spp]
      
      # Update invader population (contains the resident)
      tmat_inv = generate_tmat(inv_spp, generate_amat(s, alpha), inv, para, mu[inv_spp], a_matrix[t, inv_spp], b_matrix[t, inv_spp])
      inv[, inv_spp] = tmat_inv %*% inv[, inv_spp]
      
      # Record lambda of the invader
      lambdas[t] = sum(inv[, inv_spp])
      
      # Normalize invader population
      inv[, inv_spp] = inv[, inv_spp]/sum(inv[, inv_spp])
      
      # Update resident population to the invader
      inv[, res_spp] = res[, res_spp]
      
    }
    
    lambdas_all = c(lambdas_all, lambdas)
    
  }
  
  # Return all lambdas
  if (return_all) {return(lambdas_all)}
  
  # Return just the long-term growth rate
  return(mean(log(lambdas_all)))
}


# Invasion Analysis with Multiple Arrival Times and Season Lengths

# Wrapping inv_growth for simulating multiple initial conditions
# Input: time, length of the simulation, integer, should be multiples of season_len
#        pop, stage-specific population of the two species, matrix, where rows are stages and columns species;
#             stages 1-5 are juveniles, stage 6, adults, stage 7, dormant
#        para, baseline transition rates, in the order of fecundity (R), baseline transition rates of juvenile stages (P),
#              and adult survival (S); matrix, where columns are species
#        alpha, vector of three elements in the order of B (as in comp_coeff()), intraspecific competition of species 1 and 2;
#        mean_diff, mean difference in stages, integer, of which absolute value should not exceed s-1;
#        n_gens, season lengths, in numbers of generations, vector with integers
#        inv_spp, invader species, integer, 1 or 2
# Output: a tibble of long-term growth rates of invader population

invasion_multi = function(n_seasons, init_pop, para, alpha, mean_diff_lst, var = 0, n_gens_lst, seed = 1, periodic = F) {
  
  # Loop through season lengths
  all_out = foreach(n_gens = n_gens_lst, .combine = rbind) %dopar% {

    # Delta s longer than one generation
    # mean_diff_lst = round(c(-0.33333, -0.25, -0.16667, -0.08333, 0, 0.08333, 0.16667, 0.25, 0.33333)*n_gens*(s-1))

    out = foreach(md = mean_diff_lst, .combine = rbind) %dopar% {
      
      values = tibble()
      for (inv_spp in 1:2) {
        
        # Set seed before calculating each lambda to make sure the arrival times are the same for each pair
        # Not working -- not included in the current manuscript
        # set.seed(seed)
        
        # Calculate invasion growth of species i when species j is at equilibrium
        lambda_ij = inv_growth(n_seasons, init_pop, para, alpha, md, var, n_gens, inv_spp, periodic = periodic)
        
        # Calculate the invader's (species i's) growth rate without the resident (species j)
        pop = init_pop
        pop[7, -inv_spp] = 0
        lambda_i = inv_growth(n_seasons, pop, para, alpha, md, var, n_gens, inv_spp, periodic = periodic)
        
        values = rbind(values, tibble(inv_i = inv_spp, res_j = -inv_spp+3, 
                                      lambda_ij = lambda_ij, lambda_i = lambda_i,
                                      sensitivity = (lambda_i-lambda_ij)/lambda_i,
                                      delta_s = md, n_gens = n_gens))
      }
      
      values
    }
    
    # Use a different seed to make sure the arrival times are different for the next pair
    # Not working -- not included in the current manuscript
    # set.seed(NULL)
    
    out
    
  }
  
  return(all_out)
}


######### PLOTTING ND-RFD ##########


generate_baseplot = function(ND_min, ND_max, RFD_min, RFD_max, alpha_val = 0.25) {
  
  ND_range = seq(ND_min, ND_max, 0.001)
  nd_rfd_line = rbind(data.frame(ND = ND_range, RFD = 1/(1-ND_range), type = "upper"),
                      data.frame(ND = ND_range, RFD = 1-ND_range, type = "lower"))
  baseplot =
    ggplot() +
    
    geom_line(data = nd_rfd_line %>% filter(type == "lower"), aes(x = ND, y = log10(RFD)), color = "black") + 
    geom_line(data = nd_rfd_line %>% filter(type == "upper"), aes(x = ND, y = log10(RFD)), color = "black") + 
    
    # Coexistence
    geom_ribbon(data = nd_rfd_line %>% filter(type == "lower", ND >= 0), 
                aes(x = ND, ymin = log10(RFD), ymax = log10(1/(1-ND))), fill = "#ECB88A", alpha = alpha_val) +
    
    # Frequency-Dependent PE
    geom_ribbon(data = nd_rfd_line %>% filter(type == "upper", ND <= 0), 
                aes(x = ND, ymin = log10(RFD), ymax = log10(1-ND)), fill = "#577C8A", alpha = alpha_val) +
    
    # Species 2 wins
    geom_ribbon(data = nd_rfd_line %>% filter(type == "upper", ND >= 0), 
                aes(x = ND, ymin = log10(RFD), ymax = RFD_max), fill = "#0eb0c9", alpha = alpha_val) +
    geom_ribbon(data = nd_rfd_line %>% filter(type == "lower", ND <= 0), 
                aes(x = ND, ymin = log10(RFD), ymax = RFD_max), fill = "#0eb0c9", alpha = alpha_val) +
    
    # Species 1 wins
    geom_ribbon(data = nd_rfd_line %>% filter(type == "lower", ND >= 0), 
                aes(x = ND, ymin = RFD_min, ymax = log10(RFD)), fill = "#a61b29", alpha = alpha_val) +
    geom_ribbon(data = nd_rfd_line %>% filter(type == "upper", ND <= 0), 
                aes(x = ND, ymin = RFD_min, ymax = log10(RFD)), fill = "#a61b29", alpha = alpha_val) +
    
    # 0 lines
    geom_hline(yintercept = 0, linetype = 2, color = "#577C8A") + 
    geom_vline(xintercept = 0, linetype = 2, color = "#577C8A") + 
    theme_bw() + theme(panel.grid = element_blank()) + 
    xlim(ND_min, ND_max) + ylim(RFD_min, RFD_max)
  
  return(baseplot)
}
