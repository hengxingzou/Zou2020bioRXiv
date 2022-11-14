library(tidyverse)
library(magrittr)
library(reshape2)
library(patchwork)
library(minpack.lm)

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
# Note: relative stage difference is calculated as the arrival time of species 1 minus that of species 2,
#       i.e. the stage of species 2 minus that of species 1

generate_amat = function(s, alpha) {
  
  # Set up the matrix
  alpha_mat = matrix(0, 2*s, 2*s)
  
  # Go through rows
  for (i in 1:(2*s)) {
    
    # Go through columns
    for (j in 1:(2*s)) {
      
      # Calculate intraspecific competition (alpha11)
      if (i<(s-1) & j<(s-1)) {
        # alpha_mat[i, j] = comp_coeff(alpha[2], scal, xmid, (j-i), "v")
        alpha_mat[i, j] = comp_coeff(alpha[2], scal, xmid, (j-i), "c")
      } 
      
      # Calculate interspecific competition (alpha12)
      else if (i<(s-1) & j>s & j<(2*s-1)) {
        # alpha_mat[i, j] = comp_coeff(alpha[1], scal, xmid, (j-s-i), "v")
        alpha_mat[i, j] = comp_coeff(alpha[1], scal, xmid, (j-s-i), "c")
      } 
      
      # Calculate interspecific competition (alpha21)
      else if (i>s & i<(2*s-1) & j<(s-1)) {
        # alpha_mat[i, j] = comp_coeff(alpha[1], -scal, xmid, (i-s-j), "v")
        alpha_mat[i, j] = comp_coeff(alpha[1], -scal, xmid, (i-s-j), "c")
      } 
      
      # Calculate intraspecific competition (alpha22)
      else if (i>s & i<(2*s-1) & j>s & j<(2*s-1)){
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
  s = dim(alpha_mat)[1]/2
  
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
  para[2:(s-1), spp] = para[2:(s-1), spp]/(1 + dd[1:(s-2)])
  
  # Set up transition matrix
  tmat = matrix(0, s, s)
  
  # Set up transition matrix without density dependence
  for (i in 1:(s-2)) {
    tmat[i+1, i] = para[i+1, spp]*b
  }
  
  # Assign within-season reproduction: adults to J1
  tmat[1, s-1] = para[1, spp]*b
  
  # Assign adult survival
  tmat[s-1, s-1] = para[s, spp]*b
  
  # Assign emergence from dormancy
  tmat[1, s] = a
  
  # Assign end-of-season reproduction: adults to dormancy
  tmat[s, s-1] = para[1, spp]*(1-mu)*(1-b)
  
  # Assign remaining in dormancy 
  tmat[s, s] = 1-a
  
  return(tmat)
}


# Get Emergence Vector

# A function that decides the time of emergence from dormancy for both species
# Input: mean_diff, mean difference in stages, integer, of which absolute value should not exceed s-1;
#        var, variation from mean_diff determining the lower and upper bounds of uniform distribution from which the difference in 
#             arrival time will be drawn; when used, the sum of absolute values of mean_diff and var should not exceed s-1; integer
#        season_len, season length, integer
# Output: a sequence of time parameters a, with which 1 means emergence, of both species, in a matrix with columns being species

get_emergence = function(mean_diff, var, season_len, periodic = F) {
  
  # Get initial stage differences
  if (periodic) {delta_s = mean_diff+var}
  
  else {delta_s = round(runif(1, min = mean_diff-var, max = mean_diff+var))}
  
  # Set up "dummy" variables for emergence: 1 means emergence; columns are species
  a = matrix(0, nrow = season_len, ncol = 2)
  
  # If both species arrive at the same time
  if (delta_s == 0) {a[1, ] = 1}
  
  # Species 1 arrives first
  else if (delta_s < 0) {
    a[1, 1] = 1
    a[abs(delta_s)+1, 2] = 1} 
  
  # Species 2 arrives first
  else {
    a[1, 2] = 1
    a[abs(delta_s)+1, 1] = 1
  }
  
  return(a)
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
#        mean_diff, mean difference in stages, integer, of which absolute value should not exceed s-1;
#        var, variation from mean_diff determining the lower and upper bounds of uniform distribution from which the difference in 
#             arrival time will be drawn; when used, the sum of absolute values of mean_diff and var should not exceed s-1; integer
#        season_len, season length, integer
#        return_all, whether to return time series of both population, logical, default to TRUE;
#                    if FALSE, only returns the final population of the two species as a matrix, with columns being species
# Output: if return_all == TRUE (default), returns a list of two matrices with population dynamics: rows, stages; columns, time;
#         if return_all == FALSE, returns the final population of the two species as a matrix, with columns being species

seasons = function(time, pop, para, alpha, mean_diff, var, season_len, 
                   return_all = T, periodic = F, coexist_time = F) {
  
  # Calculate number of seasons from total time
  n_seasons = time %/% season_len
  
  # Set up matrices for recording population dynamics
  all_pop_1 = matrix(0, ncol = time, nrow = s)
  all_pop_2 = matrix(0, ncol = time, nrow = s)
  
  # Set up a vector for recording number of extinctions
  extinction_times = c()
  
  # Set up "dummy" variables for seasonality; 0 means the end of the season
  b = rep(ifelse(1:season_len %% season_len == 0, 0, 1), times = n_seasons)
  
  delta_s = mean_diff
  
  # Run simulations
  for (t in 1:time) {
    
    # Determine time step within a season
    timestep = ifelse(t %% season_len == 0, season_len, t %% season_len)
    
    # Regenerate emergence times at the beginning of the season
    if (t %% season_len == 1) {
      a_matrix = get_emergence(mean_diff, var, season_len, periodic)
      if (periodic) {var = -var}
    }
    
    # Generate transition matrices
    tmat1 = generate_tmat(1, generate_amat(s, alpha), pop, para, mu[1], a_matrix[timestep, 1], b[t])
    tmat2 = generate_tmat(2, generate_amat(s, alpha), pop, para, mu[2], a_matrix[timestep, 2], b[t])
    
    # Update population
    pop[, 1] = tmat1 %*% pop[, 1]
    pop[, 2] = tmat2 %*% pop[, 2]
    
    # Sample from previous states if one species extincts (<1e-6), then record this information
    
    if (coexist_time & (sum(pop[, 1]) < 1e-6 | sum(pop[, 2]) < 1e-6)) {
      # print(sum(pop[, 1]))
      # print(c(mean_diff, season_len))
      sampled_t = sample(1:(t-1), 1)
      # print(sampled_t)
      pop[, 1] = all_pop_1[, sampled_t]
      pop[, 2] = all_pop_2[, sampled_t]
      extinction_times = c(extinction_times, t)
    }
    
    # Record new population
    all_pop_1[, t] = pop[, 1]
    all_pop_2[, t] = pop[, 2]
  }
  
  # Return all population
  if (return_all) {return(list(all_pop_1, all_pop_2, extinction_times))}
  
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

stabledist = function(time, max_time, pop, para, alpha, mean_diff, season_len, periodic = F) {
  
  # Set up the counter
  loops = 1
  
  # Find the stable stage distribution of both species
  repeat {
    if (time*loops > max_time) {
      
      # Warning message if result still not stable after exceeding max running time
      warning(paste(c("exceeds max running time; return value may not be stable at\ ", 
                      mean_diff, "\ with season length\ ", season_len)))
      break}
    
    # Run the simulation
    new_pop = seasons(time, pop, para, alpha, mean_diff, var = 0, season_len, return_all = F, periodic = periodic)
    
    # Compare the end population of this round to the end population of the previous round
    if (sqrt(sum((c(new_pop[, 1], new_pop[, 2])-c(pop[, 1], pop[, 2]))^2) < 1e-12)) {break}
    
    # Update population
    pop = new_pop
    
    # Update counter
    loops = loops + 1
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

inv_growth = function(time, pop, para, alpha, mean_diff, var, season_len, inv_spp, return_all = F, periodic = F) {
  
  # Set up resident and invader
  res_spp = -inv_spp + 3
  
  # Calculate the stable stage distribution of the resident only
  pop[, inv_spp] = 0
  res = stabledist(time, max_time = 50*time, pop, para, alpha, mean_diff, season_len, periodic)
  
  # Set up population with the invader
  inv = res
  inv[, inv_spp] = rep(1/s, s)
  
  # Set up seasons
  lambdas = numeric(time)
  n_seasons = time %/% season_len
  b = rep(ifelse(1:season_len %% season_len == 0, 0, 1), times = n_seasons)
  
  # # Set up emergence times for the first season
  # a_matrix = get_emergence(mean_diff, var, season_len)
  
  # Set intraspecific competition of the invader to 0
  alpha[inv_spp+1] = 0
  
  for (t in 1:time) {
    
    # Determine time step within a season
    timestep = ifelse(t %% season_len == 0, season_len, t %% season_len)
    
    # Regenerate emergence times at the beginning of the season
    if (t %% season_len == 1) {
      a_matrix = get_emergence(mean_diff, var, season_len, periodic)
      if (periodic) {var = -var}
    }
    
    # Update resident population (doesn't contain the invader)
    tmat_res = generate_tmat(res_spp, generate_amat(s, alpha), res, para, mu[res_spp], a_matrix[timestep, res_spp], b[timestep])
    res[, res_spp] = tmat_res %*% res[, res_spp]
    
    # Update invader population (contains the resident)
    tmat_inv = generate_tmat(inv_spp, generate_amat(s, alpha), inv, para, mu[inv_spp], a_matrix[timestep, inv_spp], b[timestep])
    inv[, inv_spp] = tmat_inv %*% inv[, inv_spp]
    
    # Record lambda of the invader
    lambdas[t] = sum(inv[, inv_spp])
    
    # Normalize invader population
    inv[, inv_spp] = inv[, inv_spp]/sum(inv[, inv_spp])
    
    # Update resident population to the invader
    inv[, res_spp] = res[, res_spp]
  }
  
  # Return all lambdas
  if (return_all) {return(lambdas)}
  
  # Return just the long-term growth rate
  return(mean(log(lambdas)))
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
#        season_len, season length, integer
#        inv_spp, invader species, integer, 1 or 2
# Output: a tibble of long-term growth rates of invader population

invasion_multi = function(n_seasons, init_pop, para, alpha, delta_s_lst, var = 0, season_lens, seed = 1, periodic = F) {
  
  # Set up data frame for storing values
  all_out = tibble()
  
  for (ds in delta_s_lst) {
    for (len in season_lens) {
      out = foreach (inv_spp = 1:2, .combine = rbind) %dopar% {
        
        # Calculate invasion growth of species i when species j is at equilibrium
        # Set seed before calculating each lambda to make sure the arrival times are the same for each pair
        set.seed(seed)
        lambda_ij = inv_growth(n_seasons*len, init_pop, para, alpha, ds, var, len, inv_spp, periodic = periodic)
        
        # Calculate the invader's (species i's) growth rate without the resident (species j)
        pop = init_pop
        pop[7, -inv_spp] = 0
        lambda_i = inv_growth(n_seasons*len, pop, para, alpha, ds, var, len, inv_spp, periodic = periodic)
        
        values = data.frame(inv_i = inv_spp, res_j = -inv_spp+3, 
                            lambda_ij = lambda_ij, lambda_i = lambda_i,
                            sensitivity = (lambda_i-lambda_ij)/lambda_i,
                            delta_s = ds, season_len = len)
      }
      
      all_out = bind_rows(all_out, out)
      
      # Use a different seed to make sure the arrival times are different for the next pair
      set.seed(NULL)
    }
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
    
    # Numeric PE
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
