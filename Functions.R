# R Script for Simulations of Manuscript:
# Priority effects and season length shape long-term competition dynamics
# by Heng-Xing Zou and Volker H. W. Rudolf

# Packages for data wrangling

library(tidyverse)
library(magrittr)

# Packages for plotting

library(reshape2)
library(patchwork)
library(scatterpie)

# Color schemes are selected from traditional colors of China (zhongguose.com) and Japan (nipponcolors.com)

color_scheme = c("#a61b29", "#0eb0c9", "#ECB88A", "#577C8A", "#622954")
shades_xch = c("#d6979d", "#c66d76", "#b6444f", "#a61b29", "#881722", "#6a121b", "#4c0d13")
shades_xch_2 = c("#c66d76", "#6a121b")
shades_kql = c("#91dbe6", "#65ccdc", "#39bed2", "#0eb0c9", "#0c91a5", "#097180", "#07505c")
shades_kql_2 = c("#65ccdc", "#097180")
shades_063 = c("#fdd0b9", "#fcbe9d", "#fbab81", "#fb9966", "#ce7e54", "#a06241", "#73462f")


# IMPORTANT NOTES
# 1. "Cycle" and "season" are used interchangeably in this script; they both refer to "season" used in the manuscript
# 2. "deltp" in this script denotes initial stage difference (delta s) in the manuscript
# 2. Run this script before running any chunks in Analysis.Rmd
# 3. Some lines may need to be modified when running specific chunks in Analysis.Rmd; see specific instructions
# 4. Some features from this script are not implemented in actual model analysis (Analysis.Rmd)
# 5. Some functions contain hard-coded values

# Please contact the author (Heng-Xing Zou) for questions or suggestions


########## MODEL SETUP ##########


# Competition Coefficients


# A function for calculating interspecific competition coefficients that scales sigmoidally with arrival time
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
#        alpha, vector of three elements in the order of B (as in comp_coeff()), intraspecific competition of species 1 and 2
#        each competition coefficients can be either variable ("v") or constant ("c") regarding delta p
# Output: a square matrix with dimension 2*s
# Note: deltp (as defined in comp_coeff) is the arrival time of species 1 minus that of species 2,
#       i.e. the stage of species 2 minus that of species 1

generate_amat = function(s, alpha) {
  alpha_mat = matrix(0, 2*s, 2*s)
  for (i in 1:(2*s)) { # rows
    for (j in 1:(2*s)) { # columns
      if (i<s & j<s) { # intraspecific-alpha11s
        # alpha_mat[i, j] = comp_coeff(alpha[2], scal, xmid, (j-i), "v")
        alpha_mat[i, j] = comp_coeff(alpha[2], scal, xmid, (j-i), "c")
      } else if (i<s & j>s & j<2*s) { # interspecific-alpha12s
        alpha_mat[i, j] = comp_coeff(alpha[1], scal, xmid, (j-s-i), "v")
        # alpha_mat[i, j] = comp_coeff(alpha[1], scal, xmid, (j-s-i), "c")
      } else if (i>s & i<2*s & j<s) { # interspecific-alpha21s
        alpha_mat[i, j] = comp_coeff(alpha[1], -scal, xmid, (i-s-j), "v")
        # alpha_mat[i, j] = comp_coeff(alpha[1], -scal, xmid, (i-s-j), "c")
      } else if (i>s & i<2*s & j>s & j<2*s){ # intraspecific-alpha22s
        # alpha_mat[i, j] = comp_coeff(alpha[3], scal, xmid, (j-i), "v")
        alpha_mat[i, j] = comp_coeff(alpha[3], scal, xmid, (j-i), "c")
      }
    }
  }
  return(alpha_mat)
}


# Generate Transition Matrix


# A function that generates the transition matrix of populations between each time step
# Input: s, number of total stages, integer;
#        para, vector of length (s+1) with elements in the order of fecundity, baseline survival (U0) of juveniles, and adult survival
#            from juvenile stage n to n+1;
# Output: a square matrix with dimension s

generate_A = function(s, para) {
  mat = matrix(0, s, s)
  for (i in 1:(s-1)) {
    mat[i+1, i] = para[i+1] # total survival
  }
  mat[1, s] = para[1] # fecundity
  mat[s, s] = para[s+1] # adult survival
  return(mat)
}


# Scaling Function for Per-Capita Competition


# A function that returns a transition matrix scaled by intra- and interspecific competition;
#        N, a vector of population of both species, in the order of species 1 to 2, stage 1 to 6 (adult);
#        para, data frame with columns as vectors of fecundity and baseline survival
#        alpha, vector of three elements in the order of B (as in comp_coeff()), intraspecific competition of species 1 and 2
#        record, integer with default 0
# Output: by default, a transition matrix of both species with dimension 2*s * 2*s (s is the number of stages)
#         if record = 1, return interaction strengths experienced by each species in a vector: c(intra11, inter12, inter21, intra22)
#         if record = 2, return competition coefficients in a vector: c(alpha11, alpha12, alpha21, alpha22)
# Note: calculations of interaction strengths and competition coefficients do not include adults because including 0 skews weighted means

comp_sc = function(N, para, alpha, record = 0) {
  s = length(N)/2
  alpha_mat = generate_amat(s, alpha)
  intra11 = alpha_mat[1:s, 1:s] %*% N[1:s] # intraspecific 11
  inter12 = alpha_mat[1:s, (s+1):(2*s)] %*% N[(s+1):(2*s)] # interspecific 12
  inter21 = alpha_mat[(s+1):(2*s), 1:s] %*% N[1:s] # interspecific 21
  intra22 = alpha_mat[(s+1):(2*s), (s+1):(2*s)] %*% N[(s+1):(2*s)] # intraspecific 22
  para[2:(s+1), 1] = para[2:(s+1), 1]/(1 + intra11 + inter12)
  para[2:(s+1), 2] = para[2:(s+1), 2]/(1 + intra22 + inter21)
  mat1 = generate_A(s, para[[1]])
  mat2 = generate_A(s, para[[2]])
  zero_mat = matrix(0, nrow = s, ncol = s)
  if (record == 1) { # return interactions experienced by two species, weighted by stage distribution
    return(c(weighted.mean(intra11[1:(s-1)], N[1:(s-1)]), weighted.mean(inter12[1:(s-1)], N[1:(s-1)]),
             weighted.mean(inter21[1:(s-1)], N[(s+1):(2*s-1)]), weighted.mean(intra22[1:(s-1)], N[(s+1):(2*s-1)])))
  }
  if (record == 2) { # return competition coefficients weighted by stage distribution; no adults
    return(c(weighted.mean(intra11[1:(s-1)], N[1:(s-1)])/sum(N[1:(s-1)]), 
             weighted.mean(inter12[1:(s-1)], N[1:(s-1)])/sum(N[(s+1):(2*s-1)]),
             weighted.mean(inter21[1:(s-1)], N[(s+1):(2*s-1)])/sum(N[1:(s-1)]), 
             weighted.mean(intra22[1:(s-1)], N[(s+1):(2*s-1)])/sum(N[(s+1):(2*s-1)])))
  }
  return(cbind(rbind(mat1, zero_mat), rbind(zero_mat, mat2)))
}


# Get Relative Arrival Time (Initial Stage Difference)

# A function that generates a random relative arrival time
# Input: mean_time, mean relative arrival time, integer of which absolute value should not exceed s-1;
#        var: variation from mean determining the lower and upper bounds of uniform distribution from which the relative
#             arrival time will be drawn; when used, the sum of absolute values of mean_time and var should not exceed s-1
# Output: a random difference in arrival time, integer

get_deltp = function(mean_time, var) {
  return(round(runif(1, min = mean_time-var, max = mean_time+var)))
}


########## SIMULATION ##########


# Multiple Steps


# Do simulation for a specified number of time steps
# Input: time, total time steps of simulation, integer
#        para, data frame with columns as vectors of fecundity and baseline survival (see function generate_A(), para)
#        init_pop, input inital population, data frame with columns as vectors of population (see function comp_sc(), N_intra and N_inter)
#        alpha, identical as previously defined parameters
#        record, integer with default 0
# Output: by default, a time series of population of two species with specific stages, a data frame with four columns: 
#             stage, N1 (population of species 1), N2, time
#         if record = 1, returns a data frame with five columns: four interaction strengths (as defined in comp_sc()), time
#         if record = 2, returns a data frame with five columns: four competition coefficients (as defined in comp_sc()), time
#         if record = 3, returns the overall transition matrix of all time steps

multiple_steps = function(init_pop, time, para, alpha, record = 0) {
  s = nrow(init_pop)
  out = matrix(NA, nrow = 2*s, ncol = time)
  out[, 1] = c(init_pop[[1]], init_pop[[2]])
  if (record < 3 & record > 0) {
    out_alt = matrix(NA, nrow = 4, ncol = time)
    out_alt[, 1] = comp_sc(out[, 1], para, alpha, record)
  }
  all_A = diag(2*s)
  for (t in 2:time) {
    A = comp_sc(out[, t-1], para, alpha)
    new_pop = A %*% out[, t-1]
    out[, t] = new_pop
    if (record < 3 & record > 0) {
      out_alt[, t] = comp_sc(out[, t], para, alpha, record)
    } 
    if (record == 3) {all_A = A %*% all_A}
  }
  out = out[, 1:time]
  out_df = as.data.frame(t(cbind(out[1:s, ], out[(s+1):(2*s), ])))
  colnames(out_df) = 1:s
  out_df$time = rep(1:time, 2)
  out_df$species = c(rep("N1", time), rep("N2", time))
  out_df %<>% pivot_longer(cols = 1:6, names_to = "stage") %>% pivot_wider(names_from = species, values_from = value)
  if (record < 3 & record > 0) {
    out_alt = out_alt[, 1:time]
    out_alt_df = as.data.frame(t(out_alt))
    colnames(out_alt_df) = c("intra11", "inter12", "inter21", "intra22")
    out_alt_df$time = 1:time
    return(out_alt_df)
  } 
  if (record == 3) {return(all_A)}
  return(out_df)
}


# Different Arrival Time


# Do simulation with multiple time steps and different arrival time of species
# Input: as defined in multiple_steps(), except
#        deltp: difference in arrival time, integer of which absolute value should not exceed s-1; 
#               if deltp < 0, species 1 arrives earlier; if deltp > 0, species 2 arrives earlier
# Output: as described in multiple_steps()

diff_time = function(time, para, init_pop, alpha, deltp, record = 0) {
  if (deltp == 0) {
    return(multiple_steps(init_pop, time, para, alpha, record))
  }
  if (deltp < 0) { # 1 arrives earlier
    pop = data.frame("N1" = init_pop[[1]], "N2" = rep(0, s))
  } else {
    pop = data.frame("N1" = rep(0, s), "N2" = init_pop[[2]])
  }
  out_1 = multiple_steps(pop, abs(deltp)+1, para, alpha)
  if (record != 0) {
    out_alt_1 = multiple_steps(pop, abs(deltp)+1, para, alpha, record)
  }
  new_pop = out_1 %>% tail(s) %>% select(N1:N2)
  out_1 = out_1[-((nrow(out_1)-s+1):(nrow(out_1))), ]
  for (i in 1:2) {
    if (sum(new_pop[i]) == 0) {
      new_pop[i] = init_pop[i]
    }
  }
  out_2 = multiple_steps(new_pop, time-abs(deltp), para, alpha)
  out_2$time = out_2$time + abs(deltp)
  if (record < 3 & record > 0) {
    out_alt_2 = multiple_steps(new_pop, time-abs(deltp), para, alpha, record)
    out_alt_2$time = out_alt_2$time + abs(deltp)
    out_alt_1 = out_alt_1[-nrow(out_alt_1), ]
    return(rbind(out_alt_1, out_alt_2))
  }
  if (record == 3) {
    out_alt_2 = multiple_steps(new_pop, time-abs(deltp), para, alpha, record)
    return(out_alt_2 %*% out_alt_1)
  }
  return(rbind(out_1, out_2))
}


# Cycles


# Do simulation with multiple cycles (seasons), allowing difference in arrival time at the beginning of each cycle, either consistent or variable (random)
# Input: as defined in diff_time(), except
#        mean_time, mean difference in arrival time, integer of which absolute value should not exceed s-1; 
#            if mean_time < 0, species 1 arrives earlier on average; if mean_time > 0, species 2 arrives earlier on average
#        var, variation of difference in arrival time between cycles, vector with interger elements of which absolute value should not exceed s-1
#        cycle_len, length of each cycle (season), nonzero integer divisible by s
# Output: if record = 0 (default), 1 or 2, as described in diff_time(), except for an additional column of variation
#         if record = 3, returns the initial population of the next season, calculated from final population of the final season

cycles = function(time, para, init_pop, alpha, mean_time, var = c(0), cycle_len, record = 0) {
  all_N = tibble()
  if (record < 3 & record > 0) {
    all_out = tibble()
  }
  cycle_tp = seq(1, time+cycle_len, cycle_len)
  for (v in var) {
    pop = init_pop
    N = tibble()
    if (record < 3 & record > 0) {
      out_alts = tibble()
    }
    for (i in 1:(length(cycle_tp)-1)) {
      out = diff_time(cycle_len, para, pop, alpha, get_deltp(mean_time, v))
      if (record < 3 & record > 0) {
        out_alt = diff_time(cycle_len, para, pop, alpha, get_deltp(mean_time, v), record)
      }
      new_pop_1 = tail(out, s)$N1
      new_pop_2 = tail(out, s)$N2
      pop = data.frame("N1" = rep(0, s), "N2" = rep(0, s))
      pop[1, 1] = new_pop_1[s]*para[1, 1]*surv[1] # only the adult stage contributes to next cycle
      pop[1, 2] = new_pop_2[s]*para[1, 2]*surv[2]
      out$time = out$time + (cycle_len)*(i-1)
      out$variation = rep(v, length(out$time))
      old_N = rep(0, 2*s)
      if (i > 1) {old_N = c(tail(N, s)$N1, tail(N, s)$N2)}
      N = bind_rows(N, out)
      if (record < 3 & record > 0) {
        out_alt$time = out_alt$time + (cycle_len)*(i-1)
        out_alt$variation = rep(v, length(out_alt$time))
        out_alts = bind_rows(out_alts, out_alt)
      }
      new_N = c(new_pop_1, new_pop_2)
      if (record == 3 & sqrt(sum((new_N-old_N)^2)) < 1e-12) {
        return(pop) # stable initial population
        break
      }
    }
    all_N = bind_rows(all_N, N)
    if (record < 3 & record > 0) {
      all_out = bind_rows(all_out, out_alts)
    }
  }
  if (record < 3 & record > 0) return(all_out)
  if (record == 3) return(pop)
  return(all_N)
}


######### SIMULATIONS ##########


# Get Stable Distribution Numerically


# Numerically find stable condition of seasonal population. Compare the initial population (calculated from final population of the last season);
# if the current initial population is close to that of the previous season (< 1e-12), break from loop and population is considered stable
# Input: as previously defined
# Output: initial population of a season during which stable population is reached; stage 1 juveniles only

stabledist = function(para, init_pop, alpha, deltp, cycle_len) {
  repeat{
    new_pop = cycles(cycle_len, para, init_pop, alpha, deltp, cycle_len = cycle_len, record = 3)
    if (sqrt(sum((c(new_pop[, 1], new_pop[, 2])-c(init_pop[, 1], init_pop[, 2]))^2)) < 1e-12) break
    init_pop = new_pop
  }
  return(new_pop)
}


# Cycles with Invasion


# A convenient wrapper of cycles() that allows one species grow until the other species invades at a small initial population (0.1)
# Input: as defined in cycle(), except
#        inv_spp, integer, invasion species (1 or 2)
#        inv_pop, invasion population (stage 1 only); default 0.01
#        invasion_time, time until the invader species appear, integer of which absolute value should not exceed time; 
#            if invasion species is not specified, invasion_time < 0, species 1 arrives first and species 2 invades; 
#            if invasion_time > 0, species 2 arrives first and species 1 invades
# Output: if record = 0 (default), 1 or 2, as described in diff_time(), except for an additional column of variation

cycles_invade = function(time, para, init_pop, alpha, mean_time, var, cycle_len, inv_pop = 0.01, record = 0, 
                         invasion_time, inv_spp) {
  s = nrow(init_pop)
  if (missing(invasion_time)) { # start invasion from resident equilibrium
    if (missing(inv_spp)) stop ("Invading species not specified")
    init_pop[[inv_spp]] = rep(0, s)
    stable_pop = stabledist(para, init_pop, alpha, mean_time, cycle_len)
    stable_pop[1, inv_spp][1] = inv_pop
    return(cycles(time, para, stable_pop, alpha, mean_time, var, cycle_len, record))
  }
}


# Deterministic Simulations


## Invasion Analysis Using Analytical Growth Rate

# **********CONTAINS HARD CODED VALUES **********

# Calculates invasion growth rate as the dominant eigenvalue of the transition matrix of the first season
# Input: as previously defined
# Output: invasion growth rate, real part of the dominant eigenvalue calculated from transition matrix of the first season

analytic_growth = function(para, init_pop, alpha, mean_time, cycle_len, 
                           inv_pop = 0.01, inv_spp) {
  init_pop[[inv_spp]] = 0
  init_pop = stabledist(para, init_pop, alpha, mean_time, cycle_len)
  init_pop[1, inv_spp] = inv_pop
  alpha[inv_spp+1] = 0 # invader does not experience intraspecific competition
  A = diff_time(cycle_len, para, init_pop, alpha, mean_time, record = 3)
  if (inv_spp == 1) {A = A[1:6, 1:6]}
  if (inv_spp == 2) {A = A[7:12, 7:12]}
  S = matrix(0, 6, 6)
  S[1, 6] = para[1, inv_spp]*surv[inv_spp]
  A = S%*%A
  lambda = eigen(A)$values[1]
  return(Re(lambda))
}


## Invasion Analysis Using Analytical Growth Rate


# Calculates invasion growth rate as difference between invasion population and final initial population (as defined in stabledist())
# Input: as previously defined
# Output: invasion growth rate, numerically determined from simulations

numeric_growth = function(para, init_pop, alpha, mean_time, cycle_len, inv_pop = 0.01, inv_spp) {
  init_pop[[inv_spp]] = 0
  init_pop = stabledist(para, init_pop, alpha, mean_time, cycle_len)
  init_pop[1, inv_spp] = inv_pop
  final_pop = stabledist(para, init_pop, alpha, mean_time, cycle_len)
  return(final_pop[1, inv_spp]-inv_pop)
}


## Phase Diagram Using Numeric Invasion Growth

# A function that creates a matrix for phase diagram of competition results with changing difference of arrival time and cycle (season) length.
#     If one species has invasion growth rate < 0, the other species wins;
#     if both > 0, the two species could coexist; if both < 0, numerical priority effects occur and species that arrives first wins.
#     This function uses numeric growth rate calculated by numeric_growth()
# Input: as defined in invasion_rates(), except
#        deltp_lst, a vector of differences in arrival time, usually generated by seq() function
# Output: a matrix with rows as differences in arrival time, columns as baseline interspecifics competition coefficients and entries as 
#         integer 1, 2, 3 or 4, denoting species 1 or 2, coexistence, or priority effects, respectively

numeric_phase_invgrowth = function(para, init_pop, alpha, deltp_lst, cycle_lens) {
  phase = matrix(NA, nrow = length(deltp_lst), ncol = length(cycle_lens))
  for (i in 1:length(deltp_lst)) {
    for (j in 1:length(cycle_lens)) {
      rate = c(NA, NA)
      for (inv_spp in 1:2) { # invading species
        rate[inv_spp] = numeric_growth(para, init_pop, alpha, deltp_lst[i], cycle_lens[j], inv_spp = inv_spp)
      }
      if (rate[1] > 0 & rate[2] > 0) { # coexistence
        phase[i, j] = 3
      } else if (rate[1] > 0 & rate[2] < 0) { # 1 wins
        phase[i, j] = 1
      } else if (rate[1] < 0 & rate[2] > 0) { # 2 wins
        phase[i, j] = 2
      } else if (rate[1] < 0 & rate[2] < 0) { # priority effects
        phase[i, j] = 4
      }
    }
  }
  return(phase)
}


## Phase Diagram Using Analytical Invasion Growth

# A function that creates a matrix for phase diagram of competition results with changing difference of arrival time and cycle (season) length.
#     If one species has invasion growth rate < 1, the other species wins;
#     if both > 1, the two species could coexist; if both < 1, numerical priority effects occur and species that arrives first wins.
#     This function uses analytical growth rate calculated by analytic_growth()
# Input: as defined in invasion_rates(), except
#        deltp_lst, a vector of differences in arrival time, usually generated by seq() function
# Output: a matrix with rows as differences in arrival time, columns as baseline interspecifics competition coefficients and entries as 
#         integer 1, 2, 3 or 4, denoting species 1 or 2, coexistence, or priority effects, respectively

analytic_phase_invgrowth = function(para, init_pop, alpha, deltp_lst, cycle_lens) {
  phase = matrix(NA, nrow = length(deltp_lst), ncol = length(cycle_lens))
  for (i in 1:length(deltp_lst)) {
    for (j in 1:length(cycle_lens)) {
      rate = c(NA, NA)
      for (inv_spp in 1:2) { # invading species
        rate[inv_spp] = analytic_growth(para, init_pop, alpha, deltp_lst[i], cycle_lens[j], inv_pop = 0.01, inv_spp = inv_spp)
      }
      if (rate[1] > 1 & rate[2] > 1) { # coexistence
        phase[i, j] = 3
      } else if (rate[1] > 1 & rate[2] < 1) { # 1 wins
        phase[i, j] = 1
      } else if (rate[1] < 1 & rate[2] > 1) { # 2 wins
        phase[i, j] = 2
      } else if (rate[1] < 1 & rate[2] < 1) { # priority effects
        phase[i, j] = 4
      }
    }
  }
  return(phase)
}


# Stochastic Simulations


## Trajectory Analysis

# Looks at the whole population dynamics (trajectory) of stochastic simulations
# Input: as previously defined, except
#        n_cycles, integer, number of cycles (seasons) used for simulations
# Output: population dynamics of all simulations
# Note: only use single value of variation, e.g. c(0) or c(3)

traj_stoc = function(n_cycles, para, init_pop, alpha, mean_time, var, cycle_len, rept, inv_spp, inv_pop = 0.01) {
  N = tibble()
  total_time = n_cycles*cycle_len
  for (i in 1:rept) {
    init_pop[[inv_spp]] = rep(0, s)
    stable_pop = stabledist(para, init_pop, alpha, mean_time, cycle_len)
    stable_pop[1, inv_spp][1] = inv_pop
    result = cycles(total_time, para, stable_pop, alpha, mean_time, var, cycle_len, record = 0)
    result %<>% pivot_longer(cols = N1:N2, names_to = "species", values_to = "population") %>% 
      group_by(time, species, variation) %>% summarize(population = sum(population))
    result$cycle_no = rep(1:n_cycles, each = cycle_len)
    result$rept = i
    N = bind_rows(N, result)
  }
  return(N)
}


## Invasion Trajectory Analysis

# A function that is similar to traj_stoc() but takes a range of cycle lengths and relative arrival times
# Inoput: as previously defined
# Output: population dynamics of all simulations

inv_traj = function(n_cycles, para, init_pop, alpha, deltp_lst, var, cycle_lens, rept) {
  all_res = tibble()
  for (mean_time in deltp_lst) {
    bind_df = tibble()
    for (spp in 1:2) {
      bind_df_2 = tibble()
      for (cycle_len in cycle_lens) {
        result = traj_stoc(n_cycles, para, init_pop, alpha, mean_time, var, cycle_len, rept, spp)
        result$cycle_len = cycle_len
        bind_df_2 = bind_rows(bind_df_2, result)
      }
      bind_df_2$inv_spp = spp
      bind_df = bind_rows(bind_df, bind_df_2)
    }
    bind_df$deltp = mean_time
    all_res = bind_rows(all_res, bind_df)
  }
  return(all_res)
}


## Parsing Data from inv_traj()

# Parse data from inv_traj() by four possible competition outcomes and into plottable format, determined by final population (minimum population 
#   at the end of the last three cycles)
# Input: return value from inv_traj()
# Output: a tibble with frequencies of each competition outcome under each set of parameters
# Note: the minimum population for deciding competition outcome, 1e-2, is equal to the default invading population

parse_invtraj = function(result) {
  result %<>% group_by(cycle_len) %>% filter(time %in% 
                                               c(max(cycle_no)*(cycle_len), (max(cycle_no)-1)*(cycle_len), (max(cycle_no)-2)*(cycle_len))) %>% 
    group_by(species, variation, rept, deltp, cycle_len, inv_spp) %>% summarize(max = max(population), min = min(population))
  output = tibble()
  for (i in unique(result$rept)) {
    for (j in unique(result$deltp)) {
      for (v in unique(result$variation)) {
        for (k in cycle_lens) {
          tb = result %>% filter(rept == i, deltp == j, variation == v, cycle_len == k)
          if (tb[tb$species == "N1" & tb$inv_spp == 1, ]$min > 1e-2 & tb[tb$species == "N2" & tb$inv_spp == 2, ]$min > 1e-2) { # coexistence
            res = 3
          } else if (tb[tb$species == "N1" & tb$inv_spp == 1, ]$min > 1e-2 & tb[tb$species == "N2" & tb$inv_spp == 2, ]$min < 1e-2) { # 1 wins
            res = 1
          } else if (tb[tb$species == "N1" & tb$inv_spp == 1, ]$min < 1e-2 & tb[tb$species == "N2" & tb$inv_spp == 2, ]$min > 1e-2) { # 2 wins
            res = 2
          } else if (tb[tb$species == "N1" & tb$inv_spp == 1, ]$min < 1e-2 & tb[tb$species == "N2" & tb$inv_spp == 2, ]$min < 1e-2) { # numeric PE
            res = 4
          }
          output = bind_rows(output, tibble(variation = v, deltp = j, cycle_len = k, outcome = res))
        }
      }
    }
  }
  return(output)
}


######### VITAL RATES CALCULATION ##########


# Average Interaction Strengths and Adult Mortality


## Interaction Strengths, Average within Cycle

# A function that calculates the average interaction strengths WITHIN each cycle (season), looking at the effect of cycle lengths
# Input: as previously defined
# Output: a data frame with columns: species, variation in arrival time, type of intreaction, number of cycles, average interaction strengths, 
#         cycle lengths

avg_int_incycle = function(n_cycles, para, init_pop, alpha, mean_time, var, cycle_lens, inv_pop, inv_spp) {
  out = tibble()
  for (i in 1:length(cycle_lens)) {
    int = cycles_invade(n_cycles*cycle_lens[i], para, init_pop, alpha, mean_time, var, cycle_lens[i], 
                        inv_pop, record = 1, inv_spp = inv_spp) # starting at resident equilibrium
    cycle_no = rep(1:n_cycles, each = cycle_lens[i])
    int$cycle_no = cycle_no
    int$cycle_len = cycle_lens[i]
    out = bind_rows(out, int)
  }
  out %<>% pivot_longer(cols = 1:4, names_to = "type", values_to = "int") %>% group_by(variation, type, cycle_len, cycle_no) %>% 
    summarize(avg_int = mean(int, na.rm = T))
  return(out)
}


## Interaction Strengths, Overall Average

# A function that calculates the average interaction strengths of the WHOLE simulation period, looking at the effect of cycle length
# Input: as previously defined
# Output: a data frame with columns: species, variation in arrival time, type of intreaction, average interaction strengths, cycle lengths

avg_int = function(n_cycles, para, init_pop, alpha, mean_time, var, cycle_lens, inv_spp = NA) {
  avg = tibble()
  for (cycle_len in cycle_lens) {
    if (!is.na(inv_spp)) {
      init_pop[[inv_spp]] = 0
      init_pop = stabledist(para, init_pop, alpha, mean_time, cycle_len)
      init_pop[1, inv_spp] = 0.01
    }
    int = cycles(n_cycles*cycle_len, para, init_pop, alpha, mean_time, var, cycle_len, record = 1) %>% 
      pivot_longer(cols = 1:4, names_to = "type", values_to = "int") %>% group_by(variation, type) %>% 
      summarize(avg_int = mean(int, na.rm = T))
    int$cycle_len = cycle_len
    avg = bind_rows(avg, int)
  }
  return(avg)
}


## Average Adult Population in Season

# Calculates average adult population within each season; can be used as a proxy of reproduction
# Input: as previously defined
# Output: a data frame with columns: species, variation, cycle number, average adult population in cycle, relative arrival time

avg_adult_incycle = function(n_cycles, para, init_pop, alpha, mean_time, var, cycle_lens, inv_spp, inv_pop = 0.01) {
  out = tibble()
  init_pop[[inv_spp]] = rep(0, s)
  stable_pop = stabledist(para, init_pop, alpha, mean_time, cycle_len)
  stable_pop[1, inv_spp][1] = inv_pop
  for (i in 1:length(cycle_lens)) {
    adult = cycles(n_cycles*cycle_lens[i], para, stable_pop, alpha, mean_time, var, cycle_lens[i]) %>% 
      pivot_longer(cols = N1:N2, names_to = "species", values_to = "population") %>% filter(stage == 6)
    cycle_no = rep(1:n_cycles, each = cycle_lens[i]*2)
    adult$cycle_no = cycle_no
    adult %<>% group_by(species, variation, cycle_no) %>% summarize(population = mean(population)) %>% ungroup()
    adult$cycle_len = cycle_lens[i]
    out = bind_rows(out, adult)
  }
  return(out)
}

