sir_simulate_discrete <- function(t, I_0, R_0, beta, gamma) {
  S_0 <- 1 - I_0 - R_0
  
  S <- c(S_0, rep(0, t))
  I <- c(I_0, rep(0, t))
  R <- c(R_0, rep(0, t))
  
  for (i in 1:t) {
    S[i + 1] <- S[i] - beta * S[i] * I[i]
    I[i + 1] <- I[i] + beta * S[i] * I[i] - gamma * I[i]
    R[i + 1] <- R[i] + gamma * I[i]
  }
  out <- as.data.frame(cbind(S, I, R))
  return(out)
}


#from Bjornstad's "Epidemics"
sir_simulate_numerical <- function(t, I_0, R_0, beta, gamma) {
  sir_sim <- function(t, y, parms) {
    # Pull state variables from y vector
    S = y[1]
    I = y[2]
    R = y[3]
    # Pull parameter values from parms vector
    beta = parms["beta"]
    mu = parms["mu"]
    gamma = parms["gamma"]
    N = parms["N"]
    # Define equations
    dS = mu * (N - S) - beta * S * I / N
    dI = beta * S * I / N - (mu + gamma) * I
    dR = gamma * I - mu * R
    res = c(dS, dI, dR)
    # Return list of gradients
    list(res)
  }
  
  times = seq(0, t, by = 1 / 10)
  parms = c(
    mu = 0,
    N = 1,
    beta = beta,
    gamma = gamma
  )
  start = c(S = 1 - I_0 - R_0, I = I_0, R = R_0)
  out = deSolve::ode(
    y = start,
    times = times,
    func = sir_sim,
    parms = parms
  )
  out = as.data.frame(out)
  
  return(out)
}



# code from https://github.com/lilywang1988/eSIR, all credit goes to Song et al.  
song_plot_tv_transmission <- function() {
  NI_complete <- c( 41, 41, 41, 45, 62, 131, 200, 270, 375, 444, 549, 729,
                    1052, 1423, 2714, 3554, 4903, 5806, 7153, 9074, 11177,
                    13522,16678,19665,22112,24953,27100,29631,31728,33366)
  RI_complete <- c(1, 1, 7, 10, 14, 20, 25, 31, 34, 45, 55, 71, 94, 
                   121, 152, 213, 252, 345, 417, 561, 650, 811, 
                   1017, 1261, 1485, 1917, 2260, 2725,3284,3754)
  N <- 58.5e6
  R <- RI_complete / N
  Y <- NI_complete / N - R #Jan13->Feb 11
  change_time <- c("01/23/2020", "02/04/2020", "02/08/2020")
  pi0<- c(1.0, 0.9, 0.5, 0.1)
  res.step <- invisible(tvt.eSIR(Y, R, begin_str = "01/13/2020", T_fin = 200,
                                 pi0 = pi0, change_time = change_time, dic = T, 
                                 casename = "Hubei_step", save_files = T, save_mcmc=F,
                                 save_plot_data = F, M = 5e3, nburnin = 2e3))
  res.step$plot_infection
}

# code from https://github.com/lilywang1988/eSIR, all credit goes to Song et al.  
song_plot_tv_transmission <- function(with_intervention = T) {
  NI_complete <- c( 41, 41, 41, 45, 62, 131, 200, 270, 375, 444, 549, 729,
                    1052, 1423, 2714, 3554, 4903, 5806, 7153, 9074, 11177,
                    13522,16678,19665,22112,24953,27100,29631,31728,33366)
  RI_complete <- c(1, 1, 7, 10, 14, 20, 25, 31, 34, 45, 55, 71, 94, 
                   121, 152, 213, 252, 345, 417, 561, 650, 811, 
                   1017, 1261, 1485, 1917, 2260, 2725,3284,3754)
  N <- 58.5e6
  R <- RI_complete / N
  Y <- NI_complete / N - R #Jan13->Feb 11
  change_time <- c("01/23/2020", "02/04/2020", "02/08/2020")
  pi0<- c(1.0, 0.9, 0.5, 0.1)
  
  if (with_intervention) {
    res.step <- invisible(tvt.eSIR(Y, R, begin_str = "01/13/2020", T_fin = 200,
                                   pi0 = pi0, change_time = change_time, dic = T, 
                                   casename = "Hubei_step", save_files = T, save_mcmc=F,
                                   save_plot_data = F, M = 5e3, nburnin = 2e3))
    res.step$plot_infection 
  } else {
    res.nopi <- tvt.eSIR(Y, R, begin_str = "01/13/2020", death_in_R = 0.4, 
                         T_fin = 200, casename = "Hubei_nopi",
                         save_files = F,save_plot_data = F,
                         M=5e3,nburnin = 2e3)
    
    res.nopi$plot_infection
  }
}

# code from https://github.com/lilywang1988/eSIR, all credit goes to Song et al.  
song_plot_tv_quarantine <- function() {
  set.seed(20192020)
  NI_complete <- c( 41, 41, 41, 45, 62, 131, 200, 270, 375, 444,
                    549, 729, 1052, 1423, 2714, 3554, 4903, 5806,
                    7153, 9074, 11177, 13522, 16678, 19665, 22112,
                    24953, 27100, 29631, 31728, 33366)
  RI_complete <- c(1, 1, 7, 10, 14, 20, 25, 31, 34, 45, 55, 
                   71, 94, 121, 152, 213, 252, 345, 417, 561,
                   650, 811, 1017, 1261, 1485, 1917, 2260,
                   2725, 3284, 3754)
  N <- 58.5e6
  R <- RI_complete / N
  Y <- NI_complete / N - R #Jan13->Feb 11
  change_time <- c("01/23/2020", "02/04/2020", "02/08/2020")
  phi0 <- c(0.1, 0.4, 0.4)
  res.q <- qh.eSIR(Y, R, begin_str = "01/13/2020",
                   phi0 = phi0,change_time = change_time, 
                   casename = "Hubei_q", save_files = T,
                   save_mcmc = F,save_plot_data = F,
                   M = 5e3,nburnin = 2e3)
  res.q$plot_infection
}


# code from https://github.com/lilywang1988/eSIR, all credit goes to Song et al.  
song_plot_tv_antibodies <- function() {
  NI_complete <- c( 1, 2, 11, 23, 31, 76, 106, 142, 150, 220, 327,
                    421, 613, 615, 967, 1578, 3038, 5704, 8403,
                    11727, 15800, 20884, 25681, 30841, 37877,
                    44876, 52410, 59648, 66663, 75833, 83948,
                    92506, 102987, 113833, 123160, 131815,
                    139875, 151061, 161779, 172348, 181026,
                    189033, 195749, 203020, 214454, 223691,
                    230597, 237474, 243382, 248416, 253519,
                    258222, 263460, 271590, 282143, 288045,
                    291996, 295106, 299691, 304372, 308314)
  RI_complete <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 
                   20, 26, 38, 52, 200, 311, 468, 662,
                   893, 1144, 1547, 2144, 2931, 3900, 
                   5133, 6201, 7657, 9483, 11536, 
                   13915, 16276, 18781, 21110, 23424,
                   26469, 29784, 32899, 35785, 37730, 
                   39207, 40703, 42164, 43488, 44723, 
                   45887, 47473, 47686, 48769, 49572,
                   50221, 52161, 52917, 54115, 54613,
                   55473, 55816, 56809, 57265, 58525)
  N <- 8.399e6
  R <- RI_complete / N
  Y <- NI_complete / N - R 

  change_time <- c("04/29/2020")
  alpha0 <- c(0.2) # 20% of the susceptible population were found immunized
  res.antibody <- eSAIR(Y, R, begin_str = "03/01/2020",
                        alpha0 = alpha0, change_time = change_time,
                        casename = "New_York_antibody", save_files = F, 
                        save_mcmc = F, M=5e2, nburnin = 2e2)
  res.antibody$plot_infection
}