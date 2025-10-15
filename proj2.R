# Background:
# This code implements a socially structured SEIR (Susceptible–Exposed–Infectious–Recovered) epidemic model that extends the basic 
# SEIR framework by incorporating household and network-based contact structures. Each individual belongs to a household and has a 
# set of regular social contacts, making infection more likely within these groups than through random mixing. The model simulates 
# epidemic dynamics over time based on parameters governing infection, exposure, and recovery rates. It enables comparison between 
# scenarios with and without social structure to study how household clustering and contact networks influence the spread and size 
# of epidemics.

# Code Structure Overview:
# 1. Initialize population: set total size n, assign individuals to households h with random sizes.
# 2. Build contact network: use get.net() to create non-household social links based on sociability beta.
# 3. Run SEIR simulation: use nseir() to model daily transitions (S→E→I→R) through household, contact, and random mixing infections.
# 4. Visualize results: use plot_nseir() to show sociability distribution and SEIR dynamics.
# 5. Compare scenarios and interpret: analyze how household and network structures affect epidemic spread versus random mixing.

# Contribution of each group member:
# Zixuan Qiu s2777279: Implemented SEIR simulation (nseir) (~33%)
# Fengyu Shen s2798951: Built household assignments & contact networks (~33%)
# Jianru Guo s2806788: Visualization and scenario comparison (~33%)

# GitHub Link: https://github.com/JianruGuo/Firstrepogroup.git



#Assuming the population size n = 1000 and the maximum household size is 5
n<-1000

#First, we need to create a vector "h" (length of h = n) to show which household each person belongs to.
#use set.seed() function to make sure every time we run this code, we will get the same "h"
set.seed(111) 
h <- rep(1:n, times = sample(1:5, n, replace = TRUE))[1:n]
#sample() is used to randomly get 1000 house sizes which ranges from 1 to 5 based on uniform distribution
#We then use the result of sample() as the "times" parameter of rep() to assign household ID to each person
#After repeating, the result may be more than 1000 since sizes are random between 1 and 5.We only keep the first 1000 people.

#Then we move to construct a function get.net(beta,h,nc) to show the possible social contact network among individuals in the population.
#get.net(beta,h,nc) takes the "sociability" parameter vector "beta", the household belonging vector "h" and the average number of 
#contacts per person "nc" (here we set nc to 15) as inputs.
#get.net(beta,h,nc) will return a list that indicates the indices of the regular (non-household) contacts of each person.
#the probability of having a link between person i and person j is calculated by "nc*beta[i]*beta[j]/((mean of beta)^2*(n-1))".
get.net <- function(beta,h,nc = 15){
  #redefine "n" inside the function to improve robustness
  n<-length(beta)
  
  #given the "sociability" parameter vector "beta", we can calculate its mean "avg_beta".
  avg_beta <- mean(beta)
  
  #create an nxn matrix P to save the probabilities of having a link between each pair of people
  #here we initialize nxn matrix P with 3(this value is arbitrary) instead of "NA" to make following check "P[i,j] != 0" proceed smoothly.
  P<-matrix(3,nrow = n, ncol = n)
  #create an nxn empty matrix and an empty list to save final results
  R<-matrix(,nrow = n, ncol = n)
  l<-list()
  
  P[outer(h,h,FUN = "==")]<-0
  #outer() is used to apply function "FUN = "=="" to every pair (h[i], h[j]) (i,j=1,2,..n) and will return an n×n logical matrix.
  #if h[i] = h[j] (i,j=1,2,..n), it indicates that person i either is the person j (the same person) or person i and j belong to the same household.
  #In this case, we will set the corresponding values to 0 because we only care about the non-household contacts of person i.
  
  #since the network we considered here is undirected, matrix P is theoretically symmetric.
  #to reduce the time of running this function, we only calculate the upper triangle of P and set all entries in the lower triangle to 0. 
  P[lower.tri(P)]<-0
  
  #loop through the upper triangle of matrix P
  for (i in 1:n){
    for (j in i:n){ 
      #if person i and j don't belong to the same household, we will calculate the probability of having a link between them based on the formula.
      if (P[i,j] != 0){
        P[i,j]<-nc*beta[i]*beta[j]/(avg_beta^2*(n-1))      
        #note that, in some extreme cases, the calculated probability value may be larger than 1, but it won't affect the final result
      }
    }
  }
  
  for (k in 1:n){ 
    #For each row in the upper triangle of matrix P, we create a random number vector "s" of the corresponding length
    s<-runif(n+1-k, min=0, max=1) 
    #if one entry in the upper triangle of matrix P equals 0, we will also set the corresponding value in "s" to 0 to avoid comparing.
    s[P[k,c(k:n)]==0]<-0 

    #Only when the connection probability between the two people exceeds the corresponding value of "s" do we regard a link as present, 
    #i.e. the two individuals are each other's regular contacts, and we record this case as 1. Otherwise, we record as 0.
    P[k,c(k:n)][P[k,c(k:n)]>s]<-1
    P[k,c(k:n)][P[k,c(k:n)]<=s]<-0
  }
  
  #use t() to transpose the upper triangle of P and further derive the full matrix R
  R<-P+t(P)
  
  #looping through every row of R, we can find the indices of value 1 (i.e. the regular contacts of each person) by using which()
  for (m in 1:n){
    l[[m]]<-c(which(R[m,]==1))
  }
  
  #return list "l"---the ith element of which is a vector of the indices of the regular (non-household) contacts of person i; 
  # "integer(0)" means this person has no contact.
  return(l)
}


# nseir function:
# Arguments:
#   beta   - sociability parameter for each person
#   h      - household ID for each person
#   alink  - list of contact indices for each person
#   alpha  - infection probabilities [household, contacts, random mixing]
#   gamma  - daily probability E->I (default 0.2)
#   delta  - daily probability I->R (default 0.4)
#   nc     - average contacts per person (default 15)
#   nt     - days to simulate (default 100)
#   pinf   - initial infection proportion (default 0.005)
#
# Returns:
#   List of daily counts: S, E, I, R, and time vector t
#
# Process:
#   Each day: record counts, transition states (E->I, I->R), then spread
#   infection via households, contacts, and random mixing

nseir <- function(beta,h,alink,alpha=c(.1,.01,.01),
                  delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
  n <- length(beta)
  beta_bar <- mean(beta)
  
  # Extract infection probabilities
  alpha_h <- alpha[1]  # household
  alpha_c <- alpha[2]  # contacts
  alpha_r <- alpha[3]  # random mixing
  
  # Initialize states: 1=S, 2=E, 3=I, 4=R
  # All people are susceptible(S=1) at the beginning
  state <- rep(1, n)
  
  # Randomly select pinf proportion to start infected (state change to 3)
  # Use max() to ensure at least 1 infected
  n_initial <- max(1, round(n * pinf))
  initial_infected <- sample(n, n_initial)
  state[initial_infected] <- 3
  
  # Storage for daily counts of all states
  S <- numeric(nt)
  E <- numeric(nt)
  I <- numeric(nt)
  R <- numeric(nt)
  
  day <- 1
  while (day <= nt) {
    # main loop: Record current state counts
    S[day] <- sum(state == 1)
    E[day] <- sum(state == 2)
    I[day] <- sum(state == 3)
    R[day] <- sum(state == 4)
    
    u <- runif(n)
    
    # state transitions using logical indexing
    state[state == 2 & u < gamma] <- 3   # E -> I
    state[state == 3 & u < delta] <- 4   # I -> R
    
    infected_idx <- which(state == 3)
    i_idx <- 1
    
    # Process infections (S -> E)
    # Need to loop through infected individuals
    while (i_idx <= length(infected_idx)) {
      i <- infected_idx[i_idx]
      
      # 1. Household infections
      # find the infectious people who live in the same household as person i 
      household_susceptible <- which(h == h[i] & state == 1 & (1:n) != i)
      if (length(household_susceptible) > 0) {
        u_h <- runif(length(household_susceptible))
        state[household_susceptible[u_h < alpha_h]] <- 2
      }
      
      # 2. Network contact infections
      # find the contacts of person i, who are also infectious
      contacts_i <- alink[[i]]
      
      if (length(contacts_i) > 0) {
        susceptible_contacts <- contacts_i[state[contacts_i] == 1]
        if (length(susceptible_contacts) > 0) {
          u_c_sub <- runif(length(susceptible_contacts))
          state[susceptible_contacts[u_c_sub < alpha_c]] <- 2
        }
      }
      
      # 3. Random mixing infections
      susceptible_idx <- which(state == 1)
      if (length(susceptible_idx) > 0) {
        p_infect <- (alpha_r * nc * beta[i] * beta[susceptible_idx]) / 
          (beta_bar^2 * (n - 1))
        u_r <- runif(length(susceptible_idx))
        state[susceptible_idx[u_r < p_infect]] <- 2
      }
      i_idx <- i_idx + 1
    }
    day <- day + 1
  }
  return(list(S = S, E = E, I = I, R = R, t = 1:nt))
}

#Visualize the SEIR epidemic simulation :
#show the distribution of beta: indicate how socially active different people are
#show the time evolution of the epidemic: indicate how the numbers of each status change over time
#compare four versions of SEIR model 

#generate one random number between 0 and 1 for each person
beta <- runif(n, 0, 1)
# generate the list that gives the indices of ith regular contacts
alink <- get.net(beta, h, nc = 15)
# generate a list containing time series vectors for S, E, I, R, t 
epi <- nseir(beta, h, alink)

# Full model
epi_full <- nseir(beta, h, alink, alpha = c(0.1, 0.01, 0.01))

# only Random-mixing
epi_random <- nseir(beta, h, alink, alpha = c(0, 0, 0.04))

# full model with Constant beta
beta_const <- rep(mean(beta), n)
epi_const <- nseir(beta_const, h, alink, alpha = c(0.1, 0.01, 0.01))

# Constant bet + Random mixing
epi_const_random <- nseir(beta_const, h, alink, alpha = c(0, 0, 0.04))

#add the sociability vector beta into the result list to draw the histogram 
epi$beta <- beta 

epi_full$beta <- beta
epi_random$beta <- beta
epi_const$beta <- beta_const
epi_const_random$beta <- beta_const

plot_nseir<- function(epi, main = "SEIR epidemic dynamics") {
  #first set up the plotting window
  #par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  
  #left panel is the histogram of the sociability distribution beta
  #in this plot, flat histogram represents relatively mixed sociability
  #narrow histogram indicates relatively similar sociability
  hist(epi$beta,
       breaks = 20, col = "skyblue",#number of histogram bins and the color
       xlab = expression(beta),#x-axis
       main = "Distribution of sociability (β)",#title
       border = "white")#remove bar borders
  box()#draw box around the plot
  
  #right panel: SEIR time series dynamics
  #This plot shows how the epidemic evolves in the population( with the number of S, E, I, R over the time period)
  ymax <- max(epi$S, epi$E, epi$I, epi$R)
  #first plot the curve for S
  plot(epi$t, epi$S, type = "l", lwd = 2, col = "black",
       ylim = c(0, ymax), xlab = "Day", ylab = "Number of individuals",
       main = main)
  #add E, I, R into the plot with different colors
  lines(epi$t, epi$E, lwd = 2, col = "blue")
  lines(epi$t, epi$I, lwd = 2, col = "red")
  lines(epi$t, epi$R, lwd = 2, col = "darkgreen")
  
  legend("right",
         legend = c("S (Susceptible)", "E (Exposed)",
                    "I (Infectious)", "R (Recovered)"),
         col = c("black", "blue", "red", "darkgreen"),
         lwd = 2, bty = "n")
}
plot_nseir(epi)
#plot 4 scenarios next to each other
par(mfcol = c(2, 4), mar = c(4, 4, 2, 1))  
plot_nseir(epi_full, main = "1. Full model")
plot_nseir(epi_random, main = "2. Random mixing only")
plot_nseir(epi_const, main = "3. Constant β, structured")
plot_nseir(epi_const_random, main = "4. Constant β + random mixing")

#Conclusion:
#Comparing model 1 with model 2, model 3 with model 4, the inclusion of household and network structure slow down the spread of infection. 
#Households and network structure separate people in groups that limiting the opportunity for the epidemic to spread.
#For random mixing, anyone in the population can be contacted by an infectious person, so the transmission is faster, which increase the epidemic peak.
#So considering househode and network structure would make predicted epidemic grows slower with smaller epidemic peak than considering only random mixing.




