#Assuming the population size n = 1000 and the maximum household size is 5
n<-1000

#First, we need to create a vector "h" (length of h = n) to show which household each person belongs to.
set.seed(111) #use set.seed() function to make sure every time we run this code, we will get the same "h"
h <- rep(1:n, times = sample(1:5, n, replace = TRUE))[1:n]
#sample() is used to randomly get 1000 house sizes which ranges from 1 to 5 based on uniform distribution, i.e. the first element of "h" means the first household size...
#We then use the result of sample() as the "times" parameter of rep() to assign household ID to each person, i.e. repeat each household ID according to its sampled size.
#After repeating, the result may be more than 1000 since sizes are random between 1 and 5.We only keep the first 1000 people.

#Then we move to construct a function get.net(beta,h,nc) to show the possible social contact network among individuals in the population.
#get.net(beta,h,nc) takes the "sociability" parameter vector "beta" (length of "beta" should be n), the household belonging vector "h" and the average number of contacts per person "nc" (here we set nc to 15) as inputs.
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
        P[i,j]<-nc*beta[i]*beta[j]/(avg_beta^2*(n-1)) #note that, in some extreme cases, the calculated probability value may be larger than 1, but it won't affect the final result
      }
    }
  }
  
  for (k in 1:n){ 
    s<-runif(n+1-k, min=0, max=1) #For each row in the upper triangle of matrix P, we create a random number vector "s" of the corresponding length, with each entry drawn from the interval [0, 1] by using runif().
    s[P[k,c(k:n)]==0]<-0 #if one entry in the upper triangle of matrix P equals 0, we will also set the corresponding value in "s" to 0 to avoid comparing.
    #Only when the connection probability between the two people exceeds the corresponding value of "s" do we regard a link as present, i.e. the two individuals are each other's regular contacts, and we record this case as 1. Otherwise, we record as 0.
    P[k,c(k:n)][P[k,c(k:n)]>s]<-1
    P[k,c(k:n)][P[k,c(k:n)]<=s]<-0
  }
  
  #use t() to transpose the upper triangle of P and further derive the full matrix R
  R<-P+t(P)
  
  #looping through every row of R, we can find the indices of value 1 (i.e. the regular contacts of each person) by using which() and save the results to the empty list "l".
  for (m in 1:n){
    l[[m]]<-c(which(R[m,]==rep(1,n)))
  }
  
  #return list "l"---the ith element of which is a vector of the indices of the regular (non-household) contacts of person i; "integer(0)" means this person has no contact.
  return(l)
}

