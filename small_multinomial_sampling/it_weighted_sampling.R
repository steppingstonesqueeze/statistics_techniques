# Weighted sampling using the simple inverse transform for multinomial distributions


it_weighted_sampling <- function(...) {
  args <- list(...)
  if (is.null(args$data)) {
    stop("Need data to sample")
  } else {
    data <- unlist(args$data) # get the data
  }
  
  len_data <- length(data)
  
  if (len_data < 2) {
    return (data) # nothing to sample at all
  } else {
    if (is.null(args$weights)) {
      weights <- rep(1.0/len_data, len_data) # each data point gets same probability of getting sampling - equiprobable
    } else {
      weights <- args$weights
      len_weights <- length(weights)
      # stop if weights length not same as dat length
      
      if (len_weights != len_data) {
        stop("Weights if specified must be same length as the data")
      } else {
        # normalise just in case to make ita proper PMF
        weights <- weights / sum(weights)
      }
    } # check if weights specified or not
  } # check on number of data points
  
  # cumulative weights - the inverse transform
  
  cum_weights <- cumsum(weights)
  
  if (is.null(args$N)) {
    # number of samples defaults to 100
    N <- 5
  } else {
    N <- max(args$N, 5) # take care of any mischief!
    
  }
  #cat(weights, "\n")
  
  output <- numeric(length = N)
  
  for (i in 1:N) {
    r <- runif(1, 0, 1)
    poss <- ifelse(r <= cum_weights, 1, 0)
    output[i] <- min(which(poss == 1))   
     }
  
  return (output)
}

my_output <- it_weighted_sampling(data = c(1,2,3,4), weights = c(0.7, 0.2, 0.08, 0.02), N = 1000)

table(my_output)

