# alias_sampler.R
args <- commandArgs(trailingOnly = TRUE)
k <- ifelse(length(args) >= 1, as.numeric(args[1]), 10000)
n <- ifelse(length(args) >= 2, as.numeric(args[2]), 100000)

# Generate a probability vector using Dirichlet
dirichlet <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha)
  return(x / sum(x))
}
probs <- dirichlet(rep(1, k))

create_alias_table <- function(probs) {
  n <- length(probs)
  scaled_probs <- probs * n
  alias <- integer(n)
  prob <- numeric(n)
  small <- which(scaled_probs < 1)
  large <- which(scaled_probs >= 1)

  while (length(small) > 0 && length(large) > 0) {
    s <- tail(small, 1); small <- head(small, -1)
    l <- tail(large, 1); large <- head(large, -1)

    prob[s] <- scaled_probs[s]
    alias[s] <- l

    scaled_probs[l] <- scaled_probs[l] - (1 - prob[s])
    if (scaled_probs[l] < 1) {
      small <- c(small, l)
    } else {
      large <- c(large, l)
    }
  }

  for (i in c(small, large)) {
    prob[i] <- 1
    alias[i] <- i
  }

  return(list(prob = prob, alias = alias))
}

alias_sample <- function(alias, prob, n) {
  samples <- integer(n)
  K <- length(prob)
  for (i in 1:n) {
    idx <- sample.int(K, 1)
    if (runif(1) < prob[idx]) {
      samples[i] <- idx
    } else {
      samples[i] <- alias[idx] + 1
    }
  }
  return(samples)
}

start <- Sys.time()
table <- create_alias_table(probs)
build_time <- Sys.time() - start

start <- Sys.time()
samples <- alias_sample(table$alias, table$prob, n)
sample_time <- Sys.time() - start

cat(sprintf("Alias table built in %.3f secondsn", as.numeric(build_time, units="secs")))
cat(sprintf("Sampled %d values in %.3f secondsn", n, as.numeric(sample_time, units="secs")))
cat("First 10 samples: ", paste(head(samples, 10), collapse = ", "), "n")
