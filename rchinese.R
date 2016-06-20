rchinese = function(size, a) {
  stopifnot(length(a) == 1)
  samples = c()
  new.category = 1
  for (i in 1:size) {
    denom = i - 1 + a
    p.new = a / denom
    ci =
      if (runif(1) < p.new) {
        new.category = new.category + 1;
        new.category - 1
      }
    else {
      counts = unname(table(factor(samples, levels = 1:(i-1))))
      pmf = counts / sum(counts)
      sample(1:(i-1), 1, prob= pmf, replace = TRUE)
    }
    samples = c(samples, ci)
  }
  samples
}

