
# load Rcpp library
library(Rcpp)


# load combine_min_var.cpp
sourceCpp("combine_min_var.cpp")


# generate example dataset (100 observations of integer 'age' and continuous 'response')
set.seed(17); age <- sample(0:40, 100, replace = T)
set.seed(17); response <- rnorm(100)
dat <- data.frame(age, response)


# tabulate frequency of each integer age
df <- as.data.frame(table(age = dat$age), responseName = "freq")


# use combine_min_var to combine ages into ordinal groupings (group sizes 3, 6, and 9)
out1 <- combine_min_var(vec = df$freq, vec_labels = df$age, n_groups = 3)
out2 <- combine_min_var(vec = df$freq, vec_labels = df$age, n_groups = 6)
out3 <- combine_min_var(vec = df$freq, vec_labels = df$age, n_groups = 9)


# fn to get group id given integer age and group membership list
GetGroupID <- function(x, memb_list) {
  srch <- lapply(memb_list, function(y) length(which(y == as.character(x))))
  grp <- which(unlist(srch) > 0)
  return(grp)
}


# apply GetGroupID to each observation within dat, for each of the three group sizes
dat$grp1 <- sapply(dat$age, GetGroupID, memb_list = out1$group_membership)
dat$grp2 <- sapply(dat$age, GetGroupID, memb_list = out2$group_membership)
dat$grp3 <- sapply(dat$age, GetGroupID, memb_list = out3$group_membership)


# plot response vs group, for each of the three group sizes
par(mfrow=c(1, 3))
boxplot(response ~ grp1, data = dat)
boxplot(response ~ grp2, data = dat)
boxplot(response ~ grp3, data = dat)

