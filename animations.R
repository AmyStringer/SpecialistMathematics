library(tidyverse)
library(gganimate)
 
#### Playing with GGANNIMATE  
## plan is to make some animations to add to my app 
## i tried making them in app and displaying them using renderImage and
## it worked in rstudio on my local machine 
## but it didn't work once deployed 
## so here we are making gifs 


#### The Idea 
## simulate data from a few different distributions
## for each trial we will do 200 simulations 
## varying the number of samples in each trial 
## between 25, 100 and 250 
## goal: plot sample means and show that more samples == more normal 


##### FOR NORMAL 

# set up the normal parameters, these will remain constant for all samples 
mu <- 169
sigma <- 10

# number of simulations will also remain constant
n_sims <- 200

## number of samples to vary with each trial 
samp <- 25 

# set up the mesh 
sim_data <- expand.grid(sample = 1:samp, 
                        sim = 1:n_sims)

# simulate the data from a normal distribution 
sim_data <- mutate(sim_data, 
                   value = rnorm(n = samp*n_sims, 
                                 mean = mu, 
                                 sd = sigma))

# calculate the mean for each of the 200 simulations 
means1 <- sim_data %>% 
  group_by(sim) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(n = 25)

# next trial - samples = 100 
samp <- 100

# set up the grid
sim_data <- expand.grid(sample = 1:samp, 
                        sim = 1:n_sims)

# simulate from the normal distribution 
sim_data <- mutate(sim_data, 
                   value = rnorm(n = samp*n_sims, 
                                 mean = mu, 
                                 sd = sigma))

# compute the means for trial 2
means2 <- sim_data %>% 
  group_by(sim) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(n = 100)

# trial 3 - samples = 250 
samp <- 250

# set up the grid 
sim_data <- expand.grid(sample = 1:samp, 
                        sim = 1:n_sims)

# simulate samples from the normal distribution 
sim_data <- mutate(sim_data, 
                   value = rnorm(n = samp*n_sims, 
                                 mean = mu, 
                                 sd = sigma))

# compute the means for trial 3 
means3 <- sim_data %>% 
  group_by(sim) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(n = 250)

# join all the means datasets together into one mega data set to be used for plotting
bigboy <- bind_rows(means1, means2, means3)

p <- ggplot(bigboy, aes(x = mean)) + 
  geom_histogram(color = "black",
                 # hey its the theme colour 
                 fill =  rgb(24, 188, 156,
                             maxColorValue = 255)) +
  theme_bw() + 
  transition_states(as.factor(n),
                    transition_length = 1,
                    state_length = 1) +
  view_follow(fixed_x = TRUE) + 
  labs(title = "Sampling from the Normal Distribution", 
       subtitle = "Moving to n = {next_state}", 
       x = "Sample Mean")

p

anim_save("NormalOutfile.gif", animate(p)) # New


###### FOR BINOMIAL 

# number of simulations to remain constant
n_sims <- 200

# size and probability also to remain constant for all trials 
size <- 20
prob <- 0.2

# trial 1 - samples = 25 
samp <- 25 

# set up grid 
sim_data  <- expand.grid(sample = 1:samp, 
                         sim = 1:n_sims)

# simulate data from the binomial distribution 
sim_data <- mutate(sim_data, 
                   value = rbinom(n = samp*n_sims, 
                                  size = size, 
                                  prob = prob))

# compute the means for trial 1
means1 <- sim_data %>% 
  group_by(sim) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(n = 25)

# trial 2 - samples = 100
samp <- 100

# set up grid 
sim_data  <- expand.grid(sample = 1:samp, 
                         sim = 1:n_sims)

# simulate data from the binomial distribution 
sim_data <- mutate(sim_data, 
                   value = rbinom(n = samp*n_sims, 
                                  size = size, 
                                  prob = prob))

# compute the means for trial 2
means2 <- sim_data %>% 
  group_by(sim) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(n = 100)

# trial 2 - samples = 250 
samp <- 250 

# set up the grid 
sim_data  <- expand.grid(sample = 1:samp, 
                         sim = 1:n_sims)

# simulate data from the binomial distribution 
sim_data <- mutate(sim_data, 
                   value = rbinom(n = samp*n_sims, 
                                  size = size, 
                                  prob = prob))
# compute the means for trial 3 
means3 <- sim_data %>% 
  group_by(sim) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(n = 250)

# join all the means datasets together into a mega dataset 
bigboy <- bind_rows(means1, means2, means3)

p <- ggplot(bigboy, aes(x = mean)) + 
  geom_histogram(color = "black",
                 fill =  rgb(24, 188, 156,
                             maxColorValue = 255)) +
  theme_bw() + 
  transition_states(as.factor(n),
                    transition_length = 1,
                    state_length = 1) +
  view_follow(fixed_x = TRUE) + 
  labs(title = "Sampling from the Binomial Distribution", 
       subtitle = "Moving to n = {next_state}", 
       x = "Sample Mean")

p

anim_save("BinomialOutfile.gif", animate(p)) # New

##### UNIFORM DISTRIBUTION 
# number of simulation to remain constant for all three trials 
n_sims <- 200

# uniform dist parameters also to remain constant for all three trials 
min <- 0
max <- 4 

# trials 1, samples = 25 
samp <- 25 

# set up grid 
sim_data  <- expand.grid(sample = 1:samp, 
                         sim = 1:n_sims)

# simulate from the uniform distribution 
sim_data <- mutate(sim_data, 
                   value = runif(n = samp*n_sims, 
                                 min = min, 
                                 max = max))

# compute the means for the first trial 
means1 <- sim_data %>% 
  group_by(sim) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(n = 25)

# trial 2 samples = 100 
samp <- 100 

# set up the grid 
sim_data  <- expand.grid(sample = 1:samp, 
                         sim = 1:n_sims)

# simulate trial 2 data from the uniform distribution 
sim_data <- mutate(sim_data, 
                   value = runif(n = samp*n_sims, 
                                 min = min, 
                                 max = max))

# compute the trial 2 means 
means2 <- sim_data %>% 
  group_by(sim) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(n = 100)

# trial 3 samples = 250 
samp <- 250 

# set up the grid
sim_data  <- expand.grid(sample = 1:samp, 
                         sim = 1:n_sims)

# simulate trial 3 data from the uniform distriubtion 
sim_data <- mutate(sim_data, 
                   value = runif(n = samp*n_sims, 
                                 min = min, 
                                 max = max))

# compute the trial 3 means
means3 <- sim_data %>% 
  group_by(sim) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(n = 250)

# join them all together into the mega dataset
bigboy <- bind_rows(means1, means2, means3)

p <- ggplot(bigboy, aes(x = mean)) + 
  geom_histogram(color = "black",
                 fill =  rgb(24, 188, 156,
                             maxColorValue = 255)) +
  theme_bw() + 
  transition_states(as.factor(n),
                    transition_length = 1,
                    state_length = 1) +
  view_follow(fixed_x = TRUE) + 
  labs(title = "Sampling from the Uniform Distribution", 
       subtitle = "Moving to n = {next_state}", 
       x = "Sample Mean")

p

anim_save("UniformOutfile.gif", animate(p)) # New

