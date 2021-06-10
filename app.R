#### Application for Specialist Mathematics #### 
## Created by Amy Stringer 
## There is a README for this file, please 
## go and read this before attempting to edit 
## anything within this app.R file. 

## When editing, please make sure the indenting is consistent 


library(shiny)
library(tidyverse)     # for all the things 
library(scales)        # This is usually for percent, but i don't recall using this 
library(shinythemes)   # themes 
library(plotly)        # interactive plot interface 
library(glue)
library(gifski)

# Define UI for application that draws a histogram
ui <- fluidPage(theme = shinythemes::shinytheme("flatly"),
    
    # mathJax allows for maths text 
    withMathJax(), 
    # this little bit of code allows for inline math text using the usual latex $ $ 
    tags$div(HTML("<script type='text/x-mathjax-config' >
            MathJax.Hub.Config({
            tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
            });
            </script >
            ")),
    # this little bit of code prevents error messages from appearing in the user interface of the app 
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"
    ),
    navbarPage(title = "Specialist Mathematics Unit 4", 
                
               #### topic 3 #### 
               tabPanel(title = "Statistical Inference",  
                        
                        tabsetPanel(
                            tabPanel("Welcome", 
                                     # h1 for main headings 
                                     h1("Welcome to Unit 4 of Specialist Mathematics"),
                                     # h5 for paragraph text 
                                     h5("This app exists as a tool for students and teachers following the senior secondary Australian Curriculum for Mathematics. It was designed to assist in the explanation and understanding of the topics covered in Unit 4 of the Specialist Mathematics course."),
                                     h5("This unit consists of three main topics:"),
                                     # tags$ol allows you to create an ordered list 
                                     tags$ol(
                                         # tags$li creates a list item within the ordered list 
                                         tags$li("Integration and applications of integration"),
                                         tags$li("Rates if change and differential equations, and"),
                                         tags$li("Statistical inference")
                                     ),
                                     h5("Due to the nature of the techniques taught in the first two topics, this application will focus on topic 3: statistical inference. "),
                                     # h3 for sub headings 
                                     h3("Statistical Inference"), 
                                     h5("This topic consists of three tabs."),
                                     h5("Tab 1 is designed to illustrate how sampling from a distribution works, and some features of random sampling. "), 
                                     h5("Tab 2 illustrates some features of repeated random sampling for varying sample sizes, and, finally, "), 
                                     h5("Tab 3 gives a demonstration of confidence intervals for sample means. "),
                                     h3("Glossary"),
                                     h5("Many of the terms used in this unit will carry over from mathematical methods, but some definitions have been repeated here for convenience. "),
                                     # tags$b allows for bold text 
                                     h5(tags$b("Sample")," - a sample is a subset of the population which we select in order to make inferences about a population. When we talk about sampling from a distribution this refers to a random selection of values that are distributed according to the chosen distribution. For example, selecting 50 values from a normal distribution means that when we plot these in a histogram the shape will appear roughly normal $($more so if we select more values$)$. "),
                                     h5(tags$b("Mean"), " - a measure of average behaviour $($other methods of measuring average behaviour include median and mode$)$"),
                                     h5(tags$b("Sample mean"), " - denoted $\\bar{x}$, this represents the mean value of a given sample $($or a given subset of the population$)$. This is different to the population mean, $\\mu$, which refers to the mean of all possible data points. "),
                                     h5(tags$b("Confidence Interval"), " - An interval estimate for the population mean $\\mu$ is called a confidence interval for $\\mu$."),
                                     h5(tags$b("Confidence Level"), " - This value refers to exactly how confident we are that our interval contains the true population mean. In general, you will mostly be dealing 95% confidence as this is the most used in practise. "),
                                     h5(""),
                                     # horizontal rule before adding footnotes 
                                     tags$hr(), 
                                     # h6 for footnote text 
                                     h6("Created by Amy Stringer as part of an ARC Centre for Excellence in Mathematics and Statistics Initiative. For any enquiries regarding this app, or bug fixes, please email ", tags$a(href = "mailto:virtualreefdiver@qut.edu.au", "Virtual Reef Diver")), 
                                     h6("Some activities seen here have been adapted from ", tags$a(href = "https://seeing-theory.brown.edu/frequentist-inference/index.html", "Seeing theory"))), 
                            tabPanel(title = "Sampling Means", 
                                     
                                     sidebarLayout(
                                         sidebarPanel(
                                             # dropdown box for selection of distribution 
                                             # chose a range of different distributions to really get the point across 
                                             selectInput(inputId = "T3Dist1", label = "Select Distribution", 
                                                         choices = c("Normal", "Uniform", 
                                                                     "Poisson", "Binomial", 
                                                                     "Exponential", "Gamma"))
                                         ), 
                                         mainPanel(
                                             h3("Sampling Means"),
                                             h5("To the left you will see an option to select a distribution. Once you've done this the app will take 9 random samples from this distribution and plot the results below. You should notice that each looks slightly different, but the general shape of the data is similar. The purpose of this exercise is to show that despite the fact the samples were taken from the same distribution, the features of the sample are not the exactly the same."),
                                             # there is not a lot of interactivity in this tab so I've used plotlyOutput so that user can interact with the plot 
                                             plotlyOutput(outputId = "T3SimsPlot1"),
                                             h5("One feature that we are particularly intereted in is the mean. Below here you can see a table containing both the sample number and the mean value for that sample. The average behaviour is not consistent across all 9 samples, regardless of the type of distribution you select."),
                                             # datatable output looked much better here than the standard table 
                                             dataTableOutput(outputId = "T3means")
                                         )
                                     )
                                     
                                     ), 
                            tabPanel(title = "Repeated Sampling", 
                               # tab headings use h3 
                               h3("Distribution of Sample Means"),  
                               # paragraph text 
                               h5("Here we have a little demonstration illustrating some behaviour of the sample mean as a random variable."),
                               h5("Below are the results of a series of experiments whereby we simulated data from a variety of distributions, namely, the normal distribution, the binomial distribution and the uniform distribution. "),
                               h5("In each experiment, we generated three sets of data, each set consisting of 200 simulations. In the first set each simulation contained only 25 samples $(n = 25)$ of data from the chosen distribution, the second contained 100 samples $(n = 100)$ and the third, 250 samples $(n = 250)$. "),
                               h5("You will see two rows of plots below - the first row shows the shape of the distribution from which the samples were taken, and the second row shows the corresponding plot of the distribution of sample means for each set of simulated data. "),
                               h5("What you should notice while watching the gifs is that increasing the number of samples taken in the simulations results in the sample mean distribution becoming more and more normally distributed. This is a direct consequence of the", tags$b("Central Limit Theorem.")),
                               fluidRow(
                                   column(width = 4, 
                                          # these don't need to be interactive, it's just a line 
                                          plotOutput(outputId = "NormalPDF", height = "300px")
                                          ), 
                                   column(width = 4, 
                                          # as above 
                                          plotOutput(outputId = "BinomPMF", height = "300px")
                                          ), 
                                   column(width = 4, 
                                          # as above 
                                          plotOutput(outputId = "UnifPDF", height = "300px")
                                          )
                               ), 
                               fluidRow(
                                   column(width = 4, 
                                          # okay so these were really annoying 
                                          # i want to use plotly and ggannimate in the app and add them using renderImage 
                                          # it worked on my computer, but once deployed there was no image 
                                          # so I have to find a work around. 
                                          # animations contained here were created using the animations.R file, then saved and imported 
                                          img(src="NormalOutfile.gif", align = "left",height='400px',width="100%")
                                          ), 
                                   column(width = 4, 
                                          img(src = "BinomialOutfile.gif", align = "left",height='400px',width="100%")
                                          ), 
                                   column(width = 4, 
                                          img(src = "UniformOutfile.gif", align = "left",height='400px',width="100%")
                                          )
                               )
                                     
                            ),
                            tabPanel(title = "Confidence Intervals for Means", 
                                fluidRow(
                                    column(width = 12, 
                                           # for specialist i thought students might be more comfortable with the confidence level idea 
                                            sliderInput(inputId = "confLevel", 
                                                        label = "Slide me to select your confidence level!",
                                                        # default to 95% 
                                                        min = 0, max = 100, value = 95, width = '100%')
                                           )
                                ), 
                                sidebarLayout(
                                    sidebarPanel(
                                        actionButton(inputId = "Add10", label = "Add 10 Simulations"), 
                                        actionButton(inputId = "Add100", label = "Add 100 Simulations"),
                                        # use $$ for brackets because if i don't they disappear 
                                        h5("Confidence Interval $($CI$)$ Coverage"),
                                        # outputs in the sidebar! 
                                        # i know, i know, pie charts are awful, but this one will be good! 
                                        plotOutput(outputId = "ConfIntPie", height = "250px"), 
                                        # this will explain the proportions shown in the pie chart because text in pie charts is tricky 
                                        textOutput(outputId = "ratio"),
                                        h5("Proportion of CIs Containing the True Mean"),
                                        # this plots the proportion of "inside" against the simulation total
                                        plotlyOutput(outputId = "PlotProp", height = "250px")
                                    ), 
                                    mainPanel(
                                        # tab heading 
                                        h3("Plot of Sample Means and Associated Confidence Intervals"),
                                        # paragraph text 
                                        h5("In this experiment, you use the above slider to set the confidence level you wish to investigate. You then can add 10 or 100 simulations at a time, each simulation drawing 100 samples from a normal distribution with a mean of 0 and a standard deviation of 1."),
                                        h5("The app will then compute the mean value for each simulation and the associated confidence interval at the confidence level specified by the above slider. "),
                                        h5("Points and intervals are then coloured according to whether or not the interval contains the true mean of 0, indicated by the red dashed line."),
                                        # interactive so that the user can see the interval values 
                                        plotlyOutput(outputId = "PlotIntervals", height = "500px")
                                    )
                                )
                                )
                            )
                        )
                        
               )
               
               
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    #### Topic 3 - tab 1  #### 
    # set up a reactive dataframe to store the simulations for tab 1 
    # initialise with empty data frame 
    T3SimsData1 <- reactiveValues(data = data.frame())
    
    # we want something to get added to the reactive dataframe once the distribution has been selected 
    observeEvent(input$T3Dist1, {
        if (input$T3Dist1 == "Normal"){
            # mean and standard deviation chosen arbitrarily 
            mu <- 169
            sigma <- 10
            
            # I chosen 9 simulations so that when i facet by sim in the plotting 
            # I have a nice square 
            n_sims <- 9
            # choose a high number of samples so that the shape is clear 
            samp <- 250
            
            # set up a grid of samples by sims 
            sim_data <- expand.grid(sample = 1:samp, 
                                    sim = 1:n_sims)
            
            # randomly choose some normal data and update the reactive dataset 
            T3SimsData1$data <- mutate(sim_data, 
                               value = rnorm(n = samp*n_sims, 
                                             mean = mu, 
                                             sd = sigma))
        } else if (input$T3Dist1 == "Uniform"){
            # these are the same across all distributions 
            n_sims <- 9
            samp <- 250 
            
            # parameters chosen arbitrarily 
            min <- 0
            max <- 4 
            
            # set up the grid 
            sim_data  <- expand.grid(sample = 1:samp, 
                                     sim = 1:n_sims)
            
            
            # update the reactive dataset with uniform values 
            T3SimsData1$data <- mutate(sim_data, 
                               value = runif(n = samp*n_sims, 
                                              min = min, 
                                              max = max))
            
        } else if (input$T3Dist1 == "Poisson"){
            # as above 
            n_sims <- 9
            samp <- 250 
            
            # arbitrary parameters 
            lambda <- 4 
            
            # set up the grid 
            sim_data  <- expand.grid(sample = 1:samp, 
                                     sim = 1:n_sims)
            
            # update reactive dataframe with poisson data 
            T3SimsData1$data <- mutate(sim_data, 
                               value = rpois(n = samp*n_sims, 
                                             lambda = lambda))
            
        } else if (input$T3Dist1 == "Binomial"){
            # as above 
            n_sims <- 9
            samp <- 250
            
            # arbitrary parms 
            size <- 20
            prob <- 0.2
            
            # grid 
            sim_data  <- expand.grid(sample = 1:samp, 
                                     sim = 1:n_sims)
            
            # update reactive dataframe with binomial data 
            T3SimsData1$data <- mutate(sim_data, 
                               value = rbinom(n = samp*n_sims, 
                                              size = size, 
                                              prob = prob))
        } else if (input$T3Dist1 == "Exponential"){
            # as above 
            n_sims <- 9
            samp <- 250
            
            # arbitrary rate 
            rate <- 1.5 
            
            # grid 
            sim_data  <- expand.grid(sample = 1:samp, 
                                     sim = 1:n_sims) 
            
            # update reactive with exponential data 
            T3SimsData1$data <- mutate(sim_data, 
                               value = rexp(n = samp*n_sims, 
                                            rate = rate))
            
        } else if (input$T3Dist1 == "Gamma"){
            # as above 
            n_sims <- 9
            samp <- 250
            
            # arbitrary parameters 
            shape <- 2.0
            scale <- 2.0
            
            # grid 
            sim_data  <- expand.grid(sample = 1:samp, 
                                     sim = 1:n_sims) 
            
            # update reactive dataframe to contain gamma values 
            T3SimsData1$data <- mutate(sim_data, 
                               value = rgamma(n = samp*n_sims, 
                                              shape = shape, 
                                              scale = scale))
        }
    })
    
    output$T3SimsPlot1 <- renderPlotly({
        P1 <- ggplot(data = T3SimsData1$data, aes(x = value)) + 
            geom_histogram(color = "black", 
                           # this is the theme colour. this was hard to find, don't lose it
                           fill = rgb(24, 188, 156, 
                                      maxColorValue = 255)) + 
            facet_wrap(~sim) + 
            theme_bw()
        ggplotly(P1)
    })
    
    output$T3means <- renderDataTable({
        # want the data table to show the means of all 9 simulations 
        sim_data <-  T3SimsData1$data %>% 
            group_by(sim) %>% 
            summarise(mean = round(mean(value), 2)) 
        # need to print the table to render the table 
        sim_data
    })
    
    #### Topic 3 - tab 2 #### 
    
    # tab 2 just had the probability plots for three distributions and the gifs which have been 
    # added in manually 
    
    output$NormalPDF <- renderPlot({
        # set up distribution parameters 
        # chosen arbitrarily 
        mu <- 169
        sigma <- 10
        
        # using stat_function, so only need beginning and end points 
        dat <- data.frame(x = c(135, 205))
        
        ggplot(data = dat, aes(x = x)) + 
            stat_function(fun = dnorm, 
                          args = list(mean = mu, 
                                      sd = sigma)) + 
            theme_bw() + 
            labs(title = "PDF of the Normal Distribution")
        
    }, height = 300)
    
    output$BinomPMF <- renderPlot({
        # set up distribution parameters 
        # chosen arbitrarily 
        size <- 15
        prob <- 0.7
        
        # create a sequence for the x values 
        x1 <- 0:size
        # add the binomial values 
        dat <- data.frame(x = x1, 
                          y = dbinom(
                              x1, 
                              size = size, 
                              prob = prob
                          ))
        # geom_bar for discrete distribution 
        ggplot(data = dat, aes(x = x, y = y)) +
            geom_bar(stat = "identity", 
                     fill = rgb(24, 188, 156,
                                maxColorValue = 255)) +
            theme_bw() +
            labs(title = "PMF of the Binomial Distribution")
        
    }, height = 300)
    
    output$UnifPDF <- renderPlot({
        # create a dataframe to plot the distribution with 
        dat <- data.frame(x = seq(0, 4, by = 0.1))
        # make sure that dataframe contains values sampled from a uniform distribution 
        dat$y <- dunif(dat$x, min = 0, max = 4)
        
        # geom_line for continuous distribution 
        ggplot(data = dat, aes(x = x, y = y)) + 
            geom_line() + 
            theme_bw() +
            labs(title = "PDF of the Uniform Distribution")
        
    }, height = 300)
    
    
    
    #### TOPIC 3 - TAB 3 ####
    
    # set up a reactive data frame for tab 3 
    # initialise it as empty 
    T3simsdata2 <- reactiveValues(data = data.frame())
    
    # we also want to create a dataframe to contain the pie chart info 
    # this must be reactive as it depends on the other reactive dataframe 
    T3PieDat <- reactiveValues(data = data.frame())
    
    observeEvent(input$confLevel, {
        # we want the sims data to reset once the confidence level is changed 
        T3simsdata2$data <- data.frame()
    })
    
    observeEvent(input$Add10, {
        
        # set the confidence level 
        level <- input$confLevel/100
        
        # add 10 implies we want to add another 10 simulations to our dataset 
        n_sims <- 10
        # we will stick to using 100 samples 
        samp <- 100
        
        # sampling from the standard normal distribution
        mu <- 0
        sigma <- 1
        
        # setting up the grid 
        sim_data <- expand.grid(sample = 1:samp, 
                                sim = (nrow(T3simsdata2$data) + 1):(nrow(T3simsdata2$data) + n_sims))
        # extract the samples 
        sim_data <- mutate(sim_data, 
                           height = rnorm(n = samp*n_sims, 
                                          mean = mu, 
                                          sd = sigma))
        
        # compute all the bits 
        sim_means <- sim_data %>% 
            group_by(sim) %>% 
            # we need the sample mean, degree of freedom and sd for the conf int calc
            summarise(xbar = mean(height), 
                      n = n(), 
                      s = sd(height), 
                      alpha = 1-level) %>% 
            # also need to add in the t values (using the actual formula for this one)
            mutate(t.low = qt(p = alpha/2, df = n-1),
                   t.high = qt(p = 1-alpha/2, df = n-1),
                   # then sub all the things into the formula for conf.high and conf.low 
                   conf.low = xbar + t.low*(s/sqrt(n)), 
                   conf.high = xbar + t.high*(s/sqrt(n)))
        
        # add a variable the tells us whether the true mean falls within the interval 
        sim_means <- sim_means %>% 
            mutate(inside = (mu>conf.low & mu<conf.high))
        
        # join the new data with the reactive dataset 
        T3simsdata2$data <- bind_rows(T3simsdata2$data, sim_means)
        
        # use the dataset just created to develop the pie chart dataset 
        T3PieDat$data <- T3simsdata2$data %>%
            ungroup() %>%
            group_by(inside) %>% 
            tally()
        
    })
    
    observeEvent(input$Add100, {
        
        # set the confidence level 
        level <- input$confLevel/100
        
        # add 10 implies we want to add another 10 simulations to our dataset 
        n_sims <- 100
        # we will stick to using 100 samples 
        samp <- 100
        
        # sampling from the standard normal distribution
        mu <- 0
        sigma <- 1
        
        # setting up the grid 
        sim_data <- expand.grid(sample = 1:samp, 
                                sim = (nrow(T3simsdata2$data) + 1):(nrow(T3simsdata2$data) + n_sims))
        # extract the samples 
        sim_data <- mutate(sim_data, 
                           height = rnorm(n = samp*n_sims, 
                                          mean = mu, 
                                          sd = sigma))
        
        # compute all the bits 
        sim_means <- sim_data %>% 
            group_by(sim) %>% 
            summarise(xbar = mean(height), 
                      n = n(), 
                      s = sd(height), 
                      alpha = 1-level) %>% 
            # same as above 
            mutate(t.low = qt(p = alpha/2, df = n-1),
                   t.high = qt(p = 1-alpha/2, df = n-1),
                   conf.low = xbar + t.low*(s/sqrt(n)), 
                   conf.high = xbar + t.high*(s/sqrt(n)))
        
        sim_means <- sim_means %>% 
            mutate(inside = (mu>conf.low & mu<conf.high))
        
        # update sim data 
        T3simsdata2$data <- bind_rows(T3simsdata2$data, sim_means)
        
        # update pie chart data 
        T3PieDat$data <- T3simsdata2$data %>%
            ungroup() %>%
            group_by(inside) %>% 
            tally()
        
    })
    
    
    output$ConfIntPie <- renderPlot({
        
        # this is how you make a pie chart using ggplot2 
        ggplot(T3PieDat$data, aes(x="", y=n, fill=inside))+
            geom_bar(width = 1, stat = "identity") +
            # nasty 
            coord_polar("y", start=0) +
            theme_minimal() +
            theme(
                # get rid of all those markers 
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.border = element_blank(),
                panel.grid=element_blank(),
                axis.ticks = element_blank(),
                axis.text.x=element_blank(),
                plot.title=element_text(size=14, face="bold")
            ) +
            scale_fill_manual(values = c("salmon", 
                                         # there's that theme colour 
                                          rgb(24, 188, 156,
                                              maxColorValue = 255))) 
        
        
    })
    
    output$ratio <- renderText({
        
        ## Adding text underneath the pie chart is much easier than getting the 
        # proportions to show on the actual plot in the right place 
        # i is for inside 
        i <- T3PieDat$data[T3PieDat$data$inside == TRUE, 2]
        # o is for outside 
        o <- T3PieDat$data[T3PieDat$data$inside == FALSE, 2]
        
        # compute percentages 
        inside <- round(i/(i + o) * 100, 2)
        outside <- round(o/(i + o) * 100, 2)
        paste0("Percentage Contained = ", inside, "%, Percentage Missed = ", outside, "%")
        
    })
    
    output$PlotProp <- renderPlotly({
        # take the simulation dataset 
        PropDat <- T3simsdata2$data 
        
        # use it to compute a cumulative sum of the intervals which DO contain the mean
        PropDat <- PropDat %>% 
            mutate(sum = cumsum(inside), 
                   # and then work out what proportion of total sims this is 
                   prop = sum/sim)
        
        # plot the result against number of sims 
        ggplot(PropDat, aes(x = sim, y = prop)) + 
            geom_line() + 
            # we want a line showing the conf level - the graph should approach this value 
            geom_hline(yintercept = input$confLevel/100, 
                       lty = 2, 
                       color = "red") +
            theme_bw() +
            ylim(c(0, 1)) +
            labs(x = "Simulation Number", 
                 y = "Proportion")
        
    })
    
    # this bad boy will be interactive 
    output$PlotIntervals <- renderPlotly({
        
        # here we are just plotting the sample means and their intervals 
        # the point will be the means 
        ggplot(data = T3simsdata2$data, aes(y = sim, x = xbar)) + 
            # use geom_pointrange to plot the intervals around the points 
            geom_pointrange(aes(xmin = conf.low, xmax = conf.high, color = inside)) +
            theme_bw() + 
            # fix the x limits so that we can actually see the difference when changing the confidence level 
            xlim(c(-0.6, 0.6)) +
            # add a red line indicating the true mean 
            geom_vline(xintercept = 0, lty = 3, color = "red") +
            scale_color_manual(values = c("salmon", 
                                          # hey theres that theme colour again 
                                          rgb(24, 188, 156,
                                              maxColorValue = 255))) +
            labs(x = "Sample Mean", 
                 y = "Simulation Number") 
        
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
