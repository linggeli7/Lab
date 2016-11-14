library(shiny)
library(shinythemes)

shinyUI(fluidPage(
  theme = shinytheme("readable"),
  column(width = 10,
  h2('IC50 in vitro Experiment Simulation'),
  h4(em('Lingge Li')),
  h4(em('November 2016')),
  p('This web app is built for in vitro dose-response experiment simulation. In this setting,
    for a single drug, the dose-repsonse curve is usually modelled with a sigmoid function,
    known as the Hill equation. It depends on two parameters, IC50 and m (shape) where IC50
    is often the target of interest.'),
  h4('Hypothetical drug'),
  p('Before simulating an experiment, we need to come up with a hypothetical drug. E0 is the
    number of cells we expect to see in a control well without any drug. IC50 is the median
    effect concentration (an implicit assumption is that the drug will kill
    all the cells at a very high concentraion). m the shape parameter usually ranges from 0
    to 2 with higher value making the curve more sigmoidal. The default plot is what we
    expect to observe as raw data. '),
  fluidRow(
    column(width = 3,
           numericInput("E0", label = "E0 (control cell count)", value = 100),
           offset = 3),
    column(width = 3,
           numericInput('IC50', label = "IC50 (median effect concentration)", value = 50)
           )
  ),
  fluidRow(
    column(width = 3,
           numericInput("m", label = "m (shape)", value = 1),
           checkboxInput('fa', label = 'Hill curve', value=FALSE),
           offset = 3),
    column(width = 3,
           numericInput("max", label = "Maximum concentration", value = 200)
           )
  ),
  fluidRow(
    column(width = 8,
           plotOutput("distPlot"),
           offset = 2)
  ),
  h4('Experimental design'),
  p('Here we focus on the simplest experiment where only one cell line is used on one plate.
    The experimental design involves choosing a range of drug concentrations. For example,
    on a 24-well plate, it is common to have 4 controls and 5 concentrations with 4
    replicates. In practice, drug concentrations often follow 2-fold serial dilution. Note 
    that concentration 0 corresponds to control.'),
  fluidRow(
    column(width = 3,
           numericInput("conc", label = "Drug concentration", value = 0),
           offset = 3),
    column(width = 3,
           numericInput("k", label = "# of replicates", value = 4)
           ),
    column(width = 3,
           actionButton('add', label = 'Add concentration'),
           offset = 3),
    column(width = 3,
           actionButton('clear', label = 'Drop the plate')
          )
  ),
  p('Concentrations'),
  fluidRow(
    column(width = 8,
           verbatimTextOutput('wells'),
           offset = 2)
  ),
  h4('Experiment simulation'),
  p('Due to low pay and high stress, graduate students in the lab cannot always get data
    without noise. The noise is generated with',
    withMathJax(),
    helpText('$$[E(Y)]^{\\lambda}\\cdot\\epsilon\\text{ where }\\epsilon\\sim N(0,\\sigma).$$'),
    'Obviously, sigma is the standard deviation of Gaussian noise. When lambda is 0, the noise
    would be independent. When lambda is 1, the noise would be proportional to the expected cell count. 
    In other words, there would be more variation where the number of cells is large.'),
  fluidRow(
    column(width = 3,
           numericInput("lambda", label = "lambda", value = 0),
           offset = 3),
    column(width = 3,
           numericInput("sigma", label = "sigma", value = 1)
    )
  ),
  fluidRow(
    column(width = 2,
           actionButton('cook', label = 'Simulate'),
           offset = 5)
  ),
  fluidRow(
    column(width = 8,
           plotOutput("sample"),
           offset = 2)
  ),
  h4('Data analysis'),
  p('The most popular data analysis method transforms the data to log scale and fits
    a straight line. The output shown below includes the estimated parameter values and confidence intervals.
    The R-squared value measures the goodness-of-fit and should be greater than 0.95.'),
  fluidRow(
    column(width = 2,
           actionButton('crunch', label = 'Fit data'),
           offset = 5)
  ),
  fluidRow(
    column(width = 8,
           plotOutput("trans"),
           verbatimTextOutput('conf'),
           verbatimTextOutput('summary'),
           offset = 2)
  ),
  h4('Asymptotic performance'),
  p('We can perform hundreds of simulations to assess confidence interval coverage and
    asymptotic efficiency under a specific design.'),
  fluidRow(
    column(width = 2,
           actionButton('run', label = 'Simulate 500 times'),
           offset = 5)
  ),
  fluidRow(
    column(width = 8,
           plotOutput("density"),
           verbatimTextOutput('performance'),
           offset = 2)
  ),
  offset = 1)
))
