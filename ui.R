library(shiny)
library(shinythemes)

shinyUI(fluidPage(
  theme = shinytheme("readable"),
  column(
    width = 10,
    #########################################################################################
    h2('Drug Combination in vitro Experiment Simulation'),
    h4(em('Lingge Li')),
    h4(em('November 2016')),
    p(
      'This web app is built for in vitro drug combination experiment simulation. In this setting,
      for a single drug, the dose-repsonse curve is usually modelled with a sigmoid function,
      known as the Emax equation. It depends on two parameters, IC50 and m (shape) where IC50
      is often the target of interest.'
    ),
    #########################################################################################
    tabsetPanel(
      tabPanel(
        'Drug Effects',
        h4(''),
        p(
          'Before simulating an experiment, we need to come up with a hypothetical drug. E0 is the
          number of cells we expect to see in a control well without any drug. IC50 is the median
          effect concentration (an implicit assumption is that the drug will kill
          all the cells at a very high concentraion). m the shape parameter usually ranges from 0
          to 2 with higher value making the curve more sigmoidal. The default plot is what we
          expect to observe as raw data.'
        ),
        # Controls
        fluidRow(
          column(width = 4,
                 numericInput(
                   'IC50A', label = "Drug A IC50", value = 50
                 )),
          column(width = 4,
                 numericInput(
                   'IC50B', label = "Drug B IC50", value = 50
                 )),
          column(width = 4,
                 numericInput(
                   'IC50M', label = "Mixture IC50", value = 50
                 ))
        ),
        fluidRow(
          column(
            width = 4,
            numericInput("mA", label = "m (shape)", value = 1),
            checkboxInput('fa', label = 'Emax', value =
                            FALSE)
          ),
          column(width = 4,
                 numericInput("mB", label = "m (shape)", value = 1)),
          column(width = 4,
                 numericInput("mM", label = "m (shape)", value = 1))
        ),
        # Plots
        fluidRow(
          plotOutput('distPlot')),
        fluidRow(
          plotOutput("iso"))
        ),
      #########################################################################################
      tabPanel(
        'Experimental Design',
        h4(''),
        p(
          'Here we focus on the simplest experiment where only one cell line is used on one plate.
          The experimental design involves choosing a range of drug concentrations. For example,
          on a 24-well plate, it is common to have 4 controls and 5 concentrations with 4
          replicates. In practice, drug concentrations often follow 2-fold serial dilution. Note
          that concentration 0 corresponds to control.'
        ),
        fluidRow(
          column(width = 4,
                 selectInput(
                   "drug", "Drug",
                   choices = c('A', 'B', 'Mixture')
                 )),
          column(
            width = 4,
            numericInput("conc", label = "Drug concentration", value = 0)
          ),
          column(width = 4,
                 numericInput("k", label = "# of replicates", value = 4)),
          column(
            width = 3,
            actionButton('add', label = 'Add concentration'),
            offset = 3
          ),
          column(width = 3,
                 actionButton('clear', label = 'Start over'))
        ),
        fluidRow(
          verbatimTextOutput('wells1'),
          verbatimTextOutput('wells2'),
          verbatimTextOutput('wells3')
        ),
        p(
          'The noise is generated with',
          withMathJax(),
          helpText(
            '$$[E(Y)]^{\\lambda}\\cdot\\epsilon\\text{ where }\\epsilon\\sim N(0,\\sigma).$$'
          ),
          'Obviously, sigma is the standard deviation of Gaussian noise. When lambda is 0, the noise
          would be independent. When lambda is 1, the noise would be proportional to the expected cell count.
          In other words, there would be more variation where the number of cells is large.'
        ),
        fluidRow(
          column(
            width = 3,
            numericInput("lambda", label = "lambda", value = 0),
            offset = 3
          ),
          column(width = 3,
                 numericInput(
                   "sigma", label = "sigma", value = 1
                 ))
        ),
        fluidRow(column(
          width = 2,
          actionButton('cook', label = 'Simulate'),
          offset = 5
        )),
        fluidRow(
          plotOutput("sample"))
        ),
      #########################################################################################
      tabPanel(
        'Data Analysis',
        h4(''),
        p(
          'The most popular data analysis method transforms the data to log scale and fits
          a straight line. The output shown below includes the estimated parameter values and confidence intervals.'
        ),
        fluidRow(column(
          width = 2,
          actionButton('crunch', label = 'Fit models'),
          offset = 5
        )),
        fluidRow(
          plotOutput("trans")),
        fluidRow(
          column(width = 4,
                 strong('Fit A'),
                 verbatimTextOutput('conf1')),
          column(width = 4,
                 strong('Fit B'),
                 verbatimTextOutput('conf2')),
          column(width = 4,
                 strong('Fit M'),
                 verbatimTextOutput('conf3'))
        ),
        fluidRow(
          column(width = 4,
                 strong('L20'),
                 verbatimTextOutput('index1')),
          column(width = 4,
                 strong('L50'),
                 verbatimTextOutput('index2')),
          column(width = 4,
                 strong('L80'),
                 verbatimTextOutput('index3'))
        )
        #h4('Asymptotic performance'),
        #p(
        #'We can perform hundreds of simulations to assess confidence interval coverage and
        #asymptotic efficiency under a specific design.'
        #),
        #fluidRow(column(
        #width = 2,
        #actionButton('run', label = 'Simulate 500 times'),
        #offset = 5
        #)),
        #fluidRow(column(
        #width = 8,
        #plotOutput("density"),
        #verbatimTextOutput('performance'),
        #offset = 2
        #))
        )
    ),
    offset = 1
    )
  ))
