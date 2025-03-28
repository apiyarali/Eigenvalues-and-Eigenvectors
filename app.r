#Eignevalues and Eignevectors
library(shiny)
library(shinydashboard)
library(shinyWidgets)
source("jaxmat.R")   #for displaying mathematics
stylesheet <- tags$head(tags$style(HTML('
    .main-header .logo {
      font-family: "Georgia", Times, "Times New Roman", serif;
      font-weight: bold;
      font-size: 24px;
    }
  ')
))
#The user interface
header <- dashboardHeader(title = "Eignevalues and Eignevectors",
                          titleWidth = 500)
sidebar <- dashboardSidebar(disable = TRUE)
body <- dashboardBody(
  fluidRow(stylesheet,
    column(width=4,
      actionBttn("btnmake","Make a new matrix"),
      uiOutput("matA"),br(),
      actionBttn("btnpoly","Calculate Characteristic Polynomial"),
      uiOutput("poly"),br(),
      actionBttn("btneigen","Eigenvalues"),
      uiOutput("eigen"),br(),
      actionBttn("btneigenv","Eigenvector"),
      uiOutput("eigenv1"),br(),
      uiOutput("eigenv2"),br(),
      uiOutput("eigenv3")
    ),
    column(width=4,
       h4("Cayley-Hamilton Matrix"),
       uiOutput("CH"),br()
    )
  )
)
ui <- dashboardPage(header, sidebar, body, skin = "green") #other colors available

#Functions that implement the mathematics
#This file must go into the same directory as app.R
#source(".R")
#Any images, audio, or stylesheet must go into a subfolder named www

#Additional functions are OK here, but no variables


server <- function(session, input, output) {
  #Variables that are shared among server functions (use <<-)
  
  #Characteristic Polynomial a,b,c
  i <- cbind(c(1,0,0), c(0,1,0), c(0,0,1))
  
  #Initialization
  
  
  #Functions that respond to events in the input
  
  # M <- matrix(sample(0:4, 9, replace = TRUE), nrow = 3)
  # output$matA <- renderUI(jax.matrix(M, name = "M"))
  
  # M1 <- matrix(sample(0:4, 9, replace = TRUE), nrow = 3)
  # output$matA <- renderUI(jax.matrix(M1, name = "M1"))
  
  createMatrix <- function(){
    M <- matrix(sample(0:4, 9, replace = TRUE), nrow = 3)
    return(M)
  }

  observeEvent(input$btnmake,{
    M1 <<- createMatrix()
    # M1 <- diag(3)
    # M1[2,2] <- 2
    # M1[3,3] <- 3
    
    output$matA <- renderUI(jax.matrix(M1, name = "M1"))
  })
  
  
  observeEvent(input$btnpoly,{
    l <<- round(((-det(M1[-2,-2]) %%5) - (det(M1[-3,-3]) %%5) - (det(M1[-1,-1])%%5)) %%5)
    l2 <<- round((M1[1,1] + M1[2,2] + M1[3,3]) %%5)
    l3 <<- -1
    constant <<- round(det(M1) %%5)
    
    # print(l)
    # print(l2)
    # print(c)
    
    output$poly <- renderUI(h4(withMathJax(helpText("Characteristic Polynomial is: $$",l3,"\\lambda^3 + ", l2, "\\lambda^2 +",l,"\\lambda +",constant,"$$"))))
  })
  
  eigenvalues <- function(lambda, M3){
    e <- (((4 * (lambda^3)) %%5) + ((l2 * (lambda^2)) %%5) + ((l * lambda) %%5) + constant) %%5
    return(e)
  }
  
  observeEvent(input$btneigen,{
    v <<- vector()
    
    for (i in 0:4) {
      eigen <- eigenvalues(i,M1)
      if (eigen == 0){
        v <<- c(v,i)
      }
    }
    
    # print(v)
    
    if (length(v) > 0){
      output$eigen <- renderUI(h4(withMathJax(paste("Eigenvalues are: ",vecToString(v)))))  
    } else {
      output$eigen <- renderUI(h4("Characteristic Polynomial has no solution so, this matrix has no eigenvalues"))  
    }
    
  })
  
  observeEvent(input$btneigenv,{
    
    # print((round(eigen(M1)$values,0))%%5)
    # ev <- eigen(M1)
    # print(ev$values)
    
    if (length(v) != 3) {
      output$eigenv <- renderUI(h4("You need 3 eigenvalues to determine eigenvectors"))
    } else {
      # w1 <- cbind(c(1,0,0))
      # w2 <- cbind(c(0,1,0))
      # w3 <- cbind(c(0,0,1))
      
      ev1 <- (M1 - (v[2]*i))%*%(M1 - (v[3]*i))
      # print(ev1)
      
      # ev1a <- ev1[,1]*w1
      # ev1b <- ev1[,2]*w2
      # ev1c <- ev1[,3]*w3
      
      if (ev1[1,1] == 0 && ev1[2,1] == 0 && ev1[3,1] == 0){
        if (ev1[1,2] == 0 && ev1[2,2] == 0 && ev1[3,2] == 0){
          v1 <- cbind(c(ev1[1,3],ev1[2,3],ev1[3,3]))
          output$eigenv1 <- renderUI(jax.matrix(v1, name = "V1"))
        } else {
          v1 <- cbind(c(ev1[1,2],ev1[2,2],ev1[3,2]))
          output$eigenv1 <- renderUI(jax.matrix(v1, name = "V1"))  
        }  
      } else {
        v1 <- cbind(c(ev1[1,1],ev1[2,1],ev1[3,1]))
        output$eigenv1 <- renderUI(jax.matrix(v1, name = "V1"))  
      }
      
      ev2 <- (M1 - (v[1]*i))%*%(M1 - (v[3]*i))
      # print(ev2)
      
      # ev2a <- ev2[,1]*w1 
      # ev2b <- ev2[,2]*w2
      # ev2c <- ev2[,3]*w3
      
      if (ev2[1,1] == 0 && ev2[2,1] == 0 && ev2[3,1] == 0){
        if (ev2[1,2] == 0 && ev2[2,2] == 0 && ev2[3,2] == 0){
          v2 <- cbind(c(ev2[1,3],ev2[2,3],ev2[3,3]))
          output$eigenv2 <- renderUI(jax.matrix(v2, name = "V2"))
        } else {
          v2 <- cbind(c(ev2[1,2],ev2[2,2],ev2[3,2]))
          output$eigenv2 <- renderUI(jax.matrix(v2, name = "V2"))  
        }  
      } else {
        v2 <- cbind(c(ev2[1,1],ev2[2,1],ev2[3,1]))
        output$eigenv2 <- renderUI(jax.matrix(v2, name = "V2"))  
      } 
      
      ev3 <- (M1 - (v[1]*i))%*%(M1 - (v[2]*i))
      # print(ev3)
      
      # ev3a <- ev3[,1]*w1
      # ev3b <- ev3[,2]*w2
      # ev3c <- ev3[,3]*w3
      
      if (ev3[1,1] == 0 && ev3[2,1] == 0 && ev3[3,1] == 0){
        if (ev3[1,2] == 0 && ev3[2,2] == 0 && ev3[3,2] == 0){
          v3 <- cbind(c(ev3[1,3],ev3[2,3],ev3[3,3]))
          output$eigenv3 <- renderUI(jax.matrix(v3, name = "V3"))
        } else {
          v3 <- cbind(c(ev3[1,2],ev3[2,2],ev3[3,2]))
          output$eigenv3 <- renderUI(jax.matrix(v3, name = "V3"))  
        }  
      } else {
        v3 <- cbind(c(ev3[1,1],ev3[2,1],ev3[3,1]))
        output$eigenv3 <- renderUI(jax.matrix(v3, name = "V3"))  
      }
     
      ch1 <- ((-(M1%*%M1%*%M1)%%5)+((l2*M1%*%M1)%%5)+((l3*M1)%%5)+((constant*i)%%5))
      # print(ch1)
      output$CH <- renderUI(jax.matrix(ch1, name = "CH"))
      
       
    }
  })
  
  # function from Karina's section 2 button app.r file
  vecToString <- function(v) {
    paste(v, collapse = " and ")
  }
}

#Run the app
shinyApp(ui = ui, server = server)

