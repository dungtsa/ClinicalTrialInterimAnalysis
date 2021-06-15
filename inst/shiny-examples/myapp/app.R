#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
 
library(shiny)
#library(clinfun)
#library(Biobase) 
library(knitr)
library(rmarkdown) 
library(DT)
library(lubridate) 
library(janitor)
library(shinyjs)
library(ggplot2)       
library(shinythemes)
library(survival)
library(tidyverse)

#---one-sample survival analysis
power_one_sample_logrank <- function(medH0=30, medH1=42, n=120, time.accrual=24, time.followup=24 ,n.sim=10000, con.CI=.9)
{
    #  require(survminer)
    hazard.rate0<-log(2)/medH0
    med0<-median(rexp(10000,rate=hazard.rate0))
    hazard.rate1<-log(2)/medH1
    med1<-median(rexp(10000,rate=hazard.rate1))
    
    time.total<-time.accrual+time.followup
    event.rate.H0=1-1/((time.followup)*hazard.rate0)*(exp(-time.accrual*hazard.rate0)-exp(-time.total*hazard.rate0))
    event.rate.H1=1-1/((time.followup)*hazard.rate1)*(exp(-time.accrual*hazard.rate1)-exp(-time.total*hazard.rate1))
    
    p.tmp<-numeric(0)
    for(i in 1:n.sim)
    {
        x<-rexp(n,rate=hazard.rate1)
        c1=runif(n,time.accrual,time.total)
        data1<-data.frame(S_time=x,C_time=c1,Obs_time=ifelse(x>c1,c1,x),censor=ifelse(x>c1,0,1))
        #L1=(sum(data1$censor)-event.rate.H0*n)/sqrt(event.rate.H0*n)
        s1=Surv(data1$Obs_time,data1$censor)
        fit1<-survfit(s1~1,conf.int=con.CI)
        #surv.med<-summary(fit1)$table[c("median" ,    "0.9LCL" ,    "0.9UCL")]
        surv.med<-summary(fit1)$table[-(1:6)]
        s0=1-pexp(data1$Obs_time,rate=hazard.rate0)
        d11=survdiff(s1~offset(s0))$chisq
        Event.E0<-sum(-log(s0))
        Event.O0<-sum(data1$censor)
        d1<-(Event.O0-Event.E0)/sqrt(Event.E0)
        p.tmp<- rbind(p.tmp,c(Event.O0,d1,1-pchisq(d1^2,df=1),Event.E0,d11,1-pchisq(d11,df=1),surv.med))
    }
    p.tmp
}
#---one-sample survival analysis with two-stage design for futility 

power_one_sample_logrank_interim <- function(medH0=30, medH1=42, n1=60, n2=60, time.accrual=24, time.followup=24, n.sim=10000)
{  
    tt.all<-numeric(0)
    p.list<-list()
    med.list<-c(medH1,medH0)
    for(j in 1:length(med.list))
    {
        hazard.rate0<-log(2)/medH0
        med0<-median(rexp(10000,rate=hazard.rate0))
        hazard.rate1<-log(2)/med.list[j]
        med1<-median(rexp(10000,rate=hazard.rate1))
        
        n=n1+n2
        time.total<-time.accrual+time.followup
        event.rate.H0=1-1/((time.total-time.accrual)*hazard.rate0)*(exp(-time.accrual*hazard.rate0)-exp(-time.total*hazard.rate0))
        event.rate.H1=1-1/((time.total-time.accrual)*hazard.rate1)*(exp(-time.accrual*hazard.rate1)-exp(-time.total*hazard.rate1))
        
        p0.tmp<-p.tmp<-matrix(numeric(),n.sim,4)
        for(i in 1:n.sim)
        {
            k1<-ceiling(n1/(time.accrual/2))
            #--1st stage: accrual k1 (e.g., 3) subjects per time unit (per month)--
            x1<-rexp(n1,rate=hazard.rate1)
            c10<-time.accrual/2-rep((0:(ceiling(n1/k1)-1)),each=k1)
            c10<-c10[1:n1]
            data10<-data.frame(S_time=x1,C_time=c10,Obs_time=ifelse(x1>c10,c10,x1),censor=ifelse(x1>c10,0,1))
            #  s10=Surv(data10$Obs_time,data10$censor)
            s00=1-pexp(data10$Obs_time,rate=hazard.rate0)
            Event.E0<-sum(-log(s00))
            Event.O0<-sum(data10$censor)
            d1<-(Event.O0-Event.E0)/sqrt(Event.E0)
            p0.tmp[i,]<- c(Event.O0,d1,1-pchisq(d1^2,df=1),Event.E0)
            #p0.tmp<- rbind(p0.tmp,c(Event.O0,d1,1-pchisq(d1^2,df=1),Event.E0))
            
            #--2nd stage at end of study: full analysis--  
            c1_all=runif(n,time.accrual,time.total)
            x2<-rexp(n-n1,rate=hazard.rate1)
            x_all<-c(x1,x2)
            data_all<-data.frame(S_time=x_all,C_time=c1_all,Obs_time=ifelse(x_all>c1_all,c1_all,x_all),censor=ifelse(x_all>c1_all,0,1))
            #  s1_all=Surv(data_all$Obs_time,data_all$censor)
            s00_all=1-pexp(data_all$Obs_time,rate=hazard.rate0)
            Event.E<-sum(-log(s00_all))
            Event.O<-sum(data_all$censor)
            d1_all<-(Event.O-Event.E)/sqrt(Event.E)
            p.tmp[i,]<-c(Event.O,d1_all,1-pchisq(d1_all^2,df=1),Event.E)
            #p.tmp<- rbind(p.tmp,c(Event.O,d1_all,1-pchisq(d1_all^2,df=1),Event.E))
        }
        
        p.list[[j]]<-cbind(p.tmp,p0.tmp)
        tt1<-numeric(0)
        for(p1 in seq(.05,.15,by=0.01))  # p1 as p value to deterimne significance 
            for(p2 in seq(-0.1,.5,by=0.01)) # p2 as a threshold to stop at end of 1st stage (stop if > threshold)
            {
                tt0<-as.vector(table((p.tmp[,3]<p1)&(p.tmp[,2]<0),p0.tmp[,2]>p2))/n.sim
                tt1<-rbind(tt1,c(p1,p2,tt0))
            }
        
        tt.all<-cbind(tt.all,tt1)
    }
    names(p.list)<-med.list
    dimnames(tt.all)[[2]]<-c('p.2ndStage','boundary.1st','H1_NonSig_Pass','H1_Sign_Pass','H1_NonSig_Stop','H1_Sig_Stop','p.2ndStage','boundary.1st','H0_NonSig_Pass','H0_Sign_Pass','H0_NonSig_Stop','H0_Sig_Stop')
    
    return(list(p.list = p.list, summary = tt.all, n.accrual.per.time.unit = k1))
}



# Define UI for application that draws a histogram ####
ui <- shinyUI(fluidPage(

    # Application title
    titlePanel("One Sample Interim Analysis"),
   
    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            
            numericInput("medH0",
                         "Median Survival time (months) in Control group:",
                         min = 1,
                         max = 120,
                         step = 1,
                         value = 30),
            
            numericInput("medH1",
                         "Median Survival time (months) in Treatment group:",
                         min = 1,
                         max = 120,
                         step = 1,
                         value = 42),
            
            numericInput("n.sim",
                         "number of simulations:",
                         min = 1,
                         max = 10000,
                         step = 100,
                         value = 2000),
            
            numericInput("n1",
                         "number of patients enrolled initially:",
                         min = 1,
                         max = 1000,
                         step = 1,
                         value = 120),
            
            numericInput("n2",
                         "number of additional patients enrolled:",
                         min = 1,
                         max = 1000,
                         step = 1,
                         value = 120),
            
            
            numericInput("time.accrual",
                         "Accrual time (months):",
                         min = 1,
                         max = 120,
                         step = 1,
                         value = 24),
            
            numericInput("time.followup",
                         "Follow up time (months):",
                         min = 1,
                         max = 120,
                         step = 1,
                         value = 24),
            
            # radioButtons("con.CI", "Confidence interval:",
            #              choiceNames = list("95%" , "90%"  ),
            #              choiceValues  = list( 0.95  , 0.90  ) ),
            
            hr(),
            actionButton("go", "Submit"), 
            
            hr(),
            hr(),
            textInput("endpointShort", "Short Endpoint (eg. OS, DFS", "DFS"),
            
            textInput("endpointLong", "Long Endpoint (eg. Overall Sruvival, 
                      disease-free survival", "disease-free survival"),
            
            textInput("endpointDescription", "Describe the endpoint", "the time from randomization to the earliest event which includes local, 
                      regional, or distant recurrence, new lung cancer, or death"),
            
            radioButtons('format', 'Document format', c('Word', 'PDF'),
                         
                         inline = TRUE),
            
            downloadButton('downloadReport', "Download Report"),
            
 
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("scatterplot1"),   
            verbatimTextOutput("text"),
            dataTableOutput(outputId = "table2") ,  
            verbatimTextOutput("text2"),
            dataTableOutput(outputId = "table1") 
        )
    )
))

# Define server logic ####
server <- shinyServer(function(input, output) {
    
    
    # output$text <- renderText( {"Best case simulation results" })
    # output$text2 <- renderText( {"All simulation results" })
    # 
    
    text <- eventReactive(input$go, {"Best case simulation results" })
    
    text2 <- eventReactive(input$go, {"All simulation results" })
    
    output$text <- renderText({   text()   })
    
    output$text2 <- renderText({   text2()    })
    
    
    
    
    get.result  <-  eventReactive(input$go,{ 
   
    
        power_one_sample_logrank_interim(medH0 = input$medH0, 
                                         medH1 = input$medH1, 
                                         n1 = input$n1, 
                                         n2 = input$n2, 
                                         time.accrual = input$time.accrual, 
                                         time.followup = input$time.followup, 
                                         n.sim = input$n.sim)
            
    })
    
    output$table1  <- renderDataTable({
        x <- as.data.frame(get.result()$summary); 
        return(x)} , rownames = TRUE)
    
    
    
    
    somevalues <- eventReactive(input$go,{ 
        
        
        test <- as.data.frame(get.result()$summary[,c(-1,-2)] ) %>% filter(H0_Sign_Pass < 0.05) 
        type1lt05 <- as.data.frame(test ) %>% filter(H0_Sign_Pass < 0.05) %>% arrange(desc(H0_Sign_Pass), desc(boundary.1st), desc(H1_Sign_Pass)) 
        type1lt05 <-  type1lt05[1, ]
         # list( pow = round(type1lt05[ ,2], 2),
         #             typI = round(type1lt05[ ,8], 2),
         #             cut1 = round(type1lt05[ ,6], 2),
         #             cut2 = -round(sqrt(qchisq(1-type1lt05[ ,6],df=1)),2))
        return(list( pow = round(type1lt05$H1_Sign_Pass, 3),
                     typI = round(type1lt05$H0_Sign_Pass, 3),
                     cut1 = round(type1lt05$boundary.1st, 2),
                     cut2 = -round(sqrt(qchisq(1 - type1lt05$p.2ndStage, df = 1)), 2),
                     PETH0 = round(type1lt05$H0_NonSig_Stop + type1lt05$H0_Sig_Stop, 2),
                     PETH1 = round(type1lt05$H1_NonSig_Stop + type1lt05$H1_Sig_Stop, 2))
        )
          
                                               })

    
    
    
    observe({
         print("DOES THIS WORK:")
         print(somevalues())
         print(somevalues()$pow)
         print(somevalues()$typI)
         print(somevalues()$cut1)
         print(somevalues()$cut2)
        
          plotdata <- get.result()$summary[,c(-1,-2)] 
          print(str(plotdata))
        # test<-as.data.frame(plotdata ) %>% filter(H0_Sign_Pass < 0.05) 
        # print(length(test$H0_Sign_Pass[test$H0_Sign_Pass == max(test$H0_Sign_Pass)]))
        # 
    #  print( NROW(type1lt05 <- as.data.frame(plotdata) %>% filter(H0_Sign_Pass < 0.05) %>% arrange(H0_Sign_Pass, boundary.1st, H1_Sign_Pass) ))
    # print(summary(type1lt05$H1_Sign_Pass))
    # print( type1lt05)
    })
    
    output$table2  <- renderDataTable({
        
        
        plotdata <- get.result()$summary[,c(-1,-2)] 
        test <- as.data.frame(plotdata ) %>% filter(H0_Sign_Pass < 0.05) 
        type1lt05 <- as.data.frame(plotdata) %>% filter(H0_Sign_Pass < 0.05) %>% arrange(desc(H0_Sign_Pass), desc(boundary.1st), desc(H1_Sign_Pass)) 
        type1lt05 <- head(type1lt05, length(test$H0_Sign_Pass[test$H0_Sign_Pass == max(test$H0_Sign_Pass)]))
        return(type1lt05)} , rownames = TRUE)
    
    
    output$scatterplot1 <- renderPlot({ 
        x <- as.data.frame(get.result()$summary[, c(-1,-2)]);
        print(str(x)) 
        ggplot(x, aes(x = H1_Sign_Pass, y = H0_Sign_Pass)) + geom_point() + 
             geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")
        # Change the point size, and shape 
            #geom_point(size=2, shape=23)
     

    })
    
    #---output report----
    
    output$downloadReport <- downloadHandler(
        
        filename = function() {
            
            paste('InterimReport', sep = '.',
                  
                  switch(input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'))
            
        },
        
        content = function(file) {
            
            out <- rmarkdown::render(input = 'F:\\myGitRepo\\clinical_trial_design\\clinical_trial_design\\OneSampleInterimAnalysis\\InterimReport.Rmd',
                                     
                                     output_format =
                                         
                                         switch(input$format,
                                                
                                                PDF = rmarkdown::pdf_document(),
                                                
                                                HTML = rmarkdown::html_document(),
                                                
                                                Word = rmarkdown::word_document()
                                                
                                         ),
                                     
                                     params = list(set_title = input$project_title, set_author = input$author_input)
                                     
            )
            
            #file.rename(out, file)
            file.copy(out,file)
            
        }
        
    )
    
    
})


 # END O'script  
shinyApp(ui=ui,server=server)