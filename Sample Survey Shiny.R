library(shiny)
ui<-fluidPage(fluidRow(column(offset=2,width=3,sliderInput(inputId="mu", label="Population Mean", min=0, max=100,value=8)),
                       column(width=3,sliderInput(inputId="sigma", label="Population Standard Deviation", min=0, max=10,value=4)),
                       column(width=3,sliderInput(inputId="sigma1", label="Sigma_1", min=0, max=10,value=2))),
              fluidRow(column(offset=2,width=3,sliderInput(inputId="poplnsize", label="Population Size", min=0, max=100000,value=10000)),
                       column(width=3,sliderInput(inputId="samsize", label="Sample Size", min=0, max=1000,value=50)),
                       column(width=3,sliderInput(inputId="simsize", label="Simulation Size", min=0, max=10000,value=500))),
              fluidRow(verbatimTextOutput("print1")))
server<-function(input,output){
  output$print1=renderPrint({
    set.seed(100)
    
    x=rlnorm(n=input$poplnsize,meanlog=input$mu, sdlog=input$sigma)
    z=rlnorm(n=input$poplnsize,meanlog=0, sdlog=input$sigma1)
    y=x*z
    population_total=sum(y)
    
    est.population.total.SRSWR=c()
    for(i in 1:input$simsize){
      samp.srswr=sample(y,size=input$samsize,replace=T) 
      est.population.total.SRSWR[i]=input$poplnsize*mean(samp.srswr)}
    SRSWR_mean=mean(est.population.total.SRSWR)
    SRSWR_variance=var(est.population.total.SRSWR)
    
    est.population.total.SRSWOR=c()
    for(i in 1:input$simsize){
      samp.srswor=sample(y,size=input$samsize,replace=F) 
      est.population.total.SRSWOR[i]=input$poplnsize*mean(samp.srswor)}
    SRSWOR_mean=mean(est.population.total.SRSWOR)
    SRSWOR_variance=var(est.population.total.SRSWOR)
    
    q=quantile(y,p=seq(0,1,0.3))
    s1=y[which(y<=q[2])] #1st stratum
    s2=y[which(q[2]<y & y<=q[3])] #2nd stratum
    s3=y[which(q[3]<y & y<=q[4])] #3rd stratum
    s4=y[which(q[4]<=y)] #4th stratum
    Nh=c(length(s1),length(s2),length(s3),length(s4))
    nh=(Nh/input$poplnsize)*input$samsize
    est.population.total.ST=c()
    for(i in (1:input$simsize)){
      sample1=sample(s1,nh[1],replace=T) #sample selection from 1st stratum
      sample2=sample(s2,nh[2],replace=T) #sample selection from 2nd stratum
      sample3=sample(s3,nh[3],replace=T) #sample selection from 3rd stratum
      sample4=sample(s4,nh[4],replace=T) #sample selection from 4th stratum
      samp=c(sample1,sample2,sample3,sample4)
      est.population.total.ST[i]=sum(mean(samp)*Nh)}
    ST_mean=mean(est.population.total.ST)
    ST_variance=var(est.population.total.ST)
    
    data=data.frame(y,x)
    max(data$x)
    X=sum(data$x)
    m=16926905455
    est=array(0,input$simsize)
    w=1:input$poplnsize
    for(i in (1:input$simsize)){
      for(j in (1:input$samsize)){
        r=sample(w,1)
        s=runif(1,1,m)
        if(s<=data$x[r]){
          est[i]=est[i]+((data$y[r]*X)/(data$x[r]))
        }else{
          j=j-1
        }
      }
    }
    PPSWR_mean=mean(est)
    PPSWR_variance=var(est)
    
    cat("The Population Total is", population_total)
    cat("\nMean of population total estimates by SRSWR is",SRSWR_mean,"and variance of population total estimates by SRSWR is",SRSWR_variance)
    cat("\nMean of population total estimates by SRSWOR is",SRSWOR_mean,"and variance of population total estimates by SRSWOR is",SRSWOR_variance)
    cat("\nMean of population total estimates by Stratified is",ST_mean,"and variance of population total estimates by Stratified is",ST_variance)
    cat("\nMean of population total estimates by PPSWR is",PPSWR_mean,"and variance of population total estimates by PPSWR is",PPSWR_variance)
    
  })
}
shinyApp(ui,server)
