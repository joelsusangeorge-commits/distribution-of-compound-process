# app.R
# Compound Poisson-Gamma Shiny App 
# =========================================================
#  Compound Poisson + Gamma Simulator 
# =========================================================

library(shiny)
library(shinyjqui)
library(shinyWidgets)
library(ggplot2)

# -------------------  STYLESHEET (rearranged) -------------------
blue_theme <- "
body,html {margin:0;padding:0;height:100%;overflow:hidden;font-family:Inter,system-ui,Segoe UI,Roboto;}
#wrap{height:100vh;width:100vw;padding:18px;background:#E2F0FF;box-sizing:border-box;}
.card{background:#fff;padding:14px;border-radius:10px;box-shadow:0 5px 18px rgba(0,0,0,0.1);}
#control_box{position:absolute;width:330px;left:18px;top:18px;z-index:50;}
#main{position:absolute;left:360px;right:18px;top:18px;bottom:18px;}
#plotSlot{height:calc(100% - 140px);}
#summary_area{height:120px;overflow:auto;}
#drag_btn{cursor:grab;z-index:99;width:55px;height:55px;display:flex;align-items:center;justify-content:center;
          background:#fff;border-radius:30px;box-shadow:0 6px 18px rgba(0,0,0,.14);position:absolute;left:18px;top:18px;}

.download_group{display:flex;gap:10px;}
.btn_sm{padding:6px 12px;border-radius:8px;background:#0A84FF;color:#fff;border:none;}
#titlecard{color:#0A3C79;max-width:440px;margin-bottom:12px;}

@media (max-width:1100px){
  #control_box{width:300px;}
  #main{left:330px;}
}
"
# ---------------------------------------------------------------


# ========================== UI  ================================
ui <- fluidPage(
  tags$head(tags$style(HTML(blue_theme)), tags$script(src="bind/jsfile.js")),
  div(id="wrap",
      
      jqui_draggable(div(id="drag_btn", class="card",
                         span(icon("sliders"), style="font-size:20px"), title="Drag or click to show controls")),
      
      div(id="control_box", class="card",
          h4("Simulation Inputs", style="color:#0A3C79;margin:4px 0 10px;"),
          numericInput("lam_num","λ - Poisson rate",0.5,0.01,5,0.01),
          sliderInput("lam_sl","",0.01,5,0.5,0.01),
          numericInput("mu_num","μ - rate for Xi",1,0.01,5,0.01),
          sliderInput("mu_sl","",0.01,5,1,0.01),
          numericInput("time_num","Time t",100,1,10000,1),
          sliderInput("time_sl","",1,10000,100,1),
          numericInput("sims","Simulations",2000,1,200000),
          numericInput("rseed","Set seed (optional)",NA,1,1e9),
          actionButton("go","Run",class="btn_sm"),
          actionButton("reset","Reset",class="btn_sm"),
          hr(), em("Move controls using ⚙ button")
      ),
      
      div(id="main",
          
          div(id="titlecard", class="card",
              h3("Compound Poisson – Gamma Explorer"),
              p("Interactive visualisation of S(t) where N(t)~Poisson(λt) and Xi~Exponential(μ).",
                style="margin-bottom:4px;font-size:14px;")
          ),
          
          div(id="plotSlot", class="card",
              plotOutput("hist", height="100%"),
              div(style="display:flex;justify-content:space-between;align-items:center;",
                  div(class="download_group",
                      downloadButton("dl_data","CSV",class="btn_sm"),
                      downloadButton("dl_img","Image",class="btn_sm")
                  ),
                  textOutput("stamp")
              )
          ),
          
          div(id="summary_area", class="card",
              fluidRow(
                column(6, h4("Statistics",style="color:#0A3C79;"), tableOutput("stat")),
                column(6, h4("Notes",style="color:#0A3C79;"), uiOutput("explain"))
              )
          )
      )
  )
)

# ========================= SERVER ===============================
server <- function(input, output, session){
  
  # ----- sync sliders & numeric -----
  observe({ if(!identical(input$lam_num,input$lam_sl)) updateSliderInput(session,"lam_sl",value=input$lam_num) })
  observe({ if(!identical(input$lam_sl,input$lam_num)) updateNumericInput(session,"lam_num",value=input$lam_sl) })
  observe({ if(!identical(input$mu_num,input$mu_sl)) updateSliderInput(session,"mu_sl",value=input$mu_num) })
  observe({ if(!identical(input$mu_sl,input$mu_num)) updateNumericInput(session,"mu_num",value=input$mu_sl) })
  observe({ if(!identical(input$time_num,input$time_sl)) updateSliderInput(session,"time_sl",value=input$time_num) })
  observe({ if(!identical(input$time_sl,input$time_num)) updateNumericInput(session,"time_num",value=input$time_sl) })
  
  # reset all
  observeEvent(input$reset,{
    updateNumericInput(session,"lam_num",0.5); updateSliderInput(session,"lam_sl",0.5)
    updateNumericInput(session,"mu_num",1); updateSliderInput(session,"mu_sl",1)
    updateNumericInput(session,"time_num",100); updateSliderInput(session,"time_sl",100)
    updateNumericInput(session,"sims",2000); updateNumericInput(session,"rseed",NA)
  })
  
  data_sim <- reactiveVal(NULL)
  last_params <- reactiveVal(NULL)
  
  # ---------- run simulation ----------
  observeEvent(input$go,{
    n <- as.integer(input$sims)
    if(is.na(n)||n<1) return(showNotification("Invalid nsim", type="error"))
    simulate_run()
  })
  
  simulate_run <- function(){
    n <- as.integer(input$sims)
    λ <- input$lam_sl; μ <- input$mu_sl; Tval <- input$time_sl
    if(!is.na(input$rseed)) set.seed(input$rseed)
    
    withProgress(message="Running simulation", value=0,{
      block <- ifelse(n<=5000,n,2000)
      loops <- ceiling(n/block)
      S <- numeric(n); idx <- 1
      
      for(i in seq_len(loops)){
        from <- idx; to <- min(n,idx+block-1); k <- to-from+1
        Nvals <- rpois(k,λ*Tval); temp <- numeric(k)
        if(any(Nvals>0)) temp[Nvals>0] <- rgamma(sum(Nvals>0),shape=Nvals[Nvals>0],rate=μ)
        S[from:to] <- temp; idx <- to+1; incProgress(1/loops)
      }
      
      data_sim(data.frame(S=S)); last_params(list(l=λ,m=μ,t=Tval,n=n))
    })
  }
  
  # ------------------ PLOT ------------------
  output$hist <- renderPlot({
    df <- data_sim()
    if(is.null(df)) return(ggplot()+theme_void()+geom_text(aes(0,0,label="Click Run"),size=6))
    
    ggplot(df,aes(S))+
      geom_histogram(aes(y=..density..),bins=60,fill="#9ECDFF",color="#005BBB",alpha=.9)+
      geom_density(color="#003A70",size=1)+
      labs(title="Distribution of S(t)", x="S(t)", y="Density")+
      theme_minimal()
  })
  
  # ------------------ TABLE ------------------
  output$stat <- renderTable({
    df <- data_sim(); if(is.null(df)) return(NULL)
    x <- df$S; EN <- input$lam_sl*input$time_sl; ES <- EN*(1/input$mu_sl)
    
    data.frame(
      Statistic=c("Mean","Variance","Median","2.5%","97.5%","E[N]","E[S]"),
      Value=round(c(mean(x),var(x),median(x),quantile(x,0.025),quantile(x,0.975),EN,ES),6)
    )
  })
  
  output$stamp <- renderText({
    p <- last_params(); if(is.null(p)) return("")
    paste("λ=",p$l," μ=",p$m," t=",p$t)
  })
  
  output$explain <- renderUI({
    df <- data_sim(); if(is.null(df)) return("No simulation yet.")
    mean_S <- mean(df$S); EN <- input$lam_sl*input$time_sl
    
    tagList(
      tags$ul(
        tags$li(paste("E[N(t)] =",round(EN,3))),
        tags$li(paste("E[X] =",round(1/input$mu_sl,4))),
        tags$li(paste("Mean S(t) ≈",round(mean_S,4))),
        tags$li(if(EN>50) "CLT → Near normal" else "Right-skew, many zeros possible")
      )
    )
  })
  
  # DOWNLOADS
  output$dl_data <- downloadHandler(
    filename=function()paste0("sim_",Sys.Date(),".csv"),
    content=function(f) write.csv(data_sim(),f,row.names=FALSE)
  )
  
  output$dl_img <- downloadHandler(
    filename=function()paste0("plot_",Sys.Date(),".png"),
    content=function(f){
      png(f,1200,800);df<-data_sim()
      print(ggplot(df,aes(S))+geom_histogram(aes(y=..density..),bins=60,fill="#9ECDFF",color="#005BBB")+
              geom_density()+theme_minimal());dev.off()
    }
  )
}

# ====== JS Injection (rewritten formatting) ======
js <- "
document.addEventListener('DOMContentLoaded',()=>{
  const el=document.getElementById('drag_btn');
  el.addEventListener('click',()=>{
    const p=document.getElementById('control_box');
    p.style.display=p.style.display==='none'?'block':'none';
  });
});
"

onStart <- function(){
  path <- file.path(tempdir(),"bind"); dir.create(path,showWarnings=FALSE)
  fileConn <- file.path(path,"jsfile.js"); writeLines(js,fileConn)
  shiny::addResourcePath("bind",path)
}

shinyApp(ui,server,onStart)
