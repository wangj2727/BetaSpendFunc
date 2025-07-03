
library(shiny)
library(tidyverse)
library(ldbounds)
library(shinythemes)
library(shinycssloaders)
library(shinyjs)
library(bslib)
library(ggpubr)
library(roxygen2)


#' Hwang-Shih-DeCani family (gamma family) spending function
#'
#' @param bet Type I or II error rate depending on if lower or upper boundary is calculated
#' @param tt schedule of interim looks
#' @param gam gamma parameter used in Hwang-Shih-DeCani spending function
#' @return a vector of lower/upper boundaries
#' @examples hsd_spending(c(0.5, 0.75, 1),0.1, -2),
hsd_spending<-function(tt,bet,gam){
    bet*(1-exp(-gam*tt))/(1-exp(-gam))
}

#' O'brien-fleming-like spending function
#'
#' @param bet Type I or II error rate depending on if lower or upper boundary is calculated
#' @param tt schedule of interim looks
#' @return a vector of lower/upper boundaries
#' @examples obrien_fleming_spending(c(0.5, 0.75, 1), 0.1)
obrien_fleming_spending <- function(tt,bet){
    2 - 2 * pnorm(qnorm(1 - bet / 2) / sqrt(tt))
}

#' Pocock-like spending function
#'
#' @param bet Type I or II error rate depending on if lower or upper boundary is calculated
#' @param tt schedule of interim looks
#' @return a vector of lower/upper boundaries
#' @examples pocock_spending(c(0.5, 0.75, 1),0.1)
pocock_spending <- function(tt, bet){
    bet * log(1 + (exp(1) - 1) * tt)
}

#' Power family spending function
#'
#' @param bet Type I or II error rate depending on if lower or upper boundary is calculated
#' @param tt schedule of interim looks
#' @param phi phi parameter used in the power spending function
#' @return a vector of lower/upper boundaries
#' @examples power_spending(c(0.5, 0.75, 1), 0.1, 1)
power_spending<-function(tt,bet,phi){
    bet*(tt^phi)
}

#' Calculate sample size (per arm) needed to cross the boundary for continuous endpoint
#'
#' @param mean_diff difference in means between the two arms (trial setting: two arms with 1:1 allocation ratio)
#' @param sigma within group standard deviation (assume equal standard deviation between two arms)
#' @param theta the drift parameter, the expected value of Z(1) under the alternative hypothesis
#'
#' @return
#' @examples
SS_continous <- function(mean_diff, sigma1, sigma2, theta){
    ((sigma1^2 + sigma2^2)*theta^2)/mean_diff^2
    # sigma1 and sigma2 can be the same --> equal variable between arms
}


#' Calculate sample size (per arm) needed to
#'
#' @param p1 event rate in arm1
#' @param p2 event rate in arm2
#' @param theta the drift parameter, the expected value of Z(1) under the alternative hypothesis
#'
#' @return
#' @examples
SS_binary <- function(p1, p2, theta){
    ((p1*(1-p1) + p2*(1-p2))*theta^2)/(p1-p2)^2
}


#' Calculate number of events (overall) needed to
#'
#' @param hr hazard ratio
#' @param theta the drift parameter, the expected value of Z(1) under the alternative hypothesis
#' @return
#' @examples
SS_tte <- function(hr, theta){
    4*(theta/log(hr))^2
}

#### customized theme
custom_theme <- theme_bw(base_size = 16) +
    theme(
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14)
    )


ui <- fluidPage(
    # theme = bs_theme(
    #     version = 5,
    #     bootswatch = "minty",
    #     base_font = font_google("Poppins"),
    #     heading_font = font_google("Lora"),
    #     primary = "#2C3E50"
    # ),
    theme = bs_theme(
        bg = "#FFF",
        fg = "#2C3E50",
        primary = "#EE6363",
        secondary = "#0072B2",
        success = "#009E73",
        base_font = font_google("Inter"),
        code_font = font_google("JetBrains Mono")
    ),
    useShinyjs(),

    titlePanel(div(
        h1("📈 Beta Spending Boundary Explorer", style = "font-weight: bold;")
    )),
    hr(),

    sidebarLayout(
        sidebarPanel(
            h4("Spending Parameters", class = "mt-3"),
            numericInput("alpha", "Type I Error (\u03B1):", value = .025, min = .00001, max = .49999),
            numericInput("beta", "Type II Error (\u03B2):", value = .15, min = .00001, max = .9999),
            textInput("numLook", "Information Fractions (comma-separated):", value = "0.25, 0.5, 0.75, 1"),

            h4("Upper Spending Function"),
            selectInput("NumParaUpper", "# Parameters:", c("No Parameters" = "nopara", "1 Parameter" = "onepara")),
            conditionalPanel("input.NumParaUpper == 'nopara'",
                             selectInput("sf_nop_upper", "Function:", c("Pocock-like" = "pc", "O'Brien-Fleming-like" = "obf"))),
            conditionalPanel("input.NumParaUpper == 'onepara'",
                             selectInput("sf_onep_upper", "Function:", c("Hwang-Shih-DeCani" = "hsd", "Power family" = "pf")),
                             conditionalPanel("input.sf_onep_upper == 'hsd'",
                                              numericInput("gamma_upper", "Gamma (\u2260 0):", value = -2)),
                             conditionalPanel("input.sf_onep_upper == 'pf'",
                                              numericInput("phi_upper", "Phi (>0):", value = 1))),

            h4("Lower Spending Function"),
            selectInput("NumParaLower", "# Parameters:", c("No Parameters" = "nopara", "1 Parameter" = "onepara")),
            conditionalPanel("input.NumParaLower == 'nopara'",
                             selectInput("sf_nop_lower", "Function:", c("O'Brien-Fleming-like" = "obf","Pocock-like" = "pc"))),
            conditionalPanel("input.NumParaLower == 'onepara'",
                             selectInput("sf_onep_lower", "Function:", c("Hwang-Shih-DeCani" = "hsd", "Power family" = "pf")),
                             conditionalPanel("input.sf_onep_lower == 'hsd'",
                                              numericInput("gamma_lower", "Gamma (\u2260 0):", value = -2)),
                             conditionalPanel("input.sf_onep_lower == 'pf'",
                                              numericInput("phi_lower", "Phi (>0):", value = 1))),

            h4("Endpoint Assumption"),
            # selectInput("endpoint", "Endpoint type:",
            #             choices = c("Difference in Means" = "diffmean",
            #                         "Hazard Ratio" = "HR",
            #                         "Relative Risk" = "RR")),
            # conditionalPanel("input.endpoint=='diffmean'",
            #                  numericInput("meanDiff","Difference in Means:", value=0.5)),
            # conditionalPanel("input.endpoint=='HR'",
                             numericInput("HR", "Hazard Ratio:", value = 0.6, min = .001, max = .999)
            # conditionalPanel("input.endpoint=='RR'",
            #                  numericInput("RR","Relative Risk:", value=0.3))

        ),

        mainPanel(
            tabsetPanel(
                tabPanel("About",
                         h3("About This App"),
                         p("This app helps visualize and compute asymmetric group sequential boundaries using various beta spending functions."),
                         h4("How to Use"),
                         tags$ul(
                             tags$li("Enter the desired Type I and Type II error rates."),
                             tags$li("Specify the information fractions for interim looks (comma-separated)."),
                             tags$li("Choose the upper and lower spending functions (with parameters if needed)."),
                             tags$li("Enter the assumed hazard ratio for time-to-event endpoints."),
                             tags$li("Navigate to the 'Spending Function Plot' and 'Boundary Plot' tabs to view results."),
                             tags$li("Use the download buttons to export plots and tables.")
                         ),
                         h4("Notes"),
                         tags$ul(
                             tags$li("Either asymmetric or symmetric boundaries can be produced by this app."),
                             tags$li("Sample size is based on time-to-event outcome."),
                             tags$li("Spending functions control how error is spent over time.")
                         ),
                         h4("Other information"),
                         p("Built with ", tags$code("shiny"), ", using the ", tags$code("ldbounds"), " package."),
                         tags$a(href="https://github.com/wangj2727/BetaSpendFunc", "View source code on GitHub", target="_blank")

                ),

                tabPanel("Spending Function Plot",
                         selectInput("sf_plot_type", "Plot Type:",
                                     c("Cumulative spend plot" = "sfPlot_cum",
                                       "Proportion of spend plot" = "sfPlot_prop")),
                         conditionalPanel("input.sf_plot_type == 'sfPlot_cum'",
                                          div(class = "card shadow-sm p-3", plotOutput("sfPlotCum"))),
                         conditionalPanel("input.sf_plot_type == 'sfPlot_prop'",
                                          div(class = "card shadow-sm p-3", plotOutput("sfPlotProp"))),
                         downloadButton("download_spending_plot", "Download Spending Plot (PNG)")
                ),

                tabPanel("Boundary Plot",
                         verbatimTextOutput("driftText"),
                         selectInput("bd_plot_type", "Plot Type:",
                                     c("Boundary plot (Z-score)" = "bdPlot_z",
                                       "Boundary plot (HR scale)" = "bdPlot_hr")),
                         conditionalPanel("input.bd_plot_type == 'bdPlot_z'",
                                          div(class = "card shadow-sm p-3", plotOutput("boundPlot"))),
                         conditionalPanel("input.bd_plot_type == 'bdPlot_hr'",
                                          div(class = "card shadow-sm p-3", plotOutput("boundPlotHR"))),
                         downloadButton("download_boundary_plot", "Download Boundary Plot (PNG)"),
                         downloadButton("download_table", "Download Boundary Table (CSV)")
                )
            )
        )
    )
)




server <- function(input, output) {

    ### assign parameters
    looktimes<-reactive(as.numeric(unlist(str_split(input$numLook, ","))))
    alpha <- reactive(input$alpha)
    beta <- reactive(input$beta)

    sf_temp <- reactive({
        sf_use_lower <- if(input$NumParaLower == "nopara" & input$sf_nop_lower == "pc") {
            pocock_spending(looktimes(),beta())
        }else if(input$NumParaLower == "nopara" & input$sf_nop_lower == "obf") {
            obrien_fleming_spending(looktimes(),beta())
        }else if(input$NumParaLower == "onepara" & input$sf_onep_lower == "hsd") {
            hsd_spending(looktimes(),beta(), input$gamma_lower)
        }else if(input$NumParaLower == "onepara" & input$sf_onep_lower == "pf") {
            power_spending(looktimes(),beta(),input$phi_lower)
        }

        sf_use_upper <- if(input$NumParaUpper == "nopara" & input$sf_nop_upper == "pc") {
            pocock_spending(looktimes(),alpha())
        }else if(input$NumParaUpper == "nopara" & input$sf_nop_upper == "obf") {
            obrien_fleming_spending(looktimes(),alpha())
        }else if(input$NumParaUpper == "onepara" & input$sf_onep_upper == "hsd") {
            hsd_spending(looktimes(), alpha(),  input$gamma_upper)
        }else if(input$NumParaUpper == "onepara" & input$sf_onep_upper == "pf") {
            power_spending(looktimes(),alpha(),  input$phi_upper)
        }

        data.frame(tt=looktimes(),sf_beta=sf_use_lower, sf_alpha=sf_use_upper)%>%
            bind_rows(data.frame(tt=0, sf_beta=0, sf_alpha=0))%>%
            pivot_longer(cols=sf_beta:sf_alpha)%>%
            arrange(name, tt)%>%
            mutate(name = str_remove(name, "sf_"),
                   prop_spend = ifelse(name=="alpha", round(value/alpha(),2), round(value/beta(),2)))

    })

    output$sfPlotCum <- renderPlot({

        ggplot(sf_temp(),aes(x=tt, y=value, color = name))+
            geom_line(size=1)+
            geom_point(size=2, shape=19)+
            geom_text(aes(label = round(value,3)), color="black", vjust=1.5)+
            theme_bw()+
            custom_theme +
            theme(legend.title=element_blank())+
            xlab("t")+ylab("Cumulative spent")+
            ggtitle("Cumulative spending by t:")
    })

    output$sfPlotProp <- renderPlot({
        ggplot(sf_temp(),aes(x=tt, y=prop_spend, color = name))+
            geom_line(size=1, alpha=0.5)+
            geom_point(size=2, shape=19)+
            geom_text(aes(label = round(prop_spend,3)), color="black", vjust=1.5)+
            theme_bw()+
            custom_theme +
            theme(legend.title=element_blank())+
            labs(x="t", y = "Proportion of spend",
                 title = expression(bold("Proportion of spending by t")),
                 subtitle = "Note: if Pocock-like spending function is selected, the two curves overlap.")
    })


    resultData <- reactive({
            tvec <- as.numeric(unlist(strsplit(input$numLook, ",")))
            if (any(is.na(tvec)) || !is.numeric(tvec) || any(diff(tvec) <= 0) || any(tvec <= 0 | tvec > 1)) {
                showNotification("Invalid information fractions. Use increasing values between 0 and 1.", type = "error")
                return(NULL)
            }

            alpha <- input$alpha
            bet <- input$beta
            k <- length(tvec)

            ### Setup parameters for spending function for lower bound
            if(input$NumParaLower == "nopara" & input$sf_nop_lower == "obf") {
                iuse_lower <- 1
                phi_val_lower <- NULL
                spend_func_lower <- function(tt,bet){
                    2 - 2 * pnorm(qnorm(1 - bet / 2) / sqrt(tt))
                }
            } else if (input$NumParaLower == "nopara" & input$sf_nop_lower == "pc")  {
                iuse_lower <- 2
                phi_val_lower <- NULL
                spend_func_lower <- function(tt, bet){
                    bet * log(1 + (exp(1) - 1) * tt)
                }
            } else if (input$NumParaLower == "onepara" & input$sf_onep_lower == "hsd") {
                iuse_lower <- 4
                phi_val_lower <- input$gamma_lower
                spend_func_lower <- function(tt,bet,phi=phi_val_lower){
                    bet*(1-exp(-phi*tt))/(1-exp(-phi))
                }
            } else if(input$NumParaLower == "onepara" & input$sf_onep_lower == "pf") {
                iuse_lower <- 3
                phi_val_lower <- input$phi_lower
                spend_func_lower <- function(tt,bet,phi=phi_val_lower){
                    bet*(tt^phi)
                }

            } else stop("Unknown spending function")


            ### Setup parameters for spending function for upper bound
            if(input$NumParaUpper == "nopara" & input$sf_nop_upper == "obf") {
                iuse_upper <- 1
                phi_val_upper <- NULL
                spend_func_upper <- function(tt,bet){
                    2 - 2 * pnorm(qnorm(1 - bet / 2) / sqrt(tt))
                }
            } else if (input$NumParaUpper == "nopara" & input$sf_nop_upper == "pc")  {
                iuse_upper <- 2
                phi_val_upper <- NULL
                spend_func_upper <- function(tt, bet){
                    bet * log(1 + (exp(1) - 1) * tt)
                }
            } else if (input$NumParaUpper == "onepara" & input$sf_onep_upper == "hsd") {
                iuse_upper <- 4
                phi_val_upper <- input$gamma_upper
                spend_func_upper <- function(tt,bet,phi=phi_val_upper){
                    bet*(1-exp(-phi*tt))/(1-exp(-phi))
                }
            } else if(input$NumParaUpper == "onepara" & input$sf_onep_upper == "pf") {
                iuse_upper <- 3
                phi_val_upper <- input$phi_upper
                spend_func_upper <- function(tt,bet,phi=phi_val_upper){
                    bet*(tt^phi)
                }

            } else stop("Unknown spending function")

            # Upper bounds
            if (!is.null(phi_val_upper)) {
                uvec <- ldBounds(tvec, iuse = iuse_upper, alpha = alpha, sides = 1, phi = phi_val_upper)$upper.bounds
                cc <- ldBounds(tvec, iuse = iuse_upper, alpha = bet, sides = 1, phi = phi_val_upper)$upper.bounds
            } else {
                uvec <- ldBounds(tvec, iuse = iuse_upper, alpha = alpha, sides = 1)$upper.bounds
                cc <- ldBounds(tvec, iuse = iuse_upper, alpha = bet, sides = 1)$upper.bounds
            }

            lowerlim <- qnorm(1 - alpha) + qnorm(1 - bet)
            upperlim <- cc[k] + uvec[k]

            firsticrosslow <- function(li, ti, tprevvec, lprevvec, uprevvec, theta) {
                if (length(tprevvec) == 0) return("Do not apply to 1st t")
                tcurrvec <- c(tprevvec, ti)
                lcurrvec <- c(lprevvec, li)
                ucurrvec <- c(uprevvec, 12)
                ldPower(tcurrvec, za = lcurrvec - theta * sqrt(tcurrvec),
                        zb = ucurrvec - theta * sqrt(tcurrvec))$power -
                    ldPower(tprevvec, za = lprevvec - theta * sqrt(tprevvec),
                            zb = uprevvec - theta * sqrt(tprevvec))$power
            }

            diffy <- function(li, ii, ti, tprevvec, lprevvec, uprevvec, theta, bet) {
                    firsticrosslow(li, ti, tprevvec, lprevvec, uprevvec, theta) -
                        (spend_func_lower(tvec, bet)[ii] - spend_func_lower(tvec, bet)[ii - 1])
            }

            findlowerbds <- function(theta, tvec, uvec, bet) {
                tprevvec <- c(); lprevvec <- c(); uprevvec <- c(); lvec <- numeric(length(tvec))
                k <- length(tvec)
                lvec[1] <- qnorm(spend_func_lower(tvec, bet)[1], mean = theta * sqrt(tvec[1]), sd = 1)

                lvec[k] <- uvec[k]
                if (k > 2) {
                    for (i in 2:(k - 1)) {
                        ti <- tvec[i]
                        tprevvec <- c(tprevvec, tvec[i - 1])
                        lprevvec <- c(lprevvec, lvec[i - 1])
                        uprevvec <- c(uprevvec, uvec[i - 1])
                        lvec[i] <- uniroot(diffy, ii = i, ti = ti,
                                           tprevvec = tprevvec, lprevvec = lprevvec,
                                           uprevvec = uprevvec, theta = theta, bet = bet,
                                           lower = -5, upper = 5)$root
                    }
                }
                return(lvec)
            }

            findtheta <- function(theta, bet, tvec, uvec) {
                k <- length(tvec)
                lvec <- findlowerbds(theta, tvec, uvec, bet)
                tprevvec <- tvec[1:(k - 1)]
                lprevvec <- lvec[1:(k - 1)]
                uprevvec <- uvec[1:(k - 1)]
                li <- lvec[k]; ti <- tvec[k]; ii <- k
                diffy(li, ii, ti, tprevvec, lprevvec, uprevvec, theta, bet)
            }

            rightdrift <- uniroot(findtheta, bet = bet, tvec = tvec, uvec = uvec,
                                  lower = lowerlim, upper = upperlim)$root
            lvec <- findlowerbds(rightdrift, tvec, uvec, bet)

            SS <- SS_tte(input$HR, unique(rightdrift))

            df <- data.frame(Stage = 1:k, InfoFrac = tvec,
                             Z_UpperBound = round(uvec, 3), Z_LowerBound = round(lvec, 3),SampleSize=SS)%>%
                mutate(HR_LowerBound = round(exp(lvec/sqrt(SS*tvec/4)),3),
                       HR_UpperBound = round(exp(uvec/sqrt(SS*tvec/4)),3))

            dtaplot <- data.frame(lvec, uvec, tt=tvec)%>%
                mutate_all(round,2)%>%
                pivot_longer(cols = lvec:uvec)%>%
                mutate(name= case_when(name=="lvec" ~ "Lower boundary",
                                       name=="uvec" ~ "Upper boundary",
                                       TRUE ~ name))

            dtaplot_HR <- data.frame(HR_LowerBound=df$HR_LowerBound,HR_UpperBound=df$HR_UpperBound, tt=tvec)%>%
                mutate_all(round,2)%>%
                pivot_longer(cols = HR_LowerBound:HR_UpperBound)%>%
                mutate(name= case_when(name=="HR_LowerBound" ~ "Lower boundary",
                                       name=="HR_UpperBound" ~ "Upper boundary",
                                       TRUE ~ name))

            list(dtaplot = dtaplot,dtaplot_HR=dtaplot_HR,
                 df = df, theta = rightdrift, tvec = tvec, uvec = uvec, lvec = lvec)
        })

        output$driftText <- renderPrint({
            res <- resultData()
            if (is.null(res)) return()
            cat("Estimated drift (θ):", round(res$theta, 4), "\n\n")
            print(res$df)
        })

        output$boundPlot <- renderPlot({
            res <- resultData()
            if (is.null(res)) return()

            ggplot(data=res$dtaplot,aes(x=tt, y=value, color=name))+
                geom_line(size=1)+
                geom_point(size=2, shape=19)+
                geom_text(aes(label = value), color="black", vjust=1.5)+
                theme_bw()+
                custom_theme +
                theme(legend.title=element_blank())+
                labs(x="t", y="Boundaries", title="Lower/Upper Z-score boundary plot ")
        })

        output$boundPlotHR <- renderPlot({
            res <- resultData()
            if (is.null(res)) return()

            ggplot(data=res$dtaplot_HR,aes(x=tt, y=value, color=name))+
                geom_line(size=1)+
                geom_point(size=2, shape=19)+
                geom_text(aes(label = value), color="black", vjust=1.5)+
                theme_bw()+
                custom_theme +
                theme(legend.title=element_blank())+
                labs(x="t", y="Boundaries", title="Lower/Upper boundary plot (HR Scale)")
        })

        output$download_spending_plot <- downloadHandler(
            filename = function() {
                type <- input$sf_plot_type
                if (type == "sfPlot_cum") {
                    paste("cumulative_spending_plot", Sys.Date(), ".png", sep = "")
                } else {
                    paste("proportion_spending_plot", Sys.Date(), ".png", sep = "")
                }
            },
            content = function(file) {
                type <- input$sf_plot_type
                png(file, width = 800, height = 600)

                if (type == "sfPlot_cum") {
                    print(
                        ggplot(sf_temp(), aes(x = tt, y = value, color = name)) +
                            geom_line(size = 1) +
                            geom_point(size = 2, shape = 19) +
                            geom_text(aes(label = round(value, 3)), color = "black", vjust = 1.5) +
                            theme_bw() +
                            custom_theme +
                            theme(legend.title = element_blank()) +
                            xlab("t") + ylab("Cumulative spent") +
                            ggtitle("Cumulative spending by t:")
                    )
                } else if (type == "sfPlot_prop") {
                    print(
                        ggplot(sf_temp(), aes(x = tt, y = prop_spend, color = name)) +
                            geom_line(size = 1, alpha = 0.5) +
                            geom_point(size = 2, shape = 19) +
                            geom_text(aes(label = round(prop_spend, 3)), color = "black", vjust = 1.5) +
                            theme_bw() +
                            custom_theme +
                            theme(legend.title = element_blank()) +
                            labs(
                                x = "t", y = "Proportion of spend",
                                title = expression(bold("Proportion of spending by t")),
                                subtitle = "Note: if Pocock-like spending function is selected, the two curves overlap."
                            )
                    )
                }

                dev.off()
            }
        )
        output$download_boundary_plot <- downloadHandler(
            filename = function() {
                type <- input$bd_plot_type
                if (type == "bdPlot_z") {
                    paste("boundary_plot_Z", Sys.Date(), ".png", sep = "")
                } else {
                    paste("boundary_plot_HR", Sys.Date(), ".png", sep = "")
                }
            },
            content = function(file) {
                res <- resultData()
                if (is.null(res)) return()

                png(file, width = 800, height = 600)

                if (input$bd_plot_type == "bdPlot_z") {
                    print(
                        ggplot(data = res$dtaplot, aes(x = tt, y = value, color = name)) +
                            geom_line(size = 1) +
                            geom_point(size = 2, shape = 19) +
                            geom_text(aes(label = value), color = "black", vjust = 1.5) +
                            theme_bw() +
                            theme(legend.title = element_blank()) +
                            labs(x = "t", y = "Boundaries", title = "Lower/Upper Z-score boundary plot")
                    )
                } else if (input$bd_plot_type == "bdPlot_hr") {
                    print(
                        ggplot(data = res$dtaplot_HR, aes(x = tt, y = value, color = name)) +
                            geom_line(size = 1) +
                            geom_point(size = 2, shape = 19) +
                            geom_text(aes(label = value), color = "black", vjust = 1.5) +
                            theme_bw() +
                            theme(legend.title = element_blank()) +
                            labs(x = "t", y = "Boundaries", title = "Lower/Upper boundary plot (HR Scale)")
                    )
                }

                dev.off()
            }
        )


        output$download_table <- downloadHandler(
            filename = function() {
                paste("boundary_table", Sys.Date(), ".csv", sep = "")
            },
            content = function(file) {
                res <- resultData()
                if (is.null(res)) return()
                write.csv(res$df, file, row.names = FALSE)
            }
        )
}


shinyApp(ui = ui, server = server)
