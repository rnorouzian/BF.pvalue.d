
D = read.csv("https://raw.githubusercontent.com/rnorouzian/BF.pvalue.d/master/BFdPvalue.csv")

######### Histogram Panel function for pairs() ###########

panel.hist = function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}

#pairs(USJudgeRatings[1:5], panel = panel.smooth,
 #     cex = 1.5, pch = 24, bg = "light blue",
  #    diag.panel = panel.hist, cex.labels = 2, font.labels = 2)

######### Correlation Panel Function for pairs() #########

panel.cor = function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor) #cex = cex.cor*r
}
#pairs(USJudgeRatings, lower.panel = panel.smooth, upper.panel = panel.cor)


####################################################################################
## Function to economically compute BF10, Cohen's d, and p.value from the literature
####################################################################################


BF.d.pvalue = Vectorize(function(t, n1, n2 = NA, scale = sqrt(2)/2, log = FALSE, prior = "Cauchy"){
  
   options(warn = -1)  
   t = abs(t)
   N = ifelse(is.na(n2), n1, (n1*n2)/(n1+n2))
  df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
   d = t / sqrt(N)
  
  H1 = integrate(function(delta)ifelse(prior == "Cauchy", dcauchy(delta, 0, scale), dnorm(delta) )*dt(t, df, delta*sqrt(N)), -Inf, Inf)[[1]]
  H0 = dt(t, df)
  BF10 = ifelse(log, log(H1/H0), H1/H0)
  p.value = 2*(1-pt(t, df))
  
  cbind(BF10 = BF10, p.value = p.value, d = d)
  
}, vectorize.args = c("t", "n1", "n2", "scale", "log", "prior"))


b = BF.d.pvalue(t = D$t.value, n1 = D$n1, n2 = D$n2)

BF10 = b[1, ]  ;  p.value = b[2, ]   ;  cohen.d = b[3, ]

######### General Plots:

plot(p.value, cohen.d, ylim = c(0, 1.25), pch = 21, bg = 3, 
     xlab = "P-value", ylab = "Cohen's d", font = 2, font.lab = 2, las = 1)

plot(cohen.d, BF10, xlim = c(0, 1.2), ylim = c(.1, 10), pch = 21, bg = 3, 
     xlab = "Cohen's d", ylab = "BF10", font = 2, font.lab = 2, las = 1)

plot(p.value, BF10, ylim = c(0, 20), xlim = c(0, 1), pch = 21, bg = 3, 
     xlab = "P-value", ylab = "BF10", font = 2, font.lab = 2, las = 1)


## Plot with histogram:
pairs(BF10~p.value + cohen.d, data = b, gap = .2, t = "n", diag.panel = panel.hist,
      lower.panel = NULL, las = 1, panel = panel.smooth, font = 2,
      labels = c("BF10","p-value", "Cohen's d"), col = "gray70" )

## Plot with correlation:
pairs(BF10~p.value + cohen.d, data = b, gap = .2, t = "n", las = 1, 
      labels = c("BF10","p-value", "Cohen's d"), 
      lower.panel = panel.smooth, upper.panel = panel.cor, font = 2, col = 4,
      bg = 3 )
