# Compute regression of time picks
# library(sm)

# Define the files 
ifn <- 'pick_iso_iso.txt'
afn <- 'pick_iso_aniso.txt'
tfn <- 'pick_thrm_iso.txt'

# Load the data
idat <- read.delim(ifn, header = TRUE, sep = " ")
adat <- read.delim(afn, header = TRUE, sep = " ")
tdat <- read.delim(tfn, header = TRUE, sep = " ")

vals_assign <- function(val_array)
{
  x <- val_array$x 
  z <- val_array$z 
  
  t0 <- val_array$t0
  t1 <- val_array$t1
  t2 <- val_array$t2
  t3 <- val_array$t3
  t4 <- val_array$t4
  t5 <- val_array$t5
  t6 <- val_array$t6
  t7 <- val_array$t7
  t8 <- val_array$t8
  t9 <- val_array$t9
  
  t0 <- cbind(x[t0 != -12345], z[t0 != -12345], t0[t0 != -12345] )
  t1 <- cbind(x[t1 != -12345], z[t1 != -12345], t1[t1 != -12345] )
  t2 <- cbind(x[t2 != -12345], z[t2 != -12345], t2[t2 != -12345] )
  t3 <- cbind(x[t3 != -12345], z[t3 != -12345], t3[t3 != -12345] )
  t4 <- cbind(x[t4 != -12345], z[t4 != -12345], t4[t4 != -12345] )
  t5 <- cbind(x[t5 != -12345], z[t5 != -12345], t5[t5 != -12345] )
  t6 <- cbind(x[t6 != -12345], z[t6 != -12345], t6[t6 != -12345] )
  t7 <- cbind(x[t7 != -12345], z[t7 != -12345], t7[t7 != -12345] )
  t8 <- cbind(x[t8 != -12345], z[t8 != -12345], t8[t8 != -12345] )
  t9 <- cbind(x[t9 != -12345], z[t9 != -12345], t9[t9 != -12345] )
  
  # Sort according to X
  t0 <- t0[ sort(t0[,1], index.return=T)$ix, ]
  t1 <- t1[ sort(t1[,1], index.return=T)$ix, ]
  t2 <- t2[ sort(t2[,1], index.return=T)$ix, ]
  t3 <- t3[ sort(t3[,1], index.return=T)$ix, ]
  t4 <- t4[ sort(t4[,1], index.return=T)$ix, ]
  t5 <- t5[ sort(t5[,1], index.return=T)$ix, ]
  t6 <- t6[ sort(t6[,1], index.return=T)$ix, ]
  
  return(list(t0=t0, t1=t1, t2=t2, t3=t3, t4=t4, t5=t5, t6=t6, t7=t7, t8=t8, t9=t9) )
  
}

vals_regress <- function(dat, newx)
{
  x <- dat[,1] 
  t <- dat[,3]
  
  x2 <- x^2
  mod <- lm(t ~ x + x2 )
  
  # Predict the points
  newdata <- predict(mod, list(x=newx, x2=newx^2))
  
  # Get the residual
  r <- t - predict(mod, list(x = x, x2 = x^2))
  
  return(list( reg=mod, predicted=newdata, resid=r) )
  
  
}


idat_picks <- vals_assign(idat)
tdat_picks <- vals_assign(tdat)
adat_picks <- vals_assign(adat)

# Plot pick comparisons

# dev.new(width=5, height = 9, units = "in")
png("TimePickComparison_layered.png", width = 5, height = 9, units = "in", res = 100)
  par(mai = c(1, 1, 0.1, 0.1), bg = NA)
  plot( idat_picks$t0[,1], idat_picks$t0[,3], type = "l", ylim = c(6e-6, 0), cex.lab = 1.5,
        xlab = "Distance (m)", ylab = "TWTT (s)", col = "black", lty = 2, lwd = 2)

    
  lines( adat_picks$t1[,1], adat_picks$t1[,3], col = "red", lty = 1, lwd = 3)
  lines( adat_picks$t2[,1], adat_picks$t2[,3], col = "red", lty = 1, lwd = 3)
  lines( adat_picks$t3[,1], adat_picks$t3[,3], col = "red", lty = 1, lwd = 3)
  lines( adat_picks$t4[,1], adat_picks$t4[,3], col = "red", lty = 1, lwd = 3)
  lines( adat_picks$t5[,1], adat_picks$t5[,3], col = "red", lty = 1, lwd = 3)
  lines( adat_picks$t6[,1], adat_picks$t6[,3], col = "red", lty = 1, lwd = 3)
  
  lines( tdat_picks$t1[,1], tdat_picks$t1[,3], col = "blue", lty = 1, lwd = 2)
  lines( tdat_picks$t2[,1], tdat_picks$t2[,3], col = "blue", lty = 1, lwd = 2)
  lines( tdat_picks$t3[,1], tdat_picks$t3[,3], col = "blue", lty = 1, lwd = 2)
  lines( tdat_picks$t4[,1], tdat_picks$t4[,3], col = "blue", lty = 1, lwd = 2)
  lines( tdat_picks$t5[,1], tdat_picks$t5[,3], col = "blue", lty = 1, lwd = 2)
  lines( tdat_picks$t6[,1], tdat_picks$t6[,3], col = "blue", lty = 1, lwd = 2)
  
  lines( idat_picks$t1[,1], idat_picks$t1[,3], col = "black", lty = 2, lwd = 1)
  lines( idat_picks$t2[,1], idat_picks$t2[,3], col = "black", lty = 2, lwd = 1)
  lines( idat_picks$t3[,1], idat_picks$t3[,3], col = "black", lty = 2, lwd = 1)
  lines( idat_picks$t4[,1], idat_picks$t4[,3], col = "black", lty = 2, lwd = 1)
  lines( idat_picks$t5[,1], idat_picks$t5[,3], col = "black", lty = 2, lwd = 1)
  lines( idat_picks$t6[,1], idat_picks$t6[,3], col = "black", lty = 2, lwd = 1)
dev.off()

# Compute the regression
newx <- seq(20, 560, 10)
it0 <- vals_regress(idat_picks$t0, newx)
tt0 <- vals_regress(tdat_picks$t0, newx)
at0 <- vals_regress(adat_picks$t0, newx)

it1 <- vals_regress(idat_picks$t1, newx)
tt1 <- vals_regress(tdat_picks$t1, newx)
at1 <- vals_regress(adat_picks$t1, newx)

it2 <- vals_regress(idat_picks$t2, newx)
tt2 <- vals_regress(tdat_picks$t2, newx)
at2 <- vals_regress(adat_picks$t2, newx)

it3 <- vals_regress(idat_picks$t3, newx)
tt3 <- vals_regress(tdat_picks$t3, newx)
at3 <- vals_regress(adat_picks$t3, newx)

it4 <- vals_regress(idat_picks$t4, newx)
tt4 <- vals_regress(tdat_picks$t4, newx)
at4 <- vals_regress(adat_picks$t4, newx)

it5 <- vals_regress(idat_picks$t5, newx)
tt5 <- vals_regress(tdat_picks$t5, newx)
at5 <- vals_regress(adat_picks$t5, newx)

it6 <- vals_regress(idat_picks$t6, newx)
tt6 <- vals_regress(tdat_picks$t6, newx)
at6 <- vals_regress(adat_picks$t6, newx)


# Plot the results
res_comp <- function(dat1, dat2, dat3, fn)
{
  #dev.new(width = 4, height = 4, units = "in")
  png(filename = fn, width = 4, height = 4, units = "in", res = 96 )
  par( mai = c(1, 1, 0.1, 0.1), bg = NA )
  d1 <- density(dat1 - dat2)
  d2 <- density(dat1 - dat3)
  xlimits <- c( min(d1$x, d2$x), max(d1$x, d2$x))
  ylimits <- c( 0, max(d1$y, d2$y) )
  
  
  plot( d1, col = "red", type = "l", main = "", ylim = ylimits, xlim = xlimits, 
        xlab = "Time (s)", ylab = "Frequency", cex.lab = 1.5, lwd = 3)
  lines(d2, col = "blue", lwd = 3)
  dev.off()
}


res_comp(it0$predicted, at0$predicted, tt0$predicted, "time_picks0.png")
res_comp(it1$predicted, at1$predicted, tt1$predicted, "time_picks1.png")
res_comp(it2$predicted, at2$predicted, tt2$predicted, "time_picks2.png")
res_comp(it3$predicted, at3$predicted, tt3$predicted, "time_picks3.png")
res_comp(it4$predicted, at4$predicted, tt4$predicted, "time_picks4.png")
res_comp(it5$predicted, at5$predicted, tt5$predicted, "time_picks5.png")
res_comp(it6$predicted, at6$predicted, tt6$predicted, "time_picks6.png")



mu0 <- 4*pi*1e-7 
ep0 <- 8.85418e-12

vel1 <- 1/sqrt(mu0*ep0*3.15)
vel2 <- 1/sqrt(mu0*ep0*3.17)
twtt <- abs(sqrt(350^2 + 250^2)/vel1 - sqrt(350^2 + 250^2)/vel2)
dt <- 1.6678204759907602e-09
