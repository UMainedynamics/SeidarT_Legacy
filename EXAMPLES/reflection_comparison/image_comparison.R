library(RColorBrewer)


# Define the files 
ifn <- 'isothermal_isotropic/reciever_array.csv'
afn <- 'isothermal_anisotropic/reciever_array.csv'
tfn <- 'thermal_isotropic/reciever_array.csv'

# Load the data
idat <- t(as.matrix(read.delim(ifn, header = FALSE, sep = ",")))
adat <- t(as.matrix(read.delim(afn, header = FALSE, sep = ",")))
tdat <- t(as.matrix(read.delim(tfn, header = FALSE, sep = ",")))


dt <- 1.6678204759907602e-09
dx <- 10
xi <- 20
xf <- 580

m <- dim(idat)[1]
n <- dim(idat)[2]

time_vector <- seq(1,n)*dt 
x_vector <- seq(0, m-1)*dx + xi 

gain_function <- seq(1, n)^2

for (i in 1:m) 
{
    idat[i,] <- idat[i,]*gain_function
    tdat[i,] <- tdat[i,]*gain_function
    adat[i,] <- adat[i,]*gain_function
    
}

grays <- colorRampPalette(c("grey100", "grey0"))
png("Isotropic_t_x_plot_layered.png", width = 5, height = 9, units = "in", res = 100)
# dev.new(width = 5, height = 9, units = "in")
par(mai = c(1, 1, 0.1, 0.1), bg = NA)

image(x_vector, time_vector, idat, col = grays(100), ylim = c(max(time_vector), min(time_vector) ),
      ylab = "TWTT (s)", xlab = "Source-Reciever Offset (m)")

dev.off()

png("Polythermal_t_x_plot_layered.png", width = 5, height = 9, units = "in", res = 100)
# dev.new(width = 5, height = 9, units = "in")
par(mai = c(1, 1, 0.1, 0.1), bg = NA)

image(x_vector, time_vector, tdat, col = grays(100), ylim = c(max(time_vector), min(time_vector) ),
      ylab = "TWTT (s)", xlab = "Source-Reciever Offset (m)")

dev.off()

png("Anisotropic_t_x_plot_layered.png", width = 5, height = 9, units = "in", res = 100)
# dev.new(width = 5, height = 9, units = "in")
par(mai = c(1, 1, 0.1, 0.1), bg = NA)

image(x_vector, time_vector, adat, col = grays(100), ylim = c(max(time_vector), min(time_vector) ),
      ylab = "TWTT (s)", xlab = "Source-Reciever Offset (m)")

dev.off()




pal <- brewer.pal(11, "RdBu")
dia <- idat - adat
png("iso_aniso_diff.png", width = 5, height = 9, units = "in", res = 100)
par(mai = c(1, 1, 0.1, 0.1), bg = NA)
image(x_vector, time_vector, dia, col = pal, ylim = c(max(time_vector), min(time_vector) ),
      ylab = "TWTT (s)", xlab = "Source-Reciever Offset (m)")
dev.off()

dip <- idat - tdat

png("poly_aniso_diff.png", width = 5, height = 9, units = "in", res = 100)
par(mai = c(1, 1, 0.1, 0.1), bg = NA)
image(x_vector, time_vector, dip, col = pal, ylim = c(max(time_vector), min(time_vector) ),
      ylab = "TWTT (s)", xlab = "Source-Reciever Offset (m)")
dev.off()


