# Load libraries
library("phom")
library("pracma")

# Set working directory, read data
X <- read.delim("./data/Real-world-data/Vicsek__particles_300_distance_1_noise_2_v0_0.03_box_7_timestep_600_of_600.txt", sep = "", header = FALSE)
X <- data.frame(X)
colnames(X) <- c("x", "y", "angle")
# Subset
# X <- subset(X, t == thistime)

# Calculate angles
X$angle <- Arg(complex(real = X$x, imaginary = X$y))

# Calculate scaled position
lbox <- 2 * pi
X$x <- X$x / lbox * 2 * pi
X$y <- X$y / lbox * 2 * pi

ptime <- system.time({
    # Construct distance matrix
    Xsub <- X
    Xsub <- subset(Xsub, select = c(x, y, angle))
    xgrid <- meshgrid(Xsub$x, Xsub$x)
    ygrid <- meshgrid(Xsub$y, Xsub$y)
    anglegrid <- meshgrid(Xsub$angle, Xsub$angle)
    M1 <- pmin(abs(xgrid$X - xgrid$Y), lbox - abs(xgrid$X - xgrid$Y))^2
    M2 <- pmin(abs(ygrid$X - ygrid$Y), lbox - abs(ygrid$X - ygrid$Y))^2
    M3 <- pmin(abs(anglegrid$X - anglegrid$Y), 2 * pi - abs(anglegrid$X - anglegrid$Y))^2
    M <- sqrt(M1 + M2 + M3)
})[3]
write.table(M,
    file = "vicsek_dist.txt", sep = "\t",
    row.names = FALSE
)
