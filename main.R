# Team Name:               Team Members:               Date:
# Flying Monkeys           Robbert-Jan Joling          14-01-2015
#                          Damiano Luzzi

# Import packages ---------------------------------------------------------
library(raster)


# Load data ---------------------------------------------------------------
load("data/vcfGewata.rda")
load("data/GewataB1.rda")
load("data/GewataB2.rda")
load("data/GewataB3.rda")
load("data/GewataB4.rda")
load("data/GewataB5.rda")
load("data/GewataB7.rda")
load("data/trainingPoly.rda")


# Calculations ------------------------------------------------------------
gewata <- brick(GewataB1, GewataB2, GewataB3, GewataB4, GewataB5, GewataB7)
vcfGewata[vcfGewata > 100] <- NA
covs <- addLayer(gewata, vcfGewata)
names(covs) <- c("band1", "band2", "band3", "band4", "band5", "band7", "VCF")

pairs(covs)

valuetable <- getValues(covs)
valuetable <- na.omit(valuetable)
valuetable <- as.data.frame(valuetable)

lm.test <- lm(formula = vcf2000Gewata ~ gewataB1 + gewataB2 + gewataB3 + gewataB5 + gewataB7, data = valuetable)
summary(lm.test)


names(valuetable)
predVCF <- predict(covs, model=lm.test, na.rm = T)
predVCF[predVCF < 0 | predVCF > 100] <- NA

difference <- overlay(vcfGewata, predVCF, fun = function(x, y){(x-y)}, overwrite = T,
                      filename = "output/TCDiffPredictedVSActual")
diff.squared <- difference ^ 2
diff.mean <- cellStats(diff.squared, stat = 'mean')
RMSE <- sqrt(diff.mean)


trainingPoly@data$Code <- as.numeric(trainingPoly@data$Class)
classes <- rasterize(trainingPoly, predVCF, field='Code')

covs.VCF <- addLayer(predVCF, vcfGewata, classes)
names(covs.VCF) <- c("predVCF", "VCF", "class")

covs.VCF.masked <- mask(covs.VCF, classes)
difference.covs <- overlay(covs.VCF.masked$predVCF, covs.VCF.masked$VCF, fun = function(x, y){(x-y)})
diff.cov.squared <- difference.covs^2
diff.cov.mean <- zonal(diff.cov.squared, covs.VCF.masked$class, fun = 'mean')
RMSE.cov <- sqrt(diff.cov.mean[, 2])


# Visualisation -----------------------------------------------------------
opar <- par(mfrow = c(1,3))
plot(predVCF, main = "Predicted")
plot(vcfGewata, main = "Actual")
plot(difference, main = "Difference")
par(mfrow = c(1,1))



