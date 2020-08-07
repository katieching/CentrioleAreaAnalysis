###Katie Ching
###last updated 2020-08-07
###STATISTICAL TESTS FOR AREA ANALYSIS OF CENTRIOLE STRUCTURES
###---###---###---###---###---###---###---###---###---###---###

#SET UP THE WORKSPACE
##set the working directory
setwd("C:/Users/Katie/Desktop/measuring_centriole_structures")


#MAKE NUMERIC VECTORS OF DATA VALUES
##Create a data frame from the csv file with all area measurements.
##You can access a column, such as groups, using something like AllAreas$group.
AllAreas <- read.csv("C:/Users/Katie/Desktop/S1_Data.csv", header=TRUE, sep = ",", dec = ".")

##Make a list for each type of structure (eg. pairs in culture).
##Note: These are each levels in AllAreas$group
CulturedPairs <- c()
CulturedRosettes <- c()
OEPairs <- c()
OENonPairs <- c()
OES <- c()

for (id in AllAreas$id) {
  #print(AllAreas$id[id])
  if (isTRUE(AllAreas$group[id] == "pair")){
    area = AllAreas$structure.area[id]
    CulturedPairs <- c(CulturedPairs, area)
  }
  if (isTRUE(AllAreas$group[id] == "rosette")){
    area = AllAreas$structure.area[id]
    CulturedRosettes <- c(CulturedRosettes, area)
  }
  if (isTRUE(AllAreas$group[id] == "OE; mitosis pairs")){
    area = AllAreas$structure.area[id]
    OEPairs <- c(OEPairs, area)
  }
  if (isTRUE(AllAreas$group[id] == "OE; mitosis >2")){
    area = AllAreas$structure.area[id]
    OENonPairs <- c(OENonPairs, area)
  }
  if (isTRUE(AllAreas$group[id] == "OE; S phase")){
    area = AllAreas$structure.area[id]
    OES <- c(OES, area)
  }
}

###---###---###---###---###---###---###---###---###---###---###---###

#FIND VALUES NECESSARY FOR ESTIMATING DATA AS A GAUSSIAN DISTRIBUTION
#(GAUSSIAN PROBABILITY DENSITY FUNCTION)

##Find the expected value, or mean, (mu, which will be used as an estimate of b, 
##the position of curve's peak).
muCulturedPairs = mean(CulturedPairs)
muOEPairs = mean(OEPairs)

##Find the variance (sigma squared, which will be used as an estimate of c 
##squared, the standard deviation, which describes the width of the curve).
sigma2CulturedPairs = var(x = CulturedPairs, y = NULL)
sigmaCulturedPairs = sigma2CulturedPairs^(1/2)
sigma2OEPairs = var(x = OEPairs, y = NULL)
sigmaOEPairs = sigma2OEPairs^(1/2)

###---###---###---###---###---###---###---###---###---###---###---###

#ESTIMATE THE PROBABILITY DENSITY DISTRIBUTION OF THE CENTRIOLE PAIR DATA
#AS A GAUSSIAN DISTRIBUTION.

##Make a Gaussian function with parameters of pairs in culture.
g_culture <- function(x){(1/(sigmaCulturedPairs*((2*pi)^(1/2)))) * exp((-1/2)*(((x-muCulturedPairs)/sigmaCulturedPairs)^2))}
#t = function(x){x*x}
plot(g_culture, 0, 2)

##Make a Gaussian function with parameters of pairs in mitotic cells of the OE.
g_OE <- function(x){(1/(sigmaOEPairs*((2*pi)^(1/2)))) * exp((-1/2)*(((x-muOEPairs)/sigmaOEPairs)^2))}
plot(g_OE, 0, 2)

#USE THIS TO CALCULATE THE PROBABILITY THAT A GIVEN NON-PAIR MEASUREMENT COULD ALSO
#BE A CENTRIOLE PAIR.

##Try this with some specified value, v.
v = 1.5
##Find the probability, p, that measurement v could belong to the corresponding pairs data set.
###Note: This is for centrioles in cultured cells.
p = integrate(g_culture, lower = v, upper = Inf)
p
##Now, make a function for p so that you can see the probability across different values of v.
P = function(v){integrate(g_culture, lower = v, upper = Inf)}
#
P(0.7085)
###At v = 0.7085 square microns, the probability that the structure belongs to the centrioles pair 
###group in culture is less than 1%.


##Now, repeat the above for centrioles in OE.
##Find the probability, o, that measurement v could belong to the corresponding pairs data set.
o = integrate(g_OE, lower = v, upper = Inf)
o
##Now, make a function for o so that you can see the probability across different values of v.
O = function(v){integrate(g_OE, lower = v, upper = Inf)}
#
O(0.6695)
###At v = 0.6695 square microns, the probability that the structure belongs to the centrioles pair 
###group in culture is less than 1%.


###---###---###---###---###---###---###---###---###---###---###---###

#FIND OUT HOW MANY OF THE MEASUREMENTS IN THE ROSETTE GROUPS ARE ABOVE THE CUTOFF VALUE.

#Make a list of all rosette measurements above the cutoff.
#We will use the in vitro rosette cut-off, as it is more stringent.
CulturedRosettes_notpairs <- c()
for (i in CulturedRosettes) {
  if (isTRUE(i >= 0.7085)) {
    CulturedRosettes_notpairs <- c(CulturedRosettes_notpairs, i)
  }
}
FracTooBig.culture = (length(CulturedRosettes_notpairs)/length(CulturedRosettes))

OENonPairs_toobig <- c()
for (i in OENonPairs) {
  if (isTRUE(i >= 0.7085)) {
    OENonPairs_toobig <- c(OENonPairs_toobig, i)
  }
}
FracTooBig.OE = (length(OENonPairs_toobig)/length(OENonPairs))

output <- c("Number of cultured rosette measurements:", length(CulturedRosettes), 
            "Number of measurements above cutoff:", length(CulturedRosettes_notpairs),
            "Fraction of in-culture measurements above cutoff:", FracTooBig.culture)
print(output, quote = FALSE)


output <- c("Number of non-pair measurements in OE:", length(OENonPairs), 
      "Number of measurements above cutoff:", length(OENonPairs_toobig),
      "Fraction of in situ measurements above cutoff:", FracTooBig.OE)
print(output, quote = FALSE)


###---###---###---###---###---###---###---###---###---###---###---###