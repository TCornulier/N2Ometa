# R code for the meta-analysis of fertiliser-induced N2O emissions by Hillier et al (2019)
# Copyright (C) 2019  Thomas Cornulier <tomcor.abdn@gmail.com>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Major changes
# change from v20: uses corrected data from Franck & Ulrike
# change from v31: used different parametric approximation for DD effect
# change from v32: removed M15 and M15b, replaced by M15c
setwd("C:/Users/nhy577/Documents/_toma/CONSULTING/JonHillierN2O")

library(ggplot2)
library(lattice)
library(gplots)
library(gtools)
library(mgcv)
library(xlsx)
library(RNCEP)
library(lubridate)
library(rjags)
load.module("glm") ## improved samplers for GLMs often worth loading
library(coda)
library(R2jags)
library(abind)

# for use later:
NCEP.deriv.vars<- c("WetDays.prop.year", "FTCy.year", "meanPosDDays.air.year", "meanPosDDays.soil.year", "FrostDDays.soil.year", "WetDays.exp", "SumPrecip.exp", "SumPrecip.exp.halfNorm80", "DegDays.exp.halfNorm80", "DegDays.exp", "DegDays.Soil.exp", "FTCy.exp", "PosDDays.air.exp", "PosDDays.soil.exp", "FrostDDays.soil.exp", "AvDegDays.exp", "AvPrecip.exp", "PropWetDays.exp")
NCEP.deriv.vars.transf<- c("logSumPrecip.exp.compact", "logSumPrecip.exp.halfNorm80", "logAvPrecip.exp.ct", "logWetDays.exp.compact", "PropWetDays.exp.compact", "logDegDays.exp.compact", "logDegDays.exp.halfNorm80", "AvDegDays.exp", "logDegDays.Soil.exp", "PosDDays.air.exp", "PosDDays.soil.exp", "logFTCy.exp", "mlogFrostDDays.soil.exp", "WetDays.prop.year", "FTCy.year", "meanPosDDays.air.year", "meanPosDDays.soil.year", "FrostDDays.soil.year")


jags2sam<- function(x){
	lapply(x$BUGSoutput$sims.list, FUN= function(z){array(t(z), dim= c(dim(z)[2:1], 1))})
}



.simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep = "", collapse = " ")
}
Capitalize<- function(z){
	bla<- unlist(lapply(as.list(z), function(y) ifelse(is.na(y), NA, .simpleCap(y))))
	bla
}
DF.Capitalize<- function(x, cols.sub= NULL){
	if(is.null(cols.sub)) cols.sub<- rep(T, l= ncol(x))
	for(i in 1:ncol(x)){
		if(cols.sub[i] & is.factor(x[, i])){ x[, i]<- factor(Capitalize(as.character(x[, i])))}
	}
	x
}

panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y, use= "pairwise.complete.obs"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
panel.smooth2<- function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
        abline(h= 0, col= "green")
}

# pairs2<- function(x) pairs(x, lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel = panel.hist)
pairs2<- function(x) pairs(x, lower.panel= panel.smooth2, upper.panel=panel.cor, diag.panel = panel.hist)


######################## loading data ########################
dat<- read.delim("N2O_Database_TEMPLATE_working document_02Mar2018.txt", na.string= c("#N/A", "", "NA"), skip= 0, sep= "\t", strip.white= T)
# dat<- dat[, 1:83]

# tidy column names
bla<- colnames(dat)
bla
bla<- gsub(pattern= "_", replacement= ".", x= bla, fixed = T)
bla<- gsub(pattern= "[.]+$", replacement= "", x= bla)
bla<- gsub(pattern= ".of.", replacement= ".", x= bla, fixed = T)
bla<- gsub(pattern= "..", replacement= ".", x= bla, fixed = T)
bla<- gsub(pattern= "..", replacement= ".", x= bla, fixed = T)
bla<- gsub(pattern= "annual", replacement= "an", x= bla, fixed = T)
bla<- gsub(pattern= "(average)", replacement= "rec", x= bla, fixed = T)
bla<- gsub(pattern= "(avg)", replacement= "rec", x= bla, fixed = T)
bla<- gsub(pattern= "content", replacement= "cont", x= bla, fixed = T)
bla<- gsub(pattern= "Fertilizer", replacement= "Fert", x= bla, fixed = T)
bla<- gsub(pattern= "Cumulative", replacement= "Cum", x= bla, fixed = T)
bla<- gsub(pattern= "cumul", replacement= "Cum", x= bla, fixed = T)
bla<- gsub(pattern= "Carbon.cont.SOM.or.SOC", replacement= "Carbon", x= bla, fixed = T)
colnames(dat) <- bla

############## create unique index ##############
substr(bla, 1, 20)

sum(table(paste(dat$Pub.Title, dat$Fert.type.rep, dat$Main.crop, dat$Cum.N2O, dat$N.rate, dat$Mode.appl.rep))==1)/nrow(dat)
sum(table(paste(dat$Pub.Title, dat$Fert.type.rep, dat$Main.crop, dat$Cum.N2O, dat$N.rate, dat$Mode.appl.rep, dat$liming, dat$Irrigation.exp, dat$Drainage, dat$Soil.texture, dat$pH, dat$Irrigation.exp, dat$Emissionfactor.N2O, dat$Crop.yield.DM))==1)/nrow(dat)

# import conversion tables
NO3.table<- 		read.delim("NO3avail_lookup_12Mar18.txt")
precrop.4.table<- 	read.delim("Tables for regrouping_ULE_2017-7-14_sheet2precrop.txt")
Main.Crop.6.table<- read.delim("Tables for regrouping_ULE_2017-7-14_sheet3maincrop.txt")
N2O.Method.table<- 	read.delim("Tables for regrouping_ULE_2017-7-14_sheet4method.txt")
Nuptake.table<-		read.delim("Nuptake_groups_ULE.txt")

# save.image("DataImport20180306.RData")


###############################################################################


col2rowname<- function(x, col.name){
	stopifnot(length(unique(x[, col.name])) == nrow(x))
	row.names(x)<- as.character(x[, col.name])
	x
}

Main.Crop.6.table[!is.na(Main.Crop.6.table$Main.crop) & Main.Crop.6.table$Main.crop == "sugar beet", ]
Main.Crop.6.table[!is.na(Main.Crop.6.table$Main.crop) & Main.Crop.6.table$Main.crop == "potato", ] # potato and sugar beet inconsistently classified as either row crops or vegetables

NO3.table<- 		DF.Capitalize(NO3.table[, 1:3], cols.sub= c(T, T, T))          
precrop.4.table<- 	DF.Capitalize(precrop.4.table, cols.sub= c(F, T, T))  
Main.Crop.6.table<- DF.Capitalize(Main.Crop.6.table, cols.sub= c(F, T, T))

# NO3.table$NO3.avail.v2<- NO3.table$NO3.availability.proxy
# NO3.table$NO3.avail.v2[!is.na(NO3.table$NO3.availability.proxy) & NO3.table$NO3.availability.proxy %in% 2:5]<- 2
# NO3.table$NO3.avail.v2[!is.na(NO3.table$NO3.availability.proxy) & NO3.table$NO3.availability.proxy > 5]<- NO3.table$NO3.availability.proxy[!is.na(NO3.table$NO3.availability.proxy) & NO3.table$NO3.availability.proxy > 5] - 3

# NO3.table$NO3.avail.v3<- NO3.table$NO3.availability.proxy
# NO3.table$NO3.avail.v3[!is.na(NO3.table$NO3.availability.proxy) & NO3.table$NO3.availability.proxy %in% 0:2]<- 1 # remove level zero and deal with it by interaction with Fert01 in model
# NO3.table$NO3.avail.v3[!is.na(NO3.table$NO3.availability.proxy) & NO3.table$NO3.availability.proxy %in% 4:5]<- 2
# NO3.table$NO3.avail.v3[!is.na(NO3.table$NO3.availability.proxy) & NO3.table$NO3.availability.proxy > 5]<- NO3.table$NO3.availability.proxy[!is.na(NO3.table$NO3.availability.proxy) & NO3.table$NO3.availability.proxy > 5] - 2

# NO3.table$NO3.avail.f3<- "lowNO3"
# NO3.table$NO3.avail.f3[is.na(NO3.table$NO3.availability.proxy)]<- NA
# NO3.table$NO3.avail.f3[!is.na(NO3.table$NO3.availability.proxy) & NO3.table$NO3.avail.v3 >= 4]<- "highNO3"
# NO3.table$NO3.avail.f3[NO3.table$Fertilizer.type_GROUP == "Control"]<- "noFert"
# NO3.table$NO3.avail.f3<- factor(NO3.table$NO3.avail.f3)


# tmpNO3<- NO3.table[which(!duplicated(NO3.table[, 2:3])), -1]
# tmpNO3<- tmpNO3[order(tmpNO3$NO3.availability.proxy), ]
# tmpNO3

# write.xlsx(NO3.table, "NO3avail_lookup.xlsx", row.names= F)
NO3.table$Fertilizer.type_rep<- as.character(NO3.table$Fertilizer.type_rep)
NO3.table$Fertilizer.type_rep[is.na(NO3.table$Fertilizer.type_rep)]<- "NA"
NO3.table<- col2rowname(unique(NO3.table), "Fertilizer.type_rep")
precrop.4.table$precrop<- as.character(precrop.4.table$precrop)
precrop.4.table$precrop[is.na(precrop.4.table$precrop)]<- "NA"
precrop.4.table<- col2rowname(unique(precrop.4.table), "precrop")
Main.Crop.6.table$Main.Crop_GROUP<- as.character(Main.Crop.6.table$Main.Crop_GROUP)
Main.Crop.6.table$Main.Crop_GROUP[is.na(Main.Crop.6.table$Main.Crop_GROUP)]<- "NA"
Main.Crop.6.table<- col2rowname(unique(Main.Crop.6.table[, 2:3]), "Main.Crop_GROUP") # exclude first column as there are duplicates
N2O.Method.table$Method_N2O<- as.character(N2O.Method.table$Method_N2O)
N2O.Method.table$Method_N2O[is.na(N2O.Method.table$Method_N2O)]<- "NA"
N2O.Method.table<- col2rowname(unique(N2O.Method.table), "Method_N2O")
Nuptake.table$Main.Crop_GROUP<- as.character(Nuptake.table$Main.Crop_GROUP)
Nuptake.table$Main.Crop_GROUP[is.na(Nuptake.table$Main.Crop_GROUP)]<- "NA"
Nuptake.table$Main.Crop_GROUP<- Capitalize(as.character(Nuptake.table$Main.Crop_GROUP))
Nuptake.table<- col2rowname(unique(Nuptake.table), "Main.Crop_GROUP")



# bla.NO3<- NO3.table[as.character(dat$Fert.type.rep), ]
# dat$NO3.avail<- bla.NO3$NO3.availability.proxy
# table(is.na(dat$NO3.avail))
# FALSE 
#  3310

names(NO3.table)[2:3]<- c("Fert.type.GROUP", "NO3.avail")
excl<- which(names(dat) %in% "Fert.type.GROUP")
dat$Fertilizer.type_rep<- toupper(as.character(dat$Fert.type.rep))
NO3.table$Fertilizer.type_rep<- toupper(as.character(NO3.table$Fertilizer.type_rep))
dat<- merge(dat[, -excl], NO3.table[, c(1:3)], by= "Fertilizer.type_rep", sort= F, all.x= T)
dat[is.na(dat$Fert.type.rep) & (!is.na(dat$Pub.N)), c("Fert.type.GROUP", "NO3.avail")]<- NO3.table["NA", c(2:3)] # from Frank: "NA for Fertilizer.type_rep stems from two papers, which did not specify the type of mineral fertilizer used. However, since it is mineral N we have assumed a "Mineral_N_mix". I think we can keep it"

# simplify NO3.avail further
NO3convert<- read.delim("NO3avail_lookup_2ndPass.txt")
dat<- merge(dat, NO3convert[, -2], by= "Fert.type.GROUP", sort= F, all.x= T)

names(precrop.4.table)[3]<- c("precrop4")
dat<- merge(dat, precrop.4.table[, c(1,3)], by= "precrop", sort= F, all.x= T)

# bla.precrop<- precrop.4.table[as.character(dat$precrop), ]
# dat$precrop4<- bla.precrop$Precrop3cat..Perenial..Legumes..Other..Rice.
# table(is.na(dat$precrop4))

bla.maincrop<- Main.Crop.6.table[Capitalize(as.character(dat$Main.Crop.GROUP)), ]
bla.maincrop$MainCrop6cat..Legumes..Others..Cereals..Vegetables..Grasslands..Rice.<- as.character(bla.maincrop$MainCrop6cat..Legumes..Others..Cereals..Vegetables..Grasslands..Rice.)
bla.maincrop$MainCrop6cat..Legumes..Others..Cereals..Vegetables..Grasslands..Rice.[bla.maincrop$MainCrop6cat..Legumes..Others..Cereals..Vegetables..Grasslands..Rice.== "Grassland/pasture/steppe/savanah"]<- "Grasslands"
dat$maincrop6<- factor(bla.maincrop$MainCrop6cat..Legumes..Others..Cereals..Vegetables..Grasslands..Rice.)
table(is.na(dat$maincrop6))

bla.N2Omethod<- N2O.Method.table[as.character(dat$Method.N2O), ]
dat$Method.N2O.GROUP<- bla.N2Omethod$Method_N2O_GROUP
table(is.na(dat$Method.N2O.GROUP))

dat$Nuptake.group<- as.character(Nuptake.table[Capitalize(as.character(dat$Main.Crop.GROUP)), "N.uptake.group"])
dat$Nuptake.group[dat$Nuptake.group == "NA"]<- NA
dat$Nuptake.group<- factor(dat$Nuptake.group)

# check number of NAs per column
sort(unlist(lapply(dat, function(x) round(sum(!is.na(x))/nrow(dat)*100))), decreasing= T)
sort(unlist(lapply(dat, function(x) round(sum(!is.na(x))))), decreasing= T)

table(dat$Fert.type.GROUP, useNA= "always")

table(dat$Main.Crop.GROUP, useNA= "always")


dat$Method.N2O[dat$Method.N2O == "cc"]<- "closed chamber"
dat$Method.N2O[is.na(dat$Method.N2O)]<- "closed chamber" # assume closed chamber was used when method not specified (N= 56)
dat$Method.N2O<- factor(dat$Method.N2O)
table(dat$Method.N2O, useNA= "always")


dat$CF <- as.character(dat$Fert.type.GROUP)
dat$CF[!dat$CF == "Control"] <- "Fertilized"
dat$CF<- factor(dat$CF)

dat<- dat[!is.na(dat$Fert.type.GROUP), ]
dat$Fert.type.GROUP<- factor(dat$Fert.type.GROUP)

dat<- dat[!is.na(dat$Main.Crop.GROUP), ]
crop.short.names<- c('bare/fallow'= "bare.fal", 'cereal_legume rotations'= "cerealeg", 'cereals'= "cereals", 'grassland/pasture/steppe/savanah'= "grasslands", 'legumes'= "legumes", 'perennials (trees/forests/plantations)'= "perennials", 'rice'= "rice", 'row crops (sugarcane, cotton, etc.)'= "row.crops", 'vegetables'= "vegetables")
dat$Main.Crop.GROUP<- crop.short.names[as.character(dat$Main.Crop.GROUP)]
dat$Main.Crop.GROUP<- factor(dat$Main.Crop.GROUP)

dat$Multicrops<- "OneCrop"
dat$Multicrops[dat$Ncrops>1]<- "MoreThanOne"
dat$Multicrops<- factor(dat$Multicrops, levels= c("OneCrop", "MoreThanOne"))
dat$OneCrop<- as.numeric(dat$Multicrops == "OneCrop")

# clean any unused factor level
for(i in 1:ncol(dat)){if(is.factor(dat[, i])) dat[, i]<- factor(dat[, i])}

dat$Fert<- as.character(dat$Fert.type.GROUP)
dat$Fert[dat$Fert.type.GROUP == "Mix_N + NI"]<- "N.NI"
dat$Fert[dat$Fert.type.GROUP == "Mineral_N_mix"]<- "MinN"
dat$Fert[dat$Fert.type.GROUP == "Mix_N_mineral_organic"]<- "OrgMinN"
dat$Fert[dat$Fert.type.GROUP == "Ammoniumnitrate"]<- "AmmNit"
dat$Fert[dat$Fert.type.GROUP == "Urea + additives"]<- "UreaAdd"
dat$Fert[dat$Fert.type.GROUP == "organic"]<- "Org"
dat$Fert[dat$Fert.type.GROUP == "organic + NI"]<- "Nitrate"
dat$Fert[dat$Fert.type.GROUP == "Ammonium"]<- "Amm"
dat$Fert[dat$Fert.type.GROUP == "Urea + NI"]<- "UreaNI"
dat$Fert[dat$Fert.type.GROUP == "Anhydrous ammonia"]<- "AnhyAmm"
dat$Fert<- factor(dat$Fert, levels= c("Control", "Nitrate", "Urea", "OrgMinN", "N.NI", "MinN", "AmmNit", "UreaAdd", "Org", "Amm", "UreaNI", "AnhyAmm"))
table(dat$Fert, useNA= "always")
FertMat<- model.matrix(~Fert-1, dat) # create an indicator variable for each treatment
dat<- cbind(dat, FertMat)

# correct some typos in decimal separator of numeric variables
dat$Longitude<- as.numeric(gsub(",", ".", as.character(dat$Longitude)))
dat$Cum.N2O<- as.numeric(gsub(",", ".", as.character(dat$Cum.N2O)))
# correct NAs read as characters in precrop4
dat$precrop4[dat$precrop4 == "NA"]<- NA

# add a little noise to the coordinates
dat$Lat<- rnorm(nrow(dat), dat$Latitude, sd= 0.001)
dat$Lon<- rnorm(nrow(dat), dat$Longitude, sd= 0.001)

dat$N2O.net<- (as.numeric(as.character(dat$Emissionfactor.N2O)) / 100) * dat$N.rate
min(dat$N2O.net, na.rm= T) # -625 # this is Pub.N == 180 - please check as net N2O records look dodgy
dat$logN2O.net<- log(dat$N2O.net + 1 + abs(min(dat$N2O.net, na.rm= T)))
dat$transfN2O.net<- (dat$N2O.net + 0.3) ^ 0.1

dat$N2O.cum<- dat$Cum.N2O
min(dat$N2O.cum, na.rm= T) # -1.65
dat$logN2O.cum<- log(dat$N2O.cum + 1 + abs(min(dat$N2O.cum, na.rm= T)))

dat$logChamb.size<- log(dat$Chamber.size)

dat$logNdays<- log(dat$Length.exp.N2O + 5)

plot(dat$logNdays ~ dat$Length.exp.N2O)
plot(dat$Length.exp.N2O ~ dat$N2O.meas.days); abline(0, 1) # length of experiment vs. nb of days during which N2O was measured apparently...
hist(dat$logNdays, nclass= 50)

# discuss with all what to do with negative emissions:
# - censor negative N2O emissions? (not all studies have reported? different processes involved?)
# - remove them?
# - leave as negative (but have all studies reported negative emissions??)

# for now assume we'll censor everything negative
dat$logN2O.cum.pos<- log(ifelse(dat$N2O.cum >=0, dat$N2O.cum, 0) + 1)




###########################################################################
########### Add temperature and precipitation from gridded data ###########
###########################################################################

# create metadata of lat long, start and end years by site
extent<- by(dat, dat$ExpID, function(x){
	tmp<- data.frame(	ExpID= unique(as.character(x$ExpID)), 
						LAT= as.numeric(unique(x$Latitude)), 
						LON= as.numeric(unique(x$Longitude)),
						START.year= as.numeric(unique(as.character(x$Year.Exp))),
						START.month= as.numeric(unique(x$Month.year.START)),
						DURATION.days= as.numeric(unique(x$Length.exp.N2O)))
})
extent<- data.frame(do.call(rbind, extent))
head(extent)

# compute end year + safe margin for data download
extent$END.year<- extent$START.year + ceiling(extent$DURATION.days / 365)

# require(lattice)
# lats<- equal.count(extent$LAT, number=3, overlap=.1)
# xyplot(extent$START.month ~ ifelse(extent$DURATION.days > 365, 365, extent$DURATION.days)|lats)
# no consistent starting month

extent$START.day<- as.POSIXct("1099/01/01", tz= "GMT") # give dummy date to give column a POSIXct class
for(i in 1:nrow(extent)){ 
	print(i)
	try(extent$START.day[i]<- as.POSIXct(paste(extent$START.year[i], extent$START.month[i], "15", sep= "/"), tz= "GMT")) # replace value with start date, assuming experiment started the 15th day of the starting month
}
# (ignore errors)
extent$START.day[extent$START.day < as.POSIXct("1199/01/01", tz= "GMT")]<- NA # replace dummy values by NA

# compute experiment end date by adding duration to start date
extent$END.day<- extent$START.day + extent$DURATION.days * 3600 * 24 

# create a study ID sharing identical experiment characteristics (location, timing)
extent$studyID<- apply(extent, 1, FUN= function(x) paste(x[-1], collapse= ""))
dat$studyID<- apply(data.frame(as.numeric(dat$Latitude), as.numeric(dat$Longitude), as.numeric(as.character(dat$Year.Exp)), as.numeric(dat$Month.year.START), as.numeric(dat$Length.exp.N2O)), 1, FUN= function(x) paste(x, collapse= ""))

# make 2 "copies" of the data set with different lengths
extent.all<- extent # copy with all nrow= N observations (3217)
extent<- extent[!duplicated(extent[, -1]), ] # copy with unique experiment characteristics (670 rows)
rownames(extent)<- extent$studyID

#################### download data from NCEP

# NCEP.list<- list()
# for(j in 1:nrow(extent)){
# print(paste("################################", j, "/", nrow(extent), "################################\n"))
# 	NCEP.list[[as.character(extent$studyID[j])]]<- list()
# 	if(!any(is.na(extent[j, ]))) {
# 		for(i in c("air.2m", "prate.sfc", "tmp.0-10cm")){
# 			lat<- extent$LAT[j]
# 			lon<- extent$LON[j]
# 			bla<- NCEP.gather(variable= i, level= "gaussian", months.minmax= c(1, 12), years.minmax= unlist(extent[j, c("START.year", "END.year")]), lat.southnorth= lat + c(0, 0.00001), lon.westeast= lon + c(0, 0.00001), reanalysis2 = T, return.units = TRUE, status.bar= F)
# 			nearest.lat<- which.min(abs(as.numeric(dimnames(bla)[[1]])-lat))[1] # [1] in the unlikely case that they are equidistant
# 			nearest.lon<- which.min(abs(as.numeric(dimnames(bla)[[2]])-lon))[1]
# 			bla<- bla[nearest.lat, nearest.lon, ]
# 			day<- substr(names(bla), 1, 10)
# 			NCEP.list[[as.character(extent$studyID[j])]][["StudyID"]]<- as.character(extent$studyID[j])
# 			NCEP.list[[as.character(extent$studyID[j])]][["Date.POSIXct"]]<- ymd(unique(day), tz= "GMT")
# 			NCEP.list[[as.character(extent$studyID[j])]][[i]]<- tapply(bla, day, mean)
# 		}
# 		NCEP.list[[as.character(extent$studyID[j])]]<- data.frame(	StudyID= NCEP.list[[as.character(extent$studyID[j])]][["StudyID"]],
# 													Date.POSIXct= NCEP.list[[as.character(extent$studyID[j])]][["Date.POSIXct"]],
# 													air.2m= NCEP.list[[as.character(extent$studyID[j])]][["air.2m"]],
# 													prate.sfc= NCEP.list[[as.character(extent$studyID[j])]][["prate.sfc"]],
# 													#weasd.sfc= NCEP.list[[as.character(extent$studyID[j])]][["weasd.sfc"]],
# 													'tmp.0-10cm'= NCEP.list[[as.character(extent$studyID[j])]][["tmp.0-10cm"]])
# 	}
# 	else { NCEP.list[[as.character(extent$studyID[j])]]<- NA }
# }

# names(NCEP.list)<- extent$studyID
# save(NCEP.list, file= "NCEPdownload20180308.RData")
load("NCEPdownload20180308.RData")

NCEP.nonempty<- unlist(lapply(NCEP.list, FUN= function(x) is.data.frame(x)))
table(NCEP.nonempty)
# FALSE  TRUE 
#    96   574 # 14% missing
   
# NCEP.dat<- do.call(rbind, NCEP.list[NCEP.nonempty])
# NCEP.dat$air.2m<- as.numeric(NCEP.dat$air.2m) - 273.15 # Kelvin to Celcius
# NCEP.dat$prate.sfc<- NCEP.dat$prate.sfc * 24 * 60 * 60  # Kg/m^2/s to mm/day
# NCEP.dat$tmp.0.10cm<- NCEP.dat$tmp.0.10cm - 273.15 # Kelvin to Celcius

table(NCEP.nonempty[extent.all$studyID]) # 13% missing (419)

######### Compute derived variables per year over the years (+1) when experiment was run #########
# compute wet day probability for the whole year
extent$WetDays.prop.year<- NA
for(i in extent$studyID[NCEP.nonempty]){
	Nyr<- extent$END.year[extent$studyID == i] - extent$START.year[extent$studyID == i] + 1
	extent$WetDays.prop.year[extent$studyID == i]<- sum((NCEP.list[[i]]$prate.sfc * 24 * 60 * 60) > 1, na.rm= T)/(Nyr*365)
}
plot(extent$WetDays.prop.year)

# compute number of freeze-thaw cycles per year
extent$FTCy.year<- NA
for(i in extent$studyID[NCEP.nonempty]){
	Nyr<- extent$END.year[extent$studyID == i] - extent$START.year[extent$studyID == i] + 1
	extent$FTCy.year[extent$studyID == i]<- sum(running(NCEP.list[[i]]$tmp.0.10cm - 273.15, fun= function(x){
		x[1]<= -1 & x[2] > -1}, width= 2, pad= F), na.rm= T)/Nyr
}
plot(extent$FTCy.year)
table(extent$FTCy.year>0)
# FALSE  TRUE 
#   132   442

# example with ~ 35 FT cycles over 2 full years:
plot(NCEP.list[[" 39.10-117.702011 5 3652012"]]$tmp.0.10cm-273.5, type= "l"); abline(h= 0) # quite smooth, no apparent need to add more smoothness to it

# compute positive degree-days for air or soil
extent$meanPosDDays.air.year<- NA
for(i in extent$studyID[NCEP.nonempty]){
	Nyr<- extent$END.year[extent$studyID == i] - extent$START.year[extent$studyID == i] + 1
	extent$meanPosDDays.air.year[extent$studyID == i]<- sum(ifelse(NCEP.list[[i]]$air.2m - 273.15 > 0, NCEP.list[[i]]$air.2m - 273.15, 0), na.rm= T)/(Nyr*365)
}
plot(extent$meanPosDDays.air.year)

extent$meanPosDDays.soil.year<- NA
for(i in extent$studyID[NCEP.nonempty]){
	Nyr<- extent$END.year[extent$studyID == i] - extent$START.year[extent$studyID == i] + 1
	extent$meanPosDDays.soil.year[extent$studyID == i]<- sum(ifelse(NCEP.list[[i]]$tmp.0.10cm - 273.15 > 0, NCEP.list[[i]]$tmp.0.10cm - 273.15, 0), na.rm= T)/(Nyr*365)
}
plot(extent$meanPosDDays.soil.year)
table(extent$meanPosDDays.soil.year == 0)
# FALSE  TRUE 
#   521    62 # looks dodgy as suggests that 62 experimental sites had no positive surface soil temperature for the whole year
  
# compute below-frost degree-days for soil
extent$FrostDDays.soil.year<- NA
for(i in extent$studyID[NCEP.nonempty]){
	extmp<- extent[extent$studyID == i, ]
	extent$FrostDDays.soil.year[extent$studyID == i]<- sum(ifelse(NCEP.list[[i]]$tmp.0.10cm - 273.15 < -1, NCEP.list[[i]]$tmp.0.10cm - 273.15, 0), na.rm= T)/Nyr
}
plot(extent$FrostDDays.soil.year)



######### Compute derived variables over the course of the experiment only #########
# compute N wet days over course of experiment
extent$WetDays.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	extent$WetDays.exp[extent$studyID == i]<- sum((NCEP.list[[i]]$prate.sfc[exp.extent] * 24 * 60 * 60) > 1, na.rm= T)
}
plot(extent$WetDays.exp)

# compute proportion of wet days over course of experiment
extent$PropWetDays.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	extent$PropWetDays.exp[extent$studyID == i]<- sum((NCEP.list[[i]]$prate.sfc[exp.extent] * 24 * 60 * 60) > 1, na.rm= T) / sum(exp.extent)
}
plot(extent$PropWetDays.exp)

# compute Average of Precipitations over course of experiment
extent$AvPrecip.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	extent$AvPrecip.exp[extent$studyID == i]<- mean(NCEP.list[[i]]$prate.sfc[exp.extent] * 24 * 60 * 60, na.rm= T) # converting Kg/m^2/s to mm/day
}
plot(extent$AvPrecip.exp)

# compute Sum of Precipitations over course of experiment
extent$SumPrecip.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	extent$SumPrecip.exp[extent$studyID == i]<- sum(NCEP.list[[i]]$prate.sfc[exp.extent] * 24 * 60 * 60, na.rm= T) # converting Kg/m^2/s to mm/day
}
plot(extent$SumPrecip.exp)

# compute Sum of Precipitations over course of experiment (Front-weighted)
extent$SumPrecip.exp.halfNorm80<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	count.T<- sum(exp.extent) # duration of experiment
	extent$SumPrecip.exp.halfNorm80[extent$studyID == i]<- sum((NCEP.list[[i]]$prate.sfc[exp.extent] * 24 * 60 * 60) * dnorm((1:count.T) %% 365, mean = 0, sd = 80)+0.001, na.rm= T) # converting Kg/m^2/s to mm/day
	 # modulus - remainder of the integer product (looping over 365 days)
	 # front-weighting with a 0.001 + half-normal(mu= 0, sd= 80)
}
plot(extent$SumPrecip.exp.halfNorm80)
plot(extent$SumPrecip.exp, extent$SumPrecip.exp.halfNorm80)


# compute Sum of positive degree-days for air temperature over course of experiment
extent$DegDays.exp.halfNorm80<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	tmp<- NCEP.list[[i]]$air.2m[exp.extent] - 273.15 # tmp.0.10cm
	count.T<- sum(exp.extent) # duration of experiment
	extent$DegDays.exp.halfNorm80[extent$studyID == i]<- sum(ifelse(tmp<0, 0, tmp) * dnorm((1:count.T) %% 365, mean = 0, sd = 80)+0.001, na.rm= T) # converting K to celcius
	 # modulus - remainder of the integer product (looping over 365 days)
	 # front-weighting with a 0.001 + half-normal(mu= 0, sd= 80)
}
plot(extent$DegDays.exp.halfNorm80)
plot(extent$DegDays.exp, extent$DegDays.exp.halfNorm80)

# compute Sum of positive degree-days for air temperature over course of experiment
extent$DegDays.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	tmp<- NCEP.list[[i]]$air.2m[exp.extent] - 273.15 # tmp.0.10cm
	extent$DegDays.exp[extent$studyID == i]<- sum(ifelse(tmp<0, 0, tmp), na.rm= T) # converting K to celcius
}
plot(extent$DegDays.exp)

# compute Average of positive degree-days for air temperature over course of experiment
extent$AvDegDays.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	tmp<- NCEP.list[[i]]$air.2m[exp.extent] - 273.15 # tmp.0.10cm
	extent$AvDegDays.exp[extent$studyID == i]<- mean(ifelse(tmp<0, 0, tmp), na.rm= T) # converting K to celcius
}
plot(extent$AvDegDays.exp)

# compute Sum of positive degree-days for soil surface over course of experiment
extent$DegDays.Soil.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	tmp<- NCEP.list[[i]]$tmp.0.10cm[exp.extent] - 273.15
	extent$DegDays.Soil.exp[extent$studyID == i]<- sum(ifelse(tmp<0, 0, tmp), na.rm= T) # converting K to celcius
}
plot(extent$DegDays.Soil.exp)

# compute number of freeze-thaw cycles
extent$FTCy.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	extent$FTCy.exp[extent$studyID == i]<- sum(running(NCEP.list[[i]]$tmp.0.10cm[exp.extent] - 273.15, fun= function(x){
		x[1]<= -1 & x[2] > -1}, width= 2, pad= F), na.rm= T)
}
plot(extent$FTCy.exp)
table(extent$FTCy.exp>0)
# FALSE  TRUE 
#  351   232

# compute positive degree-days for air or soil
extent$PosDDays.air.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	extent$PosDDays.air.exp[extent$studyID == i]<- sum(ifelse(NCEP.list[[i]]$air.2m[exp.extent] - 273.15 > 0, NCEP.list[[i]]$air.2m[exp.extent] - 273.15, 0), na.rm= T)
}
plot(extent$PosDDays.air.exp)

extent$PosDDays.soil.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	extent$PosDDays.soil.exp[extent$studyID == i]<- sum(ifelse(NCEP.list[[i]]$tmp.0.10cm[exp.extent] - 273.15 > 0, NCEP.list[[i]]$tmp.0.10cm[exp.extent] - 273.15, 0), na.rm= T)
}
plot(extent$PosDDays.soil.exp)

# compute below-frost degree-days for soil
extent$FrostDDays.soil.exp<- NA
for(i in extent$studyID[NCEP.nonempty]){
	exp.extent<- NCEP.list[[i]]$Date.POSIXct > extent$START.day[extent$studyID == i] & NCEP.list[[i]]$Date.POSIXct < extent$END.day[extent$studyID == i]
	extent$FrostDDays.soil.exp[extent$studyID == i]<- sum(ifelse(NCEP.list[[i]]$tmp.0.10cm[exp.extent] - 273.15 < -1, NCEP.list[[i]]$tmp.0.10cm[exp.extent] - 273.15, 0), na.rm= T)
}
plot(extent$FrostDDays.soil.exp)




# plot all climate predictors
pdf("ClimatePredictors.pdf", width= 6, height= 10)
par(mfcol= c(5, 2), mar= c(3.1, 3.1, 2.1, 0.5))	
plot(extent$WetDays.exp, main= "WetDays")
plot(extent$PosDDays.air.exp, main= "PosDDays.air")
plot(extent$PosDDays.soil.exp, main= "PosDDays.soil")
plot(extent$FTCy.exp, main= "FTCy")
plot(extent$FrostDDays.soil.exp, main= "FrostDDays.soil")
hist(extent$WetDays.exp, nclass= 100)
hist(extent$PosDDays.air.exp, nclass= 100)
hist(extent$PosDDays.soil.exp, nclass= 100)
hist(extent$FTCy.exp, nclass= 100)
hist(extent$FrostDDays.soil.exp, nclass= 100)
dev.off()

# add climate predictors to all experiments in "extent.all"
extent.all<- merge(extent.all, extent[, 8:ncol(extent)], by= "studyID", all.x= T, sort= F)
# add the latest ones if done incrementally
# extent.all<- merge(extent.all, extent[, c("studyID", "PropWetDays.exp")], by= "studyID", all.x= T, sort= F)

# save(extent, extent.all, file= "Data_Clim_20180321.RData")
# load("Data_Clim_20180321.RData")

####################################################################
########### Sub model for Bulk density and soil porosity ###########
####################################################################
dat$logCarbon<- log(dat$Carbon.cont.SOC.in + 1)

# convert texture into particle proportions (using region mean from http://www.fao.org/3/a-i3794e.pdf page 192)
dat$text2clay<- NA
dat$text2silt<- NA
dat$text2sand<- NA

table(dat$Soil.texture)

dat[dat$Soil.texture %in% c("silty clay loam", "silty caly loam"), c("text2clay", "text2silt", "text2sand")]<- c(35, 55, 10)
dat[dat$Soil.texture %in% c("silty loam", "Silty loam", "loamy silt"), c("text2clay", "text2silt", "text2sand")]<- c(15, 65, 20)
dat[dat$Soil.texture %in% c("silt"), c("text2clay", "text2silt", "text2sand")]<- c(5, 88, 7)
dat[dat$Soil.texture %in% c("loam"), c("text2clay", "text2silt", "text2sand")]<- c(20, 40, 40)
dat[dat$Soil.texture %in% c("sandy loam"), c("text2clay", "text2silt", "text2sand")]<- c(10, 30, 60)
dat[dat$Soil.texture %in% c("loamy sand"), c("text2clay", "text2silt", "text2sand")]<- c(5, 15, 80)
dat[dat$Soil.texture %in% c("sand"), c("text2clay", "text2silt", "text2sand")]<- c(5, 5, 90)
dat[dat$Soil.texture %in% c("sandy clay loam"), c("text2clay", "text2silt", "text2sand")]<- c(27, 13, 60)
dat[dat$Soil.texture %in% c("sandy clay"), c("text2clay", "text2silt", "text2sand")]<- c(43, 7, 50)
dat[dat$Soil.texture %in% c("clay"), c("text2clay", "text2silt", "text2sand")]<- c(52, 22, 26)
dat[dat$Soil.texture %in% c("clay loam", "Clay loam", "loamy clay"), c("text2clay", "text2silt", "text2sand")]<- c(35, 33, 32)
dat[dat$Soil.texture %in% c("silty clay"), c("text2clay", "text2silt", "text2sand")]<- c(47, 46, 7)

dat[dat$Soil.texture %in% c("organic"), c("clay", "silt", "sand", "logCarbon")] # all "organic" soils have no particle proportions but have organic carbon measured. Assume this is a loam / clay-loam.
pairs2(dat[, c("logCarbon", "clay", "silt", "sand")])
dat[dat$Soil.texture %in% c("organic"), c("text2clay", "text2silt", "text2sand")]<- c(25, 35, 40)

apply(dat[, c("text2clay", "text2silt", "text2sand")], 2, mean, na.rm= T)
# text2clay text2silt text2sand 
# 33.29709  33.31556  33.38735

table(is.na(dat$Soil.texture)) # TRUE = 264 , FALSE = 2953
table(is.na(dat$Soil.texture), is.na(dat$clay)) # TRUE/TRUE = 99

# dat[is.na(dat$Soil.texture), c("text2clay", "text2silt", "text2sand")]<- c(33, 33, 34) # assign average in data set

dat$clay.imputed<- ifelse(is.na(dat$clay), dat$text2clay, dat$clay)
dat$silt.imputed<- ifelse(is.na(dat$silt), dat$text2silt, dat$silt)
dat$sand.imputed<- ifelse(is.na(dat$sand), dat$text2sand, dat$sand)

hist(dat$Bulk.density)

pairs2(dat[, c("clay", "silt", "logCarbon", "pH", "Bulk.density")])
pairs2(dat[, c("clay.imputed", "silt.imputed", "logCarbon", "pH", "Bulk.density")])

dat$logsilt<- log(dat$silt + 1)
dat$logclay<- log(dat$clay + 1)
pairs2(dat[, c("logclay", "logsilt", "logCarbon", "pH", "Bulk.density")])

dat$fPub.N<- factor(dat$Pub.N)

coplot(Bulk.density ~ clay | logCarbon, data= dat, number= 6)
silt.resid<- residuals(lm(Bulk.density ~ silt, data= dat, na.action= na.exclude))
coplot(silt.resid ~ clay | logCarbon, data= dat, number= 6)

Mod.Bulk.1<- gam(Bulk.density ~ logCarbon * clay + silt, data= dat)
summary(Mod.Bulk.1)
par(mfrow= c(1, 1)); vis.gam(Mod.Bulk.1)
par(mfrow= c(2, 2)); gam.check(Mod.Bulk.1)

Mod.Bulk.1.imputed<- gam(Bulk.density ~ logCarbon * clay.imputed + silt.imputed, data= dat)
summary(Mod.Bulk.1.imputed)
par(mfrow= c(1, 2)); vis.gam(Mod.Bulk.1.imputed); vis.gam(Mod.Bulk.1.imputed, view= c("silt.imputed", "clay.imputed"))
par(mfrow= c(2, 2)); gam.check(Mod.Bulk.1.imputed)
# some differences, no obvious winner

Mod.Bulk.2<- gamm(Bulk.density ~ logCarbon * clay + silt, random = list(fPub.N= ~ 1), data= dat)
summary(Mod.Bulk.2$gam)
dev.new()
par(mfrow= c(1, 2)); plot(Mod.Bulk.2$gam, all.terms= T, select= 3); vis.gam(Mod.Bulk.2$gam)
par(mfrow= c(2, 2)); gam.check(Mod.Bulk.2$gam)
# some changes in parameter estimates

# not clear one model is preferable to the other

Mod.Bulk.3<- gamm(Bulk.density ~ te(logCarbon, clay) + silt, data= dat, random = list(fPub.N= ~ 1))
summary(Mod.Bulk.3$gam)
par(mfrow= c(1, 1)); vis.gam(Mod.Bulk.3$gam, view= c("logCarbon", "clay"), theta= -30, phi= 25, too.far= 0.15)
# clearly too wiggly

Mod.Bulk.2.imputed<- gamm(Bulk.density ~ logCarbon * clay.imputed + silt.imputed, random = list(fPub.N= ~ 1), data= dat)
summary(Mod.Bulk.2.imputed$gam)
dev.new()
par(mfrow= c(1, 2)); plot(Mod.Bulk.2.imputed$gam, all.terms= T, select= 3); vis.gam(Mod.Bulk.2.imputed$gam)
par(mfrow= c(2, 2)); gam.check(Mod.Bulk.2.imputed$gam)
# use this one

lc<- ifelse(is.na(dat$logCarbon), quantile(dat$logCarbon, 0.05, na.rm= T), dat$logCarbon) # assume low carbon content if carbon value missing
predBD<- predict(Mod.Bulk.2.imputed$gam, newdata= cbind(logCarbon= lc, dat[, c("clay.imputed", "silt.imputed")], fPub.N= "1"))
dat$Bulk.density.imputed<- ifelse(is.na(dat$Bulk.density), predBD, dat$Bulk.density)

dat$Clay2Silt1ct<- ((dat$silt.imputed + dat$clay.imputed * 2) - mean((dat$silt.imputed + dat$clay.imputed * 2), na.rm= T)) / 200
dat$Clay1Silt1ct<- ((dat$silt.imputed + dat$clay.imputed) - mean((dat$silt.imputed + dat$clay.imputed), na.rm= T)) / 100
dat$Clay1ct<- (dat$clay.imputed - mean(dat$clay.imputed, na.rm= T)) / 100

dat$Clay2Silt1<- (dat$silt.imputed + dat$clay.imputed * 2) / 2 # no centering, no scaling for graphical purposes
dat$Clay1Silt1<- dat$silt.imputed + dat$clay.imputed # no centering, no scaling for graphical purposes
dat$logclay.imp<- log(dat$clay.imputed / 100 + 1)
dat$logclay.imp.ct<- log(dat$clay.imputed / 100 + 1) - 0.24

#################################### more formatting of dat ##############################################
# remove problematic studies
# dat<- dat[!dat$Pub.N %in% unique(dat[(dat$CF == "Control" & dat$N2O.net > 0), "Pub.N"]), ]
# dat<- dat[!dat$Pub.N %in% unique(dat[(dat$CF == "Control" & dat$N.rate > 0), "Pub.N"]), ]
# about 90 data points removed above

# combine dat with NCEP-derived weather indices
dat<- merge(dat, extent.all[, c("studyID", "ExpID", NCEP.deriv.vars)], by= "ExpID", all.x= T, sort= F)
# add the latest weather variables:
# dat<- dat[, names(dat) != "PropWetDays.exp"] # if adding new predictors incrementally
# dat<- merge(dat, extent.all[, c("ExpID", "PropWetDays.exp")], by= "ExpID", all.x= T, sort= F)

min(dat$N2O.net, na.rm= T) # -10.3
plot(dat$N2O.net)
dat$logN2O.net<- log(dat$N2O.net + 1 + abs(min(dat$N2O.net, na.rm= T)))
dat$vegetables<- dat$Main.Crop.GROUP == "vegetables"
dat$FertOrg.OrgMinN<- factor(dat$FertOrg==1 | dat$FertOrgMinN==1)

dat$N.appl<- dat$No.application
dat$N.appl[dat$No.application ==1]<- "1"
dat$N.appl[dat$No.application %in% c(2, 3)]<- "2-3"
dat$N.appl[dat$No.application > 3]<- "4+"
dat$N.appl1<- as.numeric(dat$N.appl == "1")
dat$N.appl2.3<- as.numeric(dat$N.appl == "2-3")
dat$N.appl4.<- as.numeric(dat$N.appl == "4+")
# standardize N.rate by log-transforming + centering
mean(log(dat$N.rate + 1)[log(dat$N.rate + 1)>0], na.rm= T) #  mean for non-zero values is 5.153328
dat$logNrate.cent<- log(dat$N.rate + 1) - 5 # centering around 5
dat$logNrate.cent.N.appl1<-   dat$logNrate.cent * dat$N.appl1 # pre-multiply. s(..., by= N.appl1:logNrate.cent) was crashing gam()
dat$logNrate.cent.N.appl2.3<- dat$logNrate.cent * dat$N.appl2.3
dat$logNrate.cent.N.appl4.<-  dat$logNrate.cent * dat$N.appl4.
# standardize logNdays by centering
mean(dat$logNdays, na.rm= T) # 5.197152
dat$logNdays.cent<- dat$logNdays - 5 # centering around 5
# standardize logNdays by centering
mean(dat$pH, na.rm= T) # 6.775199
dat$pH.cent<- dat$pH - 7 # centering around 7
dat$pH.cent.2<- dat$pH.cent^2 # squaring
dat$pH.cent.3<- dat$pH.cent^3 # cubing
dat$pH.cent.4<- dat$pH.cent^4
dat$pH.cent.5<- dat$pH.cent^5
dat$pH.cent.6<- dat$pH.cent^6

dat$pH1<- NA
dat$pH2<- NA
dat$pH3<- NA
dat$pH4<- NA
dat$pH5<- NA
dat$pH6<- NA
dat[!is.na(dat$pH.cent), paste("pH", 1:6, sep= "")]<- poly(dat$pH.cent[!is.na(dat$pH.cent)], degree= 6)

for(i in 1:ncol(dat)){if(is.factor(dat[, i])) dat[, i]<- factor(dat[, i])}

# new data set excluding controls
dat2<- dat[dat$CF == "Fertilized", ]
for(i in 1:ncol(dat2)){if(is.factor(dat2[, i])) dat2[, i]<- factor(dat2[, i])}


################# Data exploration #################
# create file "Data_DescriptiveGraphs.pdf"
source("prelim data analysis/DescriptiveGraphs_code4.txt")
source("prelim data analysis/ExploratoryGraphs_code1.txt")


################# World maps of data #################
# create file "CoordinatesCheck.pdf" and "DataMaps.pdf"
source("prelim data analysis/DescriptiveMaps_code1.txt")




dat$Pre.Leg<- dat$precrop4 == "Legumes" # 130 occurences
dat$Per2An<- dat$precrop4 == "Perennial" & dat$maincrop6 != "Grasslands" # 200 occurences

dat$Fert01<- 1; dat$Fert01[dat$CF == "Control"]<- 0

################# remove rice from analysis for now
dat<- dat[dat$Main.Crop.GROUP != "rice", ]
dat$maincrop6<- factor(dat$maincrop6)

################# remove lines with only NA values or NA in the response
dat<- dat[!(is.na(dat$fPub.N) | is.na(dat$logN2O.cum)), ]

dat$studyID.f<- factor(dat$studyID.x)

dat$logNrate.cent.v2<- log(ifelse(dat$N.rate > 0, dat$N.rate, 100) + 1) - 5 # replace zeros by arbitrary value (for plotting convenience) - to be used pre-multiplied by Fert01 to get correct estimates #  mean for non-zero values is 5.153328

dat$logFTCy.exp<- log(dat$FTCy.exp+1)
dat$logWetDays.exp<- log(dat$WetDays.exp+1)
dat$logWetDays.exp.compact<- log(dat$WetDays.exp+10) - 4 # mean(dat$logWetDays.exp.compact, na.rm= T) # 4.153867
# dat$PropWetDays.exp<- dat$WetDays.exp / dat$Length.exp.N2O
dat$PropWetDays.exp.compact<- log(dat$PropWetDays.exp+0.5) + 0.2 # mean(log(dat$PropWetDays.exp+0.5), na.rm= T) # -0.2103121
dat$logSumPrecip.exp<- log(dat$SumPrecip.exp+1)
dat$logSumPrecip.exp.compact<- log(dat$SumPrecip.exp+50) - 6 # mean(log(dat$SumPrecip.exp+50), na.rm= T) # 6.127614
dat$logDegDays.exp<- log(dat$DegDays.exp+1)
dat$logDegDays.exp.compact<- log(dat$DegDays.exp+1500) - 8 # mean(log(dat$DegDays.exp+1500), na.rm= T) # 8.323379
dat$logDegDays.Soil.exp<- log(dat$DegDays.Soil.exp+1)
dat$logSumPrecip.exp.halfNorm80<- log(dat$SumPrecip.exp.halfNorm80+1) - 0.8 # centering roughly
dat$logDegDays.exp.halfNorm80<- log(dat$DegDays.exp.halfNorm80+1) - 2 # centering roughly

dat$logAvPrecip.exp.ct<- log(dat$AvPrecip.exp+0.3) - 1 # mean(log(dat$AvPrecip.exp+0.3), na.rm= T) # 0.92
dat$AvDegDays.exp.ct<- dat$AvDegDays.exp - 15 # mean(dat$AvDegDays.exp, na.rm= T) # 14.52
dat$mlogFrostDDays.soil.exp<- -log(-dat$FrostDDays.soil.exp+1) + 4 # mean(-log(-dat$FrostDDays.soil.exp[!is.na(dat$FrostDDays.soil.exp) & dat$FrostDDays.soil.exp != 0] + 1), na.rm= T) # -4.178
dat$Frost01<- (dat$FrostDDays.soil.exp != 0) + 0


pdf("Data_explo_climate.pdf", width= 6, height= 4)
for(i in NCEP.deriv.vars.transf){
	dat2plot<- data.frame(y= dat$logN2O.cum.pos, x= dat[, i], Fertilizer= c("Control", "Fertilizer")[dat$Fert01+1], logNdays= dat$logNdays)
	p<- ggplot(dat2plot, aes(y= y, x= x, colour= logNdays)) +
	geom_point(alpha= 1, size= 0.5) +
	facet_wrap(~ Fertilizer)+
	xlab(i) + ylab("log(cN2O)") + ggtitle(i) +
	geom_smooth(col= "red", method= "loess", se= F, show.legend = F)
	print(p)
}
dev.off()


dat$maincrop4<- as.character(dat$maincrop6)
dat$maincrop4[!is.na(dat$maincrop6) & dat$maincrop6 %in% c("Legumes", "Vegetables")]<- "Veg&Legumes"
dat$maincrop4<- factor(dat$maincrop4)

dat$maincrop3<- as.character(dat$maincrop6)
dat$maincrop3[!is.na(dat$maincrop6) & dat$maincrop6 %in% c("Legumes", "Vegetables")]<- "Veg&Legumes"
dat$maincrop3[!is.na(dat$maincrop6) & dat$maincrop6 %in% c("Others", "Cereals")]<- "Cereal&Other"
dat$maincrop3<- factor(dat$maincrop3)

dat$fGrasslands<- factor(dat$maincrop6 == "Grasslands" & !is.na(dat$maincrop6))
dat$fGrasslands[is.na(dat$maincrop6)]<- NA
dat$Grasslands<- as.logical(as.character(dat$fGrasslands))

dat$Nuptake.group2<- as.character(dat$Nuptake.group)
dat$Nuptake.group2[dat$Nuptake.group %in% c("2", "3")]<- "2|3"

dat$Nuptake.group.num<- as.numeric(as.character(dat$Nuptake.group))
dat$Nuptake.group.num<- as.numeric(as.character(dat$Nuptake.group))

dat$N.appl.num<- as.numeric(factor(dat$N.appl))


# dat$NO3.avail.v2<- dat$NO3.avail
# dat$NO3.avail.v2[!is.na(dat$NO3.avail) & dat$NO3.avail %in% 2:5]<- 2
# dat$NO3.avail.v2[!is.na(dat$NO3.avail) & dat$NO3.avail > 5]<- dat$NO3.avail[!is.na(dat$NO3.avail) & dat$NO3.avail > 5] - 3

# dat$NO3.avail.v3<- dat$NO3.avail
# dat$NO3.avail.v3[!is.na(dat$NO3.avail) & dat$NO3.avail %in% 0:2]<- 1 # remove level zero and deal with it by interaction with Fert01 in model
# dat$NO3.avail.v3[!is.na(dat$NO3.avail) & dat$NO3.avail %in% 4:5]<- 2
# dat$NO3.avail.v3[!is.na(dat$NO3.avail) & dat$NO3.avail > 5]<- dat$NO3.avail[!is.na(dat$NO3.avail) & dat$NO3.avail > 5] - 2

# dat$NO3.avail.f3<- "lowNO3"
# dat$NO3.avail.f3[is.na(dat$NO3.avail)]<- NA
# dat$NO3.avail.f3[!is.na(dat$NO3.avail) & dat$NO3.avail.v3 >= 4]<- "highNO3"
# dat$NO3.avail.f3[dat$Fert01 == 0]<- "noFert"
# dat$NO3.avail.f3<- factor(dat$NO3.avail.f3)

tmp<- unique(dat[, c("NO3.avail", "NO3.avail.v2", "NO3.avail.v3", "NO3.avail.f3")])
tmp<- tmp[order(tmp$NO3.avail), ]
tmp

# make a composite of NO3 availability (from fertilizer type) and N.rate
dat$NO3.comp<- dat$N.rate * (3 - as.numeric(dat$NO3.avail.f3)) # latter codes 0 for 'noFert', 1 for 'lowFert' and 2 for 'highFert'
plot(dat$NO3.comp ~ dat$N.rate, col= as.numeric(dat$NO3.avail.f3)) # seems too ad-hoc!

# NO3 availability indicator variables
dat$highNO3<- as.numeric(dat$NO3.avail.f3 == "highNO3")
dat$lowNO3<- as.numeric(dat$NO3.avail.f3 == "lowNO3")

# pre-multiply interaction NO3.avail * N.rate
dat$logNrate.cent.lowNO3<-   dat$logNrate.cent * dat$lowNO3 # pre-multiply. s(..., by= logNrate.cent:lowNO3) was crashing gam()
dat$logNrate.cent.highNO3<-   dat$logNrate.cent * dat$highNO3

# non-centered N.rate versions starting at zero for model M15c onwards
dat$logNrate.scaled<- log(dat$N.rate + 1) / 8 # approx scaling between 0 and 1
dat$logNrate.scaled.lowNO3<-   dat$logNrate.scaled * dat$lowNO3 # pre-multiply NO3.avail (low) by N.rate. 
dat$logNrate.scaled.highNO3<-   dat$logNrate.scaled * dat$highNO3 # pre-multiply NO3.avail (high) by N.rate.

matplot(cbind(dat[, c("logNrate.cent.lowNO3", "logNrate.cent.highNO3")]))

par(mfrow= c(2, 1))
boxplot(dat$logN2O.cum.pos ~ dat$NO3.avail)
boxplot(dat$logN2O.cum.pos ~ dat$NO3.avail.v3)

# save formatted data
# save(dat, file= "Data_full_20180817.RData")
## load("Data_full_20180312.RData") # used for M9
# load("Data_full_20180817.RData")


# export list of fields as table
fieldsTable<- data.frame(Field= names(dat), Pct.missing= apply(dat, 2, FUN= function(x) round(sum(is.na(x))/nrow(dat)*100, 1)))
# write.xlsx(fieldsTable, "FieldsTableExport20180219.xlsx")

# write.xlsx(dat, "ForInfoCompleteDataExport20180329.xlsx")
write.xlsx(dat, "ForInfoCompleteDataExport20191031.xlsx")



################# World maps of predictors #################
# create file "PredictorMaps.pdf"
source("prelim data analysis/PredictorMaps_code1.txt")




################################################################################
################### Exploration of priors for missing values ###################
################################################################################

library(MASS)
fitdistr(na.omit(dat$logCarbon), "cauchy")
hist(dat$logCarbon, freq= F, nclass= 25)
mean(dat$logCarbon, na.rm= T) # 0.4852122
sd(dat$logCarbon, na.rm= T)
sd(dat$logCarbon, na.rm= T)^-2 # precision # 1.572358
lines(seq(-3, 5, l= 500), dcauchy(seq(-3, 5, l= 500), location= 0.46467060, scale= 0.43916860), col= 2)
lines(seq(-3, 5, l= 500), dnorm(seq(-3, 5, l= 500), mean= 0.4852122, sd= 0.7974882), col= 3)
# Normal prior is preferable as it fits the bulk of the data better and makes the highest values unlikely, which is what we want (unusually high carbon content should be systematically reported, so if missing it most certainly wasn't very high)

table(dat$Pre.Leg)/sum(!is.na(dat$Pre.Leg))
#      FALSE       TRUE 
# 0.93187866 0.06812134
table(dat$Per2An)/sum(!is.na(dat$Per2An))
#      FALSE       TRUE 
# 0.91616766 0.08383234

# pH
hist(dat$pH.cent, nclass= 30, freq= F, xlim= c(-4, 4))
lines(density(rbeta(100000, 1.3, 1.3)*4.5 - 2.5), col= 2)


# plot of priors
pdf("priors.pdf")
par(mfrow= c(2, 2))
# logCarbon
hist(dat$logCarbon, freq= F, nclass= 25, main= "", xlab= "logCarbon", col= grey(0.5), border= grey(0.7))
box()
axis(3, -3:4, round(exp(-3:4), 1))
mtext("Soil carbon", 3, 3)
mtext("a)", 3, 3, adj= 0)
lines(seq(-3, 5, l= 500), dnorm(seq(-3, 5, l= 500), mean= 0.4852122, sd= 0.7974882), col= 2, lwd= 2)
# pH
hist(dat$pH.cent, breaks= c(-4, -2, -1.5, -1, -0.5, 0, 0.5, 1, 2), freq= F, main= "", xlab= "pH.cent", col= grey(0.5), border= grey(0.7)) # uses classes defined from M5 onwards
box()
axis(3, -4:2, (-4:2)+7)
mtext("pH", 3, 3)
mtext("b)", 3, 3, adj= 0)
pHseq<- seq(0, 1, l= 100)
lines(pHseq*4.5 - 2.5, dbeta(pHseq, 1.3, 1.3) / 4.5, col= 2, lwd= 2)
# Pre.Leg
barplot(table(dat$Pre.Leg)/sum(!is.na(dat$Pre.Leg)), col= grey(0.5), border= grey(0.7), ylab= "Probability", xlab= "Pre.Leg")
mtext("c)", 3, 2, adj= 0)
mtext("Preceded by legumes", 3, 2)
# Per2An
barplot(table(dat$Per2An)/sum(!is.na(dat$Per2An)), col= grey(0.5), border= grey(0.7), ylab= "Probability", xlab= "Per2An")
mtext("d)", 3, 2, adj= 0)
mtext("Annual preceded by perenial", 3, 2)
dev.off()



################################################################################
#################################### MODELS ####################################
################################################################################

# exploration of data with standard gam() 
# analysis of predictors for which missing values handling returns highest value for effort
source("Rcom_exploGAM_and_MissingValuesAnalysis.txt")

######################################################################
############################# Model M1.1 #############################
######################################################################
 	
source("Rcom_M1.txt")

##################################################################
########################### Model M2 #############################
##################################################################
 	
source("Rcom_M2.txt")

##################################################################
########################### Model M3 #############################
##################################################################

# graphical exploration of more predictors from NCEP

densityplot(~ logNdays | fGrasslands, data= dat)

densityplot(~ logNrate.cent | fGrasslands, data= dat[dat$Fert01 == 1, ])

densityplot(~ N.rate | fGrasslands, data= dat[dat$Fert01 == 1, ])


pairs(dat[, c("logNdays.cent", "logSumPrecip.exp", "logDegDays.exp", "WetDays.prop.year", "meanPosDDays.air.year", "logN2O.cum.pos")], lower.panel = panel.smooth, col= rgb(0, 0, 0, 0.2)) # use temperature and precipitation over course of experiment rather than whole year

pairs(dat[, c("logSumPrecip.exp.halfNorm80", "logSumPrecip.exp", "logDegDays.exp.halfNorm80", "logDegDays.exp", "logN2O.cum.pos")], lower.panel = panel.smooth, col= rgb(0, 0, 0, 0.2)) # use temperature and precipitation over course of experiment rather than front-loaded

pairs(dat[, c("NO3.avail", "logN2O.cum.pos")], lower.panel = panel.smooth, col= rgb(0, 0, 0, 0.2))

pairs(dat[, c("logNdays.cent", "pH.cent", "logNrate.cent", "Bulk.density.imputed", "Fert01", "NO3.avail", "logSumPrecip.exp", "logDegDays.exp", "logN2O.cum.pos")], lower.panel = panel.smooth, col= rgb(0, 0, 0, 0.2))

# Model M3 setup & analysis
source("Rcom_M3.txt")



##################################################################
########################### Model M4 #############################
##################################################################

source("Rcom_M4.txt")

##################################################################
########################### Model M5 #############################
##################################################################

source("Rcom_M5.txt")

##################################################################
########################### Model M6 #############################
##################################################################

source("Rcom_M6.txt")

##################################################################
########################### Model M7 #############################
##################################################################

source("Rcom_M7.txt")

##################################################################
########################### Model M8 #############################
##################################################################

source("Rcom_M8.txt")

##################################################################
########################### Model M9 #############################
##################################################################

source("Rcom_M9.txt")

##################################################################
########################### Model M10 #############################
##################################################################

source("Rcom_M10.txt")


##################################################################
########################### Model M11 #############################
##################################################################

source("Rcom_M11.txt")

############################ explore what variance function may be appropriate ############################

source("Rcom_heteroscedasticity.txt") # using M10's residuals

# analysis very sensitive to Nslots. Not convinced any particular within-group variance heterogeneity function is required here.
# Long studies are more variable, but that's  in the control and treatment alike, therefore this is reasonably well taken care of by the studyID random effect.


##################################################################
########################### Model M12 #############################
##################################################################

source("Rcom_M12.txt")

##################################################################
########################### Model M13 #############################
##################################################################

source("Rcom_M13.txt")

##################################################################
########################### Model M14 #############################
##################################################################

source("Rcom_M14.txt")

##################################################################
########################### Model M15c #############################
##################################################################

####### Generate data sets with and without NAs (dummy imputed values for use with jagam to generate the right model matrix size)

M15c.vars<- c("logN2O.cum.pos", "pH.cent", "Grasslands", "logNrate.scaled", "logclay.imp.ct", "Fert01", "NO3.avail.v3", "studyID.f", "logDegDays.exp.compact", "NO3.avail.f3", "logSumPrecip.exp.compact", "logWetDays.exp.compact", "logCarbon", "logNrate.scaled.lowNO3", "logNrate.scaled.highNO3", "highNO3") #, "PosDDays.air.exp", "FrostDDays.soil.exp", "PosDDays.soil.exp", "Nuptake.group.num", "logFTCy.exp", "logAvPrecip.exp.ct", "AvDegDays.exp.ct"
sub.M15c<- apply(dat[, c("logN2O.cum.pos", "logNrate.cent.v2", "Fert01", "NO3.avail.v3", "studyID.f", "logWetDays.exp", "logclay.imp")], 1, FUN= function(x) !any(is.na(x))) # leave only variables for which NAs cannot be imputed, so that obs with missing data are removed # , "Bulk.density.imputed"
table(sub.M15c) 	# FALSE  TRUE 
  				#   486  2244 
dat.M15c<- dat[sub.M15c, M15c.vars]
dat.M15c.dummy<- dat.M15c # for use with jagam to generate the right model matrix size
# hard-code "Grasslands:Fert01" interaction
dat.M15c.dummy$GrassFert<- as.numeric(dat.M15c.dummy$Grasslands) * dat.M15c.dummy$Fert01


# replace NAs by any value (here median of column)
for(i in 1:ncol(dat.M15c)){
	if(is.numeric(dat.M15c.dummy[, i])){ dat.M15c.dummy[is.na(dat.M15c.dummy[, i]), i]<- median(dat.M15c.dummy[, i], na.rm= T) }
	if(is.character(dat.M15c.dummy[, i])){ dat.M15c.dummy[is.na(dat.M15c.dummy[, i]), i]<- median(dat.M15c.dummy[, i], na.rm= T) }
	if(is.factor(dat.M15c.dummy[, i])){ dat.M15c.dummy[is.na(dat.M15c.dummy[, i]), i]<- median(as.character(dat.M15c.dummy[, i]), na.rm= T) }
	if(is.logical(dat.M15c.dummy[, i])){ dat.M15c.dummy[is.na(dat.M15c.dummy[, i]), i]<- F }
}

#################### run JAGAM - Model M15c #####################
pregam.M15c<- jagam(formula= logN2O.cum.pos ~ Grasslands
    # + Fert01 #Fert01 unnecessary due to 'by= Fert01' below?
 	+ GrassFert
    # + highNO3
 	+ logCarbon
	+ s(pH.cent, bs= "cr", k= 9)
	+ te(logclay.imp.ct, logWetDays.exp.compact, k= c(5, 5))
	+ te(logclay.imp.ct, logWetDays.exp.compact, k= c(5, 5), by= logNrate.scaled.lowNO3)
	+ te(logclay.imp.ct, logWetDays.exp.compact, k= c(5, 5), by= logNrate.scaled.highNO3)
	+ s(logDegDays.exp.compact, k= 8) 
	+ s(logNrate.scaled, by= Fert01)
 	+ s(studyID.f, bs= "re"),
 	family= gaussian, data= dat.M15c.dummy, file= "./JAGS/modelM15c.txt", centred=TRUE, sp.prior = "gamma", diagonalize= FALSE)


#################### correct JAGS data #####################
# identify number of fixed effects
nfixed.M15c<- rle(colnames(pregam.M15c$jags.data$X)=="")$lengths[1] # run length encoding -> take first run of non-empty names (last is random effects)
fixed.names.M15c<- pregam.M15c$pregam$term.names[1:nfixed.M15c]

# replace dummy values in model matrix by NA where they should be (according to dat.M15c)
for(i in M15c.vars){
 	pregam.M15c$jags.data$X[is.na(dat.M15c[, i]), grep(i, pregam.M15c$pregam$term.names)]<- NA
}

# vars2impute<- apply(pregam.M15c$jags.data$X[, 1:nfixed], 2, FUN= function(x) any(is.na(x)))
vars2impute<- apply(pregam.M15c$jags.data$X, 2, FUN= function(x) any(is.na(x)))
vars2impute.colIndex<- which(vars2impute)
vars2impute.NMissing<- as.list(apply(pregam.M15c$jags.data$X[, vars2impute.colIndex], 2, FUN= function(x) sum(is.na(x))))
vars2impute.whichMissing<- list(apply(pregam.M15c$jags.data$X[, vars2impute.colIndex], 2, FUN= function(x) as.numeric(which(is.na(x)))))[[1]]
names(vars2impute.whichMissing)<- paste(gsub(":", ".", gsub("(", ".", gsub(")", "", pregam.M15c$pregam$term.names[vars2impute.colIndex], fixed= T), fixed= T)), "whichMissing", sep= ".")
names(vars2impute.NMissing)<- paste(gsub(":", ".", gsub("(", ".", gsub(")", "", pregam.M15c$pregam$term.names[vars2impute.colIndex], fixed= T), fixed= T)), "NMissing", sep= ".")

# create lookup table for pH.cent spline bases coords imputation
pH.cent.cat<- cut(dat.M15c$pH.cent, breaks= c(-4, -2, -1.5, -1, -0.5, 0, 0.5, 1, 2)) # create 8 classes of pH.cent
pH.cent.cat.prop<- round(table(pH.cent.cat)/sum(table(pH.cent.cat)), 2) # proportion of data per class
pH.cent.ref.index<- unlist(lapply(as.list(round(tapply(dat.M15c$pH.cent, pH.cent.cat, mean, na.rm= T), 1)), FUN= function(x) which.min(abs(dat.M15c$pH.cent - x)))) # index of pH.cent value closest to class mean
pH.cent.lookup<- pregam.M15c$jags.data$X[pH.cent.ref.index, grep("pH.cent", pregam.M15c$pregam$term.names)]


pregam.M15c$jags.data<- c(pregam.M15c$jags.data, vars2impute.NMissing, vars2impute.whichMissing, Grasslands.prop= list(as.vector(table(dat$Grasslands)/sum(!is.na(dat$Grasslands)))), list(pH.cent.cat.prop= as.vector(pH.cent.cat.prop)), list(pH.cent.lookup= as.matrix(pH.cent.lookup)))#,
str(pregam.M15c$jags.data)

# let JAGS auto-initialize the parametric coefficients so that each chain starts from a different point
# (in order to investigate identifiability issues with the intercept for fertilizer treatment)
pregam.M15c$jags.ini$b[1:nfixed.M15c]<- NA

#################### correct JAGS model #####################
data.frame(colIndex= vars2impute.colIndex, NMissing= unlist(vars2impute.NMissing))
#                      colIndex NMissing
# logCarbon.NMissing          4      545
# s.pH.cent.1.NMissing        5      261
# s.pH.cent.2.NMissing        6      261
# s.pH.cent.3.NMissing        7      261
# s.pH.cent.4.NMissing        8      261
# s.pH.cent.5.NMissing        9      261
# s.pH.cent.6.NMissing       10      261
# s.pH.cent.7.NMissing       11      261
# s.pH.cent.8.NMissing       12      261


#### Add to JAGS model template (NOT RUN) ####
  	# Missing data imputation with informative priors
 	for(i in 1:logCarbon.NMissing){ 
 		logC.tmp[i] ~ dnorm(0.485,1.572) # logCarbon prior
 		X[logCarbon.whichMissing[i], 4] <- logC.tmp[i]
 	}
	for(i in 1:s.pH.cent.1.NMissing){ 
		pH.index[i] ~ dcat(pH.cent.cat.prop[]) # pH.cent prior
		X[s.pH.cent.1.whichMissing[i], 5:12]  <- pH.cent.lookup[pH.index[i], ]
	}
# 	for(i in 1:GrasslandsTRUE.NMissing){ 
# 		Grasslands.tmp[i] ~ dcat(Grasslands.prop) # Grasslands.TRUE prior
# 		X[GrasslandsTRUE.whichMissing[i], 3] <- Grasslands.tmp[i] - 1
# 		X[GrasslandsTRUE.whichMissing[i], 5] <- X[GrasslandsTRUE.whichMissing[i], 3] * X[GrasslandsTRUE.whichMissing[i], 2] # Fert01:Grasslands.TRUE prior
# 	}
	
#### END NOT RUN ####

params <- c("b","rho","scale", "mu", "pH.index", "Grasslands.tmp", "tau")# , "logC.tmp"
chains <- 1
iter <- 3000
nb<- 2000

Sys.time() # 
jagsM15c.1 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
              parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
              n.chains = chains, n.iter = iter, n.burnin= nb,
              n.thin = 1, DIC= T, jags.seed= 34)
save(jagsM15c.1, file= "jagsM15c_1.RData")
Sys.time() # about 45 h for compiling + 2000 iter
jagsM15c.2 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
              parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
              n.chains = chains, n.iter = iter, n.burnin= nb,
              n.thin = 1, DIC= T, jags.seed= 83)
save(jagsM15c.2, file= "jagsM15c_2.RData")
Sys.time() # about 45 h for compiling + 2000 iter
jagsM15c.3 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
              parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
              n.chains = chains, n.iter = iter, n.burnin= nb,
              n.thin = 1, DIC= T, jags.seed= 57)
save(jagsM15c.3, file= "jagsM15c_3.RData")
Sys.time() # about 45 h for compiling + 2000 iter
jagsM15c.4 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
              parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
              n.chains = chains, n.iter = iter, n.burnin= nb,
              n.thin = 1, DIC= T, jags.seed= 389)
save(jagsM15c.4, file= "jagsM15c_4.RData")
Sys.time() # about 45 h for compiling + 2000 iter
jagsM15c.5 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
              parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
              n.chains = chains, n.iter = iter, n.burnin= nb,
              n.thin = 1, DIC= T, jags.seed= 124)
save(jagsM15c.5, file= "jagsM15c_5.RData")
Sys.time() # about 45 h for compiling + 2000 iter
jagsM15c.6 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
              parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
              n.chains = chains, n.iter = iter, n.burnin= nb,
              n.thin = 1, DIC= T, jags.seed= 720)
save(jagsM15c.6, file= "jagsM15c_6.RData")
Sys.time() # about 45 h for compiling + 2000 iter

load("jagsM15c_1.RData"); load("jagsM15c_2.RData"); load("jagsM15c_3.RData"); load("jagsM15c_4.RData"); load("jagsM15c_5.RData"); load("jagsM15c_6.RData")

jagsM15c.6c<- jagsM15c.1
jagsM15c.6c$BUGSoutput$sims.list<- mapply(rbind, jagsM15c.1$BUGSoutput$sims.list, jagsM15c.2$BUGSoutput$sims.list, jagsM15c.3$BUGSoutput$sims.list, jagsM15c.4$BUGSoutput$sims.list, jagsM15c.5$BUGSoutput$sims.list, jagsM15c.6$BUGSoutput$sims.list)
jagsM15c.6c$BUGSoutput$sims.array<- abind(	jagsM15c.1$BUGSoutput$sims.array, 
											jagsM15c.2$BUGSoutput$sims.array,
											jagsM15c.3$BUGSoutput$sims.array,
											jagsM15c.4$BUGSoutput$sims.array,
											jagsM15c.5$BUGSoutput$sims.array,
											jagsM15c.6$BUGSoutput$sims.array, along= 1)

jagsM15c.6c$n.iter<- iter*6
jagsM15c.6c$BUGSoutput$n.iter<- iter*6
jagsM15c.6c$BUGSoutput$n.burnin<- nb*6
jagsM15c.6c$BUGSoutput$n.keep<- (iter-nb)*6
jagsM15c.6c$BUGSoutput$n.sims<- (iter-nb)*6

# ##################### 12 additional re-runs of M15c with 10000 iterations each #####################
#     iter2 <- 15000
#     nb2<- 5000

#     Sys.time() # 
#     jagsM15c.11 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 2893)
#     save(jagsM15c.11, file= "jagsM15c_11.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.12 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 64)
#     save(jagsM15c.12, file= "jagsM15c_12.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.13 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 90)
#     save(jagsM15c.13, file= "jagsM15c_13.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.14 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 83)
#     save(jagsM15c.14, file= "jagsM15c_14.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.15 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 46)
#     save(jagsM15c.15, file= "jagsM15c_15.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.16 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 2378)
#     save(jagsM15c.16, file= "jagsM15c_16.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.21 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 3267)
#     save(jagsM15c.21, file= "jagsM15c_21.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.22 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 8372)
#     save(jagsM15c.22, file= "jagsM15c_22.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.23 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 645)
#     save(jagsM15c.23, file= "jagsM15c_23.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.24 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 21)
#     save(jagsM15c.24, file= "jagsM15c_24.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.25 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 12)
#     save(jagsM15c.25, file= "jagsM15c_25.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter
#     jagsM15c.26 <- jags(model.file = "./JAGS/modelM15c_NAimputation.txt", data = pregam.M15c$jags.data,
#                 parameters.to.save = params, inits = list(pregam.M15c$jags.ini),
#                 n.chains = chains, n.iter = iter2, n.burnin= nb2,
#                 n.thin = 1, DIC= T, jags.seed= 67)
#     save(jagsM15c.26, file= "jagsM15c_26.RData")
#     Sys.time() # about 45 h for compiling + 2000 iter

#     load("jagsM15c_1.RData"); load("jagsM15c_2.RData"); load("jagsM15c_3.RData"); load("jagsM15c_4.RData"); load("jagsM15c_5.RData"); load("jagsM15c_6.RData")
#     load("jagsM15c_11.RData"); load("jagsM15c_12.RData"); load("jagsM15c_13.RData"); load("jagsM15c_14.RData"); load("jagsM15c_15.RData"); load("jagsM15c_16.RData")
#     load("jagsM15c_21.RData"); load("jagsM15c_22.RData"); load("jagsM15c_23.RData"); load("jagsM15c_24.RData"); load("jagsM15c_25.RData"); load("jagsM15c_26.RData")


#     jagsM15c.6c<- jagsM15c.1
#     jagsM15c.6c$BUGSoutput$sims.list<- mapply(rbind, jagsM15c.1$BUGSoutput$sims.list, jagsM15c.2$BUGSoutput$sims.list, jagsM15c.3$BUGSoutput$sims.list, jagsM15c.4$BUGSoutput$sims.list, jagsM15c.5$BUGSoutput$sims.list, jagsM15c.6$BUGSoutput$sims.list,
#     jagsM15c.11$BUGSoutput$sims.list, jagsM15c.12$BUGSoutput$sims.list, jagsM15c.13$BUGSoutput$sims.list, jagsM15c.14$BUGSoutput$sims.list, jagsM15c.15$BUGSoutput$sims.list, jagsM15c.16$BUGSoutput$sims.list,
#     jagsM15c.21$BUGSoutput$sims.list, jagsM15c.22$BUGSoutput$sims.list, jagsM15c.23$BUGSoutput$sims.list, jagsM15c.24$BUGSoutput$sims.list, jagsM15c.25$BUGSoutput$sims.list, jagsM15c.26$BUGSoutput$sims.list)
#     jagsM15c.6c$BUGSoutput$sims.array<- abind(	jagsM15c.1$BUGSoutput$sims.array, 
#                                                 jagsM15c.2$BUGSoutput$sims.array,
#                                                 jagsM15c.3$BUGSoutput$sims.array,
#                                                 jagsM15c.4$BUGSoutput$sims.array,
#                                                 jagsM15c.5$BUGSoutput$sims.array,
#                                                 jagsM15c.6$BUGSoutput$sims.array,
#                                                 jagsM15c.11$BUGSoutput$sims.array, 
#                                                 jagsM15c.12$BUGSoutput$sims.array,
#                                                 jagsM15c.13$BUGSoutput$sims.array,
#                                                 jagsM15c.14$BUGSoutput$sims.array,
#                                                 jagsM15c.15$BUGSoutput$sims.array,
#                                                 jagsM15c.16$BUGSoutput$sims.array,
#                                                 jagsM15c.21$BUGSoutput$sims.array, 
#                                                 jagsM15c.22$BUGSoutput$sims.array,
#                                                 jagsM15c.23$BUGSoutput$sims.array,
#                                                 jagsM15c.24$BUGSoutput$sims.array,
#                                                 jagsM15c.25$BUGSoutput$sims.array,
#                                                 jagsM15c.26$BUGSoutput$sims.array, along= 1)

#     jagsM15c.6c$n.iter<- iter*6 + iter2*12
#     jagsM15c.6c$BUGSoutput$n.iter<- iter*6 + iter2*12
#     jagsM15c.6c$BUGSoutput$n.burnin<- nb*6 + nb2*12
#     jagsM15c.6c$BUGSoutput$n.keep<- (iter-nb)*6 + (iter2-nb2)*12
#     jagsM15c.6c$BUGSoutput$n.sims<- (iter-nb)*6 + (iter2-nb2)*12

###################################### end additional re-runs ######################################

samM15c<- jags2sam(jagsM15c.6c)
jamM15c<- sim2jam(samM15c,pregam.M15c$pregam)

pdf("GAM.M15coutputBayes_splines.pdf")
	plot(jamM15c, pages= 0, scale= 0, scheme= 2)
dev.off()

pdf("GAM.M15coutputBayes_splines_SEwithMean.pdf")
	plot(jamM15c, pages= 0, scale= 0, scheme= 2,seWithMean= T)
dev.off()

# pdf("GAM.M15coutputBayes_splines_CI.pdf")
# 	plot(jamM15c, pages= 0, scale= 0, scheme= 2, se= T, pers= T)
# dev.off()

fx.index.M15c<- which(dimnames(jagsM15c.6c$BUGSoutput$sims.array)[[3]] %in% paste("b[", 1:nfixed.M15c, "]", sep= ""))
dimnames(jagsM15c.6c$BUGSoutput$sims.array)[[3]][fx.index.M15c]<- fixed.names.M15c

jagsM15c.mcmc <- as.mcmc(jagsM15c.6c)

summ.M15c<- summary(jagsM15c.mcmc[[1]][, fixed.names.M15c], quantiles = c(0.025, 0.5, 0.975))
summ.M15c
test.M15c<- matrix(0 < summ.M15c$quantiles[, "2.5%"] | 0 > summ.M15c$quantiles[, "97.5%"])
data.frame(fixed.names.M15c, test.M15c)

# residual variance estimate
summary(jagsM15c.mcmc[[1]][, "tau"], quantiles = c(0.025, 0.5, 0.975))

tau.est<- summary(jagsM15c.mcmc[[1]][, "tau"], quantiles = c(0.025, 0.5, 0.975))$statistics["Mean"]
var.est<- 1/tau.est
var.est # 0.07026771

# Random Effect variance estimate - did not monitor lambda...
summary(jagsM15c.mcmc[[1]][, "rho[16]"], quantiles = c(0.025, 0.5, 0.975))

rho.RE.est<- summary(jagsM15c.mcmc[[1]][, "rho[16]"], quantiles = c(0.025, 0.5, 0.975))$statistics["Mean"]
var.RE.est<- 1/exp(rho.RE.est)
var.RE.est # 0.2482381

pdf("GAM.M15coutputBayes_traces_fixed.pdf", width= 7, height= 21)
	xyplot(jagsM15c.mcmc[[1]][, fixed.names.M15c])
dev.off()

pdf("GAM.M15coutputBayes_densityplots_fixed.pdf", width= 7, height= 15)
	densityplot(jagsM15c.mcmc[[1]][, fixed.names.M15c])
dev.off()

pdf("GAM.M15coutputBayes_traces_SplineCoefs.pdf", width= 21, height= 10)
	NsplineCoefs.M15c<- length(grep("b[", colnames(jagsM15c.mcmc[[1]]), fixed= T)) + length(fixed.names.M15c) - length(grep("studyID", pregam.M15c$pregam$term.names))
	SplCoefsNames.M15c<- paste("b[", (length(fixed.names.M15c)+1):(NsplineCoefs.M15c), "]", sep= "")
	datmp.M15c<- jagsM15c.mcmc[[1]][, SplCoefsNames.M15c]
	colnames(datmp.M15c)<- pregam.M15c$pregam$term.names[(length(fixed.names.M15c)+1):(NsplineCoefs.M15c)]
	xyplot(datmp.M15c, layout= c(5, 5, ceiling(length(SplCoefsNames.M15c)/25)))
dev.off()
rm(datmp.M15c)

pdf("GAM.M15coutputBayes_traces_MissingValues.pdf", width= 21, height= 10)
	NmuCoefs.M15c<- length(grep("studyID", pregam.M15c$pregam$term.names))
	datmp.mu.M15c<- jagsM15c.mcmc[[1]][, paste("b[", grep("studyID", pregam.M15c$pregam$term.names), "]", sep= "")]
	colnames(datmp.mu.M15c)<- pregam.M15c$pregam$term.names[grep("studyID", pregam.M15c$pregam$term.names)]
	xyplot(datmp.mu.M15c, layout= c(5, 5, ceiling(NmuCoefs.M15c/25)))
	dev.off()
rm(datmp.mu.M15c)

N.mu<- pregam.M15c$jags.data$n # number of predicted values (number of observations)
M15c.predicted.1<- apply(jagsM15c.mcmc[[1]][, paste("mu[", 1:N.mu, "]", sep= "")], 2, mean) # extract posterior mean of predicted values (on log scale)
# level 1 residuals (around random effects)
M15c.observed<- pregam.M15c$jags.data$y
M15c.residual.1<- M15c.observed - M15c.predicted.1 # residual's posterior mean
# random effects
M15c.re.coef<- apply(jagsM15c.mcmc[[1]][, paste("b[", grep("studyID", pregam.M15c$pregam$term.names), "]", sep= "")], 2, mean)
M15c.re.pred<- pregam.M15c$pregam$X[, grep("studyID", pregam.M15c$pregam$term.names)] %*% M15c.re.coef
# level 0 residuals (around fixed effects, adding studyID random effects)
M15c.predicted.0<- M15c.predicted.1 - M15c.re.pred
M15c.residual.0<- M15c.observed - M15c.predicted.0  # M15c.residual.1 + M15c.re.pred

pdf("GAM.M15c_StandardDiagnostics_Fig4.pdf", width= 7, height= 6)
	par(mfrow= c(2, 2), mar= c(4.1, 4.1, 3.1, 1.1))
	# observed ~ fitted.0
	plot(M15c.predicted.0, M15c.observed, col= rgb(0, 0, 0, alpha= 0.35), xlab= "Level-0 predictions", ylab= "Observed")
	abline(0, 1, col= 2)
	mtext("a)", 3, 0.6, adj= 0)
	cor(M15c.observed, M15c.predicted.0) # 0.6363811
	# observed ~ fitted.1
	plot(M15c.predicted.1, M15c.observed, col= rgb(0, 0, 0, alpha= 0.35), xlab= "Level-1 predictions", ylab= "Observed") # overfitting
	abline(0, 1, col= 2)
	mtext("b)", 3, 0.6, adj= 0)
	cor(M15c.observed, M15c.predicted.1) # 0.9504366
	# residual ~ fitted
	plot(M15c.residual.0 ~ M15c.predicted.0, col= rgb(0, 0, 0, alpha= 0.35), xlab= "Level-0 predictions", ylab= "Level-0 residuals")
	abline(h= 0, col= 2)
	mtext("c)", 3, 0.6, adj= 0)
	# distribution of random effects
	M15c.re<- apply(jagsM15c.mcmc[[1]][, paste("b[", grep("studyID.f", colnames(pregam.M15c$jags.data$X)), "]", sep= "")], 2, mean)
	hist(M15c.re, nclass= 25, main= "", xlab= "Experiment-level random effect")
	mtext("d)", 3, 0.6, adj= 0)
dev.off()

# comparison of level-1 and level-0 dispersion
plot(M15c.residual.1 ~ M15c.predicted.0, col= rgb(0, 0, 0, alpha= 0.35), xlab= "Level-0 predictions", ylab= "Level-1 residuals")
abline(h= 0, col= 2)
plot(M15c.residual.0 ~ M15c.predicted.0, col= rgb(0, 0, 0, alpha= 0.35), xlab= "Level-0 predictions", ylab= "Level-0 residuals")
abline(h= 0, col= 2)

source("PredictorsList9.txt")

pdf("GAM.M15cExplo_ResidualsPredictors.pdf", width= 6, height= 6)
for(i in VOI2){
	dat2plot<- data.frame(y= M15c.residual.0, x= dat[sub.M15c, i])
	p<- ggplot(dat2plot, aes(y= y, x= x)) +
	geom_point(col= rgb(0, 0, 0, 0.2)) +
	xlab(i) + ylab("Model residuals") + ggtitle(i) +
	geom_hline(yintercept= 0, col= "green") +
	geom_smooth(col= "red", method= "loess", se= F, )
	print(p)
}
for(i in groups){
	dat2plot<- data.frame(y= M15c.residual.0, x= factor(dat[sub.M15c, i]))
	p<- ggplot(dat2plot, aes(y= y, x= x)) +
	geom_boxplot(col= rgb(0, 0, 0, 0.5)) +
	geom_hline(yintercept= 0, col= "green") +
	xlab(i) + ylab("Model residuals") + ggtitle(i) +
	theme(axis.text.x=element_text(angle = -90, hjust = 0), plot.title = element_text(hjust = 0.5))
	print(p)

}
dev.off()


pdf("GAM.M15cExplo_ResidualsStructure.pdf", width= 12, height= 4.5)
dat2plot<- data.frame(M15c.residual.0= M15c.residual.0, M15c.predicted.0= M15c.predicted.0, Fertilizer= c("Control", "Fertilizer")[dat$Fert01[sub.M15c]+1], logNdays= dat$logNdays[sub.M15c], logN.rate= dat$logNrate.cent[sub.M15c], NO3avail= c("lowNO3", "highNO3")[dat$highNO3[sub.M15c]+1])
	
p<- ggplot(dat2plot, aes(y= M15c.residual.0, x= M15c.predicted.0, colour= logNdays)) +
	geom_point(alpha= 1, size= 1) +
	facet_wrap(~ Fertilizer * NO3avail) +
	xlab("Fitted values") + ylab("Model residuals")
print(p)

p<- ggplot(dat2plot, aes(y= M15c.residual.0, x= M15c.predicted.0, colour= logN.rate)) +
	geom_point(alpha= 1, size= 1) +
	facet_wrap(~ Fertilizer * NO3avail) +
	xlab("Fitted values") + ylab("Model residuals")
print(p)

dev.off()


####### Predictions from M15c

extent[extent$LAT > 48 & extent$LAT < 54 & extent$LON > 6 & extent$LON < 15 & extent$DURATION.days > 363 & extent$DURATION.days <367, ] # subset climate data for Germany / 365 days

hist(extent[extent$DURATION.days > 363 & extent$DURATION.days <367, "DegDays.exp"]) # subset climate data for 365 days # range ~ [1000, 10000]
log(c(0, 1000, 10000) + 1500) - 8 # -0.6867796 -0.1759540  1.3501023

hist(extent[, "DegDays.exp"], nclass= 500, xlim= c(0, 1000))
log(c(0, 1000, 10000) + 1500) - 8 # -0.6867796 -0.1759540  1.3501023

hist(extent[extent$DURATION.days > 363 & extent$DURATION.days <367, "WetDays.exp"]) # subset climate data for 365 days
range(extent[extent$DURATION.days > 363 & extent$DURATION.days <367, "WetDays.exp"], na.rm= T) # subset climate data for 365 days # range ~ [30, 200]


df.M15c.complete<- function(x){ # function to compute transformed variables for prediction
    x$logCarbon<- log(x$SOC + 1)
    x$pH.cent<- x$pH - 7
	x$GrassFert<- as.numeric(x$Grasslands) * x$Fert01
    x$logclay.imp.ct<- log(x$Clay / 100 + 1) - 0.24
    x$logWetDays.exp.compact<- log(x$WetDays.exp+10) - 4
    x$logDegDays.exp.compact<- log(x$DegDays.exp+1500) - 8
    x$logNrate.scaled<- log(x$N.rate + 1) / 8
    x$logNrate.scaled.lowNO3<- x$logNrate.scaled * x$Fert01 * x$lowNO3
    x$logNrate.scaled.highNO3<- x$logNrate.scaled * x$Fert01 * x$highNO3
    x
}

# !!!!!! do not compare predictions from Fert01 = 1 and N.rate = 0 with prediction from Fert01 = 0, 
# because N.rate = 0 is out of the range of observed values (which is logNrate.cent.v2 ~ {-2, 2}) when Fert01 == 1
# logNrate.cent.v2 ~ {-2, 2} equivalent to N.rate ~ {20, 1100}
# instead, restrict range of N.rate to observed values (minimum 20 kg) and compare with Fert01 = 0

## making predictions to illustrate model outputs
# Nrate only
ND.M15c.Nrate<- data.frame(Fert01= 1, lowNO3= 1, highNO3= 0, Grasslands = FALSE, 
	pH= 6.8, Clay= 25, SOC= 2, # pH @ median, Clay @ median, SOC @ 2%
	WetDays.exp= 120, DegDays.exp= 3600, # @ Germany over 1 year
    N.rate= seq(0, 1000, l= 50), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
ND.M15c.Nrate<- df.M15c.complete(ND.M15c.Nrate)
plot(exp(predict(jamM15c, newdata= ND.M15c.Nrate)) ~ ND.M15c.Nrate$N.rate)

# pH only
ND.M15c.pH<- data.frame(Fert01= 0, lowNO3= 0, highNO3= 0, Grasslands = FALSE, 
	pH= seq(5, 9, l= 25), Clay= 25, SOC= 2, # pH @ median, Clay @ median, SOC @ 2%
	WetDays.exp= 120, DegDays.exp= 3600, # @ Germany over 1 year
    N.rate= 200, studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
ND.M15c.pH<- df.M15c.complete(ND.M15c.pH)
plot(exp(predict(jamM15c, newdata= ND.M15c.pH)) ~ ND.M15c.pH$pH)

# do not do this (cf note above):
tmp<- data.frame(Fert01= c(0, 1), lowNO3= 0, highNO3= 0, Grasslands = FALSE, 
	Clay= 25, SOC= 2, # pH @ median, Clay @ median, SOC @ 2%
	WetDays.exp= 120, DegDays.exp= 3600, # @ Germany over 1 year
    pH= 7, N.rate= 0, studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
tmp<- df.M15c.complete(tmp)
predict(jamM15c, newdata= tmp, se.fit= F) # 1.1247374 0.6866758 this is problematic
# try this instead:
tmp<- data.frame(Fert01= c(0, 1, 1), lowNO3= c(0, 0, 1), highNO3= 0, Grasslands = FALSE, 
	Clay= 25, SOC= 2, # pH @ median, Clay @ median, SOC @ 2%
	WetDays.exp= 120, DegDays.exp= 3600, # @ Germany over 1 year
    pH= 7, N.rate= c(0, 20, 20), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
tmp<- df.M15c.complete(tmp) # 1.124737 1.099385 1.115259 is about okay: don't want second to be too different from first
predict(jamM15c, newdata= tmp, se.fit= F)

# pH and Grassland (baselines pH 7 and annual crops)
ND.M15c.basepH<- data.frame(Fert01= 0, lowNO3= 0, highNO3= 0, Grasslands = FALSE, 
	Clay= 25, SOC= 2, # pH @ median, Clay @ median, SOC @ 2%
	WetDays.exp= 120, DegDays.exp= 3600, # @ Germany over 1 year
    expand.grid(pH= seq(7, 7, l= 4), N.rate= rep(20, l= 50)), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
ND.M15c.basepH<- df.M15c.complete(ND.M15c.basepH)
range(predict(jamM15c, newdata= ND.M15c.basepH, se.fit= F))

ND.M15c.pH<- data.frame(Fert01= 1, lowNO3= 1, highNO3= 0, Grasslands = FALSE, 
	Clay= 25, SOC= 2, # pH @ median, Clay @ median, SOC @ 2%
	WetDays.exp= 120, DegDays.exp= 3600, # @ Germany over 1 year
    expand.grid(pH= seq(5, 8, l= 4), N.rate= seq(20, 400, l= 50)), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
ND.M15c.pH<- df.M15c.complete(ND.M15c.pH)
range(predict(jamM15c, newdata= ND.M15c.pH, se.fit= F))

ND.M15c.pHGrass<- data.frame(Fert01= 1, lowNO3= 1, highNO3= 0, Grasslands = TRUE, 
	Clay= 25, SOC= 2, # pH @ median, Clay @ median, SOC @ 2%
	WetDays.exp= 120, DegDays.exp= 3600, # @ Germany over 1 year
    expand.grid(pH= seq(5, 8, l= 4), N.rate= seq(20, 400, l= 50)), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
ND.M15c.pHGrass<- df.M15c.complete(ND.M15c.pHGrass)
range(predict(jamM15c, newdata= ND.M15c.pHGrass, se.fit= F))

ND.M15c.basepH.pred<- predict(jamM15c, newdata= ND.M15c.basepH, se.fit= F)
ND.M15c.pH$pred<- predict(jamM15c, newdata= ND.M15c.pH, se.fit= T)$fit - ND.M15c.basepH.pred # remove baseline (ignore baseline prediction error)
ND.M15c.pHGrass$pred<- predict(jamM15c, newdata= ND.M15c.pHGrass, se.fit= T)$fit - ND.M15c.basepH.pred # remove baseline (ignore baseline prediction error)
ND.M15c.pH$pred.se<- predict(jamM15c, newdata= ND.M15c.pH, se.fit= T)$se.fit
ND.M15c.pH$pred.exp<- exp(ND.M15c.pH$pred)
ND.M15c.pH$fpH<- factor(ND.M15c.pH$pH)
ND.M15c.pHGrass$pred.exp<- exp(ND.M15c.pHGrass$pred)
ND.M15c.pHGrass$fpH<- factor(ND.M15c.pHGrass$pH)

pdf("M15c_Fig_Nrate_pH_cropType.pdf", width= 6, height= 6)
pal.4<- c("#009E73", "#D55E00", "#0072B2", "#CC79A7") # c("#a6cee3", "#1f78b4", "#33a02c", "#b2df8a")
plot(ND.M15c.pH$pred.exp ~ ND.M15c.pH$N.rate, col= NULL, xlab= "N application rate (kg)", ylab= "Multiplier for (cumN2O+1)")
by(ND.M15c.pH, ND.M15c.pH$fpH, function(x) {lines(x$pred.exp ~ x$N.rate, col= pal.4[x$fpH], lwd= 3)})
by(ND.M15c.pHGrass, ND.M15c.pHGrass$fpH, function(x) {lines(x$pred.exp ~ x$N.rate, col= pal.4[x$fpH], lwd= 3, lty= 3)})
legend(400, 1.3, legend= c("pH= 5", "pH= 6", "pH= 7", "pH= 8"), fill= pal.4, border= "transparent", bty= "n", ncol= 2, yjust= 0)
legend(400, 1.3, legend= c("Annual crop", "Grassland"), lty= c(1, 3), lwd= 3, col= c(1), bty= "n", ncol= 2, yjust= 1)
dev.off()

# SOC (baseline SOC = 2%)
ND.M15c.baseSOC<- data.frame(Fert01= 0, lowNO3= 0, highNO3= 0, Grasslands = FALSE, 
	Clay= 25, pH= 7, # SOC @ median, Clay @ median, SOC @ 2%
	WetDays.exp= 120, DegDays.exp= 3600, # @ Germany over 1 year
    expand.grid(SOC= seq(2, 2, l= 5), N.rate= rep(20, l= 50)), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
ND.M15c.baseSOC<- df.M15c.complete(ND.M15c.baseSOC)
range(predict(jamM15c, newdata= ND.M15c.baseSOC, se.fit= F))

ND.M15c.SOC<- data.frame(Fert01= 1, lowNO3= 1, highNO3= 0, Grasslands = FALSE, 
	Clay= 25, pH= 7, # SOC @ median, Clay @ median, SOC @ 2%
	WetDays.exp= 120, DegDays.exp= 3600, # @ Germany over 1 year
    expand.grid(SOC= seq(2, 6, l= 5), N.rate= seq(20, 400, l= 50)), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
ND.M15c.SOC<- df.M15c.complete(ND.M15c.SOC)
range(predict(jamM15c, newdata= ND.M15c.SOC, se.fit= F))

ND.M15c.baseSOC.pred<- predict(jamM15c, newdata= ND.M15c.baseSOC, se.fit= F)
ND.M15c.SOC$pred<- as.vector(predict(jamM15c, newdata= ND.M15c.SOC, se.fit= T)$fit) - as.vector(ND.M15c.baseSOC.pred) # remove baseline (ignore baseline prediction error)
ND.M15c.SOC$pred.se<- predict(jamM15c, newdata= ND.M15c.SOC, se.fit= T)$se.fit
ND.M15c.SOC$pred.exp<- exp(ND.M15c.SOC$pred)
ND.M15c.SOC$fSOC<- factor(ND.M15c.SOC$SOC)

pdf("M15c_Fig_Nrate_pH_SOC_cropType.pdf", width= 9, height= 5)
par(mfrow= c(1, 2))
pal.4<- c("#a6cee3", "#1f78b4", "#33a02c", "#b2df8a")
plot(ND.M15c.pH$pred.exp ~ ND.M15c.pH$N.rate, col= NULL, xlab= "N application rate (kg)", ylab= "Multiplier for (cumN2O+1)", ylim= c(0.7, 3))
by(ND.M15c.pH, ND.M15c.pH$fpH, function(x) {lines(x$pred.exp ~ x$N.rate, col= pal.4[x$fpH], lwd= 3)})
# by(ND.M15c.pHGrass, ND.M15c.pHGrass$fpH, function(x) {lines(x$pred.exp ~ x$N.rate, col= pal.4[x$fpH], lwd= 3, lty= 3)})
legend(150, 1, legend= c("pH= 5", "pH= 6", "pH= 7", "pH= 8"), fill= pal.4, border= "transparent", bty= "n", ncol= 2, yjust= 0, cex= 0.75)
# legend(100, 1, legend= c("Annual crop", "Grassland"), lty= c(1, 3), lwd= 3, col= c(1), bty= "n", ncol= 2, yjust= 1, cex= 0.75)

pal.5<- c("#009E73", "#D55E00", "#0072B2", "#CC79A7", "#56B4E9") # c("#a6cee3", "#1f78b4", "#33a02c", "#b2df8a")
plot(ND.M15c.SOC$pred.exp ~ ND.M15c.SOC$N.rate, col= NULL, xlab= "N application rate (kg)", ylab= "Multiplier for (cumN2O+1)", ylim= c(0.7, 3))
by(ND.M15c.SOC, ND.M15c.SOC$fSOC, function(x) {lines(x$pred.exp ~ x$N.rate, col= pal.5[x$fSOC], lwd= 3)})
legend(150, 1.5, legend= c("SOC= 2", "SOC= 3", "SOC= 4", "SOC= 5", "SOC= 6"), fill= pal.5, border= "transparent", bty= "n", ncol= 2, yjust= 1, cex= 0.75)
dev.off()


# Soil (Clay) * Wetness (WetDays.exp) interaction # baseline is Clay= 25% and WetDays.exp= 0 days
# repeat for N.rate %in% c(50, 150, 500) kg/Ha
for(i in c(50, 150, 500)){
    N.sq<- 20
    N.rate.p<- i
    ND.M15c.baseWATERLOG<- data.frame(Fert01= 0, lowNO3= 0, highNO3= 0, Grasslands = FALSE, 
        SOC= 2, pH= 7, # SOC @ median, SOC @ 2%
        N.rate= N.rate.p, DegDays.exp= 3600, # @ Germany over 1 year
        expand.grid(Clay= seq(25, 25, l= N.sq), WetDays.exp= rep(0, l= N.sq)), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
    ND.M15c.baseWATERLOG<- df.M15c.complete(ND.M15c.baseWATERLOG)

    sq.clay<- seq(0, 60, l= N.sq)
    sq.WetDays<- seq(0, 200, l= N.sq)

    ND.M15c.WATERLOG.common<- data.frame(Fert01= 0, lowNO3= 0, highNO3= 0, Grasslands = FALSE, 
        SOC= 2, pH= 7, # SOC @ median, SOC @ 2%
        N.rate= N.rate.p, DegDays.exp= 3600, # @ Germany over 1 year
        expand.grid(Clay= sq.clay, WetDays.exp= sq.WetDays), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
    ND.M15c.WATERLOG.common<- df.M15c.complete(ND.M15c.WATERLOG.common)

    ND.M15c.WATERLOG.low<- data.frame(Fert01= 1, lowNO3= 1, highNO3= 0, Grasslands = FALSE, 
        SOC= 2, pH= 7, # SOC @ median, SOC @ 2%
        N.rate= N.rate.p, DegDays.exp= 3600, # @ Germany over 1 year
        expand.grid(Clay= sq.clay, WetDays.exp= sq.WetDays), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
    ND.M15c.WATERLOG.low<- df.M15c.complete(ND.M15c.WATERLOG.low)

    ND.M15c.WATERLOG.high<- data.frame(Fert01= 1, lowNO3= 0, highNO3= 1, Grasslands = FALSE, 
        SOC= 2, pH= 7, # SOC @ median, SOC @ 2%
        N.rate= N.rate.p, DegDays.exp= 3600, # @ Germany over 1 year
        expand.grid(Clay= sq.clay, WetDays.exp= sq.WetDays), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
    ND.M15c.WATERLOG.high<- df.M15c.complete(ND.M15c.WATERLOG.high)

    ND.M15c.baseWATERLOG.pred<- predict(jamM15c, newdata= ND.M15c.baseWATERLOG, se.fit= F)

    ND.M15c.WATERLOG.common$pred<- as.vector(predict(jamM15c, newdata= ND.M15c.WATERLOG.common, se.fit= T)$fit) - as.vector(ND.M15c.baseWATERLOG.pred) # remove baseline (ignore baseline prediction error)
    ND.M15c.WATERLOG.common$pred.exp<- exp(ND.M15c.WATERLOG.common$pred)
    ND.M15c.WATERLOG.common.mat<- matrix(ND.M15c.WATERLOG.common$pred.exp, length(sq.clay), length(sq.WetDays))

    ND.M15c.WATERLOG.low$pred<- as.vector(predict(jamM15c, newdata= ND.M15c.WATERLOG.low, se.fit= T)$fit) - as.vector(ND.M15c.baseWATERLOG.pred) # remove baseline (ignore baseline prediction error)
	ND.M15c.WATERLOG.low$pred.exp<- exp(ND.M15c.WATERLOG.low$pred)
    ND.M15c.WATERLOG.low.mat<- matrix(ND.M15c.WATERLOG.low$pred.exp, length(sq.clay), length(sq.WetDays))

    ND.M15c.WATERLOG.high$pred<- as.vector(predict(jamM15c, newdata= ND.M15c.WATERLOG.high, se.fit= T)$fit) - as.vector(ND.M15c.baseWATERLOG.pred) # remove baseline (ignore baseline prediction error)
	ND.M15c.WATERLOG.high$pred.exp<- exp(ND.M15c.WATERLOG.high$pred)
    ND.M15c.WATERLOG.high.mat<- matrix(ND.M15c.WATERLOG.high$pred.exp, length(sq.clay), length(sq.WetDays))

    rng<- c(range(c(ND.M15c.WATERLOG.common$pred.exp, ND.M15c.WATERLOG.low$pred.exp, ND.M15c.WATERLOG.high$pred.exp)))
    
    pal.3<- c("#009E73", "#D55E00", "#0072B2")
    pdf(paste("M15c_Fig_InterWATERLOG_Nrate", N.rate.p, ".pdf", sep= ""), width= 8, height= 8)
    par(mfrow= c(1, 1))
    theta.p<- -40
    phi.p<-  15
    persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.low.mat, col= NULL, border= NULL, xlab= "% Clay", ylab= "Number of wet days", zlab= "Multiplier for (cumN2O+1)", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p) -> res
    round(res, 3)
    sel<- dat$WetDays.exp < 200
    points(trans3d(jitter(dat$clay[sel], 2), jitter(dat$WetDays.exp[sel], 6), min(rng), pmat = res), col = rgb(0, 0, 0, 0.3), pch = 16, cex= 0.4)
    par(new= T)
    persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.common.mat, col= NULL, border= pal.3[1], xlab= "", ylab= "", zlab= "", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p)
    par(new= T)
    persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.high.mat, col= NULL, border= pal.3[2], xlab= "", ylab= "", zlab= "", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p)
    par(new= T)
    persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.low.mat, col= NULL, border= pal.3[3], xlab= "", ylab= "", zlab= "", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p) -> res

    theta.p<- 60
    phi.p<-  15
    persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.low.mat, col= NULL, border= NULL, xlab= "% Clay", ylab= "Number of wet days", zlab= "Multiplier for (cumN2O+1)", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p) -> res
    round(res, 3)
    sel<- dat$WetDays.exp < 200
    points(trans3d(jitter(dat$clay[sel], 2), jitter(dat$WetDays.exp[sel], 6), min(rng), pmat = res), col = rgb(0, 0, 0, 0.3), pch = 16, cex= 0.4)
    par(new= T)
    persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.common.mat, col= NULL, border= pal.3[1], xlab= "", ylab= "", zlab= "", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p)
    par(new= T)
    persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.high.mat, col= NULL, border= pal.3[2], xlab= "", ylab= "", zlab= "", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p)
    par(new= T)
    persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.low.mat, col= NULL, border= pal.3[3], xlab= "", ylab= "", zlab= "", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p) -> res

    dev.off()
}

# Soil (Clay) * Wetness (WetDays.exp) interaction # baseline is Clay= 25% and WetDays.exp= 0 days
# try different viewing angles
    N.sq<- 20
    N.rate.p<- 150
    ND.M15c.baseWATERLOG<- data.frame(Fert01= 0, lowNO3= 0, highNO3= 0, Grasslands = FALSE, 
        SOC= 2, pH= 7, # SOC @ median, SOC @ 2%
        N.rate= N.rate.p, DegDays.exp= 3600, # @ Germany over 1 year
        expand.grid(Clay= seq(25, 25, l= N.sq), WetDays.exp= rep(0, l= N.sq)), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
    ND.M15c.baseWATERLOG<- df.M15c.complete(ND.M15c.baseWATERLOG)

    sq.clay<- seq(0, 60, l= N.sq)
    sq.WetDays<- seq(0, 200, l= N.sq)

    ND.M15c.WATERLOG.common<- data.frame(Fert01= 0, lowNO3= 0, highNO3= 0, Grasslands = FALSE, 
        SOC= 2, pH= 7, # SOC @ median, SOC @ 2%
        N.rate= N.rate.p, DegDays.exp= 3600, # @ Germany over 1 year
        expand.grid(Clay= sq.clay, WetDays.exp= sq.WetDays), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
    ND.M15c.WATERLOG.common<- df.M15c.complete(ND.M15c.WATERLOG.common)

    ND.M15c.WATERLOG.low<- data.frame(Fert01= 1, lowNO3= 1, highNO3= 0, Grasslands = FALSE, 
        SOC= 2, pH= 7, # SOC @ median, SOC @ 2%
        N.rate= N.rate.p, DegDays.exp= 3600, # @ Germany over 1 year
        expand.grid(Clay= sq.clay, WetDays.exp= sq.WetDays), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
    ND.M15c.WATERLOG.low<- df.M15c.complete(ND.M15c.WATERLOG.low)

    ND.M15c.WATERLOG.high<- data.frame(Fert01= 1, lowNO3= 0, highNO3= 1, Grasslands = FALSE, 
        SOC= 2, pH= 7, # SOC @ median, SOC @ 2%
        N.rate= N.rate.p, DegDays.exp= 3600, # @ Germany over 1 year
        expand.grid(Clay= sq.clay, WetDays.exp= sq.WetDays), studyID.f= " 32.58 119.702007 8  8120082007-08-152007-11-04") # some random studyID.f (the first one)
    ND.M15c.WATERLOG.high<- df.M15c.complete(ND.M15c.WATERLOG.high)

    ND.M15c.baseWATERLOG.pred<- predict(jamM15c, newdata= ND.M15c.baseWATERLOG, se.fit= F)

    ND.M15c.WATERLOG.common$pred<- as.vector(predict(jamM15c, newdata= ND.M15c.WATERLOG.common, se.fit= T)$fit) - as.vector(ND.M15c.baseWATERLOG.pred) # remove baseline (ignore baseline prediction error)
    ND.M15c.WATERLOG.common$pred.exp<- exp(ND.M15c.WATERLOG.common$pred)
    ND.M15c.WATERLOG.common.mat<- matrix(ND.M15c.WATERLOG.common$pred.exp, length(sq.clay), length(sq.WetDays))

    ND.M15c.WATERLOG.low$pred<- as.vector(predict(jamM15c, newdata= ND.M15c.WATERLOG.low, se.fit= T)$fit) - as.vector(ND.M15c.baseWATERLOG.pred) # remove baseline (ignore baseline prediction error)
	ND.M15c.WATERLOG.low$pred.exp<- exp(ND.M15c.WATERLOG.low$pred)
    ND.M15c.WATERLOG.low.mat<- matrix(ND.M15c.WATERLOG.low$pred.exp, length(sq.clay), length(sq.WetDays))

    ND.M15c.WATERLOG.high$pred<- as.vector(predict(jamM15c, newdata= ND.M15c.WATERLOG.high, se.fit= T)$fit) - as.vector(ND.M15c.baseWATERLOG.pred) # remove baseline (ignore baseline prediction error)
	ND.M15c.WATERLOG.high$pred.exp<- exp(ND.M15c.WATERLOG.high$pred)
    ND.M15c.WATERLOG.high.mat<- matrix(ND.M15c.WATERLOG.high$pred.exp, length(sq.clay), length(sq.WetDays))

    rng<- c(range(c(ND.M15c.WATERLOG.common$pred.exp, ND.M15c.WATERLOG.low$pred.exp, ND.M15c.WATERLOG.high$pred.exp)))
    theta.p<- -40
    phi.p<-  15
    
    pal.3<- c("#009E73", "#D55E00", "#0072B2")
	pdf(paste("M15c_Fig_InterWATERLOG_Nrate_theta", ".pdf", sep= ""), width= 20, height= 26)
	par(mfrow= c(4, 3))
	for(i in seq(-180, 180, by= 30)){
		theta.p<- i
		persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.low.mat, col= NULL, border= NULL, xlab= "% Clay", ylab= "Number of wet days", zlab= "Multiplier for (cumN2O+1)", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p, main= theta.p) -> res
		round(res, 3)
		sel<- dat$WetDays.exp < 200
		points(trans3d(jitter(dat$clay[sel], 2), jitter(dat$WetDays.exp[sel], 6), min(rng), pmat = res), col = rgb(0, 0, 0, 0.3), pch = 16, cex= 0.4)
		par(new= T)
		persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.common.mat, col= NULL, border= pal.3[1], xlab= "", ylab= "", zlab= "", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p)
		par(new= T)
		persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.high.mat, col= NULL, border= pal.3[2], xlab= "", ylab= "", zlab= "", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p)
		par(new= T)
		persp(x= sq.clay, y= sq.WetDays, z= ND.M15c.WATERLOG.low.mat, col= NULL, border= pal.3[3], xlab= "", ylab= "", zlab= "", ticktype= "detailed", zlim= rng, theta= theta.p, phi= phi.p) -> res
	}
	dev.off()



# Export objects needed for prediction in Shiny App:
save(df.M15c.complete, jamM15c, samM15c, file= "ShinyApp/N2O_predict_v_1_0.RData")


####### thresholds determination for multiplier factors (determined by Jon from graph of spline interaction Clay*WetDays)
(exp(-0.1 + 0.24) - 1) * 100 # x$Clay threshold = 15%
log(dat$WetDays.exp+10) - 4 # wet days
# 0 is  ~ 45 wet days
#logNrate.scaled<- log(dat$N.rate + 1) / 8
# logNrate.scaled == 1 is N.rate= 2979.958


# lookup values for M15c - extract linear predictor component functions from model
smooths.curves<- plot(jamM15c, pages= 0, scale= 0, scheme= 2, se= 1) # set multiplier for SE to 1

# lookup values for Degree Days
DD.tr<- smooths.curves[[5]]$x
DD.orig<- exp(smooths.curves[[5]]$x + 8) - 1500

par(mfrow= c(1, 2))
matplot(DD.orig, cbind(smooths.curves[[5]]$fit, 
                    smooths.curves[[5]]$fit - smooths.curves[[5]]$se,
                    smooths.curves[[5]]$fit + smooths.curves[[5]]$se),
                    type= "l", lty= c(1, 3, 3), col= 1, lwd= c(3, 1, 1), 
                    ylab= "Degree-days effect", xlab= "Degree-days (Celcius)",
                    xlim= c(0, 10000))
DD.lm<- lm(smooths.curves[[5]]$fit ~ DD.tr + I(DD.tr^2) + I(DD.tr^3))
summary(DD.lm)
lines(DD.orig, predict(DD.lm), col= 2)

matplot(DD.tr, cbind(smooths.curves[[5]]$fit, 
                    smooths.curves[[5]]$fit - smooths.curves[[5]]$se,
                    smooths.curves[[5]]$fit + smooths.curves[[5]]$se),
                    type= "l", lty= c(1, 3, 3), col= 1, lwd= c(3, 1, 1), 
                    ylab= "Degree-days effect", xlab= "Degree-days (Celcius)")

lines(DD.tr, predict(DD.lm), col= 2)



# lookup values for pH
pH.cut<- cut(smooths.curves[[1]]$x + 7, breaks= c(0, seq(5, 8.5, by= 0.5), 20))
pH.lookup<- by(data.frame(fit= smooths.curves[[1]]$fit, var= smooths.curves[[1]]$se^2), INDICES= pH.cut, FUN= function(x){weighted.mean(x$fit, x$var)})

# lookup values for logNrate.scaled, starting at a minimum of 20 kg/ha
logNrate.scaled.pred<- smooths.curves[[6]]$fit[smooths.curves[[6]]$x > 0] # in practice, equivalent to a minimum N.rate ~ 20 kg/ha (dictated by the range of values present in the data)
logNrate.scaled.seq<- smooths.curves[[6]]$x[smooths.curves[[6]]$x > 0]
logNrate.scaled.lm<- lm(logNrate.scaled.pred ~ logNrate.scaled.seq + I(logNrate.scaled.seq^2))
summary(logNrate.scaled.lm)
plot(jamM15c, select= 6, se= 1)
lines(logNrate.scaled.seq, predict(logNrate.scaled.lm), col= 2, lty= 3)


# lookup values for wetdays*clay*NO3 interaction
clay.cut<- cut((exp(smooths.curves[[2]]$x + 0.24) - 1) * 100, breaks= c(0, 21, 100))
wetD.cut<- cut(exp(smooths.curves[[2]]$y + 4) - 10, breaks= c(0, 45, 365, 10000)) # limit predictions to 1 year max
clay.wetD.cut<- apply(expand.grid(clay.cut, wetD.cut), 1, paste, collapse= "")
noNO3.lookup<- by(data.frame(fit= smooths.curves[[2]]$fit, var= smooths.curves[[2]]$se^2), INDICES= clay.wetD.cut, FUN= function(x){weighted.mean(x$fit, x$var, na.rm= T)})[c(1, 3, 4, 6)]
lowNO3.lookup<- by(data.frame(fit= smooths.curves[[3]]$fit, var= smooths.curves[[3]]$se^2), INDICES= clay.wetD.cut, FUN= function(x){weighted.mean(x$fit, x$var, na.rm= T)})[c(1, 3, 4, 6)]
highNO3.lookup<- by(data.frame(fit= smooths.curves[[4]]$fit, var= smooths.curves[[4]]$se^2), INDICES= clay.wetD.cut, FUN= function(x){weighted.mean(x$fit, x$var, na.rm= T)})[c(1, 3, 4, 6)]

# no fertilizer
round(as.vector(noNO3.lookup), 3)
# low NO3
paste(round(logNrate.scaled.lm$coefficients[1], 3), 
        c("-", "", "+")[sign(as.vector(lowNO3.lookup) + logNrate.scaled.lm$coefficients[2])+2],
        abs(round(as.vector(lowNO3.lookup) + logNrate.scaled.lm$coefficients[2], 3)),
        "*LNR",
        c("-", "", "+")[sign(logNrate.scaled.lm$coefficients[3])+2],
        abs(round(logNrate.scaled.lm$coefficients[3], 3)),
        "*LNR^2", sep= "")
# high NO3
paste(round(logNrate.scaled.lm$coefficients[1] +
            summ.M15c$statistics[4,1], 3), # add intercept for highNO3
        c("-", "", "+")[sign(as.vector(highNO3.lookup) + logNrate.scaled.lm$coefficients[2])+2],
        abs(round(as.vector(highNO3.lookup + logNrate.scaled.lm$coefficients[2]), 3)),
        "*LNR",
        c("-", "", "+")[sign(logNrate.scaled.lm$coefficients[3])+2],
        abs(round(logNrate.scaled.lm$coefficients[3], 3)),
        "*LNR^2", sep= "")



pdf("GAM.M15coutputBayes_splines_Fig3_SE_restricted.pdf", width= 10, height= 8)
	par(mfrow= c(2, 2), mar= c(4.3, 4.3, 4.3, 2.6))
	se.par= 1
	# pH
	plot(jamM15c, pages= 0, scale= 0, scheme= 2, se= se.par, select= 1, shade= T, ylab= "", xaxt= "n", xlab= "")
	box()
	axis(1, -4:2, (-4:2)+7)
	mtext("pH", 1, 3)
	axis(3)
	mtext("pH.cent", 3, 3)
	mtext("a)", 3, 3, adj= 0)
	mtext("Effect of pH", 2, 3)
	rug(smooths.curves[[1]]$raw)

	# DegDays
	plot(jamM15c, pages= 0, scale= 0, scheme= 2, se= se.par, select= 5, shade= T, ylab= "", xaxt= "n", xlab= "")
	box()
	axis(1, log(c(100, 1000, 2500, 5000, 10000, 20000)+1500) - 8, c(100, 1000, 2500, 5000, 10000, 20000))
	mtext("Degree-days", 1, 3)
	axis(3)
	mtext("logDegDays.exp.compact", 3, 3)
	mtext("b)", 3, 3, adj= 0)
	mtext("Effect of Degree-days", 2, 3)
	rug(smooths.curves[[5]]$raw)
	
	# Nrate
	Nratmp<- plot(jamM15c, pages= 0, scale= 0, scheme= 2, se= se.par, select= 6, shade= T, ylab= "", xaxt= "n", xlab= "")
	box()
	axis(1, log(c(0, 5, 25, 100, 500, 2500) + 1) / 8, c(0, 5, 25, 100, 500, 2500))
	mtext("N application rate (Kg, log scale)", 1, 3)
	axis(3)
	mtext("logNrate.scaled", 3, 3)
	mtext("c)", 3, 3, adj= 0)
	mtext("Effect of N rate", 2, 3)
	rug(smooths.curves[[6]]$raw)

	# Nrate on untransformed scale
	Nratmp[[6]]$x<- exp(Nratmp[[6]]$x * 8) - 1
	matplot(Nratmp[[6]]$x, cbind(Nratmp[[6]]$fit - Nratmp[[6]]$se, Nratmp[[6]]$fit + Nratmp[[6]]$se), ylab= "", xaxt= "n", xlab= "", type= "l", col= NA, xlim= c(0, 600), ylim= c(-0.5, 1.1))
	polygon(x= c(Nratmp[[6]]$x, rev(Nratmp[[6]]$x)), y= c(Nratmp[[6]]$fit - Nratmp[[6]]$se, rev(Nratmp[[6]]$fit + Nratmp[[6]]$se)), col= grey(0.8), border= NA)
	lines(Nratmp[[6]]$x, Nratmp[[6]]$fit)
	box()
	mtext("N application rate (Kg, linear scale)", 1, 3)
	axis(1)
	mtext("d)", 3, 3, adj= 0)
	mtext("Effect of N rate", 2, 3)
	rug(exp(smooths.curves[[6]]$raw * 8) - 1)
dev.off()

pdf("GAM.M15coutputBayes_splines_Fig3_SE_unrestricted.pdf", width= 10, height= 8)
	par(mfrow= c(2, 2), mar= c(4.3, 4.3, 4.3, 2.6))
	se.par= 1
	# pH
	plot(jamM15c, pages= 0, scale= 0, scheme= 2, se= se.par, select= 1, shade= T, ylab= "", xaxt= "n", xlab= "")
	box()
	axis(1, -4:2, (-4:2)+7)
	mtext("pH", 1, 3)
	axis(3)
	mtext("pH.cent", 3, 3)
	mtext("a)", 3, 3, adj= 0)
	mtext("Effect of pH", 2, 3)
	rug(smooths.curves[[1]]$raw)

	# DegDays
	plot(jamM15c, pages= 0, scale= 0, scheme= 2, se= se.par, select= 5, shade= T, ylab= "", xaxt= "n", xlab= "")
	box()
	axis(1, log(c(100, 1000, 2500, 5000, 10000, 20000)+1500) - 8, c(100, 1000, 2500, 5000, 10000, 20000))
	mtext("Degree-days", 1, 3)
	axis(3)
	mtext("logDegDays.exp.compact", 3, 3)
	mtext("b)", 3, 3, adj= 0)
	mtext("Effect of Degree-days", 2, 3)
	rug(smooths.curves[[5]]$raw)
	
	# Nrate
	Nratmp<- plot(jamM15c, pages= 0, scale= 0, scheme= 2, se= se.par, select= 6, shade= T, ylab= "", xaxt= "n", xlab= "")
	box()
	axis(1, log(c(0, 5, 25, 100, 500, 2500) + 1) / 8, c(0, 5, 25, 100, 500, 2500))
	mtext("N application rate (Kg, log scale)", 1, 3)
	axis(3)
	mtext("logNrate.scaled", 3, 3)
	mtext("c)", 3, 3, adj= 0)
	mtext("Effect of N rate", 2, 3)
	rug(smooths.curves[[6]]$raw)

	# Nrate on untransformed scale
	Nratmp[[6]]$x<- exp(Nratmp[[6]]$x * 8) - 1
	matplot(Nratmp[[6]]$x, cbind(Nratmp[[6]]$fit - Nratmp[[6]]$se, Nratmp[[6]]$fit + Nratmp[[6]]$se), ylab= "", xaxt= "n", xlab= "", type= "l", col= NA)
	polygon(x= c(Nratmp[[6]]$x, rev(Nratmp[[6]]$x)), y= c(Nratmp[[6]]$fit - Nratmp[[6]]$se, rev(Nratmp[[6]]$fit + Nratmp[[6]]$se)), col= grey(0.8), border= NA)
	lines(Nratmp[[6]]$x, Nratmp[[6]]$fit)
	box()
	mtext("N application rate (Kg, linear scale)", 1, 3)
	axis(1)
	mtext("d)", 3, 3, adj= 0)
	mtext("Effect of N rate", 2, 3)
	rug(exp(smooths.curves[[6]]$raw * 8) - 1)
dev.off()


plot.mcmc.fx<- function(jags.mcmc.fx.mat){
	fx.mean <- apply(jags.mcmc.fx.mat, 2, mean)
	fx.q.025<- apply(jags.mcmc.fx.mat, 2, quantile, prob= 0.025)
	fx.q.975<- apply(jags.mcmc.fx.mat, 2, quantile, prob= 0.975)
	fx.q.25<- apply(jags.mcmc.fx.mat, 2, quantile, prob= 0.25)
	fx.q.75<- apply(jags.mcmc.fx.mat, 2, quantile, prob= 0.75)
	plot(y= 1:ncol(jags.mcmc.fx.mat), x= 1:ncol(jags.mcmc.fx.mat), xlim= c(min(fx.q.025), max(fx.q.975)), xlab= "Parameter posterior estimate", ylab= "", yaxt= "n", col= NULL)
	abline(h= 1:ncol(jags.mcmc.fx.mat), col= grey(0.75))
	abline(v= 0, col= grey(0.75))
	segments(y0= 1:ncol(jags.mcmc.fx.mat), y1= 1:ncol(jags.mcmc.fx.mat), 
		x0= fx.q.025, x1= fx.q.975, col= "blue1", lwd= 1.7)
	segments(y0= 1:ncol(jags.mcmc.fx.mat), y1= 1:ncol(jags.mcmc.fx.mat), 
		x0= fx.q.25, x1= fx.q.75, col= "blue4", lwd= 4)
	points(y= 1:ncol(jags.mcmc.fx.mat), x= fx.mean, bg= "blue1", col= "blue4", pch= 21)
	axis(2, at= 1:ncol(jags.mcmc.fx.mat), labels= colnames(jags.mcmc.fx.mat), las= 1)
}

pdf("GAM.M15coutputBayes_densityplots_fixed_Fig2.pdf", width= 6, height= 3.5)
par(mfrow= c(1, 1), mar= c(4.1, 12, 3.1, 1.1))
	fx.mcmc<- jagsM15c.mcmc[[1]][, fixed.names.M15c]
	fx.mcmc[, 2]<- fx.mcmc[, 2] #/ 10 # rescale
	# colnames(fx.mcmc)[2]<- paste(colnames(fx.mcmc)[2], "*", sep= "")
	plot.mcmc.fx(fx.mcmc)
dev.off()

# identify random effects
nranef.M15c<- rle(colnames(pregam.M15c$jags.data$X)=="")$lengths[3] # run length encoding -> 3rd run is random effects
ranef.names.M15c<- pregam.M15c$pregam$term.names[(ncol(pregam.M15c$jags.data$X)-nranef.M15c+1):ncol(pregam.M15c$jags.data$X)]

jamM15c$coefficients[grep("32.58 119.702007 8  8120082007-08-152007-11-04", colnames(jamM15c$X))] # 0.4315864

ranef.studyID.size<- 0 - jamM15c$coefficients[grep("studyID.f", colnames(jamM15c$X))]
ranef.ref<- which(abs(ranef.studyID.size) == min(abs(ranef.studyID.size)))
colnames(jamM15c$X)[grep("studyID.f", colnames(jamM15c$X))][ranef.ref]
# "studyID.f 53.42  -7.52200310 41220052003-10-152004-11-30"

jamM15c$coefficients[grep(" 53.42  -7.52200310 41220052003-10-152004-11-30", colnames(jamM15c$X))] # -0.001805962


# plot observed cumN2O ~ Nrate
plot(dat.M15c$logN2O.cum.pos ~ dat.M15c$logNrate.scaled)














