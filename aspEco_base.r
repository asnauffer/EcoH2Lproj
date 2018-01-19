# runs SnowMelt from EcoHydRology using T & P data

# load libraries and set working directory
library("EcoHydRology")
library("R.matlab")
setwd("~/PhD/EcoH2Lproj/")

# source codes
source("SnowMelt.R")
source("SnowMelt2L.R")

print("loading and initializing...")
asplatlon <- data.frame(readMat("../matfiles/ASPalllatlon.mat"))
aspname <- asplatlon$stnname
asplat <- asplatlon$stnlatlon.1  # lat used in SnowMelt
stnphysionum <- unlist(readMat("../matfiles/stnphysionum.mat"))
stnlatlonelevmat <- readMat("../matfiles/stnlatlonelev.mat")
stnID <- unlist(stnlatlonelevmat$stnID)

aspstnnums <- c(11,51,12,43,56,61,27,40) # 8 stations selected for this study

# init vectors
out_asp <- vector("list",1)
yrct <- 0

# loop through stations
for(istn in aspstnnums){
  aspregnum <- stnphysionum[aspname[istn]==stnID]
  print(paste("stn:",istn,aspname[istn],"region",aspregnum))
  # read and interpret data
  fn <- paste("../data/ASP/",aspname[istn],".csv",sep="")
  aspallna <- read.csv(fn,skip=8) # read entire file including NA dates
  # read and interpret date
  if(is.element(istn,c(2,3,38,57))){
    fmt <- "%m/%d/%Y"
  }else{
    fmt <- "%Y-%m-%d"
  }
  aspdatetest <- as.Date(aspallna$Date,format=fmt) # read dates according to format
  aspall <- aspallna[is.finite(aspdatetest),] # filter NA dates
  aspallswe <- aspall$Snow.Water.Equivalent
  aspallprec <- aspall$Precipitation
  aspalltmin <- aspall$Temp..Min.
  aspalltmax <- aspall$Temp..Max.
  aspalldate <- as.Date(aspall$Date,format=fmt)
  aspallY <- as.numeric(format(aspalldate,"%Y")) # extract year from dates
  aspYunique <- unique(aspallY) # get unique years
  
  # check for missing and invalid data
  lognatx <- is.na(aspalltmax) # Tmax is NA
  lognatn <- is.na(aspalltmin) # Tmin is NA
  lognaswe <- is.na(aspall$Snow.Water.Equivalent) # SWE is NA 
  logtxlttn <- aspalltmax < aspalltmin # Tmin>Tmax (causes error in SnowMelt)
  
  # check for increase in SWE without corresponding daily precipitation
  sweinctol <- 0.0 # SWE amount without precip above which data is thrown out
  swediff <- diff(aspallswe) # array of change in SWE
  logsweinc <- c(0,swediff) > sweinctol # logical - which days experienced SWE increase (0 added at beginning so array is the same length)
  logsweinc[is.na(logsweinc)] <- F # if NA, this means SWE didn't increase and can be set to False
  lognap <- is.na(aspallprec) # logical - precip is NA
  logzerop <- aspallprec==0 # logical - precip is 0
  lognop <- lognap | logzerop # logical - no precip when precip is 0 or NA
  logswenop <- logsweinc & lognop # logical - SWE increased but no corresponding precip 

  # create discard vector 
  logbadtp <- lognatx | lognatn | logtxlttn | lognap # bad T or P values
  logdiscard <- logbadtp | logswenop # logical - discard bad T or P or when SWE increases with no P
  
  for(iyr in aspYunique){
    dateseasoni <- paste(iyr-1,"-10-01",sep="")
    dateseasonf <- paste(iyr,"-07-01",sep="")
    logyr <- aspalldate>=dateseasoni & aspalldate<=dateseasonf # select dates for current snow year
    logvalidmodel <- logyr & !logdiscard # select valid dates for SnowMelt
    logvalidcompare <- logvalidmodel & !lognaswe # select valid dates for comparison from full date vector
    logvalidcompareyr <- logvalidcompare[logvalidmodel] # select valid dates for comparison from current year 
    swediffnop <- swediff[logyr & logswenop] # vector of SWE differences when 0 precip recorded 
    sumswenop <- sum(swediffnop[swediffnop>0],na.rm=T) 
    swediffdisc <- swediff[logyr & logdiscard] # vector of SWE differences when any inputs are invalid and must be discarded (what is really needed in place of swediffnop) 
    sumswebadtp <- sum(swediffdisc[swediffdisc>0],na.rm=T) 
    # check that some valid dates still remain, otherwise skip model run
    validdates <- aspalldate[logvalidmodel] 
    datediff <- diff(validdates)
    if(sum(logvalidmodel)==0){
      cat("no valid dates for",iyr,"\n")
      next
    }
    # if total SWE differences while inputs are valid > some threshold, skip this year
    if(sumswebadtp>100){
      cat("sumswebadtp =",sumswebadtp,"for",iyr,"\n")
      next
    }
    # run and time SnowMelt
    asp <- aspall[logvalidmodel,]
    aspswe <- asp$Snow.Water.Equivalent
    aspdate <- as.Date(asp$Date,format=fmt)
    
    # stock SnowMelt
    exectime <- tryCatch(round(system.time(
      sm_asp <- SnowMelt(Date=aspdate, precip_mm=asp$Precipitation,
                           Tmax_C=asp$Temp..Max., Tmin_C=asp$Temp..Min., lat_deg=asplat[istn])
    ),3)[3],error = function(e) e)
    if(inherits(exectime, "error")) {
      print(paste("SnowMelt ERROR: stn",istn,aspname[istn],iyr,exectime))
    } else {
      smswe <- sm_asp$SnowWaterEq_mm
      mae <- mean(abs(aspswe[logvalidcompareyr]-smswe[logvalidcompareyr]))
      smsm <- sm_asp$SnowMelt_mm
    }

    # 2L SnowMelt 
    exectime2L <- tryCatch(round(system.time(
      sm_asp2L <- SnowMelt2L(Date=aspdate, precip_mm=asp$Precipitation,
                           Tmax_C=asp$Temp..Max., Tmin_C=asp$Temp..Min., lat_deg=asplat[istn])
    ),3)[3],error = function(e) e)
    if(inherits(exectime2L, "error")) {
      print(paste("SnowMelt2L ERROR: stn",istn,aspname[istn],iyr,exectime2L))
      next
    }
    smswe2L <- sm_asp2L$SnowWaterEq_mm
    smsweu2L <- sm_asp2L$SnowWaterEqUpper_mm
    smswel2L <- sm_asp2L$SnowWaterEqLower_mm
    mae2L <- mean(abs(aspswe[logvalidcompareyr]-smswe2L[logvalidcompareyr]))
    smsm2L <- sm_asp2L$SnowMelt_mm
    smsmu2L <- sm_asp2L$SnowMeltUpper_mm
    smsml2L <- sm_asp2L$SnowMeltLower_mm
    
    # 2L SnowMelt with ground conduction G=0
    exectime2LG0 <- tryCatch(round(system.time(
      sm_asp2LG0 <- SnowMelt2L(Date=aspdate, precip_mm=asp$Precipitation,
                             Tmax_C=asp$Temp..Max., Tmin_C=asp$Temp..Min., lat_deg=asplat[istn], G=0)
    ),3)[3],error = function(e) e)
    if(inherits(exectime2LG0, "error")) {
      print(paste("SnowMelt2LG0 ERROR: stn",istn,aspname[istn],exectime2LG0))
      next
    }
    smswe2LG0 <- sm_asp2LG0$SnowWaterEq_mm
      mae2LG0 <- mean(abs(aspswe[logvalidcompareyr]-smswe2LG0[logvalidcompareyr]))

    # store stats
    yri <- data.frame(aspname[istn],
                      dateyr=iyr,
                      numdates=sum(logyr),
                      numbadtp=sum(logyr & logbadtp),
                      numswebadtp=sum(logyr & logdiscard),
                      sumswebadtp=sumswebadtp,
                      numvalidmodel=sum(logvalidmodel),
                      numvalidcompare=sum(logvalidcompare),
                      maxdategap=max(datediff,na.rm=T),
                      exectime=exectime,
                      mae=mae,
                      exectime2L=exectime2L,
                      mae2L=mae2L
    )
    yrfi <- data.frame(yri,
                      numnatn=sum(logyr & lognatn),
                      numnatx=sum(logyr & lognatx),
                      numtxlttn=sum(logyr & logtxlttn,na.rm=T),
                      numnap=sum(logyr & lognap),
                      numswenop=sum(logyr & logswenop),
                      sumswenop=sumswenop,
                      numnaswe=sum(logyr & lognaswe)
    )
    yrct <- yrct+1
    if(yrct==1){
      yearrun <- yri
      yearrunfull <- yrfi
    }
    else {
      yearrun[yrct,] <- yri
      yearrunfull[yrct,] <- yrfi
    }
    if(is.na(sumswebadtp)){
      print("sumswebadtp = NA")
    }
    
    
    #
    #
    # plot asp and snowmelt curves
    #
    #
    
    # set up and plot SWEs
    linew <- 5
    fnout <- paste("plots/SM2L_",iyr,"_",aspregnum,"_",istn,"_",aspname[istn],"_sumswebad",sumswebadtp,".jpg",sep="")
    jpeg(fnout,width=480*4,height=480*2,pointsize=24,quality=100)
    par(mfrow=c(2,1)) # setup layout
    par(xpd=TRUE) # turns off clipping so legend is visible
    marvals <- c(3.1, 4.1, 2.1, 2.1) # default c(5.1, 4.1, 4.1, 2.1)
    par(mar=marvals)
    plot(aspdate,aspswe,col="black","l",xlab="",ylab="SWE (mm)",lwd=linew,xlim=c(min(aspdate),max(aspdate)+70),ylim=c(0,max(aspswe,smswe,smswe2L,na.rm=T)), bty="L")
    lines(aspdate,smswe,col="red",lwd=linew)
    lines(aspdate,smswe2L,col="lightblue",lwd=linew+3)
    legtxt <- c("ASP measured","SnowMelt model","SnowMelt2L model")
    legcol <- c("black","red","lightblue")
    nextlegcol <- c("blue","darkblue","midnightblue")

    # plot additional SWEs with other albedo thresholds
    albedos <- c(0.3,0.05)
    for(ia in 1:2) {
      exectime2L <- tryCatch(round(system.time(
      sm_asp2L <- SnowMelt2L(Date=aspdate, precip_mm=asp$Precipitation,
                             Tmax_C=asp$Temp..Max.,Tmin_C=asp$Temp..Min.,lat_deg=asplat[istn],albedoThreshold=albedos[ia])
      ),3)[3],error = function(e) e)
      
      if(inherits(exectime2L, "error")) {
      print(paste("SnowMelt2L ERROR: stn",istn,aspname[istn],iyr,"ia=",ia,exectime2L))
      } else {
      smswe2L <- sm_asp2L$SnowWaterEq_mm
      legtxt <- c(legtxt,paste("SnowMelt2L (albedo SWE<",albedos[ia],")",sep=""))
      legcol <- c(legcol,nextlegcol[ia])
      lines(aspdate,smswe2L,col=nextlegcol[ia],lwd=linew-ia)
      }
    }
    
    # add title, letter, and legend
    title(paste("Station",aspname[istn],iyr))
    mtext("(a)",side=3,line=1,at=as.Date(paste(iyr-1,"-09-01",sep="")),cex=1)
    legend(max(aspdate),max(aspswe,smswe,smswe2L,na.rm=T),legtxt,col=legcol,lwd=linew,bty="n")
      
    # temperature plot
    tempdata <- data.frame((asp$Temp..Max.+asp$Temp..Min.)/2,sm_asp$SnowTemp,sm_asp2L$SnowTempUpper,sm_asp2L$SnowTempLower)
    plot(aspdate,(asp$Temp..Max.+asp$Temp..Min.)/2,col="green","l",xlab="",ylab="T (deg-C)",lwd=linew,xlim=c(min(aspdate),max(aspdate)+70),ylim=c(min(tempdata),max(tempdata)), bty="L")
    lines(aspdate,asp$Temp..Max.,col="lightgreen",lwd=1.5)
    lines(aspdate,asp$Temp..Min.,col="darkgreen",lwd=1.5)
    lines(aspdate,sm_asp$SnowTemp,col="black",lwd=6)
    lines(aspdate,sm_asp2L$SnowTempUpper,col="red",lwd=5)
    lines(aspdate,sm_asp2L$SnowTempLower,col="blue",lwd=2)
    mtext("(b)",side=3,line=1,at=as.Date(paste(iyr-1,"-09-01",sep="")),cex=1)
    legend(max(aspdate),max(asp$Temp..Max.,na.rm=T),c("Tav Air","Tmax Air","Tmin Air","EcoH 1L Snow T","EcoH 2L UpperL Snow T","EcoH 2L LowerL Snow T"),
           col=c("green","lightgreen","darkgreen","black","red","blue"),lwd=linew,bty="n")
    
    # finalize plot
    dev.off()
    cat(fnout,"plotted\n")
    
  } # next year
} # next station