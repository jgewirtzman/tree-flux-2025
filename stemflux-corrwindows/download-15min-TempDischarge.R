inUrl4  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-hfr/70/36/2983b2adba6675e805d144a05087924d" 
infile4 <- tempfile()
try(download.file(inUrl4,infile4,method="curl"))
if (is.na(file.size(infile4))) download.file(inUrl4,infile4,method="auto")


dt4 <-read.csv(infile4,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "datetime",     
                 "jd",     
                 "nb.stg",     
                 "f.nb.stg",     
                 "nl.stg",     
                 "f.nl.stg",     
                 "al.stg",     
                 "f.al.stg",     
                 "au.stg",     
                 "f.au.stg",     
                 "bgs.stg",     
                 "f.bgs.stg",     
                 "bvs.stg",     
                 "f.bvs.stg",     
                 "nb.dis",     
                 "f.nb.dis",     
                 "nl.dis",     
                 "f.nl.dis",     
                 "nt.dis",     
                 "f.nt.dis",     
                 "al.dis",     
                 "f.al.dis",     
                 "au.dis",     
                 "f.au.dis",     
                 "nb.wt",     
                 "f.nb.wt",     
                 "nl.wt",     
                 "f.nl.wt",     
                 "al.wt",     
                 "f.al.wt",     
                 "au.wt",     
                 "f.au.wt",     
                 "bgs.wt",     
                 "f.bgs.wt",     
                 "bvs.wt",     
                 "f.bvs.wt"    ), check.names=TRUE)

unlink(infile4)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

# attempting to convert dt4$datetime dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%Y-%m-%dT%H:%M" 
tmp4datetime<-as.POSIXct(dt4$datetime,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(nrow(dt4[dt4$datetime != "",]) == length(tmp4datetime[!is.na(tmp4datetime)])){dt4$datetime <- tmp4datetime } else {print("Date conversion failed for dt4$datetime. Please inspect the data and do the date conversion yourself.")}                                                                    

if (class(dt4$jd)=="factor") dt4$jd <-as.numeric(levels(dt4$jd))[as.integer(dt4$jd) ]               
if (class(dt4$jd)=="character") dt4$jd <-as.numeric(dt4$jd)
if (class(dt4$nb.stg)=="factor") dt4$nb.stg <-as.numeric(levels(dt4$nb.stg))[as.integer(dt4$nb.stg) ]               
if (class(dt4$nb.stg)=="character") dt4$nb.stg <-as.numeric(dt4$nb.stg)
if (class(dt4$f.nb.stg)!="factor") dt4$f.nb.stg<- as.factor(dt4$f.nb.stg)
if (class(dt4$nl.stg)=="factor") dt4$nl.stg <-as.numeric(levels(dt4$nl.stg))[as.integer(dt4$nl.stg) ]               
if (class(dt4$nl.stg)=="character") dt4$nl.stg <-as.numeric(dt4$nl.stg)
if (class(dt4$f.nl.stg)!="factor") dt4$f.nl.stg<- as.factor(dt4$f.nl.stg)
if (class(dt4$al.stg)=="factor") dt4$al.stg <-as.numeric(levels(dt4$al.stg))[as.integer(dt4$al.stg) ]               
if (class(dt4$al.stg)=="character") dt4$al.stg <-as.numeric(dt4$al.stg)
if (class(dt4$f.al.stg)!="factor") dt4$f.al.stg<- as.factor(dt4$f.al.stg)
if (class(dt4$au.stg)=="factor") dt4$au.stg <-as.numeric(levels(dt4$au.stg))[as.integer(dt4$au.stg) ]               
if (class(dt4$au.stg)=="character") dt4$au.stg <-as.numeric(dt4$au.stg)
if (class(dt4$f.au.stg)!="factor") dt4$f.au.stg<- as.factor(dt4$f.au.stg)
if (class(dt4$bgs.stg)=="factor") dt4$bgs.stg <-as.numeric(levels(dt4$bgs.stg))[as.integer(dt4$bgs.stg) ]               
if (class(dt4$bgs.stg)=="character") dt4$bgs.stg <-as.numeric(dt4$bgs.stg)
if (class(dt4$f.bgs.stg)!="factor") dt4$f.bgs.stg<- as.factor(dt4$f.bgs.stg)
if (class(dt4$bvs.stg)=="factor") dt4$bvs.stg <-as.numeric(levels(dt4$bvs.stg))[as.integer(dt4$bvs.stg) ]               
if (class(dt4$bvs.stg)=="character") dt4$bvs.stg <-as.numeric(dt4$bvs.stg)
if (class(dt4$f.bvs.stg)!="factor") dt4$f.bvs.stg<- as.factor(dt4$f.bvs.stg)
if (class(dt4$nb.dis)=="factor") dt4$nb.dis <-as.numeric(levels(dt4$nb.dis))[as.integer(dt4$nb.dis) ]               
if (class(dt4$nb.dis)=="character") dt4$nb.dis <-as.numeric(dt4$nb.dis)
if (class(dt4$f.nb.dis)!="factor") dt4$f.nb.dis<- as.factor(dt4$f.nb.dis)
if (class(dt4$nl.dis)=="factor") dt4$nl.dis <-as.numeric(levels(dt4$nl.dis))[as.integer(dt4$nl.dis) ]               
if (class(dt4$nl.dis)=="character") dt4$nl.dis <-as.numeric(dt4$nl.dis)
if (class(dt4$f.nl.dis)!="factor") dt4$f.nl.dis<- as.factor(dt4$f.nl.dis)
if (class(dt4$nt.dis)=="factor") dt4$nt.dis <-as.numeric(levels(dt4$nt.dis))[as.integer(dt4$nt.dis) ]               
if (class(dt4$nt.dis)=="character") dt4$nt.dis <-as.numeric(dt4$nt.dis)
if (class(dt4$f.nt.dis)!="factor") dt4$f.nt.dis<- as.factor(dt4$f.nt.dis)
if (class(dt4$al.dis)=="factor") dt4$al.dis <-as.numeric(levels(dt4$al.dis))[as.integer(dt4$al.dis) ]               
if (class(dt4$al.dis)=="character") dt4$al.dis <-as.numeric(dt4$al.dis)
if (class(dt4$f.al.dis)!="factor") dt4$f.al.dis<- as.factor(dt4$f.al.dis)
if (class(dt4$au.dis)=="factor") dt4$au.dis <-as.numeric(levels(dt4$au.dis))[as.integer(dt4$au.dis) ]               
if (class(dt4$au.dis)=="character") dt4$au.dis <-as.numeric(dt4$au.dis)
if (class(dt4$f.au.dis)!="factor") dt4$f.au.dis<- as.factor(dt4$f.au.dis)
if (class(dt4$nb.wt)=="factor") dt4$nb.wt <-as.numeric(levels(dt4$nb.wt))[as.integer(dt4$nb.wt) ]               
if (class(dt4$nb.wt)=="character") dt4$nb.wt <-as.numeric(dt4$nb.wt)
if (class(dt4$f.nb.wt)!="factor") dt4$f.nb.wt<- as.factor(dt4$f.nb.wt)
if (class(dt4$nl.wt)=="factor") dt4$nl.wt <-as.numeric(levels(dt4$nl.wt))[as.integer(dt4$nl.wt) ]               
if (class(dt4$nl.wt)=="character") dt4$nl.wt <-as.numeric(dt4$nl.wt)
if (class(dt4$f.nl.wt)!="factor") dt4$f.nl.wt<- as.factor(dt4$f.nl.wt)
if (class(dt4$al.wt)=="factor") dt4$al.wt <-as.numeric(levels(dt4$al.wt))[as.integer(dt4$al.wt) ]               
if (class(dt4$al.wt)=="character") dt4$al.wt <-as.numeric(dt4$al.wt)
if (class(dt4$f.al.wt)!="factor") dt4$f.al.wt<- as.factor(dt4$f.al.wt)
if (class(dt4$au.wt)=="factor") dt4$au.wt <-as.numeric(levels(dt4$au.wt))[as.integer(dt4$au.wt) ]               
if (class(dt4$au.wt)=="character") dt4$au.wt <-as.numeric(dt4$au.wt)
if (class(dt4$f.au.wt)!="factor") dt4$f.au.wt<- as.factor(dt4$f.au.wt)
if (class(dt4$bgs.wt)=="factor") dt4$bgs.wt <-as.numeric(levels(dt4$bgs.wt))[as.integer(dt4$bgs.wt) ]               
if (class(dt4$bgs.wt)=="character") dt4$bgs.wt <-as.numeric(dt4$bgs.wt)
if (class(dt4$f.bgs.wt)!="factor") dt4$f.bgs.wt<- as.factor(dt4$f.bgs.wt)
if (class(dt4$bvs.wt)=="factor") dt4$bvs.wt <-as.numeric(levels(dt4$bvs.wt))[as.integer(dt4$bvs.wt) ]               
if (class(dt4$bvs.wt)=="character") dt4$bvs.wt <-as.numeric(dt4$bvs.wt)
if (class(dt4$f.bvs.wt)!="factor") dt4$f.bvs.wt<- as.factor(dt4$f.bvs.wt)

