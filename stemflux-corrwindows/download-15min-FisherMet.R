# Aggregate air temp & pressure
# Fisher met station: 15-minute data ongoing since 2005-01-01 
inUrl10  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-hfr/1/34/0b439e8fea983c9e20bb2bfaf91931e6" 
infile10 <- tempfile()
try(download.file(inUrl10,infile10,method="curl"))
if (is.na(file.size(infile10))) download.file(inUrl10,infile10,method="auto")

dt10 <-read.csv(infile10,header=F 
                ,skip=1
                ,sep=","  
                , col.names=c(
                  "datetime",     
                  "jd",     
                  "airt",     
                  "f.airt",     
                  "rh",     
                  "f.rh",     
                  "dewp",     
                  "f.dewp",     
                  "prec",     
                  "f.prec",     
                  "slrr",     
                  "f.slrr",     
                  "parr",     
                  "f.parr",     
                  "netr",     
                  "f.netr",     
                  "bar",     
                  "f.bar",     
                  "wspd",     
                  "f.wspd",     
                  "wres",     
                  "f.wres",     
                  "wdir",     
                  "f.wdir",     
                  "wdev",     
                  "f.wdev",     
                  "gspd",     
                  "f.gspd",     
                  "s10t",     
                  "f.s10t"    ), check.names=TRUE)

unlink(infile10)