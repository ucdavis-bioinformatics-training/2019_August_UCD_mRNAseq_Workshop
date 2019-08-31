## first need to install the R packages JSON

#lets install to current working directory and then load
suppressWarnings(dir.create("r_lib"))
new_rlib = file.path(getwd(),"r_lib")

if (!any(rownames(installed.packages()) == "rjson")) install.packages("rjson", repos='http://cran.us.r-project.org', lib=new_rlib)

require(rjson, lib.loc=new_rlib)

samples <- readLines("samples.txt")
hts_dir <- "01-HTS_Preproc"

hts_json <- lapply(samples, function(s) {
  fromJSON(paste(readLines(file.path(hts_dir,s,paste0(s,"_htsStats.log"))),collapse=""))
})

if (length(hts_json) != length(samples)) stop("not all samples have json log files")
if (!all(apply(sapply(hts_json,names),1, function(x) length(unique(sub("_[0-9]+$","",x)))) == 1)) stop("not all json log files have the same htsteam sub apps")

# sub("_[0-9]+$","",x)
apps <- apply(sapply(hts_json,names),1, function(x) unique(sub("_[0-9]+$","",x)))
napps <- length(apps)
apps <- make.names(apps, unique = T)

check_apps <- c("hts_Stats", "hts_SeqScreener", "hts_SeqScreener.1", "hts_AdapterTrimmer", "hts_CutTrim", "hts_QWindowTrim", "hts_NTrimmer", "hts_CutTrim.1", "hts_Stats.1")
if (!(length(check_apps) == length(apps) && all(check_apps == apps))) stop("expected sequence of apps not found")

## Table of total reads output after each stage
#totalReadsOut <- as.data.frame(sapply(hts_json,function(x) sapply(x, "[[","totalFragmentsOutput")))
#rownames(totalReadsOut) <- apps
#colnames(totalReadsOut) <- samples

### Super Deduper Saturation Plot
# sd <- which(apps == "hts_SuperDeduper")
# if(length(sd) > 1) stop("More than one hts super deduper found, which makes little sense")
#
# if (length(sd) == 1){
#   ds <- lapply(hts_json,function(x) x[[sd]]$duplicate_saturation)
#   plot(cumsum(diff(sapply(ds[[2]],"[[",1L))),diff(sapply(ds[[2]],"[[",2L))/diff(sapply(ds[[2]],"[[",1L)),type='l')
# }

## Pre-Stats
statsTotalReadsIn <- sapply(hts_json,function(x) sapply(x[1], "[[","totalFragmentsInput"))
statsTotalReadsCG <- sapply(hts_json,function(x) sum(unlist(sapply(x[1],"[[","Base_composition")[c("C","G"),])))
statsTotalReadsN <- sapply(hts_json,function(x) unlist(sapply(x[1],"[[","Base_composition")[c("N"),]))
statsTotalReadsSE_BpLen <- sapply(hts_json,function(x) unlist(sapply(x[1],"[[","Single_end")[c("SE_bpLen"),]))
statsTotalReadsSE_BQ30 <- sapply(hts_json,function(x) unlist(sapply(x[1],"[[","Single_end")[c("SE_bQ30"),]))

## PhiX screen
seqScrTotalReadsIn <- sapply(hts_json,function(x) sapply(x[2], "[[","totalFragmentsInput"))
seqScrTotalReadsHits <- sapply(hts_json,function(x) unlist(sapply(x[2],"[[","Single_end")[c("SE_hits"),]))

## Seq Screen rRNA
seqScr2TotalReadsIn <- sapply(hts_json,function(x) sapply(x[3], "[[","totalFragmentsInput"))
seqScr2TotalReadsHits <- sapply(hts_json,function(x) unlist(sapply(x[3],"[[","Single_end")[c("SE_hits"),]))

## Adapter Trimmer
AdaptTotalReadsIn <- sapply(hts_json,function(x) sapply(x[4], "[[","totalFragmentsInput"))
AdaptTotalAdaptTrim <- sapply(hts_json,function(x) unlist(sapply(x[4],"[[","Single_end")[c("SE_adapterTrim"),]))
AdaptTotalBpTrim <- sapply(hts_json,function(x) unlist(sapply(x[4],"[[","Single_end")[c("SE_adapterBpTrim"),]))

## CutTrim, hard 20bp trim
CutTrimHardTotalReadsIn <- sapply(hts_json,function(x) sapply(x[5], "[[","totalFragmentsInput"))
CutTrimHardTotalSELT <- sapply(hts_json,function(x) unlist(sapply(x[5],"[[","Single_end")[c("SE_leftTrim"),]))
CutTrimHardTotalDiscard <- sapply(hts_json,function(x) unlist(sapply(x[5],"[[","Single_end")[c("SE_discarded"),]))

## Q windowtrim
QwinTotalReadsIn <- sapply(hts_json,function(x) sapply(x[6], "[[","totalFragmentsInput"))
QwinTotalSELT <- sapply(hts_json,function(x) unlist(sapply(x[6],"[[","Single_end")[c("SE_leftTrim"),]))
QwinTotalSERT <- sapply(hts_json,function(x) unlist(sapply(x[6],"[[","Single_end")[c("SE_rightTrim"),]))
QwinTotalDiscard <- sapply(hts_json,function(x) unlist(sapply(x[6],"[[","Single_end")[c("SE_discarded"),]))

## N trim
NwinTotalReadsIn <- sapply(hts_json,function(x) sapply(x[7], "[[","totalFragmentsInput"))
NwinTotalSELT <- sapply(hts_json,function(x) unlist(sapply(x[7],"[[","Single_end")[c("SE_leftTrim"),]))
NwinTotalSERT <- sapply(hts_json,function(x) unlist(sapply(x[7],"[[","Single_end")[c("SE_rightTrim"),]))
NwinTotalDiscard <- sapply(hts_json,function(x) unlist(sapply(x[7],"[[","Single_end")[c("SE_discarded"),]))

## CutTrim
CutTrimTotalReadsIn <- sapply(hts_json,function(x) sapply(x[8], "[[","totalFragmentsInput"))
CutTrimTotalDiscard <- sapply(hts_json,function(x) unlist(sapply(x[8],"[[","Single_end")[c("SE_discarded"),]))

## Post-Stats
stats2TotalReadsIn <- sapply(hts_json,function(x) sapply(x[9], "[[","totalFragmentsInput"))
stats2TotalReadsCG <- sapply(hts_json,function(x) sum(unlist(sapply(x[9],"[[","Base_composition")[c("C","G"),])))
stats2TotalReadsN <- sapply(hts_json,function(x) unlist(sapply(x[9],"[[","Base_composition")[c("N"),]))
stats2TotalReadsSE_BpLen <- sapply(hts_json,function(x) unlist(sapply(x[9],"[[","Single_end")[c("SE_bpLen"),]))
stats2TotalReadsSE_BQ30 <- sapply(hts_json,function(x) unlist(sapply(x[9],"[[","Single_end")[c("SE_bQ30"),]))


LongTable <- data.frame(
    Raw_Reads=statsTotalReadsIn,
    Raw_Bp=(statsTotalReadsSE_BpLen),
    Raw_SE_PercentQ30=statsTotalReadsSE_BQ30/statsTotalReadsSE_BpLen*100,
    Raw_Percent_CG=statsTotalReadsCG/(statsTotalReadsSE_BpLen)*100,
    Raw_Chars_N = statsTotalReadsN,
    PhiX_IN = seqScrTotalReadsIn,
    PhiX_Discard = seqScrTotalReadsHits,
    PhiX_Percent_Discard = seqScrTotalReadsHits/seqScrTotalReadsIn*100,
    rRNA_In = seqScr2TotalReadsIn,
    rRNA_Identified = seqScr2TotalReadsHits,
    rRNA_Percent_Identified = seqScr2TotalReadsHits/seqScr2TotalReadsIn*100,
    AdapterTrimmed_IN = AdaptTotalReadsIn,
    AdapterTrimmed_Reads = AdaptTotalAdaptTrim,
    AdapterTrimmed_Percent_Reads = AdaptTotalAdaptTrim/AdaptTotalReadsIn*100,
    AdapterTrimmed_BP = AdaptTotalBpTrim,
    HardTrim_IN = CutTrimHardTotalReadsIn,
    HardTrim_SE_LeftBpTrim = CutTrimHardTotalSELT,
    HardTrim_Discard = CutTrimHardTotalDiscard,
    QwindowTrimmed_IN = QwinTotalReadsIn,
    QwindowTrimmed_SE_LeftBpTrim = QwinTotalSELT,
    QwindowTrimmed_SE_RightBpTrim = QwinTotalSERT,
    QwindowTrimmed_Discard = QwinTotalDiscard,
    NcharTrimmed_IN = NwinTotalReadsIn,
    NcharTrimmed_SE_LeftBpTrim = NwinTotalSELT,
    NcharTrimmed_SE_RightBpTrim = NwinTotalSERT,
    NcharTrimmed_Discard = NwinTotalDiscard,
    MinLen_IN = CutTrimTotalReadsIn,
    MinLen_Discard = CutTrimTotalDiscard,
    Proc_Reads=stats2TotalReadsIn,
    Proc_Bp=(stats2TotalReadsSE_BpLen),
    Proc_SE_PercentQ30=stats2TotalReadsSE_BQ30/stats2TotalReadsSE_BpLen*100,
    Proc_Percent_CG=stats2TotalReadsCG/(stats2TotalReadsSE_BpLen)*100,
    Proc_Chars_N = stats2TotalReadsN,
    Final_Percent_Read=stats2TotalReadsIn/statsTotalReadsIn*100,
    Final_Percent_Bp=(stats2TotalReadsSE_BpLen)/(statsTotalReadsSE_BpLen)*100
)

rownames(LongTable) <- samples
write.table(round(LongTable,2),"summary_hts_SE.txt",sep="\t",row.names=T,col.names=T,quote=F)
