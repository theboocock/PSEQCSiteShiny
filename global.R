qc_m = read.table("merged_qc_new.txt")

library(dplyr)
library(ggplot2)
library(xtable)

## Load databases



fhet <- function(p) { 2*p*(1-p) }
##### In cases 37 females and 193 males
##### in controls 121 females and 86 males.
fhetX = function(p) { 2*(p)*(1-p) * (158/(278+158))}

get_tstv_total = function(){
  tstv <- length( qc_table[ ( qc_table$genotype == "A/G" |  qc_table$genotype == "G/A" | 
                            qc_table$genotype == "C/T" | qc_table$genotype == "T/C" ), 1 ] )
  return(tstv/(nrow(qc_table)-tstv))
}
get_tstv_included= function(included_sites){
  n = sum(included_sites)
  tstv <- length( qc_table[ ( qc_table$genotype == "A/G" |  qc_table$genotype == "G/A" | 
                                qc_table$genotype == "C/T" | qc_table$genotype == "T/C" ) & included_sites, 1 ] )
  return(tstv/(n-tstv))
}
get_tstv_excluded = function(included_sites){
  n = sum(!included_sites)
  tstv <- length( qc_table[ ( qc_table$genotype == "A/G" |  qc_table$genotype == "G/A" | 
                                qc_table$genotype == "C/T" | qc_table$genotype == "T/C" ) & !included_sites, 1 ] )
  return(tstv/(n-tstv))
}
qc_m = qc_m %>% mutate(HET =  ((qc_m$HETA + qc_m$HETU) / (qc_m$OBSA + qc_m$OBSU)) )
qc_m = qc_m %>% mutate(MAF_2 = ((qc_m$MAC.cases +  qc_m$MAC.controls)/ ((qc_m$OBSA + qc_m$OBSU)*2)))
qc_table = qc_m %>%
  select(
    pvalue = P,
    genotype = REF.ALT,
    variant = VAR,
    MAF_cases = MAF.cases,
    MAF_controls = MAF.controls,
    HET_cases = HET.cases,
    HET_controls = HET.controls,
    GQ50_num = GQ50,
    GQ10_num = GQ10,
    count_cases = CNTA,
    count_controls = CNTU,
    hwe = HWE,
    HET = SAMPLE_HET_AB,
    MAF = MAF,
    filter = FILTER.x,        
    coverage_mean = DPM,
    genotype_quality_mean =GQM,
    proportion_hetAB_30_70=PROP_DEVIANT_HET_AB_30_70,
    proportion_hetAB_20_80=PROP_DEVIANT_HET_AB_20_80,
    altAB=SAMPLE_HOM_AB,
    refAB=SAMPLE_REF_AB,
    # Filter indivdual sites also
    annotation = ANNOT,  
    hwe_cases = HWE.cases,
    hwe_controls = HWE.controls, 
    hetAB_cases = SAMPLE_HET_AB.cases,
    hetAB_controls = SAMPLE_HET_AB.controls,
    proportion_hetAB_30_70_cases = PROP_DEVIANT_HET_AB_30_70.cases,
    proportion_hetAB_30_70_controls =PROP_DEVIANT_HET_AB_30_70.controls,
    proportion_hetAB_20_80_controls = PROP_DEVIANT_HET_AB_20_80.controls,
    proportion_hetAB_20_80_cases = PROP_DEVIANT_HET_AB_20_80.cases,
    genotype_quality_mean_cases = GQM.cases,
    genotype_qualtiy_mean_controls= GQM.controls,
    coverage_mean_controls = DPM.controls,
    coverage_mean_cases = DPM.cases,
    refAB_cases = SAMPLE_REF_AB.cases,
    refAB_controls = SAMPLE_REF_AB.controls,
    altAB_cases = SAMPLE_HOM_AB.cases,
    altAB_controls = SAMPLE_HOM_AB.controls,
    filter = FILTER.x
    )

print(nrow(qc_table))
qc_table=qc_table[ !(qc_table$count_cases + qc_table$count_controls == 0),]
print(nrow(qc_table))


n = nrow(qc_table)
