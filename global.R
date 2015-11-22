qc_header = read.table("header.txt",stringsAsFactors =  F)
qc_filters = read.table("merged_qc_new.txt")

library(dplyr)
library(ggplot2)
library(xtable)

## Load databases



fhet <- function(p) { 2*p*(1-p) }
##### In cases 37 females and 193 males
##### in controls 121 females and 86 males.
fhetX = function(p) { 2*(p)*(1-p) * (158/(278+158))}

get_tstv_total = function(qc_table){
  tstv <- length( qc_table[ ( qc_table$genotype == "A/G" |  qc_table$genotype == "G/A" | 
                            qc_table$genotype == "C/T" | qc_table$genotype == "T/C" ), 1 ] )
  return(tstv/(nrow(qc_table)-tstv))
}
get_tstv_included= function(included_sites,qc_table){
  n = sum(included_sites)
  tstv <- length( qc_table[ ( qc_table$genotype == "A/G" |  qc_table$genotype == "G/A" | 
                                qc_table$genotype == "C/T" | qc_table$genotype == "T/C" ) & included_sites, 1 ] )
  return(tstv/(n-tstv))
}
get_tstv_excluded = function(included_sites,qc_table){
  n = sum(!included_sites)
  tstv <- length( qc_table[ ( qc_table$genotype == "A/G" |  qc_table$genotype == "G/A" | 
                                qc_table$genotype == "C/T" | qc_table$genotype == "T/C" ) & !included_sites, 1 ] )
  return(tstv/(n-tstv))
}
