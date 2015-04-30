
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)
library(gridExtra)
get_data_set = function(include){
  return(qc_table[include,])
}
get_new_ggplot = function(data_set){
  return(ggplot(data=data_set))
}
hwe_thresholds = function(x,thresh=0.01) {
  x=x  %>% mutate( hwe_threshold = ifelse(hwe <= thresh,paste0("<=",thresh),paste0(">",thresh)))
  return(x)
}

users=c("smilefreak")

shinyServer(function(input, output, session) {
  
  # Table, which contains all QC values , which have been analysed
  #output$file_input
  
  output$qc_table = renderDataTable({
    qc_table[,input$show_vars, drop=FALSE]
  }, options = list(orderClasses = TRUE))
  # Plots for the QC metrics and tables.
  
  
  inclusion_criteria = reactive({
    include = rep(T, length(qc_table$pvalue))
    if(input$conditionedPanels == 'QC Plots Urate'){
      include[ qc_table$count_cases + qc_table$count_controls == 0] = F
      include[ (qc_table$filter %in%  input$filters) ] = F 
      include[ input$max_prop_ref  < qc_table$refAB_cases | input$max_prop_ref < qc_table$refAB_controls ] = F
      include[ input$min_prop_alt > qc_table$altAB_cases | input$min_prop_alt > qc_table$altAB_controls] =F 
      include[input$max_prop_alt ] = F 
      include[ qc_table$genotype_quality_mean_cases < input$gqm | qc_table$genotype_qualtiy_mean_controls < input$gqm ] = F
      include[ qc_table$hetAB_cases < input$low_het_ab_prop | qc_table$hetAB_controls < input$low_het_ab_prop] = F
      include[ qc_table$hetAB_cases > input$high_het_ab_prop | qc_table$hetAB_controls > input$high_het_ab_prop ]= F
      include[ qc_table$proportion_hetAB_20_80_cases > input$prop_20_80 
               | qc_table$proportion_hetAB_20_80_controls > input$prop_20_80 ] = F
      include[ qc_table$proportion_hetAB_30_70_controls > input$prop_30_70
               | qc_table$proportion_hetAB_30_70_cases > input$prop_30_70] = F
      include[ abs(qc_table$genotype_quality_mean_cases - qc_table$genotype_qualtiy_mean_controls) >= input$mean_difference ] =F
      
      proportion_gq_10 = qc_table$GQ10_nums/(qc_table$GQ50_num + qc_table$GQ10_num)
      include [input$number_samples_gq_10 > qc_table$GQ10_nums ] = F
      #print(include) 
    }else if(input$conditionedPanels == 'All Site QC'){  
      include[ input$max_prop_ref  < qc_table$refAB_cases | input$max_prop_ref < qc_table$refAB_controls ] = F
      include[ input$min_prop_alt > qc_table$altAB_cases | input$min_prop_alt > qc_table$altAB_controls] = F 
      include[input$max_prop_alt ] = F 
      include[ qc_table$genotype_quality_mean_cases < input$gqm | qc_table$genotype_qualtiy_mean_controls < input$gqm ] = F
      include[ qc_table$hetAB_cases < input$low_het_ab_prop | qc_table$hetAB_controls < input$low_het_ab_prop] = F
      include[ qc_table$hetAB_cases > input$high_het_ab_prop | qc_table$hetAB_controls > input$high_het_ab_prop ]= F
      include[ qc_table$proportion_hetAB_20_80_cases > input$prop_20_80 
               | qc_table$proportion_hetAB_20_80_controls > input$prop_20_80 ] = F
      include[ qc_table$proportion_hetAB_30_70_controls > input$prop_30_70
               | qc_table$proportion_hetAB_30_70_cases > input$prop_30_70] = F  
    }
    return(include)
  })
  output$variant_summary = renderTable({
    include = inclusion_criteria()
    no_variants = sum(include)
    excluded_variants = sum(!include)
    return(data.frame(included_variants=no_variants,exc_variants=excluded_variants))
  })
  output$ts_tv_table = renderTable({
    inc = inclusion_criteria()
    tstv= get_tstv_total()
    tstv_inc = get_tstv_included(inc)
    tstv_exc = get_tstv_excluded(inc)
    coding = ! (qc_table$annotation %in% c("intronic","intergenic"))
    tstv_inc_coding = get_tstv_included(inc & coding)
    tstv_exc_coding = get_tstv_excluded(inc & !coding)
    df_ts = data.frame(total_tstv=tstv,tstv_inc=tstv_inc,tstv_exc=tstv_exc,coding_inc=tstv_inc_coding, coding_exc=tstv_exc_coding)
    return(df_ts)
  })
  output$variant_summary2 = renderTable({
    include = inclusion_criteria()
    no_variants = sum(include)
    excluded_variants = sum(!include)
    return(data.frame(included_variants=no_variants,exc_variants=excluded_variants))
  })
  output$ts_tv_table2 = renderTable({
    inc = inclusion_criteria()
    tstv= get_tstv_total()
    tstv_inc = get_tstv_included(inc)
    tstv_exc = get_tstv_excluded(inc)
    df_ts = data.frame(total_tstv=tstv,tstv_inc=tstv_inc,tstv_exc=tstv_exc)
    return(df_ts)
  })
  output$partitioning =  
    output$singleton_table = renderTable({
      inc = inclusion_criteria()
      case_singletons = length(qc_table[ qc_table$count_cases == 1 & qc_table$count_controls==0 & inc,1])
      controls_singletons = length(qc_table[ qc_table$count_cases == 0 & qc_table$count_controls == 1 & inc,1])
      df_s = data.frame(case_singletons=case_singletons,control_singleton=controls_singletons)
      return(df_s)  
    })
  output$hwe_table = renderTable({
    inc = inclusion_criteria()
    case_hwe = sum((qc_table$hwe_cases[inc]> input$hwe_thresh)/sum(inc))
    controls_hwe = sum((qc_table$hwe_controls[inc]> input$hwe_thresh)/sum(inc))
    hwe = sum((qc_table$hwe[inc]>input$hwe_thresh)/sum(inc))
    df_h = data.frame(hwe=hwe,hwe_cases=case_hwe,hwe_controls=controls_hwe)
    return(df_h)
  })
  output$hwe_cases_plot=renderPlot({
    p_df = get_data_set(inclusion_criteria())
    p_df = hwe_thresholds(p_df,thresh=input$hwe_thresh) 
    #  print(p_df$hwe_threshold)
    p_this = get_new_ggplot(p_df)
    p_this = p_this + ggtitle("Cases")
    p_this = p_this + geom_point(aes(y=HET_cases, x=MAF_cases,colour=hwe_threshold))
    p_this = p_this + labs(x="Minor Allele Frequency",y="Heterozygosity") + 
      scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limit=c(0,1))
    plot(p_this)
  })
  output$hwe_controls_plot=renderPlot({
    p_df = get_data_set(inclusion_criteria())
    p_df = hwe_thresholds(p_df,thresh=input$hwe_thresh) 
    #   print(p_df$hwe_threshold)
    p_this = get_new_ggplot(p_df)
    p_this = p_this + ggtitle("Controls")
    p_this = p_this + geom_point(aes(y=HET_controls,x=MAF_controls, colour =hwe_threshold))
    p_this = p_this + labs(x="Minor Allele Frequency",y="Heterozygosity") + 
      scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limit=c(0,1))
    plot(p_this)
  })
  output$allelic_balance_controls=renderPlot({
    p_df = get_data_set(inclusion_criteria())
   # p_this = get_new_ggplot(p_df)
    p_one = ggplot() + ggtitle("Singletons") + xlim(0,1)
    sing = p_df[p_df$count_cases == 0 & p_df$count_controls == 1,]
    p_one = p_one + geom_histogram(data=sing,aes(x=hetAB_controls),binwidth=0.01) +  xlab("Average read heterozygosity")  + ylab("Number of sites")
    p_second = ggplot() + ggtitle("Doubletons")   + xlim(0,1) 
    double = p_df[p_df$count_cases == 0 & p_df$count_controls == 2,]
    p_second= p_second + geom_histogram(data=double,aes(hetAB_controls),binwidth=0.01) + xlab("Average read heterozygosity")  + ylab("Number of sites")
    grid.arrange(p_one,p_second,ncol=2,main="Controls Heterozygosity")
    #  plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
  })
  output$allelic_balance_cases= renderPlot({
    p_df = get_data_set(inclusion_criteria())
  #  p_this = get_new_ggplot(p_df)
    p_one2 = ggplot() + ggtitle("Singletons") + xlim(0,1) 
    sing = p_df[p_df$count_cases == 1 & p_df$count_controls == 0,]
    p_one2 = p_one2 + geom_histogram(data=sing,aes(x=hetAB_cases),binwidth=0.01) +  xlab("Average read heterozygosity")  + ylab("Number of sites")
    p_second2 = ggplot() + ggtitle("Doubletons")   + xlim(0,1) 
    double = p_df[p_df$count_cases == 2 & p_df$count_controls == 0,]
    p_second2= p_second2 + geom_histogram(data=double,aes(hetAB_cases),binwidth=0.01) + xlab("Average read heterozygosity")  + ylab("Number of sites")
    grid.arrange(p_one2,p_second2,ncol=2,main="Cases Heterozygosity")
    #plot(p_one)
    #  plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
  })
  output$overall_allelic_balance = renderPlot({
    p_df = get_data_set(inclusion_criteria())
    p_overall= ggplot() + ggtitle("Overall") + xlim(0,1)
    sing = p_df[(p_df$count_cases + p_df$count_controls) <3, ]
    p_overall = p_overall + geom_histogram(data=sing,aes(x=HET),binwidth=0.01) +  xlab("Average read heterozygosity (Singletons and Doubletons)")  + ylab("Number of sites")
    p_overall
  })
  output$downloadSites = downloadHandler(
    filename = function() {
      paste0("sites-",Sys.Date(), '.txt')
    },
    content = function(con){
      p_df  = get_data_set(inclusion_criteria())
      #print(p_df)
      write.table(p_df$variant,quote= F , col.names=F , row.names = F,file=con)
    }
  )
  output$qqplot = renderPlot({
    p_df = get_data_set(inclusion_criteria())
    p_this = get_new_ggplot(p_df)
    p_this = p_this + ggtitle("QQPlot")
    observed = sort(p_df$pvalue)
    lobs = -(log10(observed))
    expected = c(1:length(observed))
    #print(expected)
    lexp = -(log10(expected / (length(expected)+1)))
    #  plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
    qq_gg = data.frame(obs=lobs,expect=lexp)
    p_this = p_this + geom_abline(intercept=0,slope=1,aes(colour='red')) 
    p_this = p_this + geom_point(data=qq_gg,aes(y=obs,x=expect)) + labs(list(x="Expected (-logP)",y="Observed (-logP)"))
    p_this = p_this + scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,10))
    plot(p_this)  
    #  print(lexp)
    # plot(lexp, lobs)  
    #  plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
    #   points(lexp, lobs, pch=23, cex=.4, bg="black") 
  })
  
  output$linear_model = renderTable({
    p_df = get_data_set(inclusion_criteria())
    ps = -(log10(p_df$pvalue))
    #cover_diff = abs(p_df$coverage_mean_cases - p_df$coverage_mean_controls)
    model = lm(ps ~ p_df$genotype_quality_mean_cases+p_df$genotype_qualtiy_mean_controls+
                 p_df$proportion_hetAB_20_80_cases+p_df$proportion_hetAB_20_80_controls+
                 p_df$GQ50_num+p_df$GQ10_num+p_df$hwe+ p_df$MAF_cases)
    return(summary(model))  
  })
  
  #### THIS IS THE STUFF FOR THE OTHER FULL FILE FILTER
  
  output$hwe_plot = renderPlot({
    p_df = get_data_set(inclusion_criteria())
    p_df = hwe_thresholds(p_df) 
    #  print(p_df$hwe_threshold)
    p_this = get_new_ggplot(p_df)
    p_this = p_this + ggtitle("ALL SITES")
    p_this = p_this + geom_point(aes(y=HET, x=MAF,colour=hwe_threshold))
    p_this = p_this + labs(x="Minor Allele Frequency",y="Heterozygosity") 
    #+ 
    #  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limit=c(0,1))
    plot(p_this)  
  })
  #
  
  
})
