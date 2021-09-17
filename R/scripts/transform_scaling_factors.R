source("/home/simulator/scripts/load_scaling_factors.R")


combine_all <-function(wanted_types=NULL, with_ercc=T){
  
  if(is.null(wanted_types)){
    all_types <- data.frame(type=unique(unlist(c(maynard_census$type, quantiseq$type, epic$type, miltenyi$type, monaco$type, vento_census$type, hao$type, trivaglini_reads$type))))
    wanted_types <- all_types
  }
  
  maynard_reads_f <- merge(maynard_reads, wanted_types, by="type", all.y = T)
  maynard_reads_f$source <- rep("Maynard (reads)", length(wanted_types$type))

  maynard_census_f <- merge(maynard_census, wanted_types, by="type", all.y = T)
  maynard_census_f$source <- rep("Maynard (census)", length(wanted_types$type))

  quantiseq_f <- merge(quantiseq, wanted_types, by="type", all.y = T)
  quantiseq_f$source <- rep("quantiseq", length(wanted_types$type))
  
  epic_f <- merge(epic, wanted_types, by="type", all.y = T)
  epic_f$source <- rep("epic", length(wanted_types$type))
  
  miltenyi_f <- merge(miltenyi, wanted_types, by="type", all.y = T)
  miltenyi_f$source <- rep("miltenyi", length(wanted_types$type))
  
  vento_tormo_reads_f <- merge(vento_reads, wanted_types, by="type", all.y = T)
  vento_tormo_reads_f$source <- rep("Vento-tormo (reads)", length(wanted_types$type))
  vento_tormo_reads_f[vento_tormo_reads_f$type %in% c("T cells CD4", "T cells CD8", "T regulatory cells")]$mRNA <- vento_tormo_reads_f[vento_tormo_reads_f$type == "T cells"]$mRNA
  
  vento_tormo_census_f <- merge(vento_census, wanted_types, by="type", all.y = T)
  vento_tormo_census_f$source <- rep("Vento-tormo (census)", length(wanted_types$type))
  vento_tormo_census_f[vento_tormo_census_f$type %in% c("T cells CD4", "T cells CD8", "T regulatory cells")]$mRNA <- vento_tormo_census_f[vento_tormo_census_f$type == "T cells"]$mRNA
  
  vento_tormo_psmb2_f <- merge(vento_psmb2, wanted_types, by="type",all.y = T)
  vento_tormo_psmb2_f$source <- rep("Vento-tormo (PSMB2)", length(wanted_types$type))
  vento_tormo_psmb2_f[vento_tormo_psmb2_f$type %in% c("T cells CD4", "T cells CD8", "T regulatory cells")]$mRNA <- vento_tormo_psmb2_f[vento_tormo_psmb2_f$type == "T cells"]$mRNA
   
  monaco_f <- merge(monaco, wanted_types, by="type", all.y = T)
  monaco_f$source <- rep("monaco", length(wanted_types$type))
  
  hao_reads_f <- merge(hao_reads, wanted_types,by="type", all.y=T)
  hao_reads_f$source <- rep("Hao (reads)", length(wanted_types$type))

  hao_census_f <- merge(hao_census, wanted_types,by="type", all.y=T)
  hao_census_f$source <- rep("Hao (census)", length(wanted_types$type))

  hao_psmb2_f <- merge(hao_psmb2, wanted_types, by="type",all.y=T)
  hao_psmb2_f$source <- rep("Hao (PSMB2)", length(wanted_types$type))

  trivaglini_reads_f <- merge(trivaglini_reads, wanted_types, by="type",all.y=T)
  trivaglini_reads_f$source <- rep("Trivaglini (reads)", length(wanted_types$type))
  trivaglini_spike_f <- merge(trivaglini_spike, wanted_types, by="type",all.y=T)
  trivaglini_spike_f$source <- rep("Trivaglini (spike-in %)", length(wanted_types$type))
  
  trivaglini_ercc_f <- merge(trivaglini_ercc, wanted_types, by="type",all.y=T)
  
  trivaglini_ercc_total_f <- merge(trivaglini_ercc_total, wanted_types, by="type",all.y=T)
  trivaglini_ercc_total_f$source <- rep("Trivaglini (spike-in total)", length(wanted_types$type))
  
  trivaglini_genes_per_ercc_f <- merge(trivaglini_genes_per_ercc, wanted_types, by="type",all.y=T)
  trivaglini_genes_per_ercc_f$source <- rep("Trivaglini (genes per ERCC)", length(wanted_types$type))
  
  data_lst <- list(maynard_reads_f, maynard_census_f, quantiseq_f, epic_f, miltenyi_f, vento_tormo_reads_f, vento_tormo_census_f, vento_tormo_psmb2_f, monaco_f,hao_reads_f,hao_census_f,hao_psmb2_f, trivaglini_reads_f, trivaglini_spike_f,trivaglini_ercc_total_f,trivaglini_genes_per_ercc_f)
  
  full <- rbindlist(data_lst)
  
  if(with_ercc){
    full <- rbind(full, gather(trivaglini_ercc_f, source, mRNA, ERCC.00002:ERCC.00171))
  }
  
  return(full)
}
