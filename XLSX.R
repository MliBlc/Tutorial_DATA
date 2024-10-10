library(openxlsx)

save_up_down_regulated_to_excel_limma <- function(results, comparison_name) {
  upregulated <- results[results$logFC > logfc_threshold & results$adj.P.Val <= fdr_threshold, ]
  
  downregulated <- results[results$logFC < -logfc_threshold & results$adj.P.Val <= fdr_threshold, ]
  
  wb_all <- createWorkbook()
  addWorksheet(wb_all, "All_Results")
  writeData(wb_all, sheet = "All_Results", results, rowNames = TRUE)
  saveWorkbook(wb_all, paste0(comparison_name, "Limma_all_results.xlsx"), overwrite = TRUE)
  
  wb_up <- createWorkbook()
  addWorksheet(wb_up, "Upregulated")
  writeData(wb_up, sheet = "Upregulated", upregulated, rowNames = TRUE)
  saveWorkbook(wb_up, paste0(comparison_name, "Limma_upregulated.xlsx"), overwrite = TRUE)
  
  wb_down <- createWorkbook()
  addWorksheet(wb_down, "Downregulated")
  writeData(wb_down, sheet = "Downregulated", downregulated, rowNames = TRUE)
  saveWorkbook(wb_down, paste0(comparison_name, "Limma_downregulated.xlsx"), overwrite = TRUE)
}


## Or create single workbook and add the results to it as worksheets
