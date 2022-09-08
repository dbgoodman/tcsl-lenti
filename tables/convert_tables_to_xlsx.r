library(data.table)
library(openxlsx)
table_dir = '/Volumes/HFSE/Box/tcsl/tcsl-lenti/tables'

load_csvs <- function(filenames) {
  
  table_list <- list()
  
  descriptions <- data.table()
  message('filename')
  
  for (filename in filenames) {
    file_id <- gsub('.*/(.*)\\.csv','\\1', filename)
    table <- fread(filename)
    header <- as.character(fread(filename, header=F, skip=0, nrows=1, sep = '|')[1,1])
    header <- gsub('# ','',header)
    descriptions <- rbind(descriptions, data.table('Workbook Name'= file_id, 'Description'=header))
    
    table_list[[file_id]] <- table
  }
  table_list <- c(list('descriptions'=descriptions), table_list)
  
  return(table_list)
}

table_list <- load_csvs(Sys.glob(paste0(table_dir,'/*.csv')))

write.xlsx(table_list, file = "supp_raw_data.xlsx")