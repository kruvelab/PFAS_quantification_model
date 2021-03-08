library(janitor)

read_excel_allsheets <- function(filename, tibble = TRUE) {
  data = tibble()
  sheets <- readxl::excel_sheets(filename)
  for (sheet in sheets) {
    data_this_sheet = readxl::read_excel(filename, sheet = sheet, skip = 0)
    data_this_sheet = data_this_sheet %>%
      mutate(across(everything(), as.character))
    #print(data_this_sheet)
    data_this_sheet = data_this_sheet %>%
      group_by(Compound, Filename) %>%
      summarise(Area = max(as.double(Area)%>%na.omit()),
             RT = mean(as.double(`Actual RT`) %>% na.omit()),
             `Theoretical Amt` = mean(as.double(str_replace(`Theoretical Amt`,pattern = ",", replacement = ".")) %>% na.omit() )) %>%
      ungroup()
    data = data %>%
      bind_rows(data_this_sheet)
  }
  return(data)
}



