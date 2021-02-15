library(janitor)

read_excel_allsheets <- function(filename, tibble = TRUE) {
  data = tibble()
  sheets <- readxl::excel_sheets(filename)
  for (sheet in sheets) {
    data_this_sheet = readxl::read_excel(filename, sheet = sheet, skip = 0)
    data_this_sheet = data_this_sheet %>%
      mutate(across(everything(), as.character))
    #print(data_this_sheet)
    data = data %>%
      bind_rows(data_this_sheet)
    data = data %>%
      select(Compound, RT, Area, Formula, `Theoretical Amt`)
    data = data %>%
      mutate(`Theoretical Amt` = str_replace(`Theoretical Amt`,pattern = ",", replacement = "."))
  }
  return(data)
}



