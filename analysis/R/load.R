library(tidyverse)
library(lubridate)
options(digits.secs = 3) # change options to see milliseconds in DateTime 
 

# folder <- choose.dir(caption = "Select main analysis folder")
folder <- "analysis/data/"
# load mastertable with file path information for each experimental night
mastertable <- read.csv2(file = paste0(folder, "metadata/mastertable.csv"), header = TRUE, sep = ";",
                         na.strings = "NA")

# load Table of conditions for each experimental night for each group and individual
conditions <- read.csv2(file = paste0(folder, "metadata/conditions.csv"), header = TRUE, sep = ";", dec = ".")

# helper function to read csv file with raw data and preprocess them
load_raw_csv <- function(path) {
  n_night <- read.csv2(file = path,  sep = ";", dec = ",", header = TRUE,
                      fileEncoding = "UTF-16LE", as.is = TRUE, row.names = NULL)
  
  n_night <- n_night %>%
    filter(!str_detect(DateTime, "#")) %>% 
    arrange(DateTime) %>% #sort chronologically
    mutate(DateTime = as.numeric(str_replace(DateTime, ",", ".")),
           DateTime = as.POSIXct(DateTime * (60 * 60 * 24),
                                 origin = "1899-12-30", tz = "UTC"))

  #find out first line of useful data
  firstrow <- n_night %>%
    mutate(rown = row_number(DateTime)) %>%
    filter(SystemMsg == "start") %>%
    summarise(firstrow = max(rown) + 1) %>%
    pull()
  if (is.na(firstrow) | is.infinite(firstrow) | firstrow > 2000) { # in rare cases the start line is missing
    # then take the first rewarded visit of an animal as the firstrow
    firstrow <- n_night %>%
      mutate(rown = row_number(DateTime)) %>%
      filter(IdLabel != "Test", reinforce1value > 0) %>%
      summarise(firstrow = min(rown)) %>%
      pull()
  }
 
  # find the last row number, larger than firstrow and preceding a row with "exp" in the
  # "unitLabel" column
  lastrow <- c()
  lastrow <- n_night %>%
    mutate(rown = row_number(DateTime)) %>%
    filter(SystemMsg == "end", rown > firstrow, str_detect(IdLabel, "T")) %>%
    summarise(lastrow = min(rown) - 1) %>%
    pull()
  if (is.na(lastrow) | is.infinite(lastrow)) {lastrow <- nrow(n_night)}

  #only select relevant data and discard the rest
  n_night <- n_night %>%
    slice(firstrow:lastrow)
}

nights <- as.list(mastertable$night)
paths <- as.list(paste0(folder, mastertable$path))

# path <- paths %>%  pluck(1)

# function for aggregating data from all nights and adding correct night column
aggregate_nights <- function(paths, nights) {
  map(paths, load_raw_csv) %>%
    set_names(nights) %>%
    enframe("night", "night_data") %>%
    unnest(night_data) %>%
    mutate(night = as.numeric(night))
}

allnights <- aggregate_nights(paths, nights) %>% 
  arrange(night, DateTime) %>% # sort by night and chronologically
  mutate(eventDuration = ifelse(is.na(eventDuration), sense1duration, eventDuration))

# remove columns with irrelevant data
allnights <- allnights %>%
  select(night, DateTime, IdLabel, IdRFID, unitLabel, eventDuration, reinforce1value, SystemMsg)

# mark when pump was in busy state
allnights <- allnights %>% 
  mutate(pumping = case_when(
    DateTime == min(DateTime) ~ 0,
    SystemMsg == "" ~ NA_real_,
    str_detect(SystemMsg, "start p") ~ 1,
    str_detect(SystemMsg, "end p") ~ 0)) %>%
  fill(pumping)

# repair missing identities
id_table <- conditions %>% 
  select(night, IdRFID, Id_actual = IdLabel) %>% 
  distinct()

allnights <- allnights %>% 
  left_join(id_table, by = c("night", "IdRFID")) %>%
  mutate(IdLabel = Id_actual) %>%
  select(-Id_actual)

# make new columns
allnights <- allnights %>%
  group_by(night) %>% 
  mutate(
    # for flexibility test, phase 1 = first pair of flowers rewarding, and phase 2 = second pair of flowers rewarding
    phase = ifelse(SystemMsg == "switch", 1, 0),
    phase = cumsum(phase) + 1,
    # numeric column for reward volume delivered, half a microliter per pump step
    vol = replace_na(reinforce1value, 0)/2,
    # numeric column for reward status 1=rewarded, 0=unrewarded, helps calculate reward rates
    rewarded = ifelse(vol > 0, 1, vol),
    # create location column from the unitLabels
    loc = as.numeric(str_extract(unitLabel, "[:digit:].*$"))
    ) %>%
  select(-reinforce1value) %>% 
  ungroup()

# merge the conditions table and the current data table to make new columns for dispenser properties
# such as maximal volume, sugar concentration, probability, etc., as well as experimental conditions
# nights to discard, etc.
conditions_table <- conditions %>% 
  select(-IdRFID) %>% 
  distinct() %>% 
  mutate(discard = as.character(discard))

add_conditions_discard_data <- function(tbl, conditions_table) {
  #discard data labeled for discarding in conditions table
  date_filters <- any(!conditions_table$discard %in% c(1, 0))
  
  tbl_w_conds <- tbl %>%
    left_join(conditions_table, by = c("night", "IdLabel")) %>%
    select(night, IdLabel, loc, everything()) %>% 
    mutate(discard = replace_na(discard, 0))
  
  if (!date_filters) {
    tbl_filtered <- tbl_w_conds %>% 
             filter(discard == 0) %>%
             select(-discard)
  } else {
    tbl_filtered <- tbl_w_conds %>% 
      group_by(night, IdLabel) %>% 
      mutate(start_date = as.POSIXct(as.Date(min(DateTime))),
             time_after = str_extract(discard, "(?<=>)[:digit:]{2}:[:digit:]{2}"),
             time_before = str_extract(discard, "(?<=<)[:digit:]{2}:[:digit:]{2}"),
             discard = as.numeric(discard == 1),
             time_after = case_when(
               is.na(time_after) & is.na(time_before) ~ NA_POSIXct_,
               is.na(time_after) & !is.na(time_before) ~ start_date,
               hm(time_after) < hm("16:00") ~ start_date + days(1) + hm(time_after),
               TRUE ~ start_date + hm(time_after)
             ),
             time_before = case_when(
               is.na(time_before) & is.na(time_after) ~ NA_POSIXct_,
               is.na(time_before) ~ start_date + days(2),
               hm(time_before) < hm("16:00") ~ start_date + days(1) + hm(time_before),
               TRUE ~ start_date + hm(time_before)
             ),
             discard_times = ifelse(DateTime %within% interval(time_after, time_before), 1, 0),
             discard_times = replace_na(discard_times, 0),
             discard = discard_times == 1 | discard == 1)  %>% 
      filter(!discard) %>%
      ungroup() %>% 
      select(-time_before, -time_after, -start_date, -discard_times, -discard)
  }
tbl_filtered
}

allnights <- allnights %>%
  add_conditions_discard_data(conditions_table) %>% 
  # left_join(flower_table, by = c("night", "IdLabel", "loc")) %>% 
  # sort chronologically again
  arrange(DateTime) %>% 
  mutate(loc = factor(loc),
         weight = as.numeric(as.character(weight)))

allnights <- allnights %>%
  filter(hour(DateTime) > 15 | hour(DateTime) < 4) %>% 
  #boolean column for whether an event was a choice at a flower (RFID signal and IR beam interrupted):
  mutate(choice = !is.na(loc) & str_detect(unitLabel, "Cond")) %>% 
  fill(group_night, group, cond)

write.table(allnights, file = paste0(folder, "EventData.csv"), sep = ";", row.names = FALSE)
