library(tidyverse)
library(dplyr)
library(lubridate)
library(chron)
library(camtrapR)
library(epitools)
library(pivottabler)
library(openxlsx)
library(purrr)

wi_data <- read.csv("wildlife-insights//sequences.csv")
str(wi_data)

for (i in 1:nrow(wi_data)) {
  out1[i] <- strsplit(wi_data$deployment_id[i], split = "_")
  wi_data$Station[i] <- paste(out1[[i]][1], "_", out1[[i]][2], sep = "")
  wi_data$Session[i] <- paste(out1[[i]][3], "_", out1[[i]][4], sep = "")
}

# rename all Sciurid, Lagomorph and Deer sequences

wi_data <- wi_data |>
  mutate_at(c('class', 'order','family', 'genus', 'species', 
              'common_name', 'Station', 'Session'), as.factor) |>
  mutate(common_name = recode(common_name, 
                              `Sciurus Species` = "Sciurid",
                              `Sciuridae Family` = "Sciurid",
                              `Western Gray Squirrel` = "Sciurid",
                              `Mexican Flying Squirrel` = "Sciurid",
                              `Eastern Gray Squirrel` = "Sciurid",
                              `Eastern Chipmunk` = "Sciurid",
                              `Eastern Cottontail` = "Lagomorph",
                              `Lagomorpha Order` = "Lagomorph",
                              `Lagomorpha Order` = "Lagomorph",
                              `Rabbit and Hare Family` = "Lagomorph",
                              `Sylvilagus Species` = "Lagomorph",
                              `Sylvilagus Species` = "Lagomorph",
                              `Cervidae Family` = "White-tailed Deer",
                              `Odocoileus Species` = "White-tailed Deer"))

levels(wi_data$common_name)

# extract Bobcat and Coyote sequences
bob_coy <- wi_data |>
  mutate_at(c('class', 'order','family', 'genus', 'species', 
              'common_name', 'Station', 'Session'), as.factor) |>
  filter(common_name == "Bobcat" | common_name == "Coyote")
str(bob_coy)
summary(bob_coy)

# create new Date column
bob_coy$newDate <- as.Date(bob_coy$start_time, format = "%m/%d/%Y")

# extract time and month and create new  columns 
bob_coy$Date.Time <- as.POSIXct(bob_coy$start_time, format = "%m/%d/%Y %H:%M") 

# time of day
bob_coy$Time <- format(bob_coy$Date.Time, format = "%H:%M") 
# month of capture
bob_coy$Month <- format(bob_coy$Date.Time, format = "%m")

# create a column of species detection from the 'group_size' information
# while there is info on group size, we need to summarize th enumber of sequences
bob_coy$Detection <- ifelse(bob_coy$group_size == 0, 0, 1)

summary(bob_coy)

# start creating a Pivot Table that summarizes
# number of sequences per month per camera for Bobcats and Coyotes separately
pt <- PivotTable$new()
pt$addData(bob_coy)
pt$addColumnDataGroups("Month")
pt$addRowDataGroups("common_name")
pt$addRowDataGroups("Station")
pt$renderPivot()
pt$defineCalculation(calculationName="Detection", 
                     summariseExpression="n()") # these are counts of events (GOOD!), 
                                                # not sum of Number.individuals
pt$renderPivot()
pt$evaluatePivot()

wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
addWorksheet(wb, "Data")
pt$writeToExcelWorksheet(wb=wb, wsName="Data", 
                         topRowNumber=1, leftMostColumnNumber=1, applyStyles=FALSE)
saveWorkbook(wb, file="WI_bob_coy_data.xlsx", overwrite = TRUE)

# extract all main prey species: squirrels, hares and deer
prey <- wi_data |>
  mutate_at(c('class', 'order','family', 'genus', 'species', 
              'common_name', 'Station', 'Session'), as.factor) |>
  filter(common_name == "Sciurid" | common_name == "Lagomorph" | 
           common_name == "White-tailed Deer")
str(prey)
summary(prey)

# create new Date column
prey$newDate <- as.Date(prey$start_time, format = "%m/%d/%Y")

# extract time and month and create new  columns 
prey$Date.Time <- as.POSIXct(prey$start_time, format = "%m/%d/%Y %H:%M") 

# time of day
prey$Time <- format(prey$Date.Time, format = "%H:%M") 
# month of capture
prey$Month <- format(prey$Date.Time, format = "%m")

# create a column of species detection from the 'group_size' information
# while there is info on group size, we need to summarize th enumber of sequences
prey$Detection <- ifelse(prey$group_size == 0, 0, 1)

summary(prey)


# start creating a Pivot Table that summarizes
# number of sequences per month per camera for Bobcats and Coyotes separately
pt <- PivotTable$new()
pt$addData(prey)
pt$addColumnDataGroups("Month")
pt$addRowDataGroups("common_name")
pt$addRowDataGroups("Station")
pt$renderPivot()
pt$defineCalculation(calculationName="Detection", 
                     summariseExpression="n()") # these are counts of events (GOOD!), 
# not sum of Number.individuals
pt$renderPivot()
pt$evaluatePivot()

wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
addWorksheet(wb, "Data")
pt$writeToExcelWorksheet(wb=wb, wsName="Data", 
                         topRowNumber=1, leftMostColumnNumber=1, applyStyles=FALSE)
saveWorkbook(wb, file="WI_prey_data.xlsx", overwrite = TRUE)



#####
setwd("C:\\Users\\viore\\Dropbox\\Columbia\\projects\\student projects\\Henry\\Marissa_data\\ohio_camera_metadata")
files <- list.files("ohio_camera_metadata")

ou_data <- plyr::ldply(list.files(), read.csv, header=TRUE)
str(ou_data)

ou_data <- ou_data |>
  mutate_at(c('Station', 'Session','Species'), as.factor) |>
  mutate(Species = recode(Species, 
                              `Canis_latrans` = "Coyote",
                              `Lynx_rufus` = "Bobcat",
                              `lynx_rufus` = "Bobcat",
                              `sciurus_sp` = "Sciurid",
                              `Sciurus_sp` = "Sciurid",
                              `Tamias_striatus` = "Sciurid",
                              `Sylvilagus_floridanus` = "Lagomorph",
                              `Glaucomys_volans` = "Sciurid",
                              `Odocoileus_virginianus` = "White-tailed Deer"))

ou_data$Number_animals <- as.integer(ou_data$Number_animals)

levels(ou_data$Species)


# extract date and create new columns for Date and Month
ou_data$newDate <- as.Date(ou_data$Date, format = "%m/%d/%Y")
ou_data$Month <- format(ou_data$newDate, format = "%m")

# extract Time and reduce to hours and minutes (no seconds)
ou_data$newTime <- chron::times(ou_data$Time)
ou_data$newTimeHM <- substr(ou_data$newTime, start = 1, stop = 5)

# create column for Detection/Non-detection
ou_data$Detection <- ifelse(ou_data$Number_animals == 0, 0, 1)
str(ou_data)

# filter data to retain only the first image in a sequence "(1).JPG"
# create sequences
out <- ou_data |>
  filter(grepl('1).JPG', filename_new))

str(out)
head(out)
summary(out)

pt <- PivotTable$new()
pt$addData(out) # all data without 
#pt$addData(out3) # data >15 min window
pt$addColumnDataGroups("Month")
pt$addRowDataGroups("Species")
pt$addRowDataGroups("Station")
pt$renderPivot()
pt$defineCalculation(calculationName="Detection", 
                     summariseExpression="n()") # these are counts of events (GOOD!), not sum of Number.individuals
pt$renderPivot()
pt$evaluatePivot()

wb <- createWorkbook(creator = Sys.getenv("USERNAME"))
addWorksheet(wb, "Data")
pt$writeToExcelWorksheet(wb=wb, wsName="Data", 
                         topRowNumber=1, leftMostColumnNumber=1, applyStyles=FALSE)
saveWorkbook(wb, file="OU_allspp.xlsx", overwrite = TRUE)



