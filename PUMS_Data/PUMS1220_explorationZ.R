################################################################
rm(list = ls())

# ---------- Download the data ----------------------------- #
# https://www2.census.gov/programs-surveys/acs/data/pums/2024/1-Year/
# Download "csv_pus.zip" and "csv_hus.zip" in the target directory, stored in "folder" below.


# ---------- Variables of Interest -----------------------------------------------------------#
# --- geographic variables --- #
# RT: Record Type; "H": "Housing Record or Group Quarters Unit", "P": "Person Record".
# SERIALNO: Housing unit/GQ person serial number; string
# DIVISION: "7": "West South Central (South Region)",
#           "2": "Middle Atlantic (Northeast region)",
#           "4": "West North Central (Midwest region)",
#           "6": "East South Central (South region)",
#           "3": "East North Central (Midwest region)",
#           "9": "Pacific (West region)",
#           "0": "Puerto Rico",
#           "5": "South Atlantic (South region)",
#           "8": "Mountain (West region)",
#           "1": "New England (Northeast region)"
# PUMA: Public use microdata area code (PUMA) with population of 100,000 or more
# STATE: State Code
#        find the state names in this link: https://www.census.gov/library/reference/code-lists/ansi/ansi-codes-for-states.html#accordion-455dd4b1cd-item-2c5946f1b7



# --- household level variables of interest, not existed in the person data --- #
# NP: Number of persons in this household
# RMSP: Number of rooms
# BDSP: Number of bedrooms
# YRBLT: When structure (house) first built:
#           "1990": "1990 to 1999",
#           "2010": "2010 to 2019",
#           "2022": "2022",
#           "1938": "N/A (GQ)",
#           "1940": "1940 to 1949",
#           "1970": "1970 to 1979",
#           "1980": "1980 to 1989",
#           "2020": "2020",
#           "2000": "2000 to 2009",
#           "2021": "2021",
#           "1939": "1939 or earlier",
#           "1950": "1950 to 1959",
#           "1960": "1960 to 1969"
# BLD: Units in structure
#               "0": "N/A (GQ)",
#               "01": "Mobile Home or Trailer",
#               "02": "One-family house detached",
#               "03": "One-family house attached",
#               "04": "2 Apartments",
#               "05": "3-4 Apartments",
#               "06": "5-9 Apartments",
#               "07": "10-19 Apartments"
#               "08": "20-49 Apartments",
#               "09": "50 or More Apartments",
#               "10": "Boat, RV, van, etc.",

# BATH: Bathtub or shower
#               "0": "N/A (GQ)"
#               "1": "Yes",
#               "2": "No",
# STOV: Stove or range
#               b .N/A (GQ)
#               1 .Yes
#               2 .No
# PLM: Complete plumbing facilities
#               "1": "Yes, has hot and cold running water, and a bathtub or shower",
#               "0": "N/A (GQ)"
#               "2": "No",
#               "9": "Case is from Puerto Rico, PLM recode not applicable"
# KIT: Complete kitchen facilities
#               "0": "N/A (GQ)",
#               "1": "Yes, has stove or range, refrigerator, and sink with a faucet",
#               "2": "No"
# SINK: Sink with a faucet
#               b .N/A (GQ)
#               1 .Yes
#               2 .No
# REFR: Refrigerator
#               b .N/A (GQ)
#               1 .Yes
#               2 .No
# RWAT: Hot and cold Running water
#               b .N/A (GQ)
#               1 .Yes
#               2 .No
#               9 .Case is from Puerto Rico, RWAT not applicable
# HFL: House heating fuel
#               b .N/A (GQ/vacant)
#               1 .Utility gas
#               2 .Bottled or tank gas (propane,butane,etc.)
#               3 .Electricity
#               4 .Fuel oil, kerosene, etc.
#               5 .Coal or coke
#               6 .Wood
#               7 .Solar energy
#               8 .Other fuel
#               9 .No fuel used

# TEN: Tenure
#               b .N/A (GQ/vacant)
#               1 .Owned with mortgage or loan (include home equity loans)
#               2 .Owned free and clear
#               3 .Rented
#               4 .Occupied without payment of rent

# VALP: Property value

# GASP: Gas cost (monthly cost)
# ELEP: Electricity cost (monthly cost)
# FULP: Fuel cost (yearly cost for fuels other than gas and electricity)
# VEH: 	Vehicles (1 ton or less) available
#           "5": "5 vehicles",
#           "2": "2 vehicles",
#           "-1": "N/A (GQ/vacant)",
#           "3": "3 vehicles",
#           "0": "No vehicles",
#           "4": "4 vehicles",
#           "1": "1 vehicle",
#           "6": "6 or more vehicles"
# HHL: Household language
#           "0": "N/A (GQ/vacant)",
#           "2": "Spanish",
#           "3": "Other Indo-European languages",
#           "1": "English Only",
#           "4": "Asian and Pacific Island languages",
#           "5": "Other Language"
# NOC: Number of own children in household (unweighted)
#           min: 1; max: 19.
# ----------------------------------------------------------------------------------------#


# --- person level variables of interest, not exist in the household data --- #
# JWMNP: Travel time to work (minutes)
# SCHL: the education level, the larger of SCHL, the higher the education
# WAGP: the salary income or wage
# PINCP: the total income, which is >= WAGP
# AGEP: age
# SEX: "1": "Male", "2": "Female".
# RAC1P: Recoded detailed race code
             # "1": "White alone",
             # "2": "Black or African American alone",
             # "3": "American Indian alone",
             # "4": "Alaska Native alone",  
             # "5": "American Indian and Alaska Native tribes specified; or American Indian or Alaska Native, not specified and no other races",
             # "6": "Asian alone",
             # "7": "Native Hawaiian and Other Pacific Islander alone",
             # "8": "Some other race alone"             
             # "9": "Two or More Races"
# MSP: marital status
             # "0": "N/A (age less than 15 years)",
             # "1": "Now Married, Spouse Present",             
             # "2": "Now Married, Spouse Absent",
             # "3": "Widowed",
             # "4": "Divorced",
             # "5": "Separated",
             # "6": "Never Married"
# WKHP: Usual hours worked per week past 12 months (8 out of 17 NA's)
             # min: 1; max: 99. 
# WKWN: Weeks worked during past 12 months (8 out of 17 NA's)
             # min: 1; max: 52. 
# ESR: Employment status recode (27 out of 173 NA's)
             # "0": "N/A (less than 16 years old)",
             # "1": "Civilian employed, at work",
             # "2": "Civilian employed, with a job but not at work",
             # "3": "Unemployed",
             # "4": "Armed Forces, At Work"
             # "5": "Armed Forces, With a Job But Not At Work",
             # "6": "Not in Labor Force".
# COW: Class of worker (7 out of 17 NA's)
             # "0": "N/A (less than 16 years old/NILF who last worked more than 5 years ago or never worked)",
             # "1": "Employee of a private for-profit co. or bus., or of an individual, for wages, salary, commissions",
             # "2": "Employee of a private not-for-profit, tax-exempt, or charitable organization",
             # "3": "Local government employee (city, county, etc.)",
             # "4": "State government employee",
             # "5": "Federal government employee",
             # "6": "Self-employed in own not incorporated business, professional practice, or farm",
             # "7": "Self-employed in own incorporated business, professional practice or farm",
             # "8": "Working without pay in family business or farm",
             # "9": "Unemployed and last worked 5 years ago or earlier or never worked"
# CIT: Citizenship status (almost no NA)
             # "1": "Born in the United States",
             # "2": "Born in Puerto Rico, Guam, the U.S. Virgin Islands, or Northern Marianas"
             # "3": "Born abroad of U.S. citizen parent or parents",
             # "4": "U.S. citizen by naturalization",
             # "5": "Not a U.S. citizen",
# DIS: Disability recode (almost no NA)
             # "1": "With a disability". "2": "Without a disability"
# INDP: Industry recode for 2023 and later based on 2022 IND codes (7 out of 17 NA's)
             # totally 266 industry categories
# OCCP: Occupation (7 out of 17 NA's)
             # 531 occupation categories
# --------------------------------------------------------------------------------------- #



# X: PINCP is personal variable
# Y: VALP: Property value, household variable






##################################################################
# ----------- Read the person data -------------------------- #
rm(list = ls())

######################################################################################################################
## ------ Please adjust the following folder to the directory 
# ------- which contains the downloaded data from https://www2.census.gov/programs-surveys/acs/data/pums/2024/1-Year/
######################################################################################################################
# folder = "D:\\Research\\20251214_PUMS"
folder = Sys.getenv("PUMS_DATA_DIR", unset = "data")

Z_vars_P <- c(
  # Demographics
  "AGEP", "SEX", "MSP", "RAC1P",
  
  # Human capital
  "SCHL", "DIS"
  
)

Z_var_H = c(
  # Language
  "HHL", 
  
  # Housing size / structure
  "RMSP", "YRBLT", "BLD"
)

Rsps = "VALP" # household
Prdct = "PINCP" # personal

# the person data is separated into psam_pusa.csv and psam_pusb.csv

psam_pusa = read.csv(file.path(folder, "csv_pus", "psam_pusa.csv"), header = T)

str(psam_pusa)
psam_pusa[1,]
c("PINCP", Z_vars_P) %in% colnames(psam_pusa)

Z_var_H %in% colnames(psam_pusa)

pusa = psam_pusa[,c("RT", "SERIALNO", "DIVISION", "PUMA", "STATE", 
                    Z_vars_P,
                    "PINCP" )]
rm(psam_pusa) # remove the large dataset
pusa = pusa[complete.cases(pusa), ] # remove the missing value NA's


pusa = pusa[(pusa$AGEP>=18), ] # only consider those people over 18 years old.
table(pusa$STATE) # pusa is the data from state: 1,2,4:6, 8:13, 15:28, same as husa
save(list = c("pusa"), file = file.path(folder, "pusa.Rdata"))


rm(list = setdiff(ls(), c("Z_var_H", "Z_vars_P", "folder")))

psam_pusb = read.csv(file.path(folder, "csv_pus", "psam_pusb.csv"),
                     header = T)
pusb = psam_pusb[,c("RT", "SERIALNO", "DIVISION", "PUMA", "STATE", 
                    Z_vars_P,
                    "PINCP" )] 
rm(psam_pusb) # remove the large dataset
sum(complete.cases(pusb))
pusb = pusb[complete.cases(pusb), ] # remove the missing value NA's
pusb = pusb[(pusb$AGEP>=18), ]
table(pusb$STATE) # pusb is the data from state: 29:42, 44:51, 53:56, same as husb
save(list = c("pusb"), file = file.path(folder, "pusb.Rdata"))




# ----------- connect pusa.Rdata and pusb.Rdata --------------------- #
rm(list = setdiff(ls(), c("Z_var_H", "Z_vars_P", "folder")))

load(file.path(folder, "pusa.Rdata"))
load(file.path(folder, "pusb.Rdata"))

pus = rbind(pusa, pusb)

save(list = "pus", file = file.path(folder, "pus.Rdata"))





#----------- Pick the representative of each house by the highest PINP ---------- #
# --- Pick the person in the household with the highest income as the representative
# --- Combined this person with the household data, such as VALP, YRBLT, GASP, ...
rm(list = setdiff(ls(), c("Z_var_H", "Z_vars_P", "folder")))

load(file.path(folder, "pus.Rdata"))

anyNA(pus) # check whether there is any NA in this data frame (FALSE)

str(pus)

dim(table(pus$STATE)) # 51 states (50 states + District of Columbia)
sort(unique(pus$STATE))


table(pus$SCHL)



# -------- the following codes check the average total income for different states
library(dplyr)

sum(duplicated(pus$SERIALNO)) # a lot of duplicates, because there are many people from the same house

representative_individuals <- pus %>%
  group_by(SERIALNO) %>%  # Group by household ID
  filter(PINCP == max(PINCP, na.rm = TRUE)) %>%  # Get the individual with the highest PINCP (total income)
  slice(1) %>%  # If there are ties, take the first one
  ungroup()  # Ungroup to avoid issues in the next operations

str(representative_individuals)

# Rename columns in 'representative_individuals' to avoid conflicts
representative_individuals <- representative_individuals %>%
  rename(RT_rep = RT, DIVISION_rep = DIVISION, ST_rep = STATE, PUMA_rep = PUMA)

str(representative_individuals)

save(list = c("representative_individuals"), file = file.path(folder, "rep_ind.Rdata"))

# ----------- representative_individuals will finally be combined with the household data ---- #






############################################################################################
# 
############################################################################################
# --------- Read the household data ------------------------------------------------------ #
# the household data is separated into psam_husa.csv and psam_husb.csv
rm(list = setdiff(ls(), c("Z_var_H", "Z_vars_P", "folder")))


psam_husa = read.csv(file.path(folder, "csv_hus", "psam_husa.csv"),
                     header = T) # read the data
c("RT", "SERIALNO", "DIVISION", "PUMA", "STATE", "VALP", Z_var_H) %in% colnames(psam_husa)

nrow(psam_husa) # 827133
sum(is.na(psam_husa$SERIALNO)) # 0
sum(is.na(psam_husa$VALP)) # 321514
sum(is.na(psam_husa$HHL))      # 143009


husa = psam_husa[,c("RT", "SERIALNO", "DIVISION",  "STATE", "PUMA",
                    "VALP",
                    Z_var_H
)] # only take these variables




rm(psam_husa) # remove the large dataset
husa = husa[complete.cases(husa), ] # remove the missing value NA's
table(husa$STATE) # husa is the data from state: 1,2,4:6, 8:13, 15:28



psam_husb = read.csv( file.path(folder, "csv_hus", "psam_husb.csv"),
                     header = T)

husb = psam_husb[,c("RT", "SERIALNO", "DIVISION",  "STATE", "PUMA",
                    "VALP",
                    Z_var_H 
)]
rm(psam_husb)
husb = husb[complete.cases(husb), ]
table(husb$STATE) # husb is the data from state: 29:42, 44:51, 53:56

hus = rbind(husa, husb)

dim(table(hus$STATE)) # 51 states (50 states + District of Columbia)

sum(duplicated(hus$SERIALNO)) # 0, so there is no duplicated SERIALNO in hus








##############################################################
# --------- combine the representative_individuals (person with highest income in each family) and hus (family-wise data)
rm(husa, husb)

load(file.path(folder, "rep_ind.Rdata"))


### the combined data frame is hpus ###
### You start with a data frame hus.
### Then you do a left join with another data frame representative_individuals.
### by = "SERIALNO" means the two datasets are matched on the column SERIALNO 
### A left join keeps all rows from hus, and adds matching columns from representative_individuals. If there’s no match, you’ll get NAs.
hpus <- hus %>%
  left_join(representative_individuals, by = "SERIALNO") 

str(hpus)

hpus = hpus[,c("SERIALNO", "DIVISION", "PUMA", "STATE", 
               "VALP", "PINCP", 
               Z_var_H, Z_vars_P)]


str(hpus)

anyNA(hpus)

hpus[1:30,]

# check some cases
representative_individuals[( representative_individuals$SERIALNO == (hpus$SERIALNO[17]) ), ]
representative_individuals[( representative_individuals$SERIALNO == (hpus$SERIALNO[32]) ), ]


hpus = hpus[complete.cases(hpus), ] 

nrow(hpus) # 988629

dim(table(hpus$STATE)) # 51
sort(unique(hpus$STATE))

anyNA(hpus)

# original SEX: "1": "Male", "2": "Female"
hpus$SEX = 2 - hpus$SEX 
# Now the SEX: "1": male; "0": female

head(hpus)

sort(unique(hpus$YRBLT)) # 1939 1940 1950 1960 1970 1980 1990 2000 2010 2020 2021 2022 2023 2024
hist(hpus$YRBLT)



#########################################################
dat_tgt = hpus[hpus$PINCP>0,] # This is the data frame we will use for regression.

dim(dat_tgt) # The sample contains 980254 people as representative of the family.


length(unique(dat_tgt$STATE)) # 51 states

sort(unique(dat_tgt$STATE))
# state numbers
# c(1,  2,  4,  5,  6 , 8 , 9 ,10 ,11, 12, 13, 15, 16, 17, 18, 19, 20,
#   21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 3,8 39,
#   40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 55, 56)


dat_tgt%>%group_by(DIVISION)%>%summarise(Number_of_States = length(unique(STATE)))
#     DIVISION Number_of_States
#        <int>            <int>
# 1        1                6
# 2        2                3
# 3        3                5
# 4        4                7
# 5        5                9
# 6        6                4
# 7        7                4
# 8        8                8
# 9        9                5

rm(hpus, hus, representative_individuals)



str(dat_tgt)

# original DIS: "1": "With disability", "2": "Without disability"
dat_tgt$DIS = 2 - dat_tgt$DIS 
sum(dat_tgt$DIS)/nrow(dat_tgt) # 0.1634597
# Now the DIS: "1": "With disability", "0": "Without disability"


# --- cat variables need to be recoded --- #
unique(dat_tgt$HHL) # 1 3 2 4 5
# Household language:
#           "1": "English Only",
#           "2": "Spanish",
#           "3": "Other Indo-European languages",
#           "4": "Asian and Pacific Island languages",
#           "5": "Other Language"# 

sort(unique(dat_tgt$YRBLT)) # 1939 1940 1950 1960 1970 1980 1990 2000 2010 2020 2021 2022 2023 2024

sort(unique(dat_tgt$BLD)) # 1  2  3  4  5  6  7  8  9 10
hist(dat_tgt$BLD)

sort(unique(dat_tgt$MSP)) # 1 2 3 4 5 6
# MSP: marital status
# "1": "Now Married, Spouse Present",             
# "2": "Now Married, Spouse Absent",
# "3": "Widowed",
# "4": "Divorced",
# "5": "Separated",
# "6": "Never Married"
hist(dat_tgt$MSP)

sort(unique(dat_tgt$RAC1P)) # 1 2 3 4 5 6 7 8 9
# race:
# "1": "White alone",
# "2": "Black or African American alone",
# "3": "American Indian alone",
# "4": "Alaska Native alone",  
# "5": "American Indian and Alaska Native tribes specified; or American Indian or Alaska Native, not specified and no other races",
# "6": "Asian alone",
# "7": "Native Hawaiian and Other Pacific Islander alone",
# "8": "Some other race alone"             
# "9": "Two or More Races"

sort(unique(dat_tgt$SCHL)) # 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
# ------------------------------------------#




head(dat_tgt)


# Categorical variables recoded
# HHL
dat_tgt$HHL_Spn = 1 * (dat_tgt$HHL == 2)
dat_tgt$HHL_Other <- 1 * (dat_tgt$HHL %in% c(3, 4, 5) )

# YRBLT
dat_tgt$YRBLT1940 = 1 * (dat_tgt$YRBLT == 1940)
dat_tgt$YRBLT1950 = 1 * (dat_tgt$YRBLT == 1950)
dat_tgt$YRBLT1960 = 1 * (dat_tgt$YRBLT == 1960)
dat_tgt$YRBLT1970 = 1 * (dat_tgt$YRBLT == 1970)
dat_tgt$YRBLT1980 = 1 * (dat_tgt$YRBLT == 1980)
dat_tgt$YRBLT1990 = 1 * (dat_tgt$YRBLT == 1990)
dat_tgt$YRBLT2000 = 1 * (dat_tgt$YRBLT == 2000)
dat_tgt$YRBLT2010 = 1 * (dat_tgt$YRBLT == 2010)
dat_tgt$YRBLT2020 = 1 * (dat_tgt$YRBLT == 2020)
dat_tgt$YRBLT2021 = 1 * (dat_tgt$YRBLT == 2021)
dat_tgt$YRBLT2022 = 1 * (dat_tgt$YRBLT == 2022)
dat_tgt$YRBLT2023 = 1 * (dat_tgt$YRBLT == 2023)
dat_tgt$YRBLT2024 = 1 * (dat_tgt$YRBLT == 2024)

# BLD
dat_tgt$BLD1U = 1 * (dat_tgt$BLD %in% c(2,3) ) # one house, attached or detached
dat_tgt$BLDMU = 1 * (dat_tgt$BLD %in% c(4,5,6,7,8,9) ) # more than one apartment but not mobile home or RV

# MSP
dat_tgt$MSP_sgl = 1 * (dat_tgt$MSP == 6)

# RAC1P
dat_tgt$RAC1P2 = 1 * (dat_tgt$RAC1P == 2)
dat_tgt$RAC1P3 = 1 * (dat_tgt$RAC1P == 3)
dat_tgt$RAC1P4 = 1 * (dat_tgt$RAC1P == 4)
dat_tgt$RAC1P5 = 1 * (dat_tgt$RAC1P == 5)
dat_tgt$RAC1P6 = 1 * (dat_tgt$RAC1P == 6)
dat_tgt$RAC1P7 = 1 * (dat_tgt$RAC1P == 7)
dat_tgt$RAC1P8 = 1 * (dat_tgt$RAC1P == 8)
dat_tgt$RAC1P9 = 1 * (dat_tgt$RAC1P == 9)


# SCHL
# dat_tgt$SCHL_LESSHS = 1 * (dat_tgt$SCHL <= 15)
dat_tgt$SCHL_HSCOMP = 1 * (dat_tgt$SCHL %in% c(16, 17, 18, 19, 20) )
dat_tgt$SCHL_COCOMP = 1 * (dat_tgt$SCHL %in% c(21, 22, 23, 24) )




hist(dat_tgt[,"VALP"])
dat_tgt$logY = log(dat_tgt[,"VALP"])
hist(dat_tgt$logY)

hist(dat_tgt$PINCP)
dat_tgt$logX = log(dat_tgt$PINCP)
hist(dat_tgt$logX)













#############################################################
# --- regression ---- #
# --- Based on dat_tgt --- #
###################################################
rm(list = setdiff(ls(), c("dat_tgt", "folder")))

 
head(dat_tgt)

str(dat_tgt)


# --- fit Y
fit.Y_fml = "logY ~ logX + HHL_Spn + HHL_Other + RMSP + YRBLT1940 + YRBLT1950 + YRBLT1960 + YRBLT1970 + YRBLT1980 + YRBLT1990 + YRBLT2000 + YRBLT2010 + YRBLT2020 + YRBLT2021 + YRBLT2022 + YRBLT2023 + YRBLT2024 + BLD1U + BLDMU + AGEP + SEX + MSP_sgl + RAC1P2 + RAC1P3 + RAC1P4 + RAC1P5 + RAC1P6 + RAC1P7 + RAC1P8 + RAC1P9 + SCHL_HSCOMP + SCHL_COCOMP + DIS"
fit.Y = lm(fit.Y_fml, data = dat_tgt)

summary(fit.Y)

# Residuals:
#   Min      1Q  Median      3Q     Max 
# -7.3498 -0.3892  0.0508  0.4978  6.6582 

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  7.237e+00  1.435e-02 504.451  < 2e-16 ***
#   logX         2.760e-01  1.173e-03 235.406  < 2e-16 ***
#   HHL_Spn      1.244e-01  4.120e-03  30.200  < 2e-16 ***
#   HHL_Other    1.822e-01  4.059e-03  44.902  < 2e-16 ***
#   RMSP         5.749e-02  4.278e-04 134.405  < 2e-16 ***
#   YRBLT1940    6.930e-02  5.355e-03  12.942  < 2e-16 ***
#   YRBLT1950    1.432e-01  4.108e-03  34.857  < 2e-16 ***
#   YRBLT1960    1.614e-01  4.136e-03  39.019  < 2e-16 ***
#   YRBLT1970    1.719e-01  3.791e-03  45.357  < 2e-16 ***
#   YRBLT1980    2.543e-01  3.893e-03  65.332  < 2e-16 ***
#   YRBLT1990    3.059e-01  3.851e-03  79.426  < 2e-16 ***
#  YRBLT2000    3.611e-01  3.776e-03  95.629  < 2e-16 ***
#   YRBLT2010    4.227e-01  4.406e-03  95.923  < 2e-16 ***
#   YRBLT2020    4.825e-01  1.042e-02  46.320  < 2e-16 ***
#   YRBLT2021    4.772e-01  1.029e-02  46.394  < 2e-16 ***
#  YRBLT2022    5.194e-01  1.098e-02  47.292  < 2e-16 ***
#  YRBLT2023    5.329e-01  1.339e-02  39.791  < 2e-16 ***
#  YRBLT2024    5.264e-01  2.519e-02  20.894  < 2e-16 ***
#  BLD1U        1.332e+00  4.192e-03 317.698  < 2e-16 ***
#  BLDMU        1.465e+00  6.104e-03 240.095  < 2e-16 ***
#  AGEP         3.567e-03  6.607e-05  53.988  < 2e-16 ***
#  SEX         -3.794e-03  1.989e-03  -1.907  0.05646 .  
# MSP_sgl     -4.266e-02  2.968e-03 -14.376  < 2e-16 ***
#   RAC1P2      -1.903e-01  3.991e-03 -47.667  < 2e-16 ***
#   RAC1P3      -2.524e-01  1.084e-02 -23.271  < 2e-16 ***
#   RAC1P4      -4.899e-01  4.771e-02 -10.270  < 2e-16 ***
#   RAC1P5      -8.456e-02  2.985e-02  -2.833  0.00461 ** 
#   RAC1P6       3.162e-01  5.263e-03  60.076  < 2e-16 ***
#   RAC1P7       3.515e-01  3.199e-02  10.988  < 2e-16 ***
#   RAC1P8       1.102e-01  6.247e-03  17.649  < 2e-16 ***
#   RAC1P9       5.695e-02  3.983e-03  14.299  < 2e-16 ***
#   SCHL_HSCOMP  8.843e-02  4.318e-03  20.478  < 2e-16 ***
#   SCHL_COCOMP  3.374e-01  4.536e-03  74.384  < 2e-16 ***
#   DIS         -8.667e-02  2.732e-03 -31.723  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.9383 on 980220 degrees of freedom
# Multiple R-squared:  0.2882,	Adjusted R-squared:  0.2881 
# F-statistic: 1.202e+04 on 33 and 980220 DF,  p-value: < 2.2e-16





head(dat_tgt)

dat_tgt$eps_hat = fit.Y$residuals

# --- fit X
fit.X_fml = "logX ~ HHL_Spn + HHL_Other + RMSP + YRBLT1940 + YRBLT1950 + YRBLT1960 + YRBLT1970 + YRBLT1980 + YRBLT1990 + YRBLT2000 + YRBLT2010 + YRBLT2020 + YRBLT2021 + YRBLT2022 + YRBLT2023 + YRBLT2024 + BLD1U + BLDMU + AGEP + SEX + MSP_sgl + RAC1P2 + RAC1P3 + RAC1P4 + RAC1P5 + RAC1P6 + RAC1P7 + RAC1P8 + RAC1P9 + SCHL_HSCOMP + SCHL_COCOMP + DIS"
fit.X = lm(fit.X_fml, data = dat_tgt)

summary(fit.X)

# Residuals:
#  Min       1Q   Median       3Q      Max 
# -11.4089  -0.3921   0.0338   0.4514   3.9940 

# Coefficients:
#   Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)  1.020e+01  6.821e-03 1495.637  < 2e-16 ***
#   HHL_Spn      9.120e-03  3.549e-03    2.570   0.0102 *  
#   HHL_Other    5.857e-02  3.495e-03   16.756  < 2e-16 ***
#   RMSP         5.859e-02  3.637e-04  161.112  < 2e-16 ***
#  YRBLT1940   -3.536e-03  4.612e-03   -0.767   0.4433    
#  YRBLT1950    4.195e-02  3.538e-03   11.856  < 2e-16 ***
#  YRBLT1960    4.222e-02  3.562e-03   11.851  < 2e-16 ***
#  YRBLT1970    3.205e-02  3.265e-03    9.815  < 2e-16 ***
#  YRBLT1980    8.449e-02  3.352e-03   25.207  < 2e-16 ***
#  YRBLT1990    1.218e-01  3.315e-03   36.735  < 2e-16 ***
#  YRBLT2000    1.605e-01  3.249e-03   49.393  < 2e-16 ***
#  YRBLT2010    2.033e-01  3.790e-03   53.640  < 2e-16 ***
#  YRBLT2020    1.897e-01  8.971e-03   21.148  < 2e-16 ***
#  YRBLT2021    2.137e-01  8.857e-03   24.131  < 2e-16 ***
#  YRBLT2022    2.105e-01  9.459e-03   22.256  < 2e-16 ***
#  YRBLT2023    2.268e-01  1.153e-02   19.660  < 2e-16 ***
#  YRBLT2024    1.952e-01  2.170e-02    8.994  < 2e-16 ***
#  BLD1U        3.300e-01  3.596e-03   91.769  < 2e-16 ***
#  BLDMU        4.599e-01  5.237e-03   87.816  < 2e-16 ***
#  AGEP        -8.931e-03  5.619e-05 -158.938  < 2e-16 ***
#  SEX          2.925e-01  1.687e-03  173.359  < 2e-16 ***
#  MSP_sgl     -2.839e-01  2.540e-03 -111.767  < 2e-16 ***
#  RAC1P2      -1.534e-01  3.434e-03  -44.651  < 2e-16 ***
#  RAC1P3      -1.605e-01  9.340e-03  -17.182  < 2e-16 ***
#  RAC1P4      -2.860e-01  4.109e-02   -6.960 3.40e-12 ***
#  RAC1P5      -1.497e-01  2.571e-02   -5.824 5.76e-09 ***
#  RAC1P6       7.843e-02  4.533e-03   17.302  < 2e-16 ***
#  RAC1P7       7.157e-03  2.755e-02    0.260   0.7951    
#  RAC1P8      -5.999e-02  5.380e-03  -11.149  < 2e-16 ***
#  RAC1P9      -2.274e-02  3.430e-03   -6.628 3.40e-11 ***
#  SCHL_HSCOMP  3.015e-01  3.707e-03   81.319  < 2e-16 ***
#  SCHL_COCOMP  8.397e-01  3.814e-03  220.161  < 2e-16 ***
#  DIS         -2.575e-01  2.339e-03 -110.085  < 2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.8082 on 980221 degrees of freedom
# Multiple R-squared:  0.2708,	Adjusted R-squared:  0.2708 
# F-statistic: 1.138e+04 on 32 and 980221 DF,  p-value: < 2.2e-16

dat_tgt$Xperp = fit.X$residuals



sum(dat_tgt$eps_hat)
sum(dat_tgt$Xperp)


# --- Use linear mixed model to estimate the spatial correlation structure --- #
library(lme4)

dat_tgt <- dat_tgt %>%
  mutate(
    STATE_f  = factor(STATE),
    PUMA_f   = factor(PUMA),
    PUMA_full = interaction(STATE_f, PUMA_f, drop = TRUE) # I create this new column, because different PUMAs can have same number in different states.
  )

length(unique(dat_tgt$PUMA_full)) # 2462 PUMAs for the whole data

print(dat_tgt%>%group_by(STATE)%>% summarise(Number_of_PUMAs = length(unique(PUMA))), n=51)

str(dat_tgt)


cov.eps <- lmer(
  eps_hat ~ 0 + (1 | DIVISION) + (1 | STATE_f) + (1 | PUMA_full),
  data = dat_tgt,
  REML = FALSE # ML instead of REML. (ML means maximum likelihood. REML means restricted ML.)
)

VarCorr(cov.eps)   # the Std.Dev. column shows the random effect's standard deviations \sigma's, not variance sigma^2
# Groups    Name        Std.Dev.
# PUMA_full (Intercept) 0.18360 
# STATE_f   (Intercept) 0.19543 
# DIVISION  (Intercept) 0.25298 
# Residual              0.86189 

# Parameter       Point.Est
# \sigma_{P,eps}    0.18360
# \sigma_{S,eps}    0.19543
# \sigma_{D,eps}    0.25298
# \sigma_{e,eps}    0.86189

sigma(cov.eps)   # same as the last row in VarCorr

# profile likelihood
pr.cov.eps <- profile(cov.eps)
confint(pr.cov.eps)
#            2.5 %    97.5 %
# .sig01 0.1781696 0.1892735
# .sig02 0.1567992 0.2498636
# .sig03 0.1550931 0.4541844
# .sigma 0.8606802 0.8630962





cov.Xperp <- lmer(
  Xperp ~ 0 + (1 | DIVISION) + (1 | STATE_f) + (1 | PUMA_full),
  data = dat_tgt,
  REML = FALSE 
)


VarCorr(cov.Xperp)   # the Std.Dev. column shows the random effect's standard deviations \sigma's, not variance sigma^2
# Groups    Name        Std.Dev.
# PUMA_full (Intercept) 0.137140
# STATE_f   (Intercept) 0.092509
# DIVISION  (Intercept) 0.078655
# Residual              0.789353

# Parameter     Point.Est
# \sigma_{P,X}   0.137140
# \sigma_{S,X}   0.092509
# \sigma_{D,X}   0.078655
# \sigma_{e,X}   0.789353 

sigma(cov.Xperp)  

# profile likelihood
pr.cov.Xperp <- profile(cov.Xperp)
confint(pr.cov.Xperp)
#            2.5 %    97.5 %
# .sig01 0.13293210 0.1415388
# .sig02 0.07187226 0.1210694
# .sig03 0.03767062 0.1503290
# .sigma 0.78824815 0.7904609


save(list = c("cov.eps", "cov.Xperp", "dat_tgt",
              "fit.X", "fit.Y", "pr.cov.eps", "pr.cov.Xperp", 
              "fit.X_fml", "fit.Y_fml"),
     file = file.path(folder, "1222_PUMA_setting.Rdata"))


