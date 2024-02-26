rm(list =ls())
library(ggplot2) 

###### STEP 1: load the data ############

filename= r"(clinicaldata.xlsx)" #MODIFY to your file
Sheetnaam="sheet1"
DATA=readxl::read_xlsx(path = filename, sheet = Sheetnaam)

###### STEP 2: exclude patients (i.e. rows) ####

# exclude patients that are missing endpoint 
DATA=DATA[!is.na(DATA$endpoint),]

# exclude patients based on an exclusion criteria
DATA=DATA[DATA$NSCLC_vs_SCLC!="Thymoma"|is.na(DATA$NSCLC_vs_SCLC),] #MODIFY, in this example--> 'NSCLC_vs_SCLC' is the variable name and 'Thymoma' is the value you want to exclude
DATA=DATA[(DATA$Flowchart %in%c("inlc")),]

((DATA[duplicated(DATA$UMCGnr),]))

###### STEP 3: exclude variables, different ways: ######

# - delete columns that hold no or a single value
count_unique=apply(DATA,2,function(x) length(unique(x[!is.na(x)]))   )
print(paste("Deleted the follow columns:",paste(colnames(DATA)[count_unique<2], collapse = ', ')))
DATA=DATA[,count_unique>1]

# - delete specific variables
# I'm going to remove all of the cols other than the once I specify
DATA = DATA[c('Gender', 'Age_StartRT', 'Upper_Middle_Lower', 
             'Comorbidity_Pulmonary', 'Smoking_status','Chemo_Sequence',
             'Pneumonitis_BeforeM6', 'LungMINGTV_doseGEM')]

###### Step 4: Transforming categoral variables ######
# Binary ones
DATA$y <- ifelse(DATA$Pneumonitis_BeforeM6 %in% c('No pneumonitis', 'Grade 1'), 0, 1)

DATA$old_age <- ifelse(DATA$Age_StartRT > 63, "Yes", "No")

DATA$old_age = factor(DATA$old_age,
                      levels = c("No", "Yes"),
                      labels = c("No", "Yes"))

DATA$smoking_prevously = factor(DATA$Smoking_status, 
                                levels = c("Former smoker < 3 months", "Former smoker >= 3 months", "Never", "Current"),
                                labels = c("Yes"                     , "Yes"                      , "No", "No"))
DATA$smoking_currently = factor(DATA$Smoking_status, 
                                levels = c("Current", "Never", "Former smoker >= 3 months", "Former smoker < 3 months"),
                                labels = c("Yes", "No", "No", "No"))

DATA$smoking_cat = factor(DATA$Smoking_status, 
                    levels = c("Never", "Former smoker >= 3 months", "Former smoker < 3 months",  "Current"),
                    labels = c("Never", "Former", "Former", "Current"))


DATA$smoking_3months = factor(DATA$Smoking_status,
                              levels = c("Never"                    ,"Former smoker >= 3 months","Former smoker < 3 months", "Current"  ),
                              labels = c("Never or > 3mo non smoker", "Never or > 3mo non smoker","current or < 3mo non smoker", "current or < 3mo non smoker" ))

DATA$tumor_location = factor(DATA$Upper_Middle_Lower,
                             levels = c("Middle/hilar", "Upper", "Bronchial","Lower"),
                             labels = c("Middle or Upper", "Middle or Upper", "Middle or Upper", "Lower")) # 1 == lower
DATA$Comorbidity_Pulmonary = factor(DATA$Comorbidity_Pulmonary,
                                    levels = c("No", "Yes"), labels = c("No", "Yes"))
DATA$Sequential_chemo = factor(DATA$Chemo_Sequence,
                               levels = c("Sequential (neo-adjuvant)", "Induction and concurrent", "Adjuvant", "No chemotherapy", "Concurrent: weekly chemotherapy", "Induction only"),
                               labels = c("Yes", "No", "No", "No", "No", "No"))
DATA$Gender = factor(DATA$Gender,
                     levels = c("Male", "Female"),
                     labels = c("Male", "Female"))

# Clean up columns that are now unnecessary\
DATA = DATA[, !(colnames(DATA) %in% c("Upper_Middle_Lower", "Smoking_status", "Chemo_Sequence", "Pneumonitis_BeforeM6"))]

###### STEP 5: Include radiomics and dose variables ####

# Radiomics
rad_path = r'(.\radiomcs.xlsx)'
rad_sheetname = 'Sheet1'
RADIOMICS = readxl::read_xlsx(path = rad_path, sheet = rad_sheetname)
colnames(RADIOMICS)[1] <- 'p_id'

# selecting only the patients that are in RADIOMICS
DATA = DATA[DATA$p_id %in% DATA$p_id, ]
allDATA = merge(DATA, RADIOMICS, by='p_id')


###### STEP 7: Save your DATA to file #####
# Save data
save(allDATA, file=r"(.\data.R)")

#--> now you can finally... continue to script step 2


