# Reproducible analysis for the 2018 Nigeria Demographic and Health Survey
# Cluster-aware Multistage Selection (CMS) analysis.

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

RUN_MODE <- "full"
RANDOM_SEED <- 20251129L

# For journal replication, run from a project directory containing a data/
# folder, or set NDHS_PROJECT_DIR, NDHS_INPUT_DIR, and NDHS_OUTPUT_DIR.
PROJECT_DIR <- normalizePath(
  Sys.getenv("NDHS_PROJECT_DIR", unset = "."),
  mustWork = TRUE
)
INPUT_DIR <- Sys.getenv(
  "NDHS_INPUT_DIR",
  unset = file.path(PROJECT_DIR, "data")
)
OUTPUT_DIR <- Sys.getenv(
  "NDHS_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "outputs")
)

RUN_SETTINGS <- list(
  step1_B = 200L,
  step2_B = 200L,
  cv_folds = 10L,
  full_tuning_grid = TRUE
)

FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
RESULT_DIR <- file.path(OUTPUT_DIR, "results")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "dplyr", "ggplot2", "glmnet", "grid", "gridExtra", "haven",
  "labelled", "naniar", "openxlsx", "patchwork", "purrr",
  "sjlabelled", "splines", "stringr", "survey", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the following required R packages before running the analysis: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmnet)
  library(grid)
  library(gridExtra)
  library(haven)
  library(labelled)
  library(naniar)
  library(openxlsx)
  library(patchwork)
  library(purrr)
  library(sjlabelled)
  library(splines)
  library(stringr)
  library(survey)
  library(tibble)
  library(tidyr)
})

set.seed(RANDOM_SEED)

########################################
## Data reading and preparing dataset ##
########################################

### Read data
input_files <- c(
  PR = "NGPR7BFL.DTA",
  KR = "NGKR7BFL.DTA",
  HR = "NGHR7BFL.DTA"
)
input_paths <- setNames(
  file.path(INPUT_DIR, unname(input_files)),
  names(input_files)
)
missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs)) {
  stop(
    "Required DHS recode file(s) not found: ",
    paste(missing_inputs, collapse = ", ")
  )
}

PRdata <- read_dta(input_paths[["PR"]])
KRdata <- read_dta(input_paths[["KR"]])
HRdata <- read_dta(input_paths[["HR"]])

### Merging data
KRdata <- KRdata %>%
  left_join(PRdata %>% 
              mutate(v001 = hv001,
                     v002 = hv002,
                     b16  = hvidx),
            by = c("v001", "v002", "b16"))

# Malaria status (by RDT or microscopy) in children 6-59 months
KRdata <- KRdata %>%
  mutate(ml_diagnosis = case_when(
    hc1>=6 & hc1<=59 & (hml35!=1 & hml32!=1) ~ 0,
    (hml35==1 | hml32==1) & hc1>=6 & hc1<=59 ~ 1),
    ml_diagnosis = set_label(ml_diagnosis, label = "Parasitemia (by RDT or microscopy) in children 6-59 months"))
# Filter dataset for those with malaria testing results only
KRdata <- filter(KRdata, !is.na(ml_diagnosis))

HRtemp <- HRdata %>%
  transmute(
    v001 = hv001,
    v002 = hv002,
    hml10_1, hml10_2, hml10_3, hml10_4, hml10_5, hml10_6, hml10_7
  )

KRdata <- KRdata %>%
  left_join(HRtemp, by = c("v001", "v002"))

### Select children of age 6-59 months
KRdata <- KRdata %>%
  mutate(age = hc1)
KRdata <- filter(KRdata, age>=6, age<=59)

# Create create age cutoff variable of 24 months
KRdata <- KRdata %>%
  mutate(age_24 = case_when(
    age < 24 ~ 0,
    age >=24 ~ 1
  ))

########################################
### Categorization for each variable ###
########################################

#Is the child female? (group by sex)
KRdata <- KRdata %>%
  mutate(sex_female = case_when(
    hv104==1 ~ 0,
    hv104==2 ~ 1)) %>%
  set_value_labels(sex_female= c("Male" = 0, "Female"=1)) %>%
  set_variable_labels(sex_female= "Sex of child")
# 0 = male, 1= female

####Malaria
# Anemia in children 6-59 months
KRdata <- KRdata %>%
  mutate(ml_anemia = case_when(
    hv103==1 & hc1>=6 & hc1<=59 & hc55==0 & hv042==1 & !(hc56<80)  ~ 0,
    hv103==1 & hc1>=6 & hc1<=59 & hc55==0 & hv042==1 & hc56<80  ~ 1),
    ml_anemia = set_label(ml_anemia, label = "Anemia in children 6-59 months"))
# 0 = > 8 g/dl or inconsistent, 1 = < 8 g/dl (anemia)

#### MDD MAD MMF
# *** Foods consumed ***
KRdata <- KRdata %>%
  # country specific foods. These can be added to the foods below based on the survey. See example for nt_root & nt_meatfish below
  mutate(food1  = case_when(v414a==1  ~ 1 , v414a!=1 ~ 0)) %>%
  mutate(food2  = case_when(v414b==1  ~ 1 , v414b!=1 ~ 0)) %>%
  mutate(food3  = case_when(v414c==1  ~ 1 , v414c!=1 ~ 0)) %>%
  mutate(food4  = case_when(v414d==1  ~ 1 , v414d!=1 ~ 0)) %>%
  mutate(nt_formula  = case_when(v411a==1  ~ 1 , v411a!=1~ 0)) %>% # Given formula
  mutate(nt_milk  = case_when(v411==1  ~ 1 , v411!=1~ 0)) %>% # Given other milk
  mutate(nt_liquids= case_when(v410==1 | v412c==1 | v413==1  ~ 1 , v410!=1 | v412c!=1 | v413!=1  ~ 0)) %>% # Given other liquids
  mutate(nt_bbyfood  = case_when(v412a==1  ~ 1 , v412a!=1~ 0)) %>% # Given fortified baby food
  mutate(nt_grains  = case_when(v412a==1 | v414e==1 ~ 1 , v412a!=1 | v414e!=1 ~ 0)) %>% # Given grains
  mutate(nt_vita = case_when(v414i==1 | v414j==1 | v414k==1 ~ 1 , v414i!=1 | v414j!=1 | v414k!=1 ~ 0)) %>% # Given Vit A rich foods
  mutate(nt_frtveg  = case_when(v414l==1  ~ 1 , v414l!=1~ 0)) %>% # Given other fruits or vegetables
  mutate(nt_root  = case_when(   # Given roots and tubers  
    (v000 == "UG7" & (v414f==1 | food1==1)) | (v000 != "UG7" & v414f==1) ~ 1, 
    (v000 == "UG7" & (v414f!=1 | food1!=1)) | (v000 != "UG7" & v414f!=1) ~ 0)) %>%
  mutate(nt_nuts  = case_when(v414o==1  ~ 1 , v414o!=1~ 0)) %>% # Given nuts or legumes
  mutate(nt_meatfish  = case_when(   # Given meat, fish, shellfish, or poultry  
    (v000 == "UG7" & (v414h==1 |v414m==1 |v414n==1| food2==1)) | (v000 != "UG7" & (v414h==1 | v414m==1 | v414n==1)) ~ 1, 
    (v000 == "UG7" & !(v414h==1 |v414m==1 |v414n==1| food2==1)) | (v000 != "UG7" & !(v414h==1 | v414m==1 | v414n==1)) ~ 0)) %>%
  mutate(nt_eggs  = case_when(v414g==1  ~ 1 , v414g!=1~ 0)) %>% # Given eggs
  mutate(nt_dairy  = case_when(v414p==1 | v414v==1 ~ 1 , v414p!=1 | v414v!=1 ~ 0)) %>% # Given dairy
  mutate(nt_solids = case_when( nt_bbyfood==1 | nt_grains==1 | nt_vita==1 | nt_frtveg==1 | nt_root==1 | nt_nuts==1 | nt_meatfish==1 | 
                                  nt_eggs==1 | nt_dairy==1 | v414s==1 ~ 1 ,
                                nt_bbyfood!=1 | nt_grains!=1 | nt_vita!=1 | nt_frtveg!=1 | nt_root!=1 | nt_nuts!=1 | nt_meatfish!=1 | 
                                  nt_eggs!=1 | nt_dairy!=1 | v414s!=1 ~ 0) ) %>%
  #add labels
  set_value_labels(nt_formula = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_formula = "Child given infant formula in day/night before survey - last-born under 2 years") %>%
  set_value_labels(nt_milk = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_milk = "Child given other milk in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_liquids = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_liquids = "Child given other liquids in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_bbyfood = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_bbyfood = "Child given fortified baby food in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_grains = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_grains = "Child given grains in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_vita = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_vita = "Child given vitamin A rich food in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_frtveg = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_frtveg = "Child given other fruits or vegetables in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_root = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_root = "Child given roots or tubers in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_nuts = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_nuts = "Child given legumes or nuts in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_meatfish = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_meatfish = "Child given meat, fish, shellfish, or poultry in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_eggs = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_eggs = "Child given eggs in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_dairy = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_dairy = "Child given cheese, yogurt, or other milk products in day/night before survey- last-born under 2 years") %>%
  set_value_labels(nt_solids = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_solids = "Child given any solid or semisolid food in day/night before survey- last-born under 2 years") 

# *** Minimum feeding indicators ***
# Totmilkf = milk + formula + yogurt
KRdata$v469e <- as.numeric(KRdata$v469e)
KRdata$v469f <- as.numeric(KRdata$v469f)
KRdata$v469x <- as.numeric(KRdata$v469x)
KRdata <- KRdata %>%
  mutate(v469e_ = v469e) %>%
  replace_with_na(replace = list(v469e_ = c(8)))%>%
  mutate(v469f_ = v469f) %>%
  replace_with_na(replace = list(v469f_ = c(8)))%>%
  mutate(v469x_ = v469x) %>%
  replace_with_na(replace = list(v469x_ = c(8)))%>%
  mutate(v469e_ = coalesce(v469e_, 0),
         v469f_ = coalesce(v469f_, 0),
         v469x_ = coalesce(v469x_, 0))%>%
  mutate(totmilkf= v469e_ + v469f_ + v469x_) %>%
  mutate(nt_fed_milk  = 
           case_when(totmilkf>=2 | m4==95 & between(age,6,23) ~ 1 , 
                     totmilkf <2 | m4!=95 & between(age,6,23) ~ 0)) %>%
  set_value_labels(nt_fed_milk = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_fed_milk = "Child given milk or milk products- last-born 6-23 months")


# //Min dietary diversity
KRdata <- KRdata %>%
  # 1. breastmilk
  mutate(group1 = case_when(m4==95  ~ 1 ,  m4!=95 ~ 0)) %>% 
  #2. infant formula, milk other than breast milk, cheese or yogurt or other milk products
  mutate(group2 = case_when(nt_formula==1 | nt_milk==1 | nt_dairy==1  ~ 1 , nt_formula!=1 | nt_milk!=1 | nt_dairy!=1 ~ 0)) %>%
  #3. foods made from grains, roots, tubers, and bananas/plantains, including porridge and fortified baby food from grains
  mutate(group3  = case_when(nt_grains==1 | nt_root==1 | nt_bbyfood==1 ~ 1 , nt_grains!=1 | nt_root!=1 | nt_bbyfood!=1 ~ 0)) %>%
  #4. vitamin A-rich fruits and vegetables
  mutate(group4  = case_when(nt_vita==1  ~ 1 , nt_vita!=1 ~ 0)) %>%
  #5. other fruits and vegetables
  mutate(group5  = case_when(nt_frtveg==1 ~ 1 , nt_frtveg!=1~ 0)) %>% 
  #6. eggs
  mutate(group6  = case_when(nt_eggs==1 ~ 1 , nt_eggs!=1~ 0)) %>% 
  #7. meat, poultry, fish, and shellfish (and organ meats)
  mutate(group7  = case_when(nt_meatfish==1 ~ 1 , nt_meatfish!=1~ 0)) %>% 
  #8. legumes and nuts
  mutate(group8  = case_when(nt_nuts==1 ~ 1 , nt_nuts!=1~ 0)) %>% 
  #add the food groups
  mutate(foodsum  = group1+group2+group3+group4+group5+group6+group7+group8) %>% 
  mutate(nt_mdd  = case_when(between(age,6,23) & foodsum<5 ~ 0 , between(age,6,23) & foodsum>=5~ 1)) %>% 
  #older surveys are 4 out of 7 food groups so the foodsum would add group2-group8 and the recode the sum for 4+ as yes
  set_value_labels(nt_mdd = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_mdd = "Child with minimum dietary diversity, 5 out of 8 food groups- last-born 6-23 months")
# 0 = No; children without MDD, 1 = Yes; children with MDD

# //Min meal frequency
KRdata <- KRdata %>%
  mutate(feedings = case_when(m39>0 & m39<8 ~ totmilkf + m39 )) %>%
  mutate(nt_mmf  = 
           if_else(between(age,6,23) & (m4==95 & between(m39,2,7) & between(age,6,8)) | (m4==95 & between(m39,3,7) & between(age,9,23)) |
                     (m4!=95 & feedings>=4 & between(age,6,23)), 1, 0 )) 
KRdata[["nt_mmf"]] <- ifelse(KRdata[["age"]]<6, NA, KRdata[["nt_mmf"]])
KRdata <- KRdata %>%
  set_value_labels(nt_mmf = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_mmf = "Child with minimum meal frequency- last-born 6-23 months")
# 0 = No; children without MMF, 1 = Yes; children with MMF

# //Min acceptable diet
KRdata <- KRdata %>%
  mutate(foodsum2 = nt_grains+nt_root+nt_nuts+nt_meatfish+nt_vita+nt_frtveg+nt_eggs) %>%
  mutate(nt_mad  = 
           if_else((m4==95 & nt_mdd==1 & nt_mmf==1) | (m4!=95 & foodsum2>=4 & nt_mmf==1 & totmilkf>=2), 1, 0 )) 
KRdata[["nt_mad"]] <- ifelse(KRdata[["age"]]<6, NA, KRdata[["nt_mad"]])
KRdata <- KRdata %>%
  set_value_labels(nt_mad = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_mad = "Child with minimum acceptable diet- last-born 6-23 months")
# 0 = No; children without MAD, 1 = Yes; children with MAD

#Iron supplement in the seven days preceding the interview
KRdata <- KRdata %>%
  mutate(iron= case_when(
    h42==1 ~ 1,
    !(h42==1) ~0 ))
#1 = Yes, 0 = No

#Vitamin A supplement in the six months preceding the interview
KRdata <- KRdata %>%
  mutate(vit_a= case_when(
    (v008a<7 | h34 ==1) ~ 1,
    !(v008a<7 | h34 ==1) ~0 ))
#1 = Yes, 0 = No

#Anemia level in mother
KRdata <- KRdata %>%
  mutate(mother_anemia= case_when(
    v457==1 ~ 1,
    v457==2 ~ 2,
    v457==3 ~ 3,
    v457==4 | v457==9~ 4))
#1 = severe anemia, 2 = moderate anemia, 3= mild anemia, 4 = no anemia

#BMI status in mother 
KRdata <- KRdata %>%
  mutate(mother_bmi= case_when(
    v445 <=1699 ~ 1,
    v445 <= 1849 & v445 >= 1700 ~ 2,
    v445 >= 1850 & v445 <= 2499 ~ 3,
    v445 >= 2500 & v445 <= 2999 ~ 4,
    v445 >= 3000 ~ 5 ))
#1 = Very thin (<17), 2 = thin (17-18.49), 3 = normal (18.5-24.9), 4 = overweight (25 - 29.9), 5 = obese (>=30)

# //Currently breastfed
KRdata <- KRdata %>%
  mutate(nt_bf_curr =
           case_when(
             m4==95  ~ 1 ,
             m4 %in% c(93,94,98,99) ~ 0)) %>%
  set_value_labels(nt_bf_curr = c("Yes" = 1, "No"=0  )) %>%
  set_variable_labels(nt_bf_curr = "Currently breastfeeding - last-born under 2 years")
# 0 = ever breastfed (not currently breastfeeding) + never breastfed + don't know, 1 = currently (still) breastfeeding

# //Height-for-age
KRdata <- KRdata %>%
  mutate(nt_ch_ht_age =
           case_when(
             hv103==1 &  hc70 < -200  ~ 1 ,
             hv103==1 &  hc70 >= -200 ~ 0,
             hc70>=9996 ~ 99)) %>%
  replace_with_na(replace = list(nt_ch_ht_age = c(99)))
# 1 = stunting, 0 = normal


# //Weight-for-height
KRdata <- KRdata %>%
  mutate(nt_ch_wt_ht =
           case_when(
             hv103==1 &  hc72< -200  ~ 2 ,
             hv103==1 &  hc72 <= 200 & hc72>= -200 ~ 1 ,
             hv103==1 &  hc72 <= 9990 & hc72>= 200 ~ 0 ,
             hc71>=9996 ~ 99)) %>%
  replace_with_na(replace = list(nt_ch_wt_ht = c(99)))
# 2 = wasting, 1 = normal, 0 = overweight

# //Weight-for-age
KRdata <- KRdata %>%
  mutate(nt_ch_wt_age =
           case_when(
             hv103==1 &  hc71 < -200 ~ 2, 
             hv103==1 &  hc71>= -200 & hc71 <= 200 ~ 1 ,
             hv103==1 &  hc71> 200 & hc71<9990 ~ 0 ,
             hc72>=9996 ~ 99)) %>%
  replace_with_na(replace = list(nt_ch_wt_age = c(99)))
# 2 = underweight for age, 1 = normal, 0 = overweight for age

# //Have electricity
KRdata <- KRdata %>% 
  mutate(ph_electric = hv206) %>%
  set_value_labels(ph_electric = c("Yes" = 1, "No"=0)) %>%
  set_variable_labels(ph_electric = "Have electricity") 
# 0 = No, 1 = Yes

# //Finished Flooring material
KRdata <- KRdata %>%
  mutate(ph_floor = hv213) %>%
  set_variable_labels(ph_floor = "Flooring material") %>%
  mutate(ph_floor = case_when(ph_floor %in% 10:22 ~ 0, ph_floor %in% 30:35 ~ 1))
# 0 = Natural, Rudimentary, 1 = Finished

# //2 or more rooms for sleeping
KRdata <- KRdata %>%
  mutate(ph_rooms_sleep =
           case_when(
             hv216==1  ~ 0 ,  
             hv216>=2  ~ 1 ,
             hv216==0 | is.na(hv216) ~ 9)) %>%
  set_value_labels(ph_rooms_sleep = c( "Missing"=9, "Two or more" =1, "One" =0)) %>%
  set_variable_labels(ph_rooms_sleep = "Rooms for sleeping")
# 1 = 1 room, 2 = 2 or more rooms

# //Cook inside
KRdata <- KRdata %>%
  mutate(ph_cook_place =
           case_when(
             hv241<9 & hv241!=3 & hv241!=6 & hv241!=2  ~ hv241 ,
             hv241==3 ~ 0,
             hv241==6 ~ 1,
             hv241==2 ~ 1,
             hv226==95  ~ 1,
             hv241>=9 | is.na(hv241) ~ 9)) %>%
  set_value_labels(ph_cook_place = c("Missing"=9, "Other" =6, "No food cooked in household" =4, "Outdoors" =0, "Inside a building"=1)) %>%
  set_variable_labels(ph_cook_place = "Place for cooking")
# 1 = In the house, inside a seperate building, no food cooked in house, 0= outside

# //Own agricultural land
KRdata <- KRdata %>%
  mutate(ph_agriland = 
           case_when(
             hv244==0 | is.na(hv244) ~ 0,
             hv244==1 ~ 1)) %>%
  set_value_labels(ph_agriland = c("Yes" = 1, "No"=0)) %>%
  set_variable_labels(ph_agriland = "Owns agricultural land")
# 1 = Yes, 0 = No

# //Own livestock
KRdata <- KRdata %>%
  mutate(ph_animals = 
           case_when(
             hv246==0 | is.na(hv246) ~ 0,
             hv246==1 ~ 1)) %>%
  set_value_labels(ph_animals = c("Yes" = 1, "No"=0)) %>%
  set_variable_labels(ph_animals = "Owns livestock or farm animals")
# 1 = Yes, 0 = No

#Mother used internet at least once a week
KRdata <- KRdata %>%
  mutate(internet_fq = case_when(
    v171b==2 | v171b ==3 ~ 1,
    !(v171b==2 | v171b ==3) ~ 0))
#1 = Yes, 0 = No

#Mother exposed to radio or news paper or TV at least once a week (exposure to mass media)
KRdata <- KRdata %>%
  mutate(mass_media = case_when(
    (v157==2 | v157==3) | (v158 == 2|v158==3) | (v159==2 | v159==3) ~ 1,
    !(v157==2 | v157==3) | (v158 == 2|v158==3) | (v159==2 | v159==3) ~ 0))
#1 = Yes, 0 = No

#Access to improved sanitation facility
KRdata <- KRdata %>%
  mutate(sanitation_imp = case_when(
    hv205==11 | hv205==12 | hv205==13 | hv205==15 | hv205==21 | hv205 ==22 | hv205 == 41 ~ 1,
    !(hv205==11 | hv205==12 | hv205==13 | hv205==15 | hv205==21 | hv205 ==22 | hv205 == 41) ~0 ))
#1 = Yes, 0 = No

#Place of handwashing
KRdata <- KRdata %>%
  mutate(handwash_place = case_when(
    hv230a == 1 ~ 1,
    hv230a == 2 ~ 2,
    !hv230a %in% c(1,2) ~ 3
  ))
#1 = Fixed, 2= Mobile, 3= Not observed/Not in dwelling

#Presence of water & soap at handwashing place
KRdata <- KRdata %>%
  mutate(soap_water = case_when(
    hv230b == 1 & hv232 == 1 ~ 1,
    TRUE ~ 0
  ))
#1 = Yes, 0 = No

#Presence of ash, mud, sand in handwashing place
KRdata <- KRdata %>%
  mutate(handwash_dirty = case_when(
    hv232b == 1 ~ 1,
    TRUE ~ 0
  ))
#1 = Yes, 0 = No

#Own source for water
KRdata <- KRdata %>%
  mutate(own_water = case_when(
    hv235 == 1 | hv235 == 2 ~ 1,
    TRUE ~ 0
  ))
#1 = In own dwelling/yard, 0 = elsewhere

#Shared toilet with other households
KRdata <- KRdata %>%
  mutate(share_toilet = case_when(
    v160 == 0 ~ 0,
    v160 %in% c(1, 7) ~ 1,
    TRUE ~ 9
  ))
#1 = Yes, 0 = No, 9= NA

#Place of delivery
KRdata <- KRdata %>%
  mutate(delivery_place = case_when(
    m15 %in% c(10, 11, 12, 96) ~ 1,
    m15 %in% c(20, 21, 22, 23, 26) ~ 2,
    m15 %in% c(30, 31, 36) ~ 3
  ))
#1 = Own or other's homes, #2= Public/Government facilities, #3 = Private facilities

#Took drugs for intestinal parasites during pregnancy
KRdata <- KRdata %>%
  mutate(preg_giparasites = case_when (
    m60 == 1 ~ 1,
    TRUE ~ 0
  ))
#1 = Yes, 0 = No

#Disposal of youngest child's stool when not using toilet
KRdata <- KRdata %>%
  mutate(stool_disposal = case_when(
    v465 %in% c(1,2) ~ 1,
    TRUE ~ 0
  ))
#1 = disposed appropriately in toilet/latrines, 0 = buried, left in open, rinsed into drain or ditch, or thrown into garbage

#Access to improved water source
KRdata <- KRdata %>%
  mutate(water_imp = case_when(
    hv201==11 | hv201==14 | hv201==21 | hv201==31 | hv201==41 | hv201==51 | hv201==61 | hv201==62 | hv201==71~ 1,
    !(hv201==11 | hv201==14 | hv201==21 | hv201==31 | hv201==41 | hv201==51 | hv201==61 | hv201==62 | hv201==71) ~0 ))
#1 = Yes, 0 = No

#Number of household members
KRdata <- KRdata %>%
  mutate(hh_members = case_when(
    hv009>=2 & hv009 <=4 ~ 1,
    hv009>=5 & hv009 <=7 ~ 2,
    hv009>=8 ~ 3))
#1 = 2~4 members, 2 = 5~7 members, 3 = 8 or more members

# *** Wealth quintile ***
KRdata <- KRdata %>%
  mutate(ph_wealth_quint =
           case_when(hv102==1 ~ hv270)) %>%
  set_variable_labels(ph_wealth_quint = "Wealth quintile - dejure population")
# 1 = poorest (low 20%), 2 = poorer, 3 = middle, 4 = richer, 5 = richest (high 20%)

# Mother education level
KRdata <- KRdata %>%
  mutate(mother_eduyr=case_when(hc61 == 0 | hc61==8 | hc61==9 ~ 0,
                                hc61 == 1 ~ 1, 
                                hc61 == 2 ~ 2,
                                hc61 == 3 ~ 3))
#0 = no education, 1 = primary, 2 = secondary, 3 = higher

# Mother had primary or higher education
KRdata <- KRdata %>%
  mutate(mother_eduyr_primary=case_when( hc1 >=6 & hc1 <= 59 & (hc61 == 0 | hc61==8 | hc61==9 )~ 0,
                                         hc1 >=6 & hc1 <= 59 & hc61 >=1 & hc61 <=3 ~ 1))
#1 = primary, secondary, higher, 0 = no education

#Father/Partner education level
KRdata <- KRdata %>%
  mutate(partner_eduyr= case_when(
    v701==0 | v701==8 | v701== 9 ~ 0,
    v701==1 ~ 1,
    v701==2 ~ 2,
    v701==3 ~ 3))
#0 = no education, 1 = primary, 2 = secondary, 3 = higher

#Mother's occupation
KRdata <- KRdata %>%
  mutate(mother_occupation= case_when(
    v717==0 | v717==98 | v717==99 | v717==6 ~ 0,
    v717==1 | v717==2 | v717 ==3 | v717==6 | v717 == 7~ 1,
    v717==4 | v717==5 | v717 ==8 | v717==9 ~ 2))
#0 = no occupation or household & domestic, 1= White-collar; professional/managerial/clerical/sales/services, 2 = Blue collar; agricultural/skilled & unskilled manual

#Husband/partner's occupation
KRdata <- KRdata %>%
  mutate(partner_occupation= case_when(
    v705==0 | v705==98 | v705==99 | v705==6 ~ 0,
    v705==1 | v705==2 | v705 ==3 | v705 == 7~ 1,
    v705==4 | v705==5 | v705 ==8 | v705==9 ~ 2))
#0 = no occupation or household & domestic, 1= White-collar; professional/managerial/clerical/sales/services, 2 = Blue collar; agricultural/skilled & unskilled manual

#Household owns ITN
KRdata <- KRdata %>%
  mutate(ml_itnhh = case_when(
    hml10_1==1 | hml10_2==1 | hml10_3==1 | hml10_4==1 | hml10_5==1 | hml10_6==1 | hml10_7 ==1   ~ 1,
    TRUE ~ 0),
    ml_itnhh = add_labels(ml_itnhh, labels = c("No"=0, "Yes"=1)),
    ml_itnhh = set_label(ml_itnhh, label = "Household owns at least one ITN"))
#1 = Yes, 0 = No

# Number of ITNs per household
KRdata <- KRdata %>%
  mutate(itnhh_01 = case_when(hml10_1==1 ~ 1,TRUE ~ 0)) %>%
  mutate(itnhh_02 = case_when(hml10_2==1 ~ 1,TRUE ~ 0)) %>%
  mutate(itnhh_03 = case_when(hml10_3==1 ~ 1,TRUE ~ 0)) %>%
  mutate(itnhh_04 = case_when(hml10_4==1 ~ 1,TRUE ~ 0)) %>%
  mutate(itnhh_05 = case_when(hml10_5==1 ~ 1,TRUE ~ 0)) %>%
  mutate(itnhh_06 = case_when(hml10_6==1 ~ 1,TRUE ~ 0)) %>%
  mutate(itnhh_07 = case_when(hml10_7==1 ~ 1,TRUE ~ 0)) %>%
  mutate(ml_numitnhh = itnhh_01 + itnhh_02 + itnhh_03 + itnhh_04 + itnhh_05 + itnhh_06 + itnhh_07,
         ml_numitnhh = set_label(ml_numitnhh, label = "Number of ITNs per household"))

# Household owns > 1 ITN per 2 members
KRdata <- KRdata %>%
  mutate(ml_potuse = ml_numitnhh*2,
         ml_potuse = set_label(ml_potuse, label = "Potential ITN users in household"))

KRdata <- KRdata %>%
  mutate(ml_hhaccess0 =ml_potuse/hv013) %>%
  mutate(ml_hhaccess = case_when(
    hv013==0 ~ 99,
    ml_hhaccess0 >= 1   ~ 1,
    TRUE   ~ 0),
    ml_hhaccess = set_label(ml_hhaccess, label = "Households with >1 ITN per 2 household members"))%>%
  replace_with_na(replace = list(ml_hhaccess = c(99)))
#1 = Yes, 0 = No

# Slept under an ITN 
KRdata <- KRdata %>%
  mutate(ml_netcat = case_when(
    hml12==0  ~ 0,
    hml12==1|hml12==2|hml12==3  ~ 1,),
    ml_netcat = set_label(ml_netcat, label = "Mosquito net categorization"))
#0= did not sleep under an ITN net, 1 = slept under an ITN net

# Geography classifications by regions of Nigeria
KRdata <- KRdata %>%
  mutate(geography_regions=case_when(
    hv024==1 ~ 1,
    hv024==2 ~ 2,
    hv024==3 ~ 3,
    hv024==4 | hv024==5 ~ 4,
    hv024==6 ~ 5)) %>%
  set_value_labels(geography_regions = c("North central"=1, "North east"=2, "North west "=3, "South east & South south"=4, "South west"=5)) %>%
  set_variable_labels(geography_regions = "Geography classifications by regions of Nigeria")
#1 = North Central, 2 = North East, 3 = North West, 4 = South East & South South, 6 = South West

# Lives in rural regions
KRdata <- KRdata %>%
  mutate(geography_rural=
           case_when(
             hv025 == 2 ~ 1,
             !(hv025==2) ~ 0)) %>%
  set_value_labels(geography_rural = c("Lives in rural regions"=1, "Does not live in rural regions"=0)) %>%
  set_variable_labels(geography_rural = "Lives in rural egions")
# 1 = Yes, 0 = No

#Visited health facility last 12 months
KRdata <- KRdata %>%
  mutate(visit_health = case_when(
    v394==1~ 1,
    !(v394==1) ~0 ))
#1 = Yes, 0 = No

#Mother has major problems in accessing healthcare
KRdata <- KRdata %>%
  mutate(access_health = case_when(
    v467a == 1 | v467b== 1 | v467c ==1 | v467d ==1 | v467f ==1 ~ 1,
    TRUE ~ 0))
#1 = Yes either major problems in know where to go/financial problems/distance problems/getting permission/not wanting to go alone, 0 = No

########################################
####### Handling missing values ########
########################################

###Replace with 1
#Wealth quintile = 127/8447 -> all into lowest wealth quintile
KRdata <- KRdata %>% 
  mutate_at(c('ph_wealth_quint'), ~replace_na(.,1))

###Replace with 0
#Anemia in child = 2/8447 -> all in to no anemia
#Partner education = 421/8447 -> all into no education
#Mother occupation = 3736/8447 -> all into no occupation
#Partner occupation = 445/8447 -> all into no occupation
#MDD = 5537/8447 -> all into no as they are above 2
#MMF = 107/8447 -> all into no
#MAD = 5/8447 -> all into no
#Stunting = 7/8447 -> all into no
#Wasting = 3/8447 -> all into no
#Overweight by height = 36/8447 -> all into no
#Underweight = 8/8447 -> all into no
#Overweight by age = 13/8447 -> all into no
#Finished flooring material = 5/8447 -> all into no

KRdata <- KRdata %>% 
  mutate_at(c('partner_eduyr', 'mother_occupation',
              'partner_occupation', 'nt_mdd', 'nt_mmf', 'nt_mad',
              'nt_ch_ht_age', 'nt_ch_wt_ht', 'nt_ch_wt_age', 'ph_floor'), ~replace_na(.,0))

###Replace with 3
#Mother BMI = 59/8447 -> all into normal

KRdata <- KRdata %>%
  mutate_at(c('mother_bmi'), ~ replace_na(.,3))

###Replace with 4
#Mother anemia = 104/8447 -> all into normal 

KRdata <- KRdata %>% 
  mutate_at(c('mother_anemia'), ~replace_na(.,4))

###################################
# Create community-level variables#
###################################

# ------------------------------------------------
# 1. Helper flags at individual level
# ------------------------------------------------

KRdata <- KRdata %>%
  mutate(
    w = v005 / 1e6,
    
    # Individual SES flags
    is_poor         = ph_wealth_quint %in% c(1, 2),     # poorest/poorer
    mother_edu_ge1  = mother_eduyr >= 1,                # ≥ primary
    partner_edu_ge1 = partner_eduyr >= 1,               # ≥ primary
    
    # Occupation flags (you can refine as needed)
    mother_occup  = mother_occupation != 0,
    partner_occup = partner_occupation != 0,
  )


# Make sure binary variables are numeric 0/1 (some are haven-labelled).
to_numeric01 <- function(x) {
  if (inherits(x, c("haven_labelled", "haven_labelled_spss"))) {
    x <- haven::zap_labels(x)
  }
  if (is.factor(x)) {
    x <- as.character(x)
  }
  as.numeric(x)
}

KRdata <- KRdata %>%
  mutate(
    ph_animals      = to_numeric01(ph_animals),
    ph_agriland     = to_numeric01(ph_agriland),
    sanitation_imp  = to_numeric01(sanitation_imp),
    water_imp       = to_numeric01(water_imp),
    ph_floor        = to_numeric01(ph_floor),
    ph_electric     = to_numeric01(ph_electric),
    mass_media      = to_numeric01(mass_media),
    geography_rural = to_numeric01(geography_rural),
    ml_itnhh        = to_numeric01(ml_itnhh)
  )

# ------------------------------------------------
# 2. Cluster-level (community) summaries (PSU level)
# ------------------------------------------------
weighted_prop <- function(x, w) {
  stats::weighted.mean(as.numeric(x), w = w, na.rm = TRUE)
}

comm_vars <- KRdata %>%
  group_by(v021) %>%
  summarise(
    # ---- SES / education ----
    comm_pov_prop        = weighted_prop(is_poor,         w),
    comm_motheredu_prop  = weighted_prop(mother_edu_ge1,  w),
    comm_partneredu_prop = weighted_prop(partner_edu_ge1, w),
    
    # ---- Occupation context (optional; keep a few) ----
    comm_mother_occupation_prop  = weighted_prop(mother_occup,  w),
    comm_partner_occupation_prop = weighted_prop(partner_occup, w),
    
    # ---- Livestock & agricultural land ----
    comm_animals_prop  = weighted_prop(ph_animals,  w),
    comm_agriland_prop = weighted_prop(ph_agriland, w),
    
    # ---- Water, sanitation, housing, electricity ----
    comm_sanitation_imp_prop = weighted_prop(sanitation_imp, w),
    comm_water_imp_prop      = weighted_prop(water_imp,      w),
    comm_floor_finished_prop = weighted_prop(ph_floor,       w),
    comm_electric_prop       = weighted_prop(ph_electric,    w),
    
    # ---- Mass media ----
    comm_massmedia_prop = weighted_prop(mass_media, w),
    # ---- ITN ownership context ----
    comm_itnhh_prop     = weighted_prop(ml_itnhh, w),
    
    .groups = "drop"
  ) %>%
  # Scale to "per 10% increase" (for modeling)
  mutate(
    comm_pov_prop_scaled        = comm_pov_prop        * 10,
    comm_motheredu_prop_scaled  = comm_motheredu_prop  * 10,
    comm_partneredu_prop_scaled = comm_partneredu_prop * 10,
    comm_mother_occupation_prop_scaled  = comm_mother_occupation_prop  * 10,
    comm_partner_occupation_prop_scaled = comm_partner_occupation_prop * 10,
    comm_animals_prop_scaled       = comm_animals_prop       * 10,
    comm_agriland_prop_scaled      = comm_agriland_prop      * 10,
    comm_sanitation_imp_prop_scaled = comm_sanitation_imp_prop * 10,
    comm_water_imp_prop_scaled      = comm_water_imp_prop      * 10,
    comm_floor_finished_prop_scaled = comm_floor_finished_prop * 10,
    comm_electric_prop_scaled       = comm_electric_prop       * 10,
    comm_massmedia_prop_scaled      = comm_massmedia_prop      * 10,
    comm_itnhh_prop_scaled          = comm_itnhh_prop          * 10
  )

# ------------------------------------------------
# 3. Merge PSU-level community vars back to individual data
# ------------------------------------------------

KRdata <- KRdata %>%
  left_join(comm_vars, by = "v021")

########################################
############ Create table 1 ############
########################################

# helper: TRUE if numeric 0/1 (ignores NA)
is_binary01_num <- function(z) is.numeric(z) && all(z[!is.na(z)] %in% c(0,1))

# helper: TRUE if 2-level factor
is_two_level_factor <- function(z) is.factor(z) && nlevels(z) == 2

# Coerce 2-level factors to 0/1 with 1 = second level (preserves NA)
to01 <- function(z) {
  if (is_binary01_num(z)) return(z)
  if (is_two_level_factor(z)) return(as.numeric(z == levels(z)[2]))
  z
}

# ---- Always-keep (theory/DAG) set (updated to include 10%-scaled community vars) ----
forced_vars_candidate <- c(
  "ph_wealth_quint",
  "mother_eduyr",
  "geography_regions",
  "comm_pov_prop_scaled",       # community poverty per 10% increase
  "comm_motheredu_prop_scaled"  # community ≥ primary maternal education per 10% increase
)

# ---- Build the subset used for Table 1 & modeling ----
# (now including 10%-scaled community-level variables for modeling)
KR_subset <- dplyr::select(
  KRdata,
  # design & outcome
  "v021", "v022", "v005",
  "ml_diagnosis", "age_24", "age",
  
  # individual-level variables
  "sex_female",
  "ph_wealth_quint",
  "mother_eduyr", "partner_eduyr", "mother_occupation",
  "partner_occupation", "nt_mdd", "nt_mmf", "nt_mad",
  "iron", "vit_a", "mother_anemia", "mother_bmi", "nt_ch_ht_age",
  "nt_ch_wt_ht","nt_ch_wt_age", "nt_bf_curr", "ml_itnhh", "ml_hhaccess", "ml_netcat",
  "ph_electric", "ph_floor", "ph_rooms_sleep", "ph_cook_place", "ph_agriland",
  "ph_animals", "internet_fq", "mass_media", "sanitation_imp",
  "water_imp", "hh_members", "geography_regions", "geography_rural",
  "visit_health", "access_health", "share_toilet", "delivery_place",
  "preg_giparasites", "stool_disposal", "handwash_place", "soap_water", 
  "handwash_dirty", "own_water",
  
  # community-level (PSU-level) variables for modeling (10%-scaled)
  "comm_pov_prop_scaled",
  "comm_motheredu_prop_scaled",
  "comm_partneredu_prop_scaled",
  "comm_mother_occupation_prop_scaled",
  "comm_partner_occupation_prop_scaled",
  "comm_animals_prop_scaled",
  "comm_agriland_prop_scaled",
  "comm_sanitation_imp_prop_scaled",
  "comm_water_imp_prop_scaled",
  "comm_floor_finished_prop_scaled",
  "comm_electric_prop_scaled",
  "comm_massmedia_prop_scaled",
  "comm_itnhh_prop_scaled"
)

# ---- Create survey design for the subset table ----
KR_subset_svy <- svydesign(
  id = ~v021,
  strata = ~v022,
  weights = ~I(v005/1e6),
  nest = TRUE,
  data = KR_subset
)

# Build an EDA/modeling copy that converts only 2-level factors to 0/1
# continuous & 10%-scaled community proportions stay continuous
KR_subset <- KR_subset %>%
  mutate(across(everything(), to01))

# =========================
# Table 1 (individual-level)
# =========================
# Regroup headings according to your requested domains by renaming
# (order of variables kept largely as before, but headings changed).

variable_names <- c(
  # ---- Socio-demographic ----
  "Socioeconomic and demographic",
  "Female sex", "Household wealth quintile", " Richest", " Richer", " Middle", " Poorer", " Poorest", 
  "Education level of the mother", " Higher", " Secondary", " Primary"," None", 
  "Education level of the father", " Higher", " Secondary", " Primary"," None",
  "Occupation of the mother", " Office work", " Manual labor", " None",
  "Occupation of the father", " Office work", " Manual labor", " None",
  "Area of residency", " North Central", " North East", " North West", " South South or South East", " South West",
  "Household lives in rural regions", 
  
  # ---- Child & Maternal nutrition ----
  "Child & Maternal nutrition",
  "Fulfills minimum dietary diversity (MDD)", "Fulfills minimum meal frequency (MMF)", "Fulfills minimum acceptable diet (MAD)",
  "Given iron supplement within 7 days", "Child given vitamin A supplement within 6 months",
  "Anemia level of the mother", " Normal", " Mild", " Moderate", " Severe", 
  "BMI status of the mother", " Normal", " Very thin", " Thin", " Overweight", " Obese",
  "Height-for-age", " With stunting", " Normal", 
  "Weight-for-height", " With wasting", " Normal", " Overweight",
  "Weight-for-age", " Underweight", " Normal", " Overweight",
  "Currently breastfed", 
  
  # ---- Household structure and physical environment ----
  "Household structure and physical environment",
  "Number of household members", " 2-4 members", " 5-7 members", " 8 or more members", 
  "Access to electricity", "House has finished flooring material",
  "2 or more rooms for sleeping", "Cooks inside", "Owns agricultural land", "Owns livestocks", 
  "Access to improved sanitation facility", "Shares toilet with others", 
  "Place of handwashing", " At home: fixed", " At home: mobile", " Not at home", 
  "Has soap and water when washing hands", "Ash, mud, sand present at handwashing place",
  "Access to improved water source",  "Owns water source", 
  
  # ---- Malaria prevention and health-related behaviors ----
  "Prevention and health-related behaviors",
  "Household owns ITN", "Household owns more than 1 ITN per 2 members",
  "Slept under an ITN last night", 
  "Mother uses the internet at least once a week", "Mother exposed to mass media at least once a week", 
  "Visited health facility in the last 12 months", "Mother had problems in accessing healthcare", "Appropriate disposal of children's stool", 
  "Mother received drugs for intestinal parasites during pregnancy", 
  "Place of delivery", " Own or other's home", " Government facility", " Private facility"
)

variable_code <- c(
  NA, 
  "sex_female", rep("ph_wealth_quint", 6), rep("mother_eduyr", 5), 
  rep("partner_eduyr", 5), rep("mother_occupation", 4), rep("partner_occupation", 4),
  rep("geography_regions", 6), "geography_rural",
  NA,
  "nt_mdd", "nt_mmf", "nt_mad", "iron", "vit_a",
  rep("mother_anemia",5), rep("mother_bmi", 6),
  "nt_ch_ht_age", "nt_ch_ht_age", "nt_ch_ht_age", 
  "nt_ch_wt_ht", "nt_ch_wt_ht", "nt_ch_wt_ht", "nt_ch_wt_ht", 
  "nt_ch_wt_age", "nt_ch_wt_age", "nt_ch_wt_age", "nt_ch_wt_age", 
  "nt_bf_curr", 
  NA,
  "hh_members", "hh_members", "hh_members", "hh_members",
  "ph_electric", "ph_floor", "ph_rooms_sleep", "ph_cook_place", "ph_agriland",  "ph_animals", 
  "sanitation_imp", "share_toilet", 
  rep("handwash_place", 4), "soap_water", "handwash_dirty",
  "water_imp", "own_water",
  NA,
  "ml_itnhh", "ml_hhaccess", "ml_netcat",
  "internet_fq", "mass_media", "visit_health", "access_health", "stool_disposal",
  "preg_giparasites", rep("delivery_place", 4)
)

response_list <- c(
  NA,
  1, NA, 5, 4, 3, 2, 1, NA, 3, 2, 1, 0, NA, 3, 2, 1, 0, NA, 1, 2, 0, NA, 1, 2, 0,
  NA, 1, 2, 3, 4, 5, 1, 
  NA,
  1, 1, 1, 1, 1, 
  NA, 4, 3, 2, 1, NA, 3, 1, 2, 4, 5, 
  NA, 1, 0, NA, 2, 1, 0, NA, 2, 1, 0, 
  1,
  NA,
  NA, 1, 2, 3,
  1, 1, 1, 1, 1, 1, 
  1, 1, 
  NA, 1, 2, 3,
  1, 1, 1, 1, 
  NA,
  1, 1,  1,
  1, 1, 1, 1, 1,
  1, NA, 1, 2, 3
)

reference_list <- c(
  NA,
  0, NA, 5, 5, 5, 5, 5, NA, 3, 3, 3, 3, NA, 3, 3, 3, 3, NA, 1, 1, 1, NA, 1, 1, 1,
  NA, 1, 1, 1, 1, 1, 0, 
  NA,
  0, 0, 0, 0, 0, 
  NA, 4, 4, 4, 4, NA, 3, 3, 3, 3, 3, 
  NA, 0, 0, NA, 1, 1, 1, NA, 1, 1, 1, 
  0,
  NA,
  NA, 1, 1, 1,
  0, 0, 0, 0, 0, 0, 
  0, 0, 
  NA, 1, 1, 1,
  0, 0, 0, 0, 
  NA,
  0, 0, 0,
  0, 0, 0, 0, 0,
  0, NA, 1, 1, 1
)

# -------------------------------------------
# Build Table 1 structure
# -------------------------------------------
table1 <- data.frame(matrix(NA, ncol = 10, nrow = length(variable_names)))
colnames(table1) <- c(
  'Variable', 'variable_code', 'response', 'reference',
  'Younger children N (%)', 'Younger children PR (95% CI)',
  'p-value', 'Older children N (%)', 'Older children PR (95% CI)', 'p-value '
)
table1$Variable      <- variable_names
table1$variable_code <- variable_code
table1$response      <- response_list
table1$reference     <- reference_list

# -------------------------------------------
# Set reference levels per row (for categorical variables)
# -------------------------------------------

for (i in seq_len(nrow(table1))) {
  var  <- table1$variable_code[i]
  resp <- table1$response[i]
  if (is.na(var) | is.na(resp)) next
  ref_level <- table1$reference[i]
  KR_subset[[var]] <- relevel(factor(KR_subset[[var]]), ref = as.character(ref_level))
}

# Rebuild design (factors now updated)
KR_subset_svy <- svydesign(
  id = ~v021,
  strata = ~v022,
  weights = ~I(v005/1e6),
  nest = TRUE,
  data = KR_subset
)

# -------------------------------------------
# Helper to compute weighted N(%) by age_24 band
# -------------------------------------------
weighted_npct <- function(design, var, level, band_value) {
  dsub <- subset(design, age_24 == band_value)
  tab  <- svytable(reformulate(c(var), response = NULL), dsub)
  prp  <- prop.table(tab)
  lvl_char <- as.character(level)
  n   <- as.numeric(tab[lvl_char])
  pct <- as.numeric(prp[lvl_char]) * 100
  if (is.na(n)) return("—")
  paste0(round(n, 0), " (", round(pct, 1), "%)")
}

# -------------------------------------------
# MAIN LOOP: build N(%) and age-adjusted PRs
# -------------------------------------------
df_age <- 4

for (i in seq_len(nrow(table1))) {
  var  <- table1$variable_code[i]
  resp <- table1$response[i]
  ref  <- table1$reference[i]
  if (is.na(var) | is.na(resp)) next
  
  # Weighted N(%) by bands
  table1$`Younger children N (%)`[i] <- weighted_npct(KR_subset_svy, var, resp, band_value = 0)
  table1$`Older children N (%)`[i]   <- weighted_npct(KR_subset_svy, var, resp, band_value = 1)
  
  # If response equals reference, mark as Reference for PR columns
  if (identical(resp, ref)) {
    table1$`Younger children PR (95% CI)`[i] <- 'Reference'
    table1$`p-value`[i]                      <- NA
    table1$`Older children PR (95% CI)`[i]   <- 'Reference'
    table1$`p-value `[i]                     <- NA
    next
  }
  
  # Coefficient name (for 0/1 and coded factors this matches paste0(var, level))
  coef_name <- paste0(var, resp)
  
  # Younger band (age_24 == 0): survey-weighted Poisson (log link) + ns(age)
  d_y  <- subset(KR_subset_svy, age_24 == 0)
  ok_y <- try(svytable(as.formula(paste0("~", var)), d_y), silent = TRUE)
  if (!inherits(ok_y, "try-error")) {
    fit_y <- try(
      svyglm(
        as.formula(paste0("ml_diagnosis ~ ns(age, ", df_age, ") + ", var)),
        design = d_y,
        family = quasipoisson(link = "log")
      ),
      silent = TRUE
    )
    if (!inherits(fit_y, "try-error") && !is.na(coef(fit_y)[coef_name])) {
      est  <- coef(fit_y)[coef_name]
      se   <- sqrt(vcov(fit_y)[coef_name, coef_name])
      PR   <- exp(est)
      LCL  <- exp(est - 1.96 * se)
      UCL  <- exp(est + 1.96 * se)
      pval <- coef(summary(fit_y))[coef_name, "Pr(>|t|)"]
      table1$`Younger children PR (95% CI)`[i] <- sprintf('%.2f (%.2f–%.2f)', PR, LCL, UCL)
      table1$`p-value`[i] <- sprintf('%.3f', pval)
    }
  }
  
  # Older band (age_24 == 1)
  d_o  <- subset(KR_subset_svy, age_24 == 1)
  ok_o <- try(svytable(as.formula(paste0("~", var)), d_o), silent = TRUE)
  if (!inherits(ok_o, "try-error")) {
    fit_o <- try(
      svyglm(
        as.formula(paste0("ml_diagnosis ~ ns(age, ", df_age, ") + ", var)),
        design = d_o,
        family = quasipoisson("log")
      ),
      silent = TRUE
    )
    if (!inherits(fit_o, "try-error") && !is.na(coef(fit_o)[coef_name])) {
      est  <- coef(fit_o)[coef_name]
      se   <- sqrt(vcov(fit_o)[coef_name, coef_name])
      PR   <- exp(est)
      LCL  <- exp(est - 1.96 * se)
      UCL  <- exp(est + 1.96 * se)
      pval <- coef(summary(fit_o))[coef_name, "Pr(>|t|)"]
      table1$`Older children PR (95% CI)`[i] <- sprintf('%.2f (%.2f–%.2f)', PR, LCL, UCL)
      table1$`p-value `[i] <- sprintf('%.3f', pval)
    }
  }
}

# ---- Keep the always-keep list handy for later steps (using 10%-scaled vars) ----
always_keep <- forced_vars_candidate

# ============================================
# Add Community-level block to Table 1
# (mean (25th–75th percentile), PSU-level, unweighted)
# ============================================

# Variables at PSU level (from comm_vars)
community_vars <- c(
  "comm_pov_prop",
  "comm_motheredu_prop",
  "comm_partneredu_prop",
  "comm_mother_occupation_prop",
  "comm_partner_occupation_prop",
  "comm_animals_prop",
  "comm_agriland_prop",
  "comm_sanitation_imp_prop",
  "comm_water_imp_prop",
  "comm_floor_finished_prop",
  "comm_electric_prop",
  "comm_massmedia_prop",
  "comm_itnhh_prop"
)

community_labels <- c(
  "Household is poor or poorer",
  "Mother received primary education or higher",
  "Father received primary education or higher",
  "Mother is employed",
  "Father is employed",
  "Household owns livestock",
  "Household owns agricultural land",
  "Access to improved sanitation facility",
  "Access to improved water source",
  "House has finished flooring material",
  "Access to electricity",
  "Mother exposed to mass media at least once a week",
  "Household owns ITN"
)
comm_mean_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return("—")
  q <- 100 * quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  m <- 100 * median(x, na.rm = TRUE)
  sprintf("%.1f%% (%.1f%%–%.1f%%)", m, q[1], q[2])
}

comm_stats <- vapply(
  community_vars,
  function(v) comm_mean_iqr(comm_vars[[v]]),
  FUN.VALUE = character(1)
)

community_variable_names <- c(
  # ---- Community socio-demographic and economic factors ----
  "Socioeconomic and demographic",
  "Household is poor or poorer",
  "Mother received primary education or higher",
  "Father received primary education or higher",
  "Mother is employed",
  "Father is employed",
  
  # ---- Community household assets and physical environment ----
  "Household structure and physical environment",
  "Household owns livestock",
  "Household owns agricultural land",
  "House has finished flooring material",
  "Access to electricity",
  "Access to improved sanitation facility",
  "Access to improved water source",
  
  # ---- Community information and malaria-prevention resources ----
  "Prevention and health-related behaviors",
  "Mother exposed to mass media at least once a week",
  "Household owns ITN"
)

community_variable_code <- c(
  # heading
  NA,
  # socio-demographic/economic
  "comm_pov_prop",
  "comm_motheredu_prop",
  "comm_partneredu_prop",
  "comm_mother_occupation_prop",
  "comm_partner_occupation_prop",
  # heading
  NA,
  # assets & physical env
  "comm_animals_prop",
  "comm_agriland_prop",
  "comm_floor_finished_prop",
  "comm_electric_prop",
  "comm_sanitation_imp_prop",
  "comm_water_imp_prop",
  # heading
  NA,
  # info & malaria-prevention
  "comm_massmedia_prop",
  "comm_itnhh_prop"
)

# map comm_stats into the non-heading rows
comm_stats_full <- rep(NA_character_, length(community_variable_names))
idx_non_heading <- which(!is.na(community_variable_code))
comm_stats_full[idx_non_heading] <- comm_stats[match(
  community_variable_code[idx_non_heading],
  community_vars
)]

community_block <- data.frame(
  Variable                      = community_variable_names,
  variable_code                 = community_variable_code,
  response                      = NA,
  reference                     = NA,
  `Median proportion (25%-75%)`   = comm_stats_full,
  stringsAsFactors              = FALSE,
  check.names                   = FALSE
)

# You can then print table1 with kable/flextable as before.
# Example:
# kable(table1, align = 'l', booktabs = TRUE) %>%
#   kable_styling(full_width = FALSE, font_size = 11)

###########################################################
## 0. Setup: age-band and analysis objects
###########################################################

# Assumed existing: KR_subset with
# - ml_diagnosis (0/1)
# - age (months)
# - v021 (PSU), v022 (strata), v005 (sampling weight)
# - covariates (including community-level variables)

# ------------------------------------------------
# Age band: 6–23.9 vs 24–59 months
# ------------------------------------------------
KR_subset <- KR_subset %>%
  mutate(
    age_band = cut(
      age,
      breaks = c(6, 24, 60),   # [6,24) and [24,60)
      right  = FALSE,
      labels = c("6-23", "24-59")
    )
  )

# Relevel so 6–23 is reference
KR_subset$age_band <- relevel(KR_subset$age_band, ref = "6-23")

# Basic pieces
y_var         <- "ml_diagnosis"
age_band_var  <- "age_band"
id_var        <- "v021"
strata_var    <- "v022"
weight_var    <- "v005"

w   <- KR_subset[[weight_var]] / 1e6
psu <- KR_subset[[id_var]]

## --------------------------------------------------------
## Global caches for efficiency
## --------------------------------------------------------
step1_cache_env <- new.env(parent = emptyenv())
step2_cache_env <- new.env(parent = emptyenv())

run_3step_with_thresholds <- function(
    KR_subset,
    ## Table-1 metadata
    variable_names,
    variable_code,
    response_list,
    community_variable_names,
    community_variable_code,
    community_labels = NULL,
    
    ## Core variable names
    y_var        = "ml_diagnosis",
    age_band_var = "age_band",
    id_var       = "v021",
    strata_var   = "v022",
    weight_var   = "v005",
    
    ## Always-keep variables in Step 1
    always_keep = c(
      "ph_wealth_quint",
      "mother_eduyr",
      "geography_regions",
      "comm_pov_prop_scaled",
      "comm_motheredu_prop_scaled"
    ),
    
    ## Family choice for Step 2 & Step 3
    family_type = "quasipoisson",
    
    ## STEP 1 thresholds
    step1_stab_cut          = 0.60,  # stability cut for lasso selection
    step1_B                 = 200,   # number of lasso CV repeats
    step1_min_psu_per_fold  = 1,     # m_min used to define K
    step1_lambda_type       = c("lambda.min", "lambda.1se"),
    
    ## STEP 2 thresholds
    step2_p_max    = 0.20,
    step2_stab_cut = 0.60,
    step2_B        = 200,
    random_seed    = 20251129L,
    
    ## STEP 3 thresholds
    step3_alpha_main   = 0.05,
    step3_alpha_int    = 0.05,
    step3_qic_tol      = 0,
    step3_rel_mean_th  = 0.10,
    step3_rel_max_th   = 0.20
) {
  step1_lambda_type <- match.arg(step1_lambda_type)
  family_type       <- match.arg(family_type)
  
  ## Family object used in Step 2 (glm) and Step 3 (svyglm)
  fam_glm <- switch(
    family_type,
    "quasipoisson"  = quasipoisson("log"),
    "quasibinomial" = quasibinomial("log")
  )
  
  ## ---------------------------------------------------------
  ## Helper: build a readable config name for sensitivity runs
  ## ---------------------------------------------------------
  make_config_name <- function() {
    sprintf(
      paste0(
        "step1_stab_%.2f",
        "_lambda_%s",
        "_step2_p_%.2f_stab_%.2f",
        "_step3_p_%.2f_qic_%g_meanpr_%.2f_maxpr_%.2f",
        "_fam_%s"
      ),
      step1_stab_cut,
      step1_lambda_type,
      step2_p_max, step2_stab_cut,
      step3_alpha_main,
      step3_qic_tol,
      step3_rel_mean_th,
      step3_rel_max_th,
      family_type
    )
  }
  
  ## ---------------------------------------------------------
  ## Shared helper functions
  ## ---------------------------------------------------------
  
  collapse_to_parent <- function(x) {
    x <- unique(x)
    x <- gsub("([0-9]+|_[0-9]+)$", "", x)
    x <- gsub("[0-9]+$", "", x)
    x <- gsub("(\\.[0-9]+)$", "", x)
    x <- gsub("([12]$)", "", x)
    x <- gsub("([a-z]+)[0-9]+$", "\\1", x)
    x <- gsub("\\.[A-Za-z0-9]+$", "", x)
    x <- gsub("^(.*?)[0-9]+$", "\\1", x)
    unique(x)
  }
  
  make_psu_folds_stratified <- function(psu, strata, K = 5) {
    foldid <- rep(NA_integer_, length(psu))
    strata_levels <- unique(strata)
    for (s in strata_levels) {
      idx_s   <- which(strata == s)
      psu_s   <- psu[idx_s]
      psu_ids <- unique(psu_s)
      psu_folds <- sample(rep(1:K, length.out = length(psu_ids)))
      names(psu_folds) <- psu_ids
      foldid[idx_s] <- psu_folds[as.character(psu_s)]
    }
    foldid
  }

  build_boot_idx_stratified <- function(id, strata) {
    idx_by_psu <- split(seq_along(id), id)
    psu_by_stratum <- split(id, strata)
    psu_by_stratum <- lapply(psu_by_stratum, unique)
    sampled_rows <- integer(0)

    for (s in names(psu_by_stratum)) {
      psus_s <- psu_by_stratum[[s]]
      sampled_psus <- sample(psus_s, size = length(psus_s), replace = TRUE)
      sampled_rows <- c(
        sampled_rows,
        unlist(idx_by_psu[as.character(sampled_psus)], use.names = FALSE)
      )
    }

    sampled_rows
  }
  
  qic_like <- function(fit) fit$deviance + 2 * length(coef(fit))
  
  wald_block <- function(fit, pattern) {
    beta <- coef(fit)
    V    <- vcov(fit)
    
    keep <- grep(pattern, names(beta), value = TRUE)
    if (!length(keep)) return(NA_real_)
    
    keep <- keep[!is.na(beta[keep])]
    if (!length(keep)) return(NA_real_)
    
    b  <- beta[keep]
    Vb <- V[keep, keep, drop = FALSE]
    
    if (length(keep) == 1L) {
      W  <- as.numeric((b^2) / Vb)
      df <- 1L
    } else {
      W  <- as.numeric(t(b) %*% solve(Vb, b))
      df <- qr(Vb)$rank
    }
    pchisq(W, df = df, lower.tail = FALSE)
  }
  
  build_int_pattern <- function(v, age_band_var) {
    paste0(
      "(",
      "^", age_band_var, ".*:", v,
      ")|(",
      "^", v, ":", age_band_var, ".*",
      ")"
    )
  }
  
  ## ---------------------------------------------------------
  ## Weights & PSU
  ## ---------------------------------------------------------
  w   <- KR_subset[[weight_var]] / 1e6
  psu <- KR_subset[[id_var]]
  
  ## =========================================================
  ## STEP 1 — Weighted Lasso with age_band×X (cached by lambda)
  ## =========================================================
  
  all_candidates <- setdiff(
    names(KR_subset),
    c(id_var, strata_var, weight_var, y_var, age_band_var)
  )
  
  ## Key for Step 1 cache: lambda, replicates, fold rule, and seed.
  step1_key <- paste(
    step1_lambda_type,
    step1_B,
    step1_min_psu_per_fold,
    random_seed,
    sep = "|"
  )
  
  if (exists(step1_key, envir = step1_cache_env, inherits = FALSE)) {
    base1 <- get(step1_key, envir = step1_cache_env, inherits = FALSE)
    
    freq_tbl              <- base1$freq_tbl
    freq_ageint           <- base1$freq_ageint
    step1_all_candidates  <- base1$step1_all_candidates
    colnames_mm           <- base1$colnames_mm
    col_parent            <- base1$col_parent
    col_type              <- base1$col_type
    penalty               <- base1$penalty
    
  } else {
    ## --- Build design matrix and run stability LASSO once for this lambda ---
    form_main <- as.formula(
      paste(
        y_var, "~",
        "age",
        "+", age_band_var,
        "+",
        paste(always_keep, collapse = " + "),
        "+",
        paste(setdiff(all_candidates, c(always_keep, "age")),
              collapse = " + ")
      )
    )
    
    mm_main <- model.matrix(form_main, data = KR_subset)[, -1, drop = FALSE]
    y       <- KR_subset[[y_var]]
    
    tm          <- attr(mm_main, "terms")
    assign_vec  <- attr(mm_main, "assign")
    term_labels <- attr(tm, "term.labels")
    
    col_parent_main_raw <- vapply(assign_vec, function(k) {
      if (k == 0) "(Intercept)" else term_labels[k]
    }, character(1))
    col_parent_main <- col_parent_main_raw
    
    candidate_for_int <- setdiff(
      all_candidates,
      c("age", age_band_var, always_keep)
    )
    
    build_age_band_int_matrix <- function(parent_var) {
      f_int <- as.formula(paste("~", age_band_var, ":", parent_var))
      model.matrix(f_int, data = KR_subset)[, -1, drop = FALSE]
    }
    
    mm_int_list <- lapply(candidate_for_int, build_age_band_int_matrix)
    mm_int      <- if (length(mm_int_list)) do.call(cbind, mm_int_list) else NULL
    
    if (is.null(mm_int)) {
      mm <- mm_main
    } else {
      mm <- cbind(mm_main, mm_int)
    }
    colnames_mm <- colnames(mm)
    
    ## Map columns to parents & type
    col_parent <- character(length(colnames_mm))
    col_type   <- character(length(colnames_mm))
    
    col_to_parent_main <- function(col) {
      base <- strsplit(col, ":", fixed = TRUE)[[1]][1]
      base
    }
    
    for (j in seq_along(colnames_mm)) {
      cname <- colnames_mm[j]
      if (grepl(paste0("^", age_band_var), cname) && grepl(":", cname)) {
        tmp <- sub(paste0("^", age_band_var, "[^:]*:"), "", cname)
        tmp <- collapse_to_parent(tmp)
        col_parent[j] <- tmp
        col_type[j]   <- "age_int"
      } else {
        tmp <- col_to_parent_main(cname)
        tmp <- collapse_to_parent(tmp)
        col_parent[j] <- tmp
        col_type[j]   <- "main"
      }
    }
    
    penalty <- rep(1, ncol(mm))
    names(penalty) <- colnames_mm
    
    penalty[col_parent %in% c("age", age_band_var) & col_type == "main"] <- 0
    for (ak in always_keep) {
      idx_ak_main <- which(col_parent == ak & col_type == "main")
      penalty[idx_ak_main] <- 0
    }
    
    run_one_bootstrap_lasso <- function(seed, K = 5) {
      set.seed(seed)
      rows <- build_boot_idx_stratified(
        id = KR_subset[[id_var]],
        strata = KR_subset[[strata_var]]
      )

      mm_b <- mm[rows, , drop = FALSE]
      y_b <- y[rows]
      w_b <- w[rows]
      psu_b <- KR_subset[[id_var]][rows]
      strata_b <- KR_subset[[strata_var]][rows]

      foldid <- make_psu_folds_stratified(
        psu    = psu_b,
        strata = strata_b,
        K      = K
      )
      
      cvfit <- cv.glmnet(
        x = mm_b, y = y_b,
        family = "binomial",
        weights = w_b,
        foldid  = foldid,
        alpha   = 1,
        standardize    = TRUE,
        penalty.factor = penalty,
        type.measure   = "deviance",
        parallel       = FALSE
      )
      
      beta <- coef(cvfit, s = step1_lambda_type)
      sel  <- rownames(beta)[as.numeric(beta) != 0]
      sel  <- setdiff(sel, "(Intercept)")
      
      list(selected_cols = sel, cvfit = cvfit)
    }
    
    ## Stability selection: K based on PSU-stratum structure
    tab_psu_stratum <- KR_subset |>
      dplyr::distinct(!!rlang::sym(strata_var), !!rlang::sym(id_var)) |>
      dplyr::count(!!rlang::sym(strata_var), name = "n_psu")
    
    mean_psu_per_stratum <- mean(tab_psu_stratum$n_psu)
    m_min <- step1_min_psu_per_fold
    K_max_design <- min(floor(tab_psu_stratum$n_psu / m_min))
    K <- min(5L, K_max_design)
    if (K < 2L) {
      stop("At least two PSUs per stratum are required for Step 1 CV folds.")
    }

    set.seed(random_seed + 1000L)
    bootstrap_seeds <- sample.int(.Machine$integer.max, step1_B)
    stab_list <- vector("list", step1_B)
    for (b in seq_len(step1_B)) {
      stab_list[[b]] <- run_one_bootstrap_lasso(
        seed = bootstrap_seeds[[b]],
        K = K
      )
    }
    
    selected_cols_all <- lapply(stab_list, `[[`, "selected_cols")
    
    selected_parent_each <- lapply(selected_cols_all, function(cols) {
      if (length(cols) == 0) return(NULL)
      idx <- match(cols, colnames_mm)
      tibble::tibble(
        column = cols,
        parent = col_parent[idx],
        type   = col_type[idx]
      )
    })
    
    var_selected_each <- lapply(selected_parent_each, function(df) {
      if (is.null(df) || nrow(df) == 0) return(character(0))
      unique(df$parent)
    })
    
    freq_tbl <- sort(table(unlist(var_selected_each)) / step1_B, decreasing = TRUE)
    
    ageint_selected_each <- lapply(selected_parent_each, function(df) {
      if (is.null(df) || nrow(df) == 0) return(character(0))
      unique(df$parent[df$type == "age_int"])
    })
    freq_ageint <- sort(table(unlist(ageint_selected_each)) / step1_B, decreasing = TRUE)
    
    step1_all_candidates <- all_candidates
    
    base1 <- list(
      freq_tbl             = freq_tbl,
      freq_ageint          = freq_ageint,
      step1_all_candidates = step1_all_candidates,
      colnames_mm          = colnames_mm,
      col_parent           = col_parent,
      col_type             = col_type,
      penalty              = penalty
    )
    assign(step1_key, base1, envir = step1_cache_env)
  }
  
  ## Now apply step1_stab_cut to cached frequencies
  stability_cut <- step1_stab_cut
  stable_vars   <- names(freq_tbl)[freq_tbl >= stability_cut]
  step1_keepers <- union(always_keep, stable_vars)
  
  age_int_candidates_stable_step1 <- names(freq_ageint)[freq_ageint >= stability_cut]
  
  step1_keepers_collapsed <- collapse_to_parent(step1_keepers)
  
  step1_results <- list(
    always_keep               = always_keep,
    stable_vars               = stable_vars,
    selection_freq            = as.data.frame(freq_tbl),
    selection_freq_ageint     = as.data.frame(freq_ageint),
    step1_keepers             = step1_keepers,
    age_int_candidates_stable = age_int_candidates_stable_step1,
    step1_keepers_collapsed   = step1_keepers_collapsed,
    column_map = data.frame(
      column  = colnames_mm,
      parent  = col_parent,
      type    = col_type,
      penalty = penalty,
      stringsAsFactors = FALSE
    ),
    all_candidates = step1_all_candidates,
    stability_cut  = stability_cut,
    B              = step1_B
  )
  
  ## =========================================================
  ## STEP 2 — Age-common vs age-varying via PSU bootstrap
  ## =========================================================
  
  vars_to_factor <- intersect(step1_keepers_collapsed, names(KR_subset))
  KR_subset[vars_to_factor] <- lapply(KR_subset[vars_to_factor], function(x) {
    if (is.character(x) || is.factor(x)) factor(x) else x
  })
  
  fit_step2_model <- function(dat, vars) {
    vars_norage <- setdiff(vars, c("age", "age_yrs", "age_24", age_band_var))
    rhs_main <- c(age_band_var, "age", vars_norage)
    rhs_int  <- unlist(lapply(vars_norage, function(v) paste0(age_band_var, ":", v)))
    form <- reformulate(c(rhs_main, rhs_int), response = y_var)
    glm(form,
        family  = fam_glm,
        data    = dat,
        weights = dat[[weight_var]] / 1e6)
  }
  
  extract_agevary_vars_step2 <- function(fit, vars, age_band_var, p_max = 0.20) {
    out <- character(0)
    for (v in vars) {
      pat <- paste0("^", age_band_var, "[^:]*:", v)
      p_block <- wald_block(fit, pattern = pat)
      if (!is.na(p_block) && p_block < p_max) {
        out <- c(out, v)
      }
    }
    unique(out)
  }
  
  step2_vars <- step1_keepers_collapsed
  
  ## Step2 cache key: includes family_type
  step2_key <- paste(
    step1_lambda_type,
    sprintf("stab1=%.2f", step1_stab_cut),
    sprintf("pmax=%.2f", step2_p_max),
    sprintf("B=%d", step2_B),
    sprintf("fam=%s", family_type),
    sprintf("seed=%d", random_seed),
    sep = "|"
  )
  
  if (exists(step2_key, envir = step2_cache_env, inherits = FALSE)) {
    base2 <- get(step2_key, envir = step2_cache_env, inherits = FALSE)
    fit2_full          <- base2$fit2_full
    agevary_boot_list  <- base2$agevary_boot_list
    step2_B            <- base2$B_BOOT
  } else {
    ## Full model once for this (step1, p_max, family) combo
    fit2_full  <- fit_step2_model(KR_subset, step2_vars)
    
    set.seed(random_seed + 2000L)
    agevary_boot_list <- vector("list", step2_B)
    
    for (b in seq_len(step2_B)) {
      rows <- build_boot_idx_stratified(
        id     = KR_subset[[id_var]],
        strata = KR_subset[[strata_var]]
      )
      boot_d <- KR_subset[rows, , drop = FALSE]
      boot_d <- droplevels(boot_d)
      
      fit_b <- try(fit_step2_model(boot_d, step2_vars), silent = TRUE)
      if (inherits(fit_b, "try-error")) {
        agevary_boot_list[[b]] <- character(0)
      } else {
        agevary_boot_list[[b]] <- extract_agevary_vars_step2(
          fit_b, step2_vars, age_band_var, p_max = step2_p_max
        )
      }
    }
    
    base2 <- list(
      fit2_full         = fit2_full,
      agevary_boot_list = agevary_boot_list,
      B_BOOT            = step2_B
    )
    assign(step2_key, base2, envir = step2_cache_env)
  }
  
  ## Now apply step2_stab_cut to cached bootstrap results
  agevary_full <- extract_agevary_vars_step2(
    fit2_full, step2_vars, age_band_var, p_max = step2_p_max
  )
  
  all_sel   <- unlist(agevary_boot_list, use.names = FALSE)
  stab_freq <- sort(table(all_sel) / step2_B, decreasing = TRUE)
  stab_tbl  <- data.frame(
    variable  = names(stab_freq),
    stability = as.numeric(stab_freq),
    row.names = NULL
  )
  
  stab_cut <- step2_stab_cut
  age_varying_candidates_stable <- stab_tbl %>%
    dplyr::filter(stability >= stab_cut) %>%
    dplyr::arrange(dplyr::desc(stability)) %>%
    dplyr::pull(variable) %>%
    unique()
  
  age_common_candidates <- setdiff(
    step1_keepers_collapsed,
    age_varying_candidates_stable
  )
  
  step2_results <- list(
    model_full                    = fit2_full,
    stability_table               = stab_tbl,
    age_varying_candidates_stable = age_varying_candidates_stable,
    age_common_candidates         = age_common_candidates,
    p_max                         = step2_p_max,
    stab_cut                      = stab_cut,
    B_BOOT                        = step2_B
  )
  
  ## =========================================================
  ## STEP 3 — Survey-weighted purposeful selection
  ## =========================================================
  
  KR_subset_svy <- svydesign(
    ids     = as.formula(paste0("~", id_var)),
    strata  = as.formula(paste0("~", strata_var)),
    weights = as.formula(paste0("~I(", weight_var, "/1e6)")),
    data    = KR_subset,
    nest    = TRUE
  )
  
  check_confounding_prediction <- function(fit_full,
                                           fit_red,
                                           design,
                                           rel_thresh_mean = 0.10,
                                           rel_thresh_max  = 0.20) {
    dat <- design$variables
    
    X_full <- model.matrix(formula(fit_full), data = dat)
    X_red  <- model.matrix(formula(fit_red),  data = dat)
    
    b_full <- coef(fit_full)
    b_red  <- coef(fit_red)
    
    common_full <- intersect(colnames(X_full), names(b_full))
    X_full      <- X_full[, common_full, drop = FALSE]
    b_full      <- b_full[common_full]
    
    common_red <- intersect(colnames(X_red), names(b_red))
    X_red      <- X_red[, common_red, drop = FALSE]
    b_red      <- b_red[common_red]
    
    mu_full <- exp(as.numeric(X_full %*% b_full))
    mu_red  <- exp(as.numeric(X_red  %*% b_red))
    
    rel_diff <- abs(mu_full - mu_red) / pmax(mu_full, 1e-8)
    
    mean_rel <- mean(rel_diff, na.rm = TRUE)
    max_rel  <- max(rel_diff, na.rm = TRUE)
    
    list(
      mean_rel = mean_rel,
      max_rel  = max_rel,
      flag     = (mean_rel > rel_thresh_mean) || (max_rel > rel_thresh_max)
    )
  }
  
  refit_with_sets <- function(main_vars, interaction_vars, design) {
    rhs_main <- c(age_band_var, "age", main_vars)
    rhs_int  <- if (length(interaction_vars)) {
      unlist(lapply(interaction_vars, function(v) paste0(age_band_var, ":", v)))
    } else character(0)
    form <- reformulate(c(rhs_main, rhs_int), response = y_var)
    fit  <- svyglm(form, design = design, family = fam_glm)
    list(fit = fit, qic = qic_like(fit), form = form)
  }
  
  age_like <- c("age", age_band_var, "age_yrs", "age_24")
  all_main_vars_step3 <- setdiff(step1_keepers_collapsed, age_like)
  
  age_varying_step3 <- intersect(step2_results$age_varying_candidates_stable, all_main_vars_step3)
  age_common_step3  <- setdiff(all_main_vars_step3, age_varying_step3)
  
  main_terms_current       <- all_main_vars_step3
  interaction_vars_current <- age_varying_step3
  protected_vars           <- always_keep
  
  init_fit <- refit_with_sets(
    main_vars        = main_terms_current,
    interaction_vars = interaction_vars_current,
    design           = KR_subset_svy
  )
  
  fit_current <- init_fit$fit
  qic_current <- init_fit$qic
  
  step3_log <- tibble::tibble(
    var       = character(),
    phase     = character(),
    component = character(),
    p         = numeric(),
    dQIC      = numeric(),
    mean_rel  = numeric(),
    max_rel   = numeric(),
    conf_flag = logical()
  )
  
  alpha_int   <- step3_alpha_int
  alpha_main  <- step3_alpha_main
  qic_tol     <- step3_qic_tol
  rel_mean_th <- step3_rel_mean_th
  rel_max_th  <- step3_rel_max_th
  
  ## ----------------------
  ## Backward selection
  ## ----------------------
  
  repeat {
    drop_happened <- FALSE
    
    ## (a) Interactions
    if (length(interaction_vars_current) > 0L) {
      cand_int <- tibble::tibble(
        var       = character(),
        p         = numeric(),
        dQIC      = numeric(),
        mean_rel  = numeric(),
        max_rel   = numeric(),
        conf_flag = logical()
      )
      
      for (v in interaction_vars_current) {
        pat <- build_int_pattern(v, age_band_var)
        p_block <- wald_block(fit_current, pattern = pat)
        
        int_keep <- setdiff(interaction_vars_current, v)
        fit_red_obj <- refit_with_sets(
          main_vars        = main_terms_current,
          interaction_vars = int_keep,
          design           = KR_subset_svy
        )
        
        qic_red <- fit_red_obj$qic
        dQIC    <- qic_red - qic_current
        
        conf_obj <- check_confounding_prediction(
          fit_full        = fit_current,
          fit_red         = fit_red_obj$fit,
          design          = KR_subset_svy,
          rel_thresh_mean = rel_mean_th,
          rel_thresh_max  = rel_max_th
        )
        conf_flag <- conf_obj$flag
        
        step3_log <- dplyr::bind_rows(
          step3_log,
          tibble::tibble(
            var       = v,
            phase     = "backward",
            component = "interaction",
            p         = p_block,
            dQIC      = dQIC,
            mean_rel  = conf_obj$mean_rel,
            max_rel   = conf_obj$max_rel,
            conf_flag = conf_flag
          )
        )
        
        cand_int <- dplyr::bind_rows(
          cand_int,
          tibble::tibble(
            var       = v,
            p         = p_block,
            dQIC      = dQIC,
            mean_rel  = conf_obj$mean_rel,
            max_rel   = conf_obj$max_rel,
            conf_flag = conf_flag
          )
        )
      }
      
      ## Backward rule for interaction:
      ## p > alpha_int AND dQIC <= qic_tol AND not confounder
      eligible_int <- cand_int %>%
        dplyr::filter(
          !is.na(p),
          p > alpha_int,
          dQIC <= qic_tol,
          !conf_flag
        )
      
      if (nrow(eligible_int) > 0L) {
        ## 1) least QIC worsening (smallest dQIC)
        ## 2) if ties, largest p
        ## 3) if still ties, alphabetical
        chosen_int <- eligible_int %>%
          dplyr::arrange(dQIC, dplyr::desc(p), var) %>%
          dplyr::slice(1L)
        v_drop <- chosen_int$var
        
        int_keep <- setdiff(interaction_vars_current, v_drop)
        fit_red_obj <- refit_with_sets(
          main_vars        = main_terms_current,
          interaction_vars = int_keep,
          design           = KR_subset_svy
        )
        
        interaction_vars_current <- int_keep
        fit_current              <- fit_red_obj$fit
        qic_current              <- fit_red_obj$qic
        drop_happened            <- TRUE
      }
    }
    
    ## (b) Main effects
    main_candidates <- setdiff(
      main_terms_current,
      union(protected_vars, interaction_vars_current)
    )
    
    if (length(main_candidates) > 0L) {
      cand_main <- tibble::tibble(
        var       = character(),
        p         = numeric(),
        dQIC      = numeric(),
        mean_rel  = numeric(),
        max_rel   = numeric(),
        conf_flag = logical()
      )
      
      for (v in main_candidates) {
        p_v <- try(
          regTermTest(fit_current, as.formula(paste("~", v)))$p,
          silent = TRUE
        )
        p_v <- if (inherits(p_v, "try-error") || !length(p_v)) {
          NA_real_
        } else {
          as.numeric(p_v)[1]
        }
        
        mains_red <- setdiff(main_terms_current, v)
        fit_red_obj <- refit_with_sets(
          main_vars        = mains_red,
          interaction_vars = interaction_vars_current,
          design           = KR_subset_svy
        )
        
        qic_red <- fit_red_obj$qic
        dQIC    <- qic_red - qic_current
        
        conf_obj <- check_confounding_prediction(
          fit_full        = fit_current,
          fit_red         = fit_red_obj$fit,
          design          = KR_subset_svy,
          rel_thresh_mean = rel_mean_th,
          rel_thresh_max  = rel_max_th
        )
        conf_flag <- conf_obj$flag
        
        step3_log <- dplyr::bind_rows(
          step3_log,
          tibble::tibble(
            var       = v,
            phase     = "backward",
            component = "main",
            p         = p_v,
            dQIC      = dQIC,
            mean_rel  = conf_obj$mean_rel,
            max_rel   = conf_obj$max_rel,
            conf_flag = conf_flag
          )
        )
        
        cand_main <- dplyr::bind_rows(
          cand_main,
          tibble::tibble(
            var       = v,
            p         = p_v,
            dQIC      = dQIC,
            mean_rel  = conf_obj$mean_rel,
            max_rel   = conf_obj$max_rel,
            conf_flag = conf_flag
          )
        )
      }
      
      ## Backward rule for main:
      ## p > alpha_main AND dQIC <= qic_tol AND not confounder
      eligible_main <- cand_main %>%
        dplyr::filter(
          !is.na(p),
          p > alpha_main,
          dQIC <= qic_tol,
          !conf_flag
        )
      
      if (nrow(eligible_main) > 0L) {
        chosen_main <- eligible_main %>%
          dplyr::arrange(dQIC, dplyr::desc(p), var) %>%
          dplyr::slice(1L)
        v_drop <- chosen_main$var
        
        mains_red <- setdiff(main_terms_current, v_drop)
        fit_red_obj <- refit_with_sets(
          main_vars        = mains_red,
          interaction_vars = interaction_vars_current,
          design           = KR_subset_svy
        )
        
        main_terms_current <- mains_red
        fit_current        <- fit_red_obj$fit
        qic_current        <- fit_red_obj$qic
        drop_happened      <- TRUE
      }
    }
    
    if (!drop_happened) break
  }
  
  ## ----------------------
  ## Forward Step 1: re-test dropped vars
  ## ----------------------
  
  step1_keepers_collapsed <- step1_results$step1_keepers_collapsed
  step1_all_candidates    <- step1_results$all_candidates
  
  step1_dropped_vars <- setdiff(step1_all_candidates, c(step1_keepers_collapsed, "age_24"))
  forward_pool <- step1_dropped_vars
  
  while (length(forward_pool)) {
    cand_info  <- tibble::tibble(
      var          = character(),
      conf_global  = logical(),
      dQIC_global  = numeric(),
      p_rank       = numeric(),
      p_main       = numeric(),
      p_int        = numeric()
    )
    cand_models <- list()
    
    for (v in forward_pool) {
      main_try <- union(main_terms_current, v)
      int_try  <- union(interaction_vars_current, v)
      
      ## Fit augmented model with v (main + interaction)
      fit_aug_obj <- refit_with_sets(
        main_vars        = main_try,
        interaction_vars = int_try,
        design           = KR_subset_svy
      )
      fit_aug <- fit_aug_obj$fit
      qic_aug <- fit_aug_obj$qic
      
      ## -----------------------------
      ## Internal pruning: interaction
      ## -----------------------------
      p_block <- NA_real_
      if (v %in% int_try) {
        pat      <- build_int_pattern(v, age_band_var)
        p_block  <- wald_block(fit_aug, pattern = pat)
        
        int_keep_local    <- setdiff(int_try, v)
        fit_red_local_obj <- refit_with_sets(
          main_vars        = main_try,
          interaction_vars = int_keep_local,
          design           = KR_subset_svy
        )
        qic_red_local <- fit_red_local_obj$qic
        dQIC_local    <- qic_red_local - qic_aug
        
        conf_obj_local <- check_confounding_prediction(
          fit_full        = fit_aug,
          fit_red         = fit_red_local_obj$fit,
          design          = KR_subset_svy,
          rel_thresh_mean = rel_mean_th,
          rel_thresh_max  = rel_max_th
        )
        conf_flag_local <- conf_obj_local$flag
        
        step3_log <- dplyr::bind_rows(
          step3_log,
          tibble::tibble(
            var       = v,
            phase     = "forward_step1",
            component = "interaction",
            p         = p_block,
            dQIC      = dQIC_local,
            mean_rel  = conf_obj_local$mean_rel,
            max_rel   = conf_obj_local$max_rel,
            conf_flag = conf_flag_local
          )
        )
        
        ## Drop interaction locally if:
        ## non-significant, no QIC gain, not confounder
        if (!is.na(p_block) && p_block > alpha_int &&
            dQIC_local <= qic_tol &&
            !conf_flag_local) {
          int_try <- int_keep_local
          fit_aug <- fit_red_local_obj$fit
          qic_aug <- fit_red_local_obj$qic
        }
      }
      
      ## -----------------------------
      ## Internal pruning: main effect
      ## -----------------------------
      p_main_rank <- NA_real_
      if (v %in% main_try && !(v %in% always_keep) && !(v %in% int_try)) {
        p_main_rank <- try(
          regTermTest(fit_aug, as.formula(paste("~", v)))$p,
          silent = TRUE
        )
        p_main_rank <- if (inherits(p_main_rank, "try-error") || !length(p_main_rank)) {
          NA_real_
        } else {
          as.numeric(p_main_rank)[1]
        }
        
        mains_red_local <- setdiff(main_try, v)
        fit_red_local_obj <- refit_with_sets(
          main_vars        = mains_red_local,
          interaction_vars = int_try,
          design           = KR_subset_svy
        )
        
        qic_red_local <- fit_red_local_obj$qic
        dQIC_local    <- qic_red_local - qic_aug
        
        conf_obj_local <- check_confounding_prediction(
          fit_full        = fit_aug,
          fit_red         = fit_red_local_obj$fit,
          design          = KR_subset_svy,
          rel_thresh_mean = rel_mean_th,
          rel_thresh_max  = rel_max_th
        )
        conf_flag_local <- conf_obj_local$flag
        
        step3_log <- dplyr::bind_rows(
          step3_log,
          tibble::tibble(
            var       = v,
            phase     = "forward_step1",
            component = "main",
            p         = p_main_rank,
            dQIC      = dQIC_local,
            mean_rel  = conf_obj_local$mean_rel,
            max_rel   = conf_obj_local$max_rel,
            conf_flag = conf_flag_local
          )
        )
        
        ## Drop main locally if:
        ## non-significant, no QIC gain, not confounder
        if (!is.na(p_main_rank) && p_main_rank > alpha_main &&
            dQIC_local <= qic_tol &&
            !conf_flag_local) {
          main_try <- mains_red_local
          fit_aug  <- fit_red_local_obj$fit
          qic_aug  <- fit_red_local_obj$qic
        }
      }
      
      ## If both main and interaction are gone, skip this candidate
      v_survives <- (v %in% main_try) || (v %in% int_try)
      if (!v_survives) {
        cand_info <- dplyr::bind_rows(
          cand_info,
          tibble::tibble(
            var          = v,
            conf_global  = FALSE,
            dQIC_global  = NA_real_,
            p_rank       = NA_real_,
            p_main       = NA_real_,
            p_int        = NA_real_
          )
        )
        next
      }
      
      ## -----------------------------
      ## Global confounding & p-ranks
      ## -----------------------------
      conf_obj_global <- check_confounding_prediction(
        fit_full        = fit_aug,
        fit_red         = fit_current,
        design          = KR_subset_svy,
        rel_thresh_mean = rel_mean_th,
        rel_thresh_max  = rel_max_th
      )
      conf_flag_global <- conf_obj_global$flag
      
      ## p-values for main & interaction (for ranking & “significant” flag)
      if (v %in% main_try && is.na(p_main_rank)) {
        p_main_rank <- try(
          regTermTest(fit_aug, as.formula(paste("~", v)))$p,
          silent = TRUE
        )
        p_main_rank <- if (inherits(p_main_rank, "try-error") || !length(p_main_rank)) {
          NA_real_
        } else {
          as.numeric(p_main_rank)[1]
        }
      }
      
      p_int_rank <- NA_real_
      if (v %in% int_try) {
        pat_rank   <- build_int_pattern(v, age_band_var)
        p_int_rank <- wald_block(fit_aug, pattern = pat_rank)
      }
      
      if (all(is.na(c(p_main_rank, p_int_rank)))) {
        p_rank <- NA_real_
      } else {
        p_rank <- min(c(p_main_rank, p_int_rank), na.rm = TRUE)
      }
      
      dQIC_global <- qic_aug - qic_current
      
      cand_info <- dplyr::bind_rows(
        cand_info,
        tibble::tibble(
          var          = v,
          conf_global  = conf_flag_global,
          dQIC_global  = dQIC_global,
          p_rank       = p_rank,
          p_main       = p_main_rank,
          p_int        = p_int_rank
        )
      )
      
      cand_models[[v]] <- list(
        main_try = main_try,
        int_try  = int_try,
        fit_aug  = fit_aug,
        qic_aug  = qic_aug
      )
    } 
    if (!nrow(cand_info)) break
    
    cand_info <- cand_info %>%
      dplyr::mutate(
        sig_main    = !is.na(p_main) & (p_main <= alpha_main),
        sig_int     = !is.na(p_int)  & (p_int  <= alpha_int),
        sig_any     = sig_main | sig_int,
        qic_improve = !is.na(dQIC_global) & (dQIC_global < -qic_tol),
        eligible    = sig_any | qic_improve | conf_global
      )
    
    eligible <- cand_info %>% dplyr::filter(eligible)
    
    if (!nrow(eligible)) break
    
    ## Rank among eligible: 
    ## 1) most QIC improvement (smallest dQIC_global)
    ## 2) then lowest p_rank
    ## 3) then alphabetical
    rank_df <- eligible %>%
      dplyr::mutate(p_rank_tmp = ifelse(is.na(p_rank), 1, p_rank)) %>%
      dplyr::arrange(dQIC_global, p_rank_tmp, var)
    
    chosen   <- rank_df %>% dplyr::slice(1L)
    v_star   <- chosen$var
    model_star <- cand_models[[v_star]]
    
    ## Adopt chosen augmented model
    main_terms_current       <- model_star$main_try
    interaction_vars_current <- model_star$int_try
    fit_current              <- model_star$fit_aug
    qic_current              <- model_star$qic_aug
    
    forward_pool <- setdiff(forward_pool, v_star)
  }
  
  ## ----------------------
  ## Forward: ageX for age-common vars
  ## ----------------------
  
  ageX_candidates_all <- setdiff(age_common_step3, interaction_vars_current)
  ageX_pool <- ageX_candidates_all
  
  while (length(ageX_pool)) {
    cand_info  <- tibble::tibble(
      var          = character(),
      conf_global  = logical(),
      dQIC_global  = numeric(),
      p_rank       = numeric()
    )
    cand_models <- list()
    
    for (v in ageX_pool) {
      main_try <- main_terms_current
      int_try  <- union(interaction_vars_current, v)
      
      fit_aug_obj <- refit_with_sets(
        main_vars        = main_try,
        interaction_vars = int_try,
        design           = KR_subset_svy
      )
      fit_aug <- fit_aug_obj$fit
      qic_aug <- fit_aug_obj$qic
      
      ## Interaction block p & local pruning
      pat     <- build_int_pattern(v, age_band_var)
      p_block <- wald_block(fit_aug, pattern = pat)
      
      int_keep_local    <- setdiff(int_try, v)
      fit_red_local_obj <- refit_with_sets(
        main_vars        = main_try,
        interaction_vars = int_keep_local,
        design           = KR_subset_svy
      )
      qic_red_local <- fit_red_local_obj$qic
      dQIC_local    <- qic_red_local - qic_aug
      
      conf_obj_local <- check_confounding_prediction(
        fit_full        = fit_aug,
        fit_red         = fit_red_local_obj$fit,
        design          = KR_subset_svy,
        rel_thresh_mean = rel_mean_th,
        rel_thresh_max  = rel_max_th
      )
      conf_flag_local <- conf_obj_local$flag
      
      step3_log <- dplyr::bind_rows(
        step3_log,
        tibble::tibble(
          var       = v,
          phase     = "forward_ageX",
          component = "interaction",
          p         = p_block,
          dQIC      = dQIC_local,
          mean_rel  = conf_obj_local$mean_rel,
          max_rel   = conf_obj_local$max_rel,
          conf_flag = conf_flag_local
        )
      )
      
      ## Local pruning for interaction
      if (!is.na(p_block) && p_block > alpha_int &&
          dQIC_local <= qic_tol &&
          !conf_flag_local) {
        int_try <- int_keep_local
        fit_aug <- fit_red_local_obj$fit
        qic_aug <- fit_red_local_obj$qic
      }
      
      ## If interaction dropped locally, skip this candidate
      v_survives <- v %in% int_try
      if (!v_survives) {
        cand_info <- dplyr::bind_rows(
          cand_info,
          tibble::tibble(
            var          = v,
            conf_global  = FALSE,
            dQIC_global  = NA_real_,
            p_rank       = NA_real_
          )
        )
        next
      }
      
      ## Global confounding for interaction candidate
      conf_obj_global <- check_confounding_prediction(
        fit_full        = fit_aug,
        fit_red         = fit_current,
        design          = KR_subset_svy,
        rel_thresh_mean = rel_mean_th,
        rel_thresh_max  = rel_max_th
      )
      conf_flag_global <- conf_obj_global$flag
      
      p_rank <- wald_block(fit_aug, pattern = pat)
      dQIC_global <- qic_aug - qic_current
      
      cand_info <- dplyr::bind_rows(
        cand_info,
        tibble::tibble(
          var          = v,
          conf_global  = conf_flag_global,
          dQIC_global  = dQIC_global,
          p_rank       = p_rank
        )
      )
      
      cand_models[[v]] <- list(
        main_try = main_try,
        int_try  = int_try,
        fit_aug  = fit_aug,
        qic_aug  = qic_aug
      )
    } # end loop over ageX_pool
    
    if (!nrow(cand_info)) break
    
    cand_info <- cand_info %>%
      dplyr::mutate(
        sig_any     = !is.na(p_rank) & (p_rank <= alpha_int),
        qic_improve = !is.na(dQIC_global) & (dQIC_global < -qic_tol),
        eligible    = sig_any | qic_improve | conf_global
      )
    
    eligible <- cand_info %>% dplyr::filter(eligible)
    if (!nrow(eligible)) break
    
    rank_df <- eligible %>%
      dplyr::mutate(p_rank_tmp = ifelse(is.na(p_rank), 1, p_rank)) %>%
      dplyr::arrange(dQIC_global, p_rank_tmp, var)
    
    chosen    <- rank_df %>% dplyr::slice(1L)
    v_star    <- chosen$var
    model_star <- cand_models[[v_star]]
    
    interaction_vars_current <- model_star$int_try
    main_terms_current       <- model_star$main_try
    fit_current              <- model_star$fit_aug
    qic_current              <- model_star$qic_aug
    
    ageX_pool <- setdiff(ageX_pool, v_star)
  }
  
  ## Final model after backward + forward
  fit_final  <- fit_current
  form_final <- formula(fit_final)
  
  step3_results <- list(
    formula              = form_final,
    fit                  = fit_final,
    main_terms           = main_terms_current,
    interaction_vars     = interaction_vars_current,
    qic                  = qic_current,
    age_band_var         = age_band_var,
    step1_keepers        = step1_keepers_collapsed,
    age_varying_step3    = age_varying_step3,
    age_common_step3     = age_common_step3,
    alpha_main           = step3_alpha_main,
    alpha_int            = step3_alpha_int,
    qic_tol              = step3_qic_tol,
    rel_mean_th          = step3_rel_mean_th,
    rel_max_th           = step3_rel_max_th,
    family_type          = family_type
  )
  
  ## =========================================================
  ## FINAL PUBLICATION TABLE (integrated)
  ## (your table-building code, unchanged)
  ## =========================================================
  
  build_final_results_table <- function(step3_results,
                                        variable_names,
                                        variable_code,
                                        response_list,
                                        community_variable_names,
                                        community_variable_code,
                                        community_labels) {
    final_fit   <- step3_results$fit
    final_coef  <- coef(final_fit)
    final_vcov  <- vcov(final_fit)
    coef_names  <- names(final_coef)
    main_terms_current       <- step3_results$main_terms
    interaction_vars_current <- step3_results$interaction_vars
    age_band_var             <- step3_results$age_band_var
    step1_keepers_collapsed  <- step3_results$step1_keepers
    age_varying_step3        <- step3_results$age_varying_step3
    age_common_step3         <- step3_results$age_common_step3
    
    get_main_coef_names <- function(v, coef_names) {
      nms <- grep(paste0("^", v), coef_names, value = TRUE)
      nms[!grepl(":", nms)]
    }
    
    get_int_coef_names <- function(v, coef_names, age_band_var) {
      grep(paste0("^", age_band_var, ".*:", v), coef_names, value = TRUE)
    }
    
    get_coef_stats <- function(term_name) {
      if (!term_name %in% coef_names) {
        return(list(PR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_))
      }
      est <- final_coef[term_name]
      se  <- sqrt(final_vcov[term_name, term_name])
      PR  <- exp(est)
      LCL <- exp(est - 1.96 * se)
      UCL <- exp(est + 1.96 * se)
      z   <- est / se
      p   <- 2 * pnorm(-abs(z))
      list(PR = PR, LCL = LCL, UCL = UCL, p = p)
    }
    
    build_rows_for_var <- function(v) {
      main_coefs <- get_main_coef_names(v, coef_names)
      int_coefs  <- get_int_coef_names(v, coef_names, age_band_var)
      
      passed_step1 <- ifelse(v %in% step1_keepers_collapsed, "Yes", "No")
      in_final     <- v %in% union(main_terms_current, interaction_vars_current)
      has_int      <- v %in% interaction_vars_current
      
      age_type <- dplyr::case_when(
        !in_final               ~ "Not retained in final model",
        has_int                 ~ "Age-specific",
        in_final & !has_int     ~ "Age-common",
        TRUE                    ~ NA_character_
      )
      
      rows <- list()
      
      if (length(main_coefs) <= 1L) {
        level_label <- if (length(main_coefs) == 1L) main_coefs[1] else v
        
        main_stats <- if (length(main_coefs) == 1L) {
          get_coef_stats(main_coefs[1])
        } else {
          list(PR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_)
        }
        
        if (length(int_coefs) == 1L) {
          int_stats <- get_coef_stats(int_coefs[1])
        } else {
          int_stats <- list(PR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_)
        }
        
        rows[[1]] <- tibble::tibble(
          parent_variable = v,
          level_term      = level_label,
          main_PR         = main_stats$PR,
          main_LCL        = main_stats$LCL,
          main_UCL        = main_stats$UCL,
          main_p          = main_stats$p,
          int_PR          = int_stats$PR,
          int_LCL         = int_stats$LCL,
          int_UCL         = int_stats$UCL,
          int_p           = int_stats$p,
          passed_step1    = passed_step1,
          age_pattern     = age_type
        )
      } else {
        for (cm in main_coefs) {
          main_stats <- get_coef_stats(cm)
          int_name <- grep(paste0(":", cm, "$"), int_coefs, value = TRUE)
          if (length(int_name) == 1L) {
            int_stats <- get_coef_stats(int_name)
          } else {
            int_stats <- list(PR = NA_real_, LCL = NA_real_, UCL = NA_real_, p = NA_real_)
          }
          rows[[length(rows) + 1L]] <- tibble::tibble(
            parent_variable = v,
            level_term      = cm,
            main_PR         = main_stats$PR,
            main_LCL        = main_stats$LCL,
            main_UCL        = main_stats$UCL,
            main_p          = main_stats$p,
            int_PR          = int_stats$PR,
            int_LCL         = int_stats$LCL,
            int_UCL         = int_stats$UCL,
            int_p           = int_stats$p,
            passed_step1    = passed_step1,
            age_pattern     = age_type
          )
        }
      }
      dplyr::bind_rows(rows)
    }
    
    vars_final <- union(main_terms_current, interaction_vars_current)
    
    var_level_summary <- purrr::map_dfr(vars_final, build_rows_for_var)
    
    var_level_summary2 <- var_level_summary %>%
      dplyr::mutate(
        variable_code = gsub("_scaled$", "", as.character(parent_variable)),
        response      = suppressWarnings(
          as.numeric(gsub("^.*?([0-9]+)$", "\\1", as.character(level_term)))
        ),
        factor24     = ifelse(is.na(int_PR), 1, int_PR),
        factor24_LCL = ifelse(is.na(int_LCL), 1, int_LCL),
        factor24_UCL = ifelse(is.na(int_UCL), 1, int_UCL),
        PR_6_23   = main_PR,
        LCL_6_23  = main_LCL,
        UCL_6_23  = main_UCL,
        PR_24_59  = main_PR * factor24,
        LCL_24_59 = main_LCL * factor24_LCL,
        UCL_24_59 = main_UCL * factor24_UCL,
        PR_6_23_str = ifelse(
          is.na(PR_6_23),
          NA_character_,
          sprintf("%.3f (%.3f–%.3f)", PR_6_23, LCL_6_23, UCL_6_23)
        ),
        PR_24_59_str = ifelse(
          is.na(PR_24_59),
          NA_character_,
          sprintf("%.3f (%.3f–%.3f)", PR_24_59, LCL_24_59, UCL_24_59)
        )
      )
    
    age_pattern_by_parent <- var_level_summary2 %>%
      dplyr::select(parent_variable, age_pattern) %>%
      dplyr::distinct() %>%
      tibble::deframe()
    
    ## --------------------------
    ## Community block
    ## --------------------------
    
    parents_final <- unique(var_level_summary2$parent_variable)
    community_parents_order <- paste0(
      community_variable_code[!is.na(community_variable_code)],
      "_scaled"
    )
    community_parents_final <- intersect(community_parents_order, parents_final)
    
    comm_domain_for_parent <- character()
    current_dom <- NA_character_
    
    for (i in seq_along(community_variable_names)) {
      comm_label <- community_variable_names[i]
      comm_code  <- community_variable_code[i]
      
      if (is.na(comm_code)) {
        current_dom <- comm_label
      } else {
        parent_scaled <- paste0(comm_code, "_scaled")
        if (!is.na(current_dom)) {
          comm_domain_for_parent[parent_scaled] <- current_dom
        }
      }
    }
    
    community_domain_order <- unique(
      community_variable_names[is.na(community_variable_code)]
    )
    
    final_rows <- list()
    
    if (length(community_parents_final) > 0L) {
      final_rows[[length(final_rows) + 1L]] <- tibble::tibble(
        Variable             = "Community-level factors (per 10%-point increase)",
        `Age pattern`        = NA_character_,
        `PR (6–23.9 months)` = NA_character_,
        `PR (24–59 months)`  = NA_character_
      )
      
      for (dom in community_domain_order) {
        dom_parents <- community_parents_final[
          comm_domain_for_parent[community_parents_final] == dom
        ]
        dom_parents <- unique(dom_parents)
        if (!length(dom_parents)) next
        
        final_rows[[length(final_rows) + 1L]] <- tibble::tibble(
          Variable             = dom,
          `Age pattern`        = NA_character_,
          `PR (6–23.9 months)` = NA_character_,
          `PR (24–59 months)`  = NA_character_
        )
        
        dom_parents_ordered <- community_parents_order[
          community_parents_order %in% dom_parents
        ]
        
        for (cp in dom_parents_ordered) {
          vrow <- var_level_summary2 %>%
            dplyr::filter(parent_variable == cp)
          if (!nrow(vrow)) next
          
          base_code <- gsub("_scaled$", "", cp)
          
          if (!is.null(community_labels)) {
            lbl <- community_block$Variable[
              community_block$variable_code == base_code & 
                !is.na(community_block$variable_code)
            ]
            base_label <- lbl
          } else {
            base_label <- base_code
          }
          
          age_pat <- age_pattern_by_parent[[cp]]
          
          final_rows[[length(final_rows) + 1L]] <- tibble::tibble(
            Variable             = base_label,
            `Age pattern`        = age_pat,
            `PR (6–23.9 months)` = vrow$PR_6_23_str[1],
            `PR (24–59 months)`  = vrow$PR_24_59_str[1]
          )
        }
      }
    }
    
    ## --------------------------
    ## Individual-level block
    ## --------------------------
    
    parent_to_code <- var_level_summary2 %>%
      dplyr::select(parent_variable, variable_code) %>%
      dplyr::distinct() %>%
      tibble::deframe()
    
    header_codes <- variable_code[!is.na(variable_code) & is.na(response_list)] |> unique()
    
    labels_tbl <- table1 %>%   # assumes table1 exists in outer scope
      dplyr::mutate(order_index = dplyr::row_number()) %>%
      dplyr::filter(!is.na(variable_code), !is.na(response)) %>%
      dplyr::transmute(
        variable_code,
        response    = as.numeric(response),
        label       = Variable,
        is_ref      = response == reference,
        order_index = order_index
      ) %>%
      dplyr::distinct()
    
    final_rows[[length(final_rows) + 1L]] <- tibble::tibble(
      Variable             = "Individual-level factors",
      `Age pattern`        = NA_character_,
      `PR (6–23.9 months)` = NA_character_,
      `PR (24–59 months)`  = NA_character_
    )
    
    for (i in seq_along(variable_names)) {
      code  <- variable_code[i]
      resp  <- response_list[i]
      label <- variable_names[i]
      
      if (is.na(code)) {
        final_rows[[length(final_rows) + 1L]] <- tibble::tibble(
          Variable             = label,
          `Age pattern`        = NA_character_,
          `PR (6–23.9 months)` = NA_character_,
          `PR (24–59 months)`  = NA_character_
        )
        next
      }
      
      if (!is.na(code) && is.na(resp)) {
        parent <- names(parent_to_code)[parent_to_code == code][1]
        if (is.na(parent) || !(parent %in% vars_final)) next
        
        age_pat <- age_pattern_by_parent[[parent]]
        final_rows[[length(final_rows) + 1L]] <- tibble::tibble(
          Variable             = label,
          `Age pattern`        = age_pat,
          `PR (6–23.9 months)` = NA_character_,
          `PR (24–59 months)`  = NA_character_
        )
        
        labs <- labels_tbl %>%
          dplyr::filter(variable_code == code) %>%
          dplyr::arrange(order_index)
        
        for (j in seq_len(nrow(labs))) {
          lj <- labs[j, ]
          if (isTRUE(lj$is_ref)) {
            final_rows[[length(final_rows) + 1L]] <- tibble::tibble(
              Variable             = lj$label,
              `Age pattern`        = NA_character_,
              `PR (6–23.9 months)` = "Reference",
              `PR (24–59 months)`  = "Reference"
            )
          } else {
            vrow <- var_level_summary2 %>%
              dplyr::filter(
                parent_variable == parent,
                variable_code   == code,
                response        == lj$response
              )
            if (!nrow(vrow)) next
            
            final_rows[[length(final_rows) + 1L]] <- tibble::tibble(
              Variable             = lj$label,
              `Age pattern`        = NA_character_,
              `PR (6–23.9 months)` = vrow$PR_6_23_str[1],
              `PR (24–59 months)`  = vrow$PR_24_59_str[1]
            )
          }
        }
        next
      }
      
      if (!is.na(code) && !is.na(resp) && !(code %in% header_codes)) {
        parent <- names(parent_to_code)[parent_to_code == code][1]
        if (is.na(parent) || !(parent %in% vars_final)) next
        
        age_pat <- age_pattern_by_parent[[parent]]
        vrow <- var_level_summary2 %>%
          dplyr::filter(parent_variable == parent)
        
        if ("response" %in% names(vrow) && any(!is.na(vrow$response))) {
          vrow <- vrow %>% dplyr::filter(response == resp)
        }
        if (!nrow(vrow)) next
        
        final_rows[[length(final_rows) + 1L]] <- tibble::tibble(
          Variable             = label,
          `Age pattern`        = age_pat,
          `PR (6–23.9 months)` = vrow$PR_6_23_str[1],
          `PR (24–59 months)`  = vrow$PR_24_59_str[1]
        )
      }
    }
    
    dplyr::bind_rows(final_rows)
  }
  
  final_results_table <- build_final_results_table(
    step3_results             = step3_results,
    variable_names            = variable_names,
    variable_code             = variable_code,
    response_list             = response_list,
    community_variable_names  = community_variable_names,
    community_variable_code   = community_variable_code,
    community_labels          = community_labels
  )
  
  ## =========================================================
  ## Supplementary selection table
  ## =========================================================
  
  header_codes <- variable_code[!is.na(variable_code) & is.na(response_list)] |> unique()
  
  summarise_parent_selection <- function(v) {
    in_final_main <- v %in% step3_results$main_terms
    in_final_int  <- v %in% step3_results$interaction_vars
    
    age_pattern <- dplyr::case_when(
      in_final_int & in_final_main ~ "Age-varying (main + age_band×X in final model)",
      !in_final_int & in_final_main ~ "Age-common (no age_band×X in final model)",
      TRUE ~ "Not retained in final model"
    )
    
    bm <- step3_log %>%
      dplyr::filter(var == v, phase == "backward", component == "main")
    bi <- step3_log %>%
      dplyr::filter(var == v, phase == "backward", component == "interaction")
    fm <- step3_log %>%
      dplyr::filter(var == v, phase == "forward_step1", component == "main")
    fi <- step3_log %>%
      dplyr::filter(var == v & phase %in% c("forward_step1", "forward_ageX"),
                    component == "interaction")
    
    backward_main_result_step3 <- dplyr::case_when(
      in_final_main & nrow(bm) > 0 ~ "Main effect retained after backward selection",
      !in_final_main & nrow(bm) > 0 ~ "Main effect removed in backward selection",
      TRUE ~ NA_character_
    )
    backward_int_result_step3 <- dplyr::case_when(
      in_final_int & nrow(bi) > 0 ~ "Interaction retained after backward selection",
      !in_final_int & nrow(bi) > 0 ~ "Interaction removed in backward selection",
      TRUE ~ NA_character_
    )
    forward_main_result_step3 <- dplyr::case_when(
      in_final_main & nrow(fm) > 0 ~ "Main effect added in forward selection",
      TRUE ~ NA_character_
    )
    forward_int_result_step3 <- dplyr::case_when(
      in_final_int & nrow(fi) > 0 ~ "Interaction added in forward selection",
      TRUE ~ NA_character_
    )
    
    sel_freq <- NA_real_
    sf <- step1_results$selection_freq
    if (!is.null(sf) && nrow(sf) && all(c("Var1", "Freq") %in% names(sf))) {
      hit <- sf$Freq[sf$Var1 == v]
      if (length(hit)) sel_freq <- as.numeric(hit[1])
    }
    if (is.na(sel_freq) && v %in% step1_results$all_candidates) {
      sel_freq <- 0
    }
    
    if (v %in% always_keep) {
      step1_status <- "Always-keep (DAG)"
    } else if (v %in% step1_results$step1_keepers_collapsed) {
      if (!is.na(sel_freq)) {
        step1_status <- sprintf("Retained in Step 1 (stability = %.2f)", sel_freq)
      } else {
        step1_status <- "Retained in Step 1"
      }
    } else if (v %in% step1_results$all_candidates) {
      if (!is.na(sel_freq)) {
        step1_status <- sprintf("Dropped in Step 1 (stability = %.2f)", sel_freq)
      } else {
        step1_status <- "Dropped in Step 1"
      }
    } else {
      step1_status <- "Not in Step 1 candidate set"
    }
    
    stab_tbl <- step2_results$stability_table
    stab_val <- NA_real_
    if (!is.null(stab_tbl) && nrow(stab_tbl)) {
      hit2 <- stab_tbl$stability[stab_tbl$variable == v]
      if (length(hit2)) stab_val <- as.numeric(hit2[1])
    }
    if (is.na(stab_val) && (v %in% step2_results$age_varying_candidates_stable ||
                            v %in% step2_results$age_common_candidates)) {
      stab_val <- 0
    }
    
    step2_status <- dplyr::case_when(
      v %in% step2_results$age_varying_candidates_stable & !is.na(stab_val) ~
        sprintf("Age-varying candidate (stability = %.2f)", stab_val),
      v %in% step2_results$age_varying_candidates_stable ~
        "Age-varying candidate",
      v %in% step2_results$age_common_candidates & !is.na(stab_val) ~
        sprintf("Age-common candidate (stability = %.2f)", stab_val),
      v %in% step2_results$age_common_candidates ~
        "Age-common candidate",
      v %in% step1_results$step1_keepers_collapsed ~
        "Not evaluated in Step 2",
      TRUE ~ NA_character_
    )
    
    tibble::tibble(
      parent_variable             = v,
      age_pattern                 = age_pattern,
      backward_main_result_step3  = backward_main_result_step3,
      backward_int_result_step3   = backward_int_result_step3,
      forward_main_result_step3   = forward_main_result_step3,
      forward_int_result_step3    = forward_int_result_step3,
      step1_status                = step1_status,
      step2_status                = step2_status
    )
  }
  
  selection_rows <- list()
  
  selection_rows[[length(selection_rows) + 1L]] <- tibble::tibble(
    Variable_label              = "Community-level factors (per 10%-point increase)",
    parent_variable             = NA_character_,
    age_pattern                 = NA_character_,
    step1_status                = NA_character_,
    step2_status                = NA_character_,
    backward_main_result_step3  = NA_character_,
    backward_int_result_step3   = NA_character_,
    forward_main_result_step3   = NA_character_,
    forward_int_result_step3    = NA_character_
  )
  
  for (i in seq_along(community_variable_names)) {
    comm_label <- community_variable_names[i]
    comm_code  <- community_variable_code[i]
    
    if (is.na(comm_code)) {
      selection_rows[[length(selection_rows) + 1L]] <- tibble::tibble(
        Variable_label              = comm_label,
        parent_variable             = NA_character_,
        age_pattern                 = NA_character_,
        step1_status                = NA_character_,
        step2_status                = NA_character_,
        backward_main_result_step3  = NA_character_,
        backward_int_result_step3   = NA_character_,
        forward_main_result_step3   = NA_character_,
        forward_int_result_step3    = NA_character_
      )
      next
    }
    
    parent_scaled <- paste0(comm_code, "_scaled")
    sel_info <- summarise_parent_selection(parent_scaled)
    
    base_label <- if (!is.null(community_labels)) {
      lbl <- community_block$Variable[
        community_block$variable_code == comm_code & 
          !is.na(community_block$variable_code)
      ]
      if (is.null(lbl)) comm_label else lbl
    } else {
      comm_label
    }
    
    selection_rows[[length(selection_rows) + 1L]] <- sel_info %>%
      dplyr::mutate(Variable_label = base_label) %>%
      dplyr::select(
        Variable_label, parent_variable,
        age_pattern, step1_status, step2_status,
        backward_main_result_step3,
        backward_int_result_step3,
        forward_main_result_step3,
        forward_int_result_step3
      )
  }
  
  selection_rows[[length(selection_rows) + 1L]] <- tibble::tibble(
    Variable_label              = "Individual-level factors",
    parent_variable             = NA_character_,
    age_pattern                 = NA_character_,
    step1_status                = NA_character_,
    step2_status                = NA_character_,
    backward_main_result_step3  = NA_character_,
    backward_int_result_step3   = NA_character_,
    forward_main_result_step3   = NA_character_,
    forward_int_result_step3    = NA_character_
  )
  
  for (i in seq_along(variable_names)) {
    code  <- variable_code[i]
    resp  <- response_list[i]
    label <- variable_names[i]
    
    if (is.na(code)) {
      selection_rows[[length(selection_rows) + 1L]] <- tibble::tibble(
        Variable_label              = label,
        parent_variable             = NA_character_,
        age_pattern                 = NA_character_,
        step1_status                = NA_character_,
        step2_status                = NA_character_,
        backward_main_result_step3  = NA_character_,
        backward_int_result_step3   = NA_character_,
        forward_main_result_step3   = NA_character_,
        forward_int_result_step3    = NA_character_
      )
      next
    }
    
    if (!is.na(code) && is.na(resp)) {
      v <- code
      sel_info <- summarise_parent_selection(v)
      selection_rows[[length(selection_rows) + 1L]] <- sel_info %>%
        dplyr::mutate(Variable_label = label) %>%
        dplyr::select(
          Variable_label, parent_variable,
          age_pattern, step1_status, step2_status,
          backward_main_result_step3,
          backward_int_result_step3,
          forward_main_result_step3,
          forward_int_result_step3
        )
      next
    }
    
    if (!is.na(code) && !is.na(resp) && !(code %in% header_codes)) {
      v <- code
      sel_info <- summarise_parent_selection(v)
      selection_rows[[length(selection_rows) + 1L]] <- sel_info %>%
        dplyr::mutate(Variable_label = label) %>%
        dplyr::select(
          Variable_label, parent_variable,
          age_pattern, step1_status, step2_status,
          backward_main_result_step3,
          backward_int_result_step3,
          forward_main_result_step3,
          forward_int_result_step3
        )
    }
  }
  
  var_selection_table_all <- dplyr::bind_rows(selection_rows)
  
  ## =========================================================
  ## Final return object
  ## =========================================================
  
  res <- list(
    step1_results           = step1_results,
    step2_results           = step2_results,
    step3_results           = step3_results,
    var_selection_table_all = var_selection_table_all,
    final_results_table     = final_results_table,
    step3_log               = step3_log
  )
  attr(res, "config_name") <- make_config_name()
  res
}

## --------------------------------------------------------
## Runner and config grid (now also tuning over family_type)
## --------------------------------------------------------

run_configs <- function(config_list) {
  out <- vector("list", length(config_list))
  for (i in seq_along(config_list)) {
    th <- config_list[[i]]
    res <- do.call(
      run_3step_with_thresholds,
      c(
        list(
          KR_subset                = KR_subset,
          variable_names           = variable_names,
          variable_code            = variable_code,
          response_list            = response_list,
          community_variable_names = community_variable_names,
          community_variable_code  = community_variable_code,
          community_labels         = community_labels
        ),
        th
      )
    )
    nm <- attr(res, "config_name")
    out[[i]] <- res
    names(out)[i] <- nm
  }
  out
}

## --------------------------------------------------------
## Build the configs under your rules (× family_type)
## --------------------------------------------------------

if (RUN_SETTINGS$full_tuning_grid) {
  step1_stab_vals <- c(0.30, 0.60)
  lambda_types <- c("lambda.min", "lambda.1se")
  qic_vals <- c(-2, 0, 2)
  conf_pairs <- list(c(0.10, 0.20), c(0.20, 0.40))

  ## (step2_p_max, step3_alpha_main = step3_alpha_int)
  p_tuples <- list(
    list(step2_p_max = 0.10, alpha = 0.05),
    list(step2_p_max = 0.15, alpha = 0.10),
    list(step2_p_max = 0.15, alpha = 0.05),
    list(step2_p_max = 0.20, alpha = 0.15),
    list(step2_p_max = 0.20, alpha = 0.10),
    list(step2_p_max = 0.20, alpha = 0.05)
  )
} else {
  # Minimal grid used only by the local smoke-test script.
  step1_stab_vals <- 0.30
  lambda_types <- "lambda.min"
  qic_vals <- 2
  conf_pairs <- list(c(0.20, 0.40))
  p_tuples <- list(list(step2_p_max = 0.20, alpha = 0.05))
}

family_types <- c("quasipoisson")

configs <- list()
idx <- 1L

for (stab in step1_stab_vals) {
  for (lam in lambda_types) {
    for (pt in p_tuples) {
      for (qic in qic_vals) {
        for (cp in conf_pairs) {
          for (fam in family_types) {
            configs[[idx]] <- list(
              step1_stab_cut     = stab,
              step2_stab_cut     = stab,                 # in tandem with step1
              step1_lambda_type  = lam,
              step2_p_max        = pt$step2_p_max,
              step3_alpha_main   = pt$alpha,
              step3_alpha_int    = pt$alpha,             # always equal in a given run
              step3_qic_tol      = qic,
              step3_rel_mean_th  = cp[1],
              step3_rel_max_th   = cp[2],
              family_type        = fam,
              step1_B            = RUN_SETTINGS$step1_B,
              step2_B            = RUN_SETTINGS$step2_B,
              random_seed        = RANDOM_SEED
            )
            idx <- idx + 1L
          }
        }
      }
    }
  }
}

if (RUN_SETTINGS$full_tuning_grid && length(configs) != 144L) {
  stop("The full tuning grid must contain 144 configurations.")
}

sens_results <- run_configs(configs)
sens_results_save_light <- lapply(
  sens_results,
  function(x) {
    list(
      var_selection_table_all = x$var_selection_table_all,
      final_results_table     = x$final_results_table
    )
  }
)

saveRDS(
  sens_results_save_light,
  file = file.path(RESULT_DIR, "cms_tuning_results_light.rds"),
  compress = "gzip"
)
saveRDS(
  sens_results,
  file = file.path(RESULT_DIR, "cms_tuning_results_full.rds"),
  compress = "gzip"
)
## -------------------------------
## 1) PSU-stratified CV folds
## -------------------------------
make_psu_folds_stratified <- function(psu, strata, K = 10) {
  foldid <- rep(NA_integer_, length(psu))
  strata_levels <- unique(strata)
  
  for (s in strata_levels) {
    idx_s   <- which(strata == s)
    psu_s   <- psu[idx_s]
    psu_ids <- unique(psu_s)
    
    if (length(psu_ids) == 1L) {
      # only 1 PSU in this stratum: put everything in fold 1
      folds_s <- 1L
      names(folds_s) <- psu_ids
    } else {
      K_s <- min(K, length(psu_ids))
      folds_s <- sample(rep(1:K_s, length.out = length(psu_ids)))
      names(folds_s) <- psu_ids
    }
    
    foldid[idx_s] <- folds_s[as.character(psu_s)]
  }
  foldid
}

set.seed(RANDOM_SEED + 3000L)
KR_subset$cv_fold <- make_psu_folds_stratified(
  psu    = KR_subset[[id_var]],
  strata = KR_subset[[strata_var]],
  K      = RUN_SETTINGS$cv_folds
)

## -------------------------------
## 2) Weighted Brier score
## -------------------------------
brier_weighted <- function(p_hat, y, w) {
  num <- sum(w * (p_hat - y)^2, na.rm = TRUE)
  den <- sum(w, na.rm = TRUE)
  if (den == 0) return(NA_real_)
  num / den
}

## -------------------------------
## 3) CV Brier for one final model
##    (uses step3_results$formula)
## -------------------------------

cv_brier_for_result <- function(res,
                                KR_data = KR_subset,
                                K       = RUN_SETTINGS$cv_folds,
                                verbose = FALSE) {
  step3 <- res$step3_results
  
  form_final  <- step3$formula
  family_type <- step3$family_type
  
  fam_glm <- switch(
    family_type,
    "quasipoisson"  = quasipoisson("log"),
    "quasibinomial" = quasibinomial("log"),
    stop("Unknown family_type: ", family_type)
  )
  
  folds <- sort(unique(KR_data$cv_fold))
  folds <- folds[folds %in% seq_len(K)]
  
  brier_per_fold <- numeric(length(folds))
  brier_per_fold[] <- NA_real_
  
  for (i in seq_along(folds)) {
    k <- folds[i]
    if (verbose) message("  Fold ", k, " of ", max(folds))
    
    train_idx <- which(KR_data$cv_fold != k)
    test_idx  <- which(KR_data$cv_fold == k)
    
    train_dat <- KR_data[train_idx, , drop = FALSE]
    test_dat  <- KR_data[test_idx,  , drop = FALSE]
    
    ## 1) survey design on training data
    design_train <- tryCatch(
      svydesign(
        ids     = as.formula(paste0("~", id_var)),
        strata  = as.formula(paste0("~", strata_var)),
        weights = as.formula(paste0("~I(", weight_var, "/1e6)")),
        data    = train_dat,
        nest    = TRUE
      ),
      error = function(e) {
        warning("svydesign failed in fold ", k, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(design_train)) next
    
    ## 2) fit same formula on training data
    fit_k <- tryCatch(
      svyglm(form_final, design = design_train, family = fam_glm),
      error = function(e) {
        warning("svyglm failed in fold ", k, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(fit_k)) next
    
    ## ---------------------------
    ## 3) Build aligned test design matrix
    ## ---------------------------
    
    # Add row id so we can subset y & w consistently
    test_dat2 <- test_dat
    test_dat2$.row_id <- seq_len(nrow(test_dat2))
    
    # Model frame with .row_id so we only keep complete cases
    form_with_id <- update(form_final, . ~ . + .row_id)
    mf_test <- tryCatch(
      model.frame(form_with_id, data = test_dat2, na.action = na.omit),
      error = function(e) {
        warning("model.frame failed in fold ", k, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(mf_test) || !nrow(mf_test)) next
    
    used_ids  <- mf_test$.row_id
    test_used <- test_dat2[match(used_ids, test_dat2$.row_id), , drop = FALSE]
    
    # Design matrix for predictors only
    X <- tryCatch(
      model.matrix(delete.response(terms(fit_k)), data = mf_test),
      error = function(e) {
        warning("model.matrix failed in fold ", k, ": ", conditionMessage(e))
        return(NULL)
      }
    )
    if (is.null(X) || !nrow(X)) next
    
    beta <- coef(fit_k)
    
    # Add any missing columns (all zeros)
    missing_cols <- setdiff(names(beta), colnames(X))
    if (length(missing_cols) > 0L) {
      X_missing <- matrix(0, nrow = nrow(X), ncol = length(missing_cols))
      colnames(X_missing) <- missing_cols
      X <- cbind(X, X_missing)
    }
    
    # Drop any extra columns not in beta
    extra_cols <- setdiff(colnames(X), names(beta))
    if (length(extra_cols) > 0L) {
      X <- X[, !(colnames(X) %in% extra_cols), drop = FALSE]
    }
    
    # Ensure same order as coefficients
    X <- X[, names(beta), drop = FALSE]
    
    # Linear predictor & inverse link
    eta   <- as.numeric(X %*% beta)
    p_hat <- fit_k$family$linkinv(eta)
    
    ## ---------------------------
    ## 4) Brier on aligned rows
    ## ---------------------------
    y <- test_used[[y_var]]
    w <- test_used[[weight_var]] / 1e6
    
    if (!(length(p_hat) == length(y) && length(y) == length(w))) {
      warning(
        "Length mismatch in fold ", k,
        ": p=", length(p_hat),
        " y=", length(y),
        " w=", length(w),
        ". Brier set to NA for this fold."
      )
      next
    }
    
    brier_per_fold[i] <- brier_weighted(p_hat, y, w)
  }
  
  mean(brier_per_fold, na.rm = TRUE)
}

## -------------------------------
## 4) Loop over sens_results
## -------------------------------

cv_brier_results <- map_dfr(seq_along(sens_results), function(i) {
  res_i <- sens_results[[i]]
  
  config_name <- attr(res_i, "config_name")
  if (is.null(config_name) && !is.null(names(sens_results)[i])) {
    config_name <- names(sens_results)[i]
  }
  
  if (is.null(config_name)) {
    config_name <- paste0("config_", i)
  }
  
  brier_cv <- cv_brier_for_result(
    res    = res_i,
    KR_data = KR_subset,
    K      = RUN_SETTINGS$cv_folds,
    verbose = FALSE
  )
  
  tibble(
    config_index = i,
    config_name  = config_name,
    cv_brier     = brier_cv
  )
})

## Rank: smaller Brier = better
cv_brier_results <- cv_brier_results %>%
  arrange(cv_brier)

best_brier <- min(cv_brier_results$cv_brier, na.rm = TRUE)

rel_tol <- 0.001  # first 0.1% Brier-score bin above the minimum
cv_brier_results_filter <- cv_brier_results %>%
  filter(cv_brier <= best_brier * (1 + rel_tol))

# 1) Start from your filtered CV results (within ΔBrier window)
configs_ok <- cv_brier_results_filter$config_name

## --------------------
## Tuning model summary
## --------------------

compute_model_metrics <- function(res, name = NA_character_) {
  step3    <- res$step3_results
  fit      <- step3$fit
  qic_val  <- step3$qic
  
  ## ---- Pearson dispersion ----
  pearson_res  <- residuals(fit, type = "pearson")
  df_res       <- fit$df.residual
  pearson_disp <- sum(pearson_res^2, na.rm = TRUE) / df_res
  
  ## ---- Calibration correlation (by deciles of fitted) ----
  p_hat <- as.numeric(fitted(fit))
  mf    <- model.frame(fit)
  y     <- model.response(mf)
  w     <- stats::weights(fit)
  if (is.null(w)) w <- rep(1, length(y))
  
  q_breaks <- quantile(p_hat, probs = seq(0, 1, length.out = 11), na.rm = TRUE)
  q_breaks <- unique(q_breaks)
  
  if (length(q_breaks) < 3L) {
    calib_corr <- NA_real_
  } else {
    dec   <- cut(p_hat, breaks = q_breaks, include.lowest = TRUE, labels = FALSE)
    valid <- !is.na(dec) & !is.na(y) & !is.na(p_hat) & !is.na(w)
    
    if (!any(valid)) {
      calib_corr <- NA_real_
    } else {
      dec_u    <- sort(unique(dec[valid]))
      mean_obs <- mean_fit <- numeric(length(dec_u))
      
      for (i in seq_along(dec_u)) {
        g   <- dec_u[i]
        idx <- valid & dec == g
        ww  <- w[idx]
        yy  <- y[idx]
        pp  <- p_hat[idx]
        if (sum(ww) > 0) {
          mean_obs[i] <- sum(ww * yy) / sum(ww)
          mean_fit[i] <- sum(ww * pp) / sum(ww)
        } else {
          mean_obs[i] <- NA_real_
          mean_fit[i] <- NA_real_
        }
      }
      
      keep <- !is.na(mean_obs) & !is.na(mean_fit)
      if (sum(keep) >= 3L) {
        calib_corr <- suppressWarnings(cor(mean_obs[keep], mean_fit[keep]))
      } else {
        calib_corr <- NA_real_
      }
    }
  }
  
  ## ---- Survey-design Rao–Scott F via anova(null, fit, test = "F") ----
  rs_stat <- NA_real_  # F statistic
  rs_df   <- NA_real_  # denominator df
  rs_p    <- NA_real_  # numeric-ish p (stored as character here)
  
  # Try to recover the survey design from the fit object
  design_obj <- NULL
  if (!is.null(fit$survey.design)) {
    design_obj <- fit$survey.design
  } else if (!is.null(fit$design)) {
    design_obj <- fit$design
  }
  
  if (is.null(design_obj)) stop("No survey design stored in fit.")
  
  fam    <- fit$family
  y_name <- all.vars(formula(fit))[1]
  
  form_null <- as.formula(paste(y_name, "~ 1"))
  
  null <- svyglm(
    formula = form_null,
    design  = design_obj,
    family  = fam
  )
  
  # test="Chisq" or "F" both print "Working 2logLR" on your system
  aov_txt <- capture.output(anova(null, fit, test = "F"))
  
  line_LR <- aov_txt[grepl("Working 2logLR", aov_txt)]
  line_df <- aov_txt[grepl("denominator df", aov_txt)]
  
  if (length(line_LR) == 0L) stop("No 'Working 2logLR' line found.")
  
  # Parse Working 2logLR
  LR_val <- suppressWarnings(
    as.numeric(sub(".*Working 2logLR *= *([0-9.eE+-]+).*", "\\1", line_LR))
  )
  
  # Parse denominator df
  dfden <- NA_real_
  if (length(line_df) > 0L) {
    dfden <- suppressWarnings(
      as.numeric(sub(".*denominator df *= *([0-9]+).*", "\\1", line_df))
    )
  }
  
  # Parse p-value string (may be "<2.22e-16")
  p_str <- sub(".*p *= *", "", line_LR)
  p_str <- gsub("[ \t]", "", p_str)
  p_str_num <- as.numeric(gsub("<", "", p_str))
  # store as character but clipped for very small values
  if (!is.na(p_str_num) && p_str_num < 10^-6) {
    p_val <- "0.000001"
  } else {
    p_val <- sprintf("%.6f", p_str_num)
  }
  
  rs_stat <- LR_val
  rs_df   <- dfden
  rs_p    <- p_val
  
  ## ---- Number of terms (for parsimony) ----
  # You can choose whether to count intercept or not; here we exclude intercept.
  n_terms <- length(coef(fit)) - 1L
  
  ## ---- Safe config_name handling ----
  cfg_attr <- attr(res, "config_name")
  cfg_name <- if (!missing(name) && length(name) == 1L) {
    as.character(name)
  } else if (!is.null(cfg_attr)) {
    as.character(cfg_attr)
  } else {
    NA_character_
  }
  
  tibble::tibble(
    config_name     = cfg_name,
    pearson_disp    = pearson_disp,
    qic             = qic_val,
    calib_corr      = calib_corr,
    rao_scott_stat  = rs_stat,  # F or 2logLR-like stat
    rao_scott_df    = rs_df,    # denominator df
    rao_scott_p     = rs_p,     # p-value (character)
    n_terms         = n_terms   # for parsimony filter
  )
}

sens_metrics <- purrr::imap_dfr(
  sens_results,
  ~ compute_model_metrics(.x, name = .y)
)

# 2) Join to sens_metrics to get n_terms, QIC, dispersion, etc.
parsity_table <- cv_brier_results_filter %>%
  inner_join(sens_metrics, by = "config_name") %>%
  # keep only the configs in your Brier window (should already be, but explicit)
  filter(config_name %in% configs_ok) %>%
  # arrange by parsimony first, then Brier, then QIC as tie-breakers
  arrange(n_terms, cv_brier, qic)

best_res <- sens_results[[parsity_table$config_name[1]]]
step1_results <- best_res$step1_results
step2_results <- best_res$step2_results
step3_results <- best_res$step3_results
step3_log      <- best_res$step3_log

final_results_table <- best_res$final_results_table
var_selection_table_all <- best_res$var_selection_table_all


###########################################################
## Final publication table (exportable)
###########################################################


###########################################################
## Rebuild var_selection_table_all + Table 2S
##  - Uses final chosen model only (best_res)
##  - Restores detailed Step 3 reasons and forward logs
###########################################################

## 0. Pull final model components from best_res
best_res       <- sens_results[[parsity_table$config_name[1]]]
step1_results  <- best_res$step1_results
step2_results  <- best_res$step2_results
step3_results  <- best_res$step3_results
step3_log      <- best_res$step3_log

final_results_table      <- best_res$final_results_table
var_selection_table_all  <- best_res$var_selection_table_all  # will be overwritten

## Core fit objects
final_fit   <- step3_results$fit
final_coef  <- coef(final_fit)
final_vcov  <- vcov(final_fit)
coef_names  <- names(final_coef)

main_terms_current       <- step3_results$main_terms
interaction_vars_current <- step3_results$interaction_vars
age_band_var             <- step3_results$age_band_var
step1_keepers_collapsed  <- step3_results$step1_keepers
age_varying_step3        <- step3_results$age_varying_step3
age_common_step3         <- step3_results$age_common_step3

## Variables that behaved like "age"
age_like <- c("age", age_band_var, "age_yrs", "age_24")

## Always-keep (DAG) variables
always_keep <- c(
  "ph_wealth_quint",
  "mother_eduyr",
  "geography_regions",
  "comm_pov_prop_scaled",
  "comm_motheredu_prop_scaled"
)
protected_vars <- always_keep

## Initial Step 3 pools (as earlier)
initial_main_terms <- setdiff(step1_results$step1_keepers_collapsed, age_like)
initial_interaction_vars <- intersect(
  step2_results$age_varying_candidates_stable,
  initial_main_terms
)

## Recover Step 3 logs (as before)
back_main_log <- step3_log %>%
  dplyr::filter(phase == "backward", component == "main")

back_int_log  <- step3_log %>%
  dplyr::filter(phase == "backward", component == "interaction")

fwd_main_log  <- step3_log %>%
  dplyr::filter(grepl("^forward", phase), component == "main")

fwd_int_log   <- step3_log %>%
  dplyr::filter(grepl("^forward", phase), component == "interaction")

## ------------------------------------------------------
## Helper: Step 1 / Step 2 candidate sets
## ------------------------------------------------------

# 1) All Step 1 candidates
step1_all_candidates <- if ("all_candidates" %in% names(step1_results)) {
  step1_results$all_candidates
} else if ("selection_freq" %in% names(step1_results)) {
  unique(step1_results$selection_freq$Var1)
} else {
  character(0)
}

# 2) Step 2 age-varying & age-common candidates
age_varying_candidates_stable <- if ("age_varying_candidates_stable" %in% names(step2_results)) {
  step2_results$age_varying_candidates_stable
} else if ("age_varying_candidates" %in% names(step2_results)) {
  step2_results$age_varying_candidates
} else {
  character(0)
}

age_common_candidates <- if ("age_common_candidates" %in% names(step2_results)) {
  step2_results$age_common_candidates
} else {
  character(0)
}

## ------------------------------------------------------
## Helper: detailed Step 3 reason string
## ------------------------------------------------------

build_reason_detailed <- function(sub_log,
                                  alpha   = 0.05,
                                  qic_tol = 0) {
  if (!nrow(sub_log)) return(NA_character_)
  
  # p-value criterion
  has_p <- any(!is.na(sub_log$p) & sub_log$p <= alpha)
  p_min <- if (has_p) min(sub_log$p, na.rm = TRUE) else NA_real_
  
  # QIC criterion: dQIC > 0 means dropping worsens QIC
  has_qic  <- any(!is.na(sub_log$dQIC) & sub_log$dQIC > qic_tol)
  dQIC_max <- if (has_qic) max(sub_log$dQIC, na.rm = TRUE) else NA_real_
  
  # Confounding criterion
  has_conf  <- any(sub_log$conf_flag %in% TRUE)
  meanR_max <- if (has_conf) max(sub_log$mean_rel, na.rm = TRUE) else NA_real_
  maxR_max  <- if (has_conf) max(sub_log$max_rel, na.rm = TRUE) else NA_real_
  
  pieces <- character(0)
  if (has_p && is.finite(p_min)) {
    pieces <- c(pieces, sprintf("p-value (p=%.3f)", p_min))
  }
  if (has_qic && is.finite(dQIC_max)) {
    pieces <- c(pieces, sprintf("QIC (dQIC=%.2f)", dQIC_max))
  }
  if (has_conf && is.finite(meanR_max) && is.finite(maxR_max)) {
    pieces <- c(
      pieces,
      sprintf(
        "confounding (mean dPR=%.1f%%, max dPR=%.1f%%)",
        100 * meanR_max, 100 * maxR_max
      )
    )
  }
  
  if (!length(pieces)) {
    return("Selected")
  }
  paste("Selected:", paste(pieces, collapse = "/"))
}

## ------------------------------------------------------
## Minimal build_rows_for_var:
##  - Only used to pull Step 3 reasons + age_pattern
## ------------------------------------------------------

build_rows_for_var <- function(v,
                               alpha_main = 0.05,
                               alpha_int  = 0.05,
                               qic_tol    = 0) {
  passed_step1 <- ifelse(v %in% step1_keepers_collapsed, "Yes", "No")
  in_final     <- v %in% union(main_terms_current, interaction_vars_current)
  has_int      <- v %in% interaction_vars_current
  
  age_type <- dplyr::case_when(
    !in_final               ~ "Not retained in final model",
    has_int                 ~ "Age-specific",
    in_final & !has_int     ~ "Age-common",
    TRUE                    ~ NA_character_
  )
  
  ## Backward main (Step 3)
  if (v %in% protected_vars) {
    backward_main_result <- "Selected: DAG identified variable"
  } else if (v %in% interaction_vars_current) {
    backward_main_result <- "Selected: kept due to interaction term"
  } else {
    back_main_sub <- back_main_log %>% dplyr::filter(var == v)
    if (nrow(back_main_sub) && v %in% c(initial_main_terms, initial_interaction_vars)) {
      main_kept <- v %in% main_terms_current
      if (main_kept) {
        backward_main_result <- build_reason_detailed(back_main_sub, alpha_main, qic_tol)
      } else {
        backward_main_result <- "Not selected"
      }
    } else {
      backward_main_result <- NA_character_
    }
  }
  
  ## Backward interaction (Step 3)
  if (!(v %in% age_varying_step3)) {
    if (v %in% age_common_candidates) {
      backward_int_result <- "Not tested due to Step 2 results"
    } else {
      backward_int_result <- NA_character_
    }
  } else {
    back_int_sub <- back_int_log %>% dplyr::filter(var == v)
    if (nrow(back_int_sub) && v %in% c(initial_main_terms, initial_interaction_vars)) {
      int_kept <- v %in% interaction_vars_current
      if (int_kept) {
        backward_int_result <- build_reason_detailed(back_int_sub, alpha_int, qic_tol)
      } else {
        backward_int_result <- "Not selected"
      }
    } else {
      backward_int_result <- NA_character_
    }
  }
  
  ## Forward main (Step 3)
  fwd_main_sub <- fwd_main_log %>% dplyr::filter(var == v)
  if (nrow(fwd_main_sub)) {
    main_kept <- v %in% main_terms_current
    if (main_kept) {
      forward_main_result <- build_reason_detailed(fwd_main_sub, alpha_main, qic_tol)
    } else {
      forward_main_result <- "Tested but not selected"
    }
  } else {
    forward_main_result <- NA_character_
  }
  
  ## Forward interaction (Step 3)
  fwd_int_sub <- fwd_int_log %>% dplyr::filter(var == v)
  if (nrow(fwd_int_sub)) {
    int_kept <- v %in% interaction_vars_current
    if (int_kept) {
      forward_int_result <- build_reason_detailed(fwd_int_sub, alpha_int, qic_tol)
    } else {
      forward_int_result <- "Tested but not selected"
    }
  } else {
    forward_int_result <- NA_character_
  }
  
  tibble::tibble(
    parent_variable            = v,
    age_pattern                = age_type,
    backward_main_result_step3 = backward_main_result,
    backward_int_result_step3  = backward_int_result,
    forward_main_result_step3  = forward_main_result,
    forward_int_result_step3   = forward_int_result
  )
}

## ------------------------------------------------------
## Summarise selection for each parent variable
##  - Step 1 + Step 2 status (with stability)
##  - Step 3 reasons (backward + forward)
## ------------------------------------------------------

summarise_parent_selection <- function(v) {
  # Step 3 detailed row
  tmp <- build_rows_for_var(v)
  
  ## Step 1 selection stability
  sel_freq <- NA_real_
  if (!is.null(step1_results$selection_freq)) {
    sf <- step1_results$selection_freq
    if (nrow(sf) && "Var1" %in% names(sf) && "Freq" %in% names(sf)) {
      hit <- sf$Freq[sf$Var1 == v]
      if (length(hit)) sel_freq <- as.numeric(hit[1])
    }
  }
  
  ## If candidate but never selected, treat stability = 0
  if (is.na(sel_freq) && v %in% step1_all_candidates) {
    sel_freq <- 0
  }
  
  ## Step 1 status (always-keep / retained / dropped / not candidate)
  if (v %in% always_keep) {
    step1_status <- "Always-keep (DAG)"
  } else if (v %in% step1_keepers_collapsed) {
    if (!is.na(sel_freq)) {
      step1_status <- sprintf("Retained in Step 1 (stability = %.2f)", sel_freq)
    } else {
      step1_status <- "Retained in Step 1"
    }
  } else if (v %in% step1_all_candidates) {
    if (!is.na(sel_freq)) {
      step1_status <- sprintf("Dropped in Step 1 (stability = %.2f)", sel_freq)
    } else {
      step1_status <- "Dropped in Step 1"
    }
  } else {
    step1_status <- "Not in Step 1 candidate set"
  }
  
  ## If appears in Step 2 candidate sets but no stability, set 0
  if (is.na(sel_freq) &&
      (v %in% age_varying_candidates_stable || v %in% age_common_candidates)) {
    sel_freq <- 0
  }
  
  ## Step 2 status
  step2_status <- dplyr::case_when(
    v %in% age_varying_candidates_stable & !is.na(sel_freq) ~
      sprintf("Age-varying candidate (stability = %.2f)", sel_freq),
    v %in% age_varying_candidates_stable ~
      "Age-varying candidate",
    v %in% age_common_candidates & !is.na(sel_freq) ~
      sprintf("Age-common candidate (stability = %.2f)", sel_freq),
    v %in% age_common_candidates ~
      "Age-common candidate",
    v %in% step1_keepers_collapsed ~
      "Not evaluated in Step 2",
    TRUE ~ NA_character_
  )
  
  if (!nrow(tmp)) {
    tmp <- tibble::tibble(
      parent_variable            = v,
      age_pattern                = "Not retained in final model",
      backward_main_result_step3 = NA_character_,
      backward_int_result_step3  = NA_character_,
      forward_main_result_step3  = NA_character_,
      forward_int_result_step3   = NA_character_
    )
  } else if (!(v %in% union(main_terms_current, interaction_vars_current))) {
    tmp$age_pattern <- "Not retained in final model"
  }
  
  tmp$step1_status <- step1_status
  tmp$step2_status <- step2_status
  
  tmp
}

###########################################################
## Build var_selection_table_all:
##  - Community-level block first
##  - Then individual-level in Table 1 order
###########################################################

selection_rows <- list()

## ---- 0) Build a safe lookup for community_labels (if provided) ----

comm_label_lookup <- NULL
if (exists("community_labels")) {
  # codes for rows that actually have a code
  comm_codes_non_na  <- community_variable_code[!is.na(community_variable_code)]
  # labels should correspond 1:1 to those coded rows
  if (length(community_labels) == length(comm_codes_non_na)) {
    comm_label_lookup <- setNames(community_labels, comm_codes_non_na)
  } else {
    warning("Length of community_labels does not match non-NA community_variable_code; using default labels.")
  }
}

## 1) Community-level variables (Table1_comm order)
selection_rows[[length(selection_rows) + 1L]] <- tibble::tibble(
  Variable_label              = "Community-level factors (per 10%-point increase)",
  parent_variable             = NA_character_,
  age_pattern                 = NA_character_,
  step1_status                = NA_character_,
  step2_status                = NA_character_,
  backward_main_result_step3  = NA_character_,
  backward_int_result_step3   = NA_character_,
  forward_main_result_step3   = NA_character_,
  forward_int_result_step3    = NA_character_
)

for (i in seq_along(community_variable_names)) {
  comm_label <- community_variable_names[i]
  comm_code  <- community_variable_code[i]  # unscaled code, e.g., "comm_pov_prop"
  
  # Domain subtitles (no code)
  if (is.na(comm_code)) {
    selection_rows[[length(selection_rows) + 1L]] <- tibble::tibble(
      Variable_label              = comm_label,
      parent_variable             = NA_character_,
      age_pattern                 = NA_character_,
      step1_status                = NA_character_,
      step2_status                = NA_character_,
      backward_main_result_step3  = NA_character_,
      backward_int_result_step3   = NA_character_,
      forward_main_result_step3   = NA_character_,
      forward_int_result_step3    = NA_character_
    )
    next
  }
  
  parent_scaled <- paste0(comm_code, "_scaled")
  sel_info      <- summarise_parent_selection(parent_scaled)
  
  # Use custom label if lookup exists, otherwise default to comm_label
  if (!is.null(comm_label_lookup)) {
    lbl <- unname(comm_label_lookup[comm_code])  # will be NA if not found
    base_label <- ifelse(is.na(lbl), comm_label, lbl)
  } else {
    base_label <- comm_label
  }
  
  selection_rows[[length(selection_rows) + 1L]] <- sel_info %>%
    dplyr::mutate(Variable_label = base_label) %>%
    dplyr::select(
      Variable_label, parent_variable,
      age_pattern, step1_status, step2_status,
      backward_main_result_step3,
      backward_int_result_step3,
      forward_main_result_step3,
      forward_int_result_step3
    )
}
## 2) Individual-level variables (Table 1 order)

selection_rows[[length(selection_rows) + 1L]] <- tibble::tibble(
  Variable_label              = "Individual-level factors",
  parent_variable             = NA_character_,
  age_pattern                 = NA_character_,
  step1_status                = NA_character_,
  step2_status                = NA_character_,
  backward_main_result_step3  = NA_character_,
  backward_int_result_step3   = NA_character_,
  forward_main_result_step3   = NA_character_,
  forward_int_result_step3    = NA_character_
)

header_codes <- variable_code[!is.na(variable_code) & is.na(response_list)] |> unique()

for (i in seq_along(variable_names)) {
  code  <- variable_code[i]
  resp  <- response_list[i]
  label <- variable_names[i]
  
  # Domain subtitles
  if (is.na(code)) {
    selection_rows[[length(selection_rows) + 1L]] <- tibble::tibble(
      Variable_label              = label,
      parent_variable             = NA_character_,
      age_pattern                 = NA_character_,
      step1_status                = NA_character_,
      step2_status                = NA_character_,
      backward_main_result_step3  = NA_character_,
      backward_int_result_step3   = NA_character_,
      forward_main_result_step3   = NA_character_,
      forward_int_result_step3    = NA_character_
    )
    next
  }
  
  # Multi-level variables: header row is the parent
  if (!is.na(code) && is.na(resp)) {
    v <- code
    sel_info <- summarise_parent_selection(v)
    
    selection_rows[[length(selection_rows) + 1L]] <- sel_info %>%
      dplyr::mutate(Variable_label = label) %>%
      dplyr::select(
        Variable_label, parent_variable,
        age_pattern, step1_status, step2_status,
        backward_main_result_step3,
        backward_int_result_step3,
        forward_main_result_step3,
        forward_int_result_step3
      )
    next
  }
  
  # Binary / single-level variables
  if (!is.na(code) && !is.na(resp) && !(code %in% header_codes)) {
    v <- code
    sel_info <- summarise_parent_selection(v)
    
    selection_rows[[length(selection_rows) + 1L]] <- sel_info %>%
      dplyr::mutate(Variable_label = label) %>%
      dplyr::select(
        Variable_label, parent_variable,
        age_pattern, step1_status, step2_status,
        backward_main_result_step3,
        backward_int_result_step3,
        forward_main_result_step3,
        forward_int_result_step3
      )
  }
}

var_selection_table_all <- dplyr::bind_rows(selection_rows)

############################################################
## NOTE: final_results_table and var_selection_table_all
##       are already computed in best_res. We reuse them
##       here and do NOT recompute them.
############################################################

# Just ensure we have the final fitted model handy
final_fit  <- step3_results$fit
final_coef <- coef(final_fit)
final_vcov <- vcov(final_fit)

## ----------------------------------------
## Full regression output from final_fit
## ----------------------------------------

est <- as.numeric(final_coef)
se  <- sqrt(diag(final_vcov))
z   <- est / se
p   <- 2 * pnorm(-abs(z))

# 95% CI on log scale
lcl_log <- est - 1.96 * se
ucl_log <- est + 1.96 * se

# Exponentiate to get PRs and CI for PRs
PR      <- exp(est)
PR_lcl  <- exp(lcl_log)
PR_ucl  <- exp(ucl_log)

regression_table_raw <- data.frame(
  term          = names(final_coef),
  estimate_log  = est,
  std_error     = se,
  z_value       = z,
  p_value       = p,
  conf_low_log  = lcl_log,
  conf_high_log = ucl_log,
  PR            = PR,
  PR_low        = PR_lcl,
  PR_high       = PR_ucl,
  stringsAsFactors = FALSE
)

## ----------------------------------------
## 1. Parent-variable order & labels
##    from var_selection_table_all
## ----------------------------------------

parent_order_tbl <- var_selection_table_all %>%
  dplyr::filter(!is.na(parent_variable)) %>%
  dplyr::select(parent_variable, Variable_label) %>%
  dplyr::distinct()

parent_vars <- parent_order_tbl$parent_variable

# header_codes: multi-level variables (from Table 1)
if (!exists("header_codes")) {
  header_codes <- variable_code[!is.na(variable_code) & is.na(response_list)] |> unique()
}

## ----------------------------------------
## 2. Map each term to its parent_variable
## ----------------------------------------

map_term_to_parent <- function(term) {
  if (term == "(Intercept)") return("(Intercept)")
  hits <- parent_vars[
    vapply(parent_vars, function(p) {
      grepl(paste0("(^|:)", p, "($|[0-9\\.:])"), term)
    }, logical(1))
  ]
  if (length(hits) == 0L) return(NA_character_)
  hits[1]
}

regression_table_raw$parent_variable <- vapply(
  regression_table_raw$term,
  map_term_to_parent,
  character(1)
)

## ----------------------------------------
## 3. Basic ordering by parent
## ----------------------------------------

parent_order <- unique(parent_order_tbl$parent_variable)

extra_parents <- setdiff(
  unique(regression_table_raw$parent_variable),
  c("(Intercept)", parent_order)
)

regression_table_ordered <- regression_table_raw %>%
  dplyr::mutate(
    parent_variable = factor(
      parent_variable,
      levels = c("(Intercept)", parent_order, extra_parents)
    )
  ) %>%
  dplyr::arrange(parent_variable, term)

## ----------------------------------------
## 4. Add labels:
##    - parent label from selection table
##    - dummy-level label from Table 1
## ----------------------------------------

# Parent-level label
regression_table_ordered <- regression_table_ordered %>%
  dplyr::left_join(parent_order_tbl, by = "parent_variable") %>%
  dplyr::rename(variable_name = Variable_label) %>%
  dplyr::relocate(variable_name, .before = term)

# Lookup for dummy-level labels from Table 1
labels_tbl <- table1 %>%
  dplyr::mutate(order_index = dplyr::row_number()) %>%  # keep Table 1 row order
  dplyr::filter(!is.na(variable_code), !is.na(response)) %>%
  dplyr::transmute(
    variable_code,
    response    = as.numeric(response),
    label       = Variable,
    is_ref      = response == reference,
    order_index = order_index
  ) %>%
  dplyr::distinct()

get_label_for_term <- function(term, parent_variable, default_label) {
  # Intercept
  if (term == "(Intercept)") {
    return(ifelse(is.na(default_label), "(Intercept)", default_label))
  }
  
  # Work with last piece (after any ":")
  last_piece <- sub("^.*:", "", term)
  
  # Parse "var123" -> var + response
  m <- regexec("^(.+?)([0-9]+)$", last_piece)
  reg <- regmatches(last_piece, m)[[1]]
  if (length(reg) == 3) {
    base_var <- reg[2]
    resp_val <- suppressWarnings(as.numeric(reg[3]))
    hit <- which(
      labels_tbl$variable_code == base_var &
        labels_tbl$response      == resp_val
    )
    if (length(hit) >= 1L) {
      return(labels_tbl$label[hit[1]])
    }
  }
  
  default_label
}

regression_table_ordered$variable_name <- vapply(
  seq_len(nrow(regression_table_ordered)),
  function(i) {
    get_label_for_term(
      term            = regression_table_ordered$term[i],
      parent_variable = regression_table_ordered$parent_variable[i],
      default_label   = regression_table_ordered$variable_name[i]
    )
  },
  character(1)
)

## ----------------------------------------
## 5. Round, format p-values, and tag age-band interactions
## ----------------------------------------

regression_table <- regression_table_ordered %>%
  dplyr::mutate(
    estimate_log  = round(estimate_log, 3),
    std_error     = round(std_error, 3),
    z_value       = round(z_value, 3),
    conf_low_log  = round(conf_low_log, 3),
    conf_high_log = round(conf_high_log, 3),
    PR            = round(PR, 3),
    PR_low        = round(PR_low, 3),
    PR_high       = round(PR_high, 3),
    PR_CI         = dplyr::if_else(
      is.na(PR),
      NA_character_,
      sprintf("%.3f (%.3f–%.3f)", PR, PR_low, PR_high)
    ),
    p_value_fmt   = dplyr::case_when(
      is.na(p_value)      ~ NA_character_,
      p_value < 0.001     ~ "<0.001",
      TRUE                ~ sprintf("%.3f", round(p_value, 3))
    ),
    # Flag age-band interaction terms
    is_age_int    = grepl(paste0("^", age_band_var, ":"), term)
  ) %>%
  dplyr::mutate(
    variable_name = if_else(
      is_age_int & !is.na(variable_name) & term != "(Intercept)",
      paste0("Age × ", variable_name),
      variable_name
    )
  )

## ----------------------------------------
## 6. Inject headings and order levels:
##    - no duplicated rows for community/binary variables
##    - multi-level vars get a parent header row
##    - within parent: main, then age interaction for each level
##    - drop subtitles that have no covariates
##    - add "Individual-level factors" block header
## ----------------------------------------

# Separate intercept
intercept_row <- regression_table %>% dplyr::filter(term == "(Intercept)")
coef_rows     <- regression_table %>% dplyr::filter(term != "(Intercept)")

# Helper: empty heading row
empty_heading_row <- function(label) {
  tibble(
    variable_name   = label,
    term            = NA_character_,
    estimate_log    = NA_real_,
    std_error       = NA_real_,
    z_value         = NA_real_,
    p_value         = NA_real_,
    conf_low_log    = NA_real_,
    conf_high_log   = NA_real_,
    PR              = NA_real_,
    PR_low          = NA_real_,
    PR_high         = NA_real_,
    parent_variable = NA_character_,
    PR_CI           = NA_character_,
    p_value_fmt     = NA_character_,
    is_age_int      = NA
  )
}

# Helper: reference row for multi-level categorical variables
make_reference_row <- function(label, parent) {
  tibble(
    variable_name   = label,
    term            = NA_character_,
    estimate_log    = NA_real_,
    std_error       = NA_real_,
    z_value         = NA_real_,
    p_value         = NA_real_,
    conf_low_log    = NA_real_,
    conf_high_log   = NA_real_,
    PR              = NA_real_,
    PR_low          = NA_real_,
    PR_high         = NA_real_,
    parent_variable = parent,
    PR_CI           = "Reference",
    p_value_fmt     = NA_character_,
    is_age_int      = FALSE
  )
}

# Helper: within a parent, order rows as
# level main -> level interaction -> next level main -> ...
order_within_parent <- function(df_parent) {
  if (!nrow(df_parent)) return(df_parent)
  
  int_flag <- grepl(paste0("^", age_band_var, ":"), df_parent$term)
  main_rows <- df_parent[!int_flag, , drop = FALSE]
  int_rows  <- df_parent[int_flag, , drop = FALSE]
  
  # If no interactions for this parent, just sort by term
  if (!nrow(int_rows)) {
    return(main_rows[order(main_rows$term), , drop = FALSE])
  }
  
  # Map interaction term name -> row index
  int_index <- seq_len(nrow(int_rows))
  names(int_index) <- int_rows$term
  
  ordered <- list()
  
  # Sort mains by term to get stable level order
  main_rows <- main_rows[order(main_rows$term), , drop = FALSE]
  
  for (i in seq_len(nrow(main_rows))) {
    tm <- main_rows$term[i]
    ordered[[length(ordered) + 1L]] <- main_rows[i, ]
    
    int_term <- paste0(age_band_var, ":", tm)
    if (int_term %in% names(int_index)) {
      ordered[[length(ordered) + 1L]] <- int_rows[int_index[int_term], , drop = FALSE]
    }
  }
  
  # Any interaction rows without a main effect (should be rare) go at the end
  leftover <- setdiff(int_rows$term, paste0(age_band_var, ":", main_rows$term))
  if (length(leftover)) {
    ordered[[length(ordered) + 1L]] <- int_rows[int_rows$term %in% leftover, , drop = FALSE]
  }
  
  dplyr::bind_rows(ordered)
}

final_rows <- list()

## 6A. Intercept row
if (nrow(intercept_row) == 1L) {
  final_rows[[length(final_rows) + 1L]] <- intercept_row
}

## 6B. Walk through var_selection_table_all in order
for (i in seq_len(nrow(var_selection_table_all))) {
  this_label  <- var_selection_table_all$Variable_label[i]
  this_parent <- var_selection_table_all$parent_variable[i]
  
  if (is.na(this_label) && is.na(this_parent)) next
  
  # Headings: domain titles, community block title, etc.
  if (is.na(this_parent) && !is.na(this_label)) {
    final_rows[[length(final_rows) + 1L]] <- empty_heading_row(this_label)
    next
  }
  
  if (!is.na(this_parent)) {
    base_var  <- gsub("_scaled$", "", this_parent)
    is_comm   <- grepl("^comm_", base_var)
    is_multi  <- base_var %in% header_codes   # multi-level categoricals
    
    rows_pv <- coef_rows %>%
      dplyr::filter(parent_variable == this_parent)
    if (!nrow(rows_pv)) next
    
    if (is_multi && !is_comm && !is_na(this_label)) {
      # Parent header (e.g., "Household wealth quintile")
      final_rows[[length(final_rows) + 1L]] <- empty_heading_row(this_label)
      
      # rows_pv: all non-reference dummies (and their age interactions) for this parent
      # Add helper columns to group main & interaction for each level
      rows_pv$is_interaction <- grepl(paste0("^", age_band_var, ":"), rows_pv$term)
      rows_pv$base_dummy     <- sub("^.*:", "", rows_pv$term)  # e.g., "ph_wealth_quint1"
      rows_pv$effect_order   <- ifelse(rows_pv$is_interaction, 2L, 1L)  # main(1) then interaction(2)
      
      # Table 1 order for this variable
      labs <- labels_tbl %>%
        dplyr::filter(variable_code == base_var) %>%
        dplyr::arrange(order_index)
      
      for (j in seq_len(nrow(labs))) {
        lj <- labs[j, ]
        
        if (isTRUE(lj$is_ref)) {
          # Reference level row (e.g., "   Richest")
          final_rows[[length(final_rows) + 1L]] <-
            make_reference_row(lj$label, this_parent)
        } else {
          # Non-reference level: get main + interaction rows for this dummy
          dummy_name <- paste0(base_var, lj$response)  # e.g., "ph_wealth_quint1"
          sub_rows <- rows_pv %>%
            dplyr::filter(base_dummy == dummy_name) %>%
            dplyr::arrange(effect_order)
          
          if (nrow(sub_rows)) {
            final_rows[[length(final_rows) + 1L]] <- sub_rows
          }
        }
      }
      
    } else {
      # Original behavior for non-multi or community variables
      rows_pv_ordered <- order_within_parent(rows_pv)
      final_rows[[length(final_rows) + 1L]] <- rows_pv_ordered
    }
  }
}

# Bind to get the raw headed table
regression_table_with_headings <- dplyr::bind_rows(final_rows)

## 6C. Insert "Individual-level factors" subheadline
rt <- regression_table_with_headings

# Only add if not already present and if we have any non-community terms
if (!any(rt$variable_name == "Individual-level factors", na.rm = TRUE)) {
  comm_rows <- which(!is.na(rt$term) & grepl("^comm_", rt$term))
  if (length(comm_rows)) {
    last_comm <- max(comm_rows)
    
    has_indiv <- any(
      !is.na(rt$term[(last_comm + 1L):nrow(rt)]) &
        !grepl("^comm_", rt$term[(last_comm + 1L):nrow(rt)])
    )
    
    if (has_indiv) {
      indiv_header <- empty_heading_row("Individual-level factors")
      
      rt <- dplyr::bind_rows(
        rt[seq_len(last_comm), ],
        indiv_header,
        rt[(last_comm + 1L):nrow(rt), ]
      )
    }
  }
}

## 6D. Drop subtitle headings that have no covariates beneath them
##     (look only up to the next *top-level* heading)
is_heading <- is.na(rt$term) & !is.na(rt$variable_name)

# Top-level headings that we always keep
is_top_heading <- is_heading & rt$variable_name %in% c(
  "Community-level factors (per 10%-point increase)",
  "Individual-level factors"
)

drop <- logical(nrow(rt))
heading_idx <- which(is_heading & !is_top_heading)

for (h in heading_idx) {
  # Find the next *top-level* heading below this one
  next_top <- which(is_top_heading & seq_len(nrow(rt)) > h)
  next_top <- if (length(next_top)) min(next_top) else nrow(rt) + 1L
  
  # Rows between this heading and the next top-level heading (or end)
  if (next_top > h + 1L) {
    block_seq <- (h + 1L):(next_top - 1L)
  } else {
    block_seq <- integer(0)
  }
  
  # If there is no non-heading row (i.e., no actual coefficient rows)
  # in this block, drop this subtitle heading
  if (!length(block_seq) || !any(!is.na(rt$term[block_seq]))) {
    drop[h] <- TRUE
  }
}

regression_table_with_headings <- rt[!drop, , drop = FALSE]

## =========================================================
## Post-process regression_table_with_headings
## - Label age-band interactions as "Age*..."
## - Reorder within each parent: main, then interaction
## =========================================================

df <- regression_table_with_headings

# Keep original row index (for safety if needed later)
df$orig_row <- seq_len(nrow(df))

## 1) Identify age-band interaction terms and relabel variable_name
is_age_int2 <- !is.na(df$term) & grepl("age_band24-59:", df$term)

idx  <- which(is_age_int2 & !is.na(df$variable_name))
orig <- df$variable_name[idx]

# does the original label start with a space?
has_leading_space <- grepl("^\\s", orig)

# strip all leading spaces for the base text
base <- sub("^\\s*", "", orig)

# if it had a leading space → "  Age*X"
# if not → " Age*X"
df$variable_name[idx] <- ifelse(
  has_leading_space,
  paste0("  Age*", base),
  paste0(" Age*", base)
)

## 2) Reorder rows *within each parent_variable* block

get_base_dummy <- function(term) {
  ifelse(grepl(":", term), sub("^.*:", "", term), term)
}

for (pv in unique(na.omit(df$parent_variable))) {
  idx_p <- which(df$parent_variable == pv & !is.na(df$term))
  if (length(idx_p) <= 1L) next
  
  block <- df[idx_p, , drop = FALSE]
  
  block$is_interaction <- grepl("age_band24-59:", block$term)
  block$base_dummy     <- get_base_dummy(block$term)
  
  dummy_levels <- unique(block$base_dummy)
  block$base_dummy <- factor(block$base_dummy, levels = dummy_levels)
  
  block$effect_order <- ifelse(block$is_interaction, 2L, 1L)
  
  block <- block[order(block$base_dummy, block$effect_order), ]
  
  block$is_interaction <- NULL
  block$base_dummy     <- NULL
  block$effect_order   <- NULL
  
  df[idx_p, ] <- block
}

df$orig_row <- NULL

regression_table_with_headings <- df

## --------------------------
## Tuning model ranking
## --------------------------

# Start from cv_brier_results (config_index, config_name, cv_brier)
model_summary <- cv_brier_results %>%
  # Join diagnostics from sens_metrics
  left_join(sens_metrics, by = "config_name") %>%
  # Parse tuning parameters from config_name
  mutate(
    # Step 1: lambda calculation
    step1_lambda = str_match(config_name, "lambda_([^_]+)_")[, 2],
    
    # Step 1 & 2: bootstrap stability threshold
    stab_cut     = as.numeric(str_match(config_name, "step1_stab_([0-9.]+)_")[, 2]),
    
    # Step 2: p-value threshold (p_max)
    step2_p_max  = as.numeric(str_match(config_name, "step2_p_([0-9.]+)_stab_")[, 2]),
    
    # Step 3: p-value threshold
    step3_p      = as.numeric(str_match(config_name, "step3_p_([0-9.]+)_qic_")[, 2]),
    
    # Step 3: dQIC threshold (can be -2 / 0 / 2)
    step3_dqic   = as.numeric(str_match(config_name, "_qic_(-?[0-9.]+)_meanpr_")[, 2]),
    
    # Step 3: mean and max dPR thresholds
    mean_dPR_th  = as.numeric(str_match(config_name, "meanpr_([0-9.]+)_maxpr_")[, 2]),
    max_dPR_th   = as.numeric(str_match(config_name, "maxpr_([0-9.]+)_fam_")[, 2])
  )

best_brier <- min(model_summary$cv_brier, na.rm = TRUE)

model_summary <- model_summary %>%
  mutate(
    # Relative difference from best (0 = best)
    rel_diff = (cv_brier - best_brier) / best_brier,
    rel_diff = pmax(rel_diff, 0),    # guard against tiny negatives
    
    # 0.1% = 0.001 relative → bin 1 = [0, 0.001), bin 2 = [0.001, 0.002), ...
    brier_bin = floor(rel_diff / 0.001) + 1L
  ) %>%
  # Arrange by your rule: bin → parsimony → Brier → QIC
  arrange(brier_bin, n_terms, cv_brier, qic) 

model_summary <- model_summary %>%
  mutate(
    # group defining "essentially the same model"
    rank_group = interaction(brier_bin, n_terms, drop = TRUE),
    
    # rank is the order of unique (bin, n_terms) combos
    Rank = as.integer(factor(rank_group, levels = unique(rank_group)))
  )

model_selection_table <- model_summary %>%
  transmute(
    Model = config_index,
    `Step 1: Penalty λ selection`             = step1_lambda,
    `Step 1 & 2: bootstrap threshold`         = stab_cut,
    `Step 2: p-value threshold`               = step2_p_max,
    `Step 3: p-value threshold`               = step3_p,
    `Step 3: dQIC threshold`                  = step3_dqic,
    `Step 3: mean / max dPR thresholds`       = sprintf("%.2f / %.2f", mean_dPR_th, max_dPR_th),
    `Pearson dispersion (φ)`                  = pearson_disp,
    `QIC`                                     = qic,
    `Calibration coefficient (r)`             = calib_corr,
    `Rao–Scott 2logLR`                        = rao_scott_stat,
    `CV Brier score`                          = cv_brier,
    `ΔBrier bin (≤0.1% from best)`           = brier_bin,
    `Number of terms`                         = n_terms,
    Rank                                      = Rank
  )
model_selection_table$Model <- 1:nrow(model_selection_table)
model_selection_table$`Step 1: Penalty λ selection` <- ifelse(model_selection_table$`Step 1: Penalty λ selection` == 'lambda.1se', '1-SE',
                                                           ifelse(model_selection_table$`Step 1: Penalty λ selection` == 'lambda.min', 'Minimum-error', NA))


###########################################################
## Minimal GOF quantities for plots (no table)
###########################################################

# final_fit should already exist, e.g.
# final_fit <- best_result$step3_results$fit

## 1) Core pieces from final_fit
fam      <- family(final_fit)
fam_str  <- fam$family
link_str <- fam$link

# Response variable name
y_name <- all.vars(formula(final_fit))[1]

# Survey design object used in final_fit
design_obj <- final_fit$survey.design

# Fitted values and Pearson residuals
fit_vals <- fitted(final_fit)
res_pear <- residuals(final_fit, type = "pearson")

## 2) Calibration by deciles (to get pred_vec, obs_vec)

design_cal <- update(design_obj, pred = fit_vals)

brks <- quantile(
  design_cal$variables$pred,
  probs = seq(0, 1, 0.1),
  na.rm = TRUE
)

design_cal <- update(
  design_cal,
  decile = cut(pred, breaks = brks, include.lowest = TRUE)
)

form_y <- as.formula(paste("~", y_name))

cal <- svyby(
  form_y,
  ~ decile,
  design = design_cal,
  FUN = svymean,
  na.rm = TRUE
)

pred_tab <- svyby(
  ~ pred,
  ~ decile,
  design = design_cal,
  FUN = svymean,
  na.rm = TRUE
)

obs_vec  <- cal[[y_name]]
pred_vec <- pred_tab$pred

## -----------------------------
## 6. Combined plots: residuals + calibration (1 × 2)
## -----------------------------

df_res <- data.frame(
  Fitted            = fit_vals,
  Pearson_residuals = res_pear
)

p_res <- ggplot(df_res, aes(x = Fitted, y = Pearson_residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.4, size = 0.6) +
  labs(
    x = "Fitted prevalence",
    y = "Pearson residual",
    title    = "Residuals versus fitted values",
    subtitle = "Quasi-Poisson log-link with survey weights"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9)
  )

df_cal <- data.frame(
  Predicted = pred_vec,
  Observed  = obs_vec
)

p_cal <- ggplot(df_cal, aes(x = Predicted, y = Observed)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(size = 1.8) +
  coord_equal() +
  labs(
    x = "Mean fitted probability (10% intervals)",
    y = "Mean observed prevalence",
    title    = "Calibration by deciles of fitted probability",
    subtitle = "Points represent deciles; dashed line indicates perfect calibration"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(size = 9)
  )

combined_plot <- p_res + p_cal + plot_layout(nrow = 1)

fig2s_path <- file.path(FIGURE_DIR, "Figure_2S_residual_calibration.png")
ggsave(
  filename = fig2s_path,
  plot     = combined_plot,
  width    = 10,   # adjust as you like
  height   = 6,    # adjust as you like
  dpi      = 300
)

## =========================================================
## 1. Encode category (incl. "not_selected")
## =========================================================
summarise_parent_effect <- function(res, parent_variable, older_band = FALSE) {
  fit <- res$step3_results$fit
  beta <- coef(fit)
  vcov_mat <- vcov(fit)
  coef_names <- names(beta)

  main_names <- grep(paste0("^", parent_variable), coef_names, value = TRUE)
  main_names <- main_names[!grepl(":", main_names)]
  in_final <- parent_variable %in% union(
    res$step3_results$main_terms,
    res$step3_results$interaction_vars
  )

  if (!in_final || !length(main_names)) {
    return(tibble(in_final = FALSE, dir = "none", age_specific = FALSE))
  }

  estimates <- p_values <- numeric(length(main_names))
  for (i in seq_along(main_names)) {
    contrast <- setNames(rep(0, length(beta)), coef_names)
    contrast[main_names[[i]]] <- 1

    if (older_band) {
      int_name <- coef_names[
        startsWith(coef_names, "age_band") &
          endsWith(coef_names, paste0(":", main_names[[i]]))
      ]
      if (length(int_name)) contrast[int_name[[1]]] <- 1
    }

    estimates[[i]] <- sum(contrast * beta)
    se <- sqrt(as.numeric(t(contrast) %*% vcov_mat %*% contrast))
    p_values[[i]] <- if (is.finite(se) && se > 0) {
      2 * pnorm(-abs(estimates[[i]] / se))
    } else {
      NA_real_
    }
  }

  significant <- !is.na(p_values) & p_values < 0.05
  direction <- if (!any(significant)) {
    "none"
  } else if (all(estimates[significant] > 0)) {
    "pos"
  } else if (all(estimates[significant] < 0)) {
    "neg"
  } else {
    "none"
  }

  tibble(
    in_final = TRUE,
    dir = direction,
    age_specific = parent_variable %in% res$step3_results$interaction_vars
  )
}

model_order <- model_summary$config_name
heat_df2 <- purrr::map_dfr(seq_along(model_order), function(model_rank) {
  config_name <- model_order[[model_rank]]
  res <- sens_results[[config_name]]
  parent_variables <- unique(na.omit(res$var_selection_table_all$parent_variable))

  purrr::map_dfr(parent_variables, function(parent_variable) {
    bind_rows(
      summarise_parent_effect(res, parent_variable, older_band = FALSE) %>%
        mutate(age_band_group = "6-23.9"),
      summarise_parent_effect(res, parent_variable, older_band = TRUE) %>%
        mutate(age_band_group = "24-59")
    ) %>%
      mutate(
        model_rank = model_rank,
        Model = factor(model_rank),
        config_name = config_name,
        parent_variable = parent_variable,
        variable_code = gsub("_scaled$", "", parent_variable)
      )
  })
}) %>%
  mutate(Model = factor(model_rank, levels = seq_along(model_order)))

heat_df2 <- heat_df2 %>%
  mutate(
    cat = case_when(
      !in_final                    ~ "not_selected",
      dir == "none"                ~ "none",
      dir == "pos" & !age_specific ~ "pos_common",
      dir == "neg" & !age_specific ~ "neg_common",
      dir == "pos" &  age_specific ~ "pos_specific",
      dir == "neg" &  age_specific ~ "neg_specific",
      TRUE                         ~ "not_selected"
    )
  )

## =========================================================
## 2. Parent label table & ordered var_id
##    (community block first, then individual-level Table 1)
## =========================================================

parent_label_df <- bind_rows(
  community_block %>% dplyr::select(Variable, variable_code),
  table1          %>% dplyr::select(Variable, variable_code)
) %>%
  filter(
    !is.na(variable_code),
    !grepl("^\\s", Variable)   # drop indented child rows
  ) %>%
  distinct(variable_code, .keep_all = TRUE)

parent_label_df$Variable <- ifelse(
  grepl("comm", parent_label_df$variable_code),
  paste0("Community: ", parent_label_df$Variable),
  parent_label_df$Variable
)

parent_order_codes <- parent_label_df$variable_code

# Rebuild var_index_df in that order
parent_vars <- unique(heat_df2$parent_variable)

var_index_df <- tibble(parent_variable = parent_vars) %>%
  mutate(
    base_code = gsub("_scaled$", "", parent_variable),
    ord       = match(base_code, parent_order_codes)
  ) %>%
  arrange(ord) %>%
  mutate(var_id = row_number()) %>%
  dplyr::select(parent_variable, base_code, var_id)

# Attach ordered var_id back to heat_df2
heat_df2 <- heat_df2 %>%
  dplyr::select(-any_of("var_id")) %>%
  left_join(var_index_df %>% dplyr::select(parent_variable, var_id),
            by = "parent_variable")

var_levels <- sort(unique(heat_df2$var_id))
heat_df2$var_id <- factor(heat_df2$var_id, levels = rev(var_levels))

## =========================================================
## 3. Colour palette & labels (used by heatmaps + legend)
## =========================================================

fill_values <- c(
  not_selected = "white",
  none         = "grey80",
  pos_common   = "rosybrown2",
  neg_common   = "lightblue",
  pos_specific = "red",
  neg_specific = "dodgerblue2"
)

fill_labels <- c(
  not_selected = "Not selected in final model",
  none         = "No significant association",
  pos_common   = "↑ PR, age-common",
  neg_common   = "↓ PR, age-common",
  pos_specific = "↑ PR, age-specific",
  neg_specific = "↓ PR, age-specific"
)

## =========================================================
## 4. Base heatmap function
## =========================================================
common_y_levels <- levels(heat_df2$var_id)
heat_base <- function(df, title) {
  ggplot(df, aes(x = Model, y = var_id, fill = cat)) +
    geom_tile(color = NA) +
    geom_text(
      data = df %>% dplyr::filter(cat == "not_selected"),
      aes(label = "X"),
      size      = 2,
      colour    = "black"
    ) +
    scale_fill_manual(
      values = fill_values,
      breaks = names(fill_labels),
      labels = fill_labels,
      drop   = FALSE
    ) +
    scale_y_discrete(
      name   = "Covariates",
      limits = common_y_levels,  # <--- force same levels & order
      drop   = FALSE,
      expand = c(0, 0)
    ) +
    scale_x_discrete(
      name   = "Models",
      breaks = function(x) {
        nx   <- as.numeric(as.character(x))
        max_model <- max(nx, na.rm = TRUE)
        keep <- if (max_model >= 10) c(1, seq(10, max_model, by = 10)) else 1
        as.character(keep)
      },
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") + 
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x   = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      panel.grid    = element_blank(),
      plot.title    = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(
      title    = title
    )
}

## =========================================================
## 5. Two age-band heatmaps (no legends)
## =========================================================

p_heat_6_23 <- heat_base(
  heat_df2 %>% filter(age_band_group == "6-23.9"),
  title = "Younger children (6–23.9 months)"
) + theme(legend.position = "none")

p_heat_24_59 <- heat_base(
  heat_df2 %>% filter(age_band_group == "24-59"),
  title = "Older children (24–59 months)"
) + theme(legend.position = "none")

## =========================================================
## 6. Left-side numbered covariate list, aligned to var_id
## =========================================================

# Rebuild var_key with var_id so it matches heat_df2
var_key <- var_index_df %>%
  left_join(
    parent_label_df %>% dplyr::rename(label = Variable),
    by = c("base_code" = "variable_code")
  ) %>%
  mutate(
    label          = if_else(is.na(label), parent_variable, label),
    idx            = row_number(),
    label_numbered = sprintf("%d. %s", idx, label),
    # use same factor + levels as in heat_df2 for vertical alignment
    var_id         = factor(var_id, levels = levels(heat_df2$var_id))
  )

# Covariate list with empty x-axis baseline
var_key$x_fake <- 0
p_var_key <- ggplot(var_key, aes(x = x_fake, y = var_id)) +
  geom_text(
    aes(label = label_numbered),
    hjust = 0, size = 3
  ) +
  scale_y_discrete(
    limits = common_y_levels,  # <--- force same levels & order
    drop   = FALSE,
    expand = c(0, 0)
  ) +
  scale_x_discrete(
    name   = "",
    breaks = "",
    limits = factor(0, 1),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid   = element_blank(),
    # x-axis: keep line, hide labels & ticks (invisible axis)
    axis.text.x   = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.title.x  = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
    axis.ticks.x = element_line(),
    axis.line.x  = element_line(),
    # y-axis: we use our own labels, so hide default
    axis.title.y  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title   = element_text(size = 11, face = "bold"),
    plot.margin  = margin(t = 10, r = 0, b = 25, l = 0)
  ) +
  labs(title = "Covariate (parent variable)")
## =========================================================
## 7. Compact manual legend
## =========================================================

legend_df <- tibble(
  cat = factor(
    c("not_selected", "none", "pos_common", "neg_common",
      "pos_specific", "neg_specific"),
    levels = c("not_selected", "none", "pos_common",
               "neg_common", "pos_specific", "neg_specific")
  ),
  label = c(
    "Not selected in final model",
    "No significant association",
    "\u2191 PR, age-common",
    "\u2193 PR, age-common",
    "\u2191 PR, age-specific",
    "\u2193 PR, age-specific"
  )
)

legend_df_compact <- legend_df %>%
  mutate(y = seq(from = 1, to = 0, length.out = n()))

min_y   <- min(legend_df_compact$y)
max_y   <- max(legend_df_compact$y)
title_y <- max_y + 0.2

p_legend <- ggplot(legend_df_compact, aes(x = 1, y = y)) +
  # title just above first entry
  geom_text(
    data = data.frame(x = 1, y = title_y, label = "Direction & age pattern"),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0,
    size  = 3,
    fontface = "bold"
  ) +
  geom_tile(
    aes(fill = cat),
    width  = 0.2,
    height = 0.15,
    colour = "black"
  ) +
  geom_text(
    data = subset(legend_df_compact, cat == "not_selected"),
    aes(label = "X"),
    color = "black",
    size  = 3
  ) +
  geom_text(
    aes(x = 1.15, label = label),
    hjust = 0,
    size  = 2.6
  ) +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_x_continuous(limits = c(0.8, 2), breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(limits = c(min_y - 1, max_y + 1), breaks = NULL, expand = c(0, 0)) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title.x  = element_blank(),
    axis.title.y  = element_blank(),
    axis.text     = element_blank(),
    axis.ticks    = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )
## =========================================================
## 8. Final combined figure (4 columns)
## =========================================================

## 2. Turn these into grobs (use the *tight* versions)
g_var <- ggplotGrob(p_var_key)
g_h6  <- ggplotGrob(p_heat_6_23)
g_h24 <- ggplotGrob(p_heat_24_59)

# White filler grob
white_panel <- rectGrob(gp = gpar(col = NA, fill = "white"))
g_legend <- arrangeGrob(
  grobs = list(
    white_panel,  # top white band
    p_legend,     # legend in the middle
    white_panel   # bottom white band
  ),
  ncol    = 1,
  heights = unit(c(1, 2, 1), "null")   # same relative heights as before
)

## 4. Assemble the final layout without drawing to an interactive device.
combined_model_plot <- arrangeGrob(
  g_var, g_h6, g_h24, g_legend,
  ncol   = 4,
  widths = unit(c(2.5, 4, 4 , 1.5), "null")
)

fig1s_path <- file.path(FIGURE_DIR, "Figure_1S_model_stability.png")

ggsave(
  filename = fig1s_path,
  plot     = combined_model_plot,
  width    = 20,   # adjust as you like
  height   = 8,    # adjust as you like
  dpi      = 300
)

## --------------------------------------------------
## 0. Helper: base skeleton from table1 + community
## --------------------------------------------------

# Individual-level skeleton
supp_indiv_core <- table1 %>%
  transmute(
    covariate_name   = Variable,
    recoded_variable = variable_code,
    response         = ifelse(is.na(response), NA_character_, as.character(response))
  )

supp_indiv <- bind_rows(
  tibble(
    covariate_name   = "Individual-level factors",
    recoded_variable = NA_character_,
    response         = NA_character_
  ),
  supp_indiv_core
)

# Community-level skeleton
supp_comm_core <- community_block %>%
  transmute(
    covariate_name   = Variable,
    recoded_variable = variable_code,
    response         = NA_character_
  )

supp_comm <- bind_rows(
  tibble(
    covariate_name   = "Community-level factors",
    recoded_variable = NA_character_,
    response         = NA_character_
  ),
  supp_comm_core
)

supp_skeleton <- bind_rows(supp_comm, supp_indiv)
## --------------------------------------------------
## 1. Definitions / recoding text
##    (fill as many as you like; others can be added later)
##    Key: recoded_variable + response
##    For parent rows, use response = NA_character_
## --------------------------------------------------
recode_defs <- tribble(
  ~recoded_variable, ~response, ~recoding,
  
  # =====================================================
  # Sex of child
  # =====================================================
  "sex_female", NA_character_,
  "Derived from hv104 (sex of child).",
  
  "sex_female", "1",
  "sex_female = 1 (Female) if hv104 = 2 (female); sex_female = 0 (Male) if hv104 = 1 (male).",
  
  # =====================================================
  # Household wealth quintile (ph_wealth_quint)
  # =====================================================
  "ph_wealth_quint", NA_character_,
  "Derived from hv270 (household wealth index quintile) among de jure household members (hv102 = 1).",
  
  "ph_wealth_quint", "5",
  "ph_wealth_quint = 5 (Richest) if hv270 = 5 (highest household wealth index quintile) AND hv102 = 1 (de jure household member).",
  
  "ph_wealth_quint", "4",
  "ph_wealth_quint = 4 (Richer) if hv270 = 4 (second-highest household wealth index quintile) AND hv102 = 1 (de jure household member).",
  
  "ph_wealth_quint", "3",
  "ph_wealth_quint = 3 (Middle) if hv270 = 3 (middle household wealth index quintile) AND hv102 = 1 (de jure household member).",
  
  "ph_wealth_quint", "2",
  "ph_wealth_quint = 2 (Poorer) if hv270 = 2 (second-lowest household wealth index quintile) AND hv102 = 1 (de jure household member).",
  
  "ph_wealth_quint", "1",
  "ph_wealth_quint = 1 (Poorest) if hv270 = 1 (lowest household wealth index quintile) AND hv102 = 1 (de jure household member).",
  
  # =====================================================
  # Maternal education (mother_eduyr, hc61)
  # =====================================================
  "mother_eduyr", NA_character_,
  "Recode of hc61 (mother’s highest education level).",
  
  "mother_eduyr", "3",
  "mother_eduyr = 3 (Higher) if hc61 = 3 (higher education).",
  
  "mother_eduyr", "2",
  "mother_eduyr = 2 (Secondary) if hc61 = 2 (secondary education).",
  
  "mother_eduyr", "1",
  "mother_eduyr = 1 (Primary) if hc61 = 1 (primary education).",
  
  "mother_eduyr", "0",
  "mother_eduyr = 0 (None) if hc61 ∈ {0, 8, 9} (no education, don’t know, or missing recoded as no education).",
  
  # =====================================================
  # Paternal/partner education (partner_eduyr, v701)
  # =====================================================
  "partner_eduyr", NA_character_,
  "Recode of v701 (partner’s highest education level).",
  
  "partner_eduyr", "3",
  "partner_eduyr = 3 (Higher) if v701 = 3 (higher education).",
  
  "partner_eduyr", "2",
  "partner_eduyr = 2 (Secondary) if v701 = 2 (secondary education).",
  
  "partner_eduyr", "1",
  "partner_eduyr = 1 (Primary) if v701 = 1 (primary education).",
  
  "partner_eduyr", "0",
  "partner_eduyr = 0 (None) if v701 ∈ {0, 8, 9} (no education, don’t know, or missing recoded as no education).",
  
  # =====================================================
  # Maternal occupation (mother_occupation, v717)
  # =====================================================
  "mother_occupation", NA_character_,
  "Recode of v717 (mother’s current occupation) into three groups.",
  
  "mother_occupation", "1",
  "mother_occupation = 1 (Office work) if v717 ∈ {1, 2, 3, 6, 7} (professional/technical/managerial, clerical, sales, services, related white-collar jobs).",
  
  "mother_occupation", "2",
  "mother_occupation = 2 (Manual labor) if v717 ∈ {4, 5, 8, 9} (agricultural, skilled and unskilled manual, or other blue-collar work).",
  
  "mother_occupation", "0",
  "mother_occupation = 0 (None) if v717 ∈ {0, 98, 99} or v717 indicates no occupation/household & domestic only; missing values recoded as 0.",
  
  # =====================================================
  # Partner occupation (partner_occupation, v705)
  # =====================================================
  "partner_occupation", NA_character_,
  "Recode of v705 (partner’s current occupation) into three groups.",
  
  "partner_occupation", "1",
  "partner_occupation = 1 (Office work) if v705 ∈ {1, 2, 3, 7} (professional/technical/managerial, clerical, sales, services).",
  
  "partner_occupation", "2",
  "partner_occupation = 2 (Manual labor) if v705 ∈ {4, 5, 8, 9} (agricultural, skilled and unskilled manual, or other blue-collar work).",
  
  "partner_occupation", "0",
  "partner_occupation = 0 (None) if v705 ∈ {0, 6, 98, 99} (no occupation, household/domestic only, don’t know or missing).",
  
  # =====================================================
  # Region & rural residence
  # =====================================================
  "geography_regions", NA_character_,
  "Recode of hv024 (region of residence) into five macro-regions of Nigeria.",
  
  "geography_regions", "1",
  "geography_regions = 1 (North Central) if hv024 = 1.",
  
  "geography_regions", "2",
  "geography_regions = 2 (North East) if hv024 = 2.",
  
  "geography_regions", "3",
  "geography_regions = 3 (North West) if hv024 = 3.",
  
  "geography_regions", "4",
  "geography_regions = 4 (South South or South East) if hv024 ∈ {4, 5}.",
  
  "geography_regions", "5",
  "geography_regions = 5 (South West) if hv024 = 6.",
  
  "geography_rural", NA_character_,
  "Recode of hv025 (type of place of residence) into rural vs non-rural.",
  
  "geography_rural", "1",
  "geography_rural = 1 (Lives in rural regions) if hv025 = 2 (rural); geography_rural = 0 (Does not live in rural regions) if hv025 = 1 (urban).",
  
  # =====================================================
  # MDD / MMF / MAD (child feeding indicators)
  # =====================================================
  "nt_mdd", NA_character_,
  "Minimum dietary diversity (MDD) indicator for last-born children aged 6–23 months, derived from v4xx feeding variables and WHO 8-group definition.",
  
  "nt_mdd", "1",
  "nt_mdd = 1 (Yes) if child aged 6–23 months consumed foods from ≥5 of 8 WHO food groups in the previous day/night based on nt_* food group indicators; nt_mdd = 0 (No) otherwise.",
  
  "nt_mmf", NA_character_,
  "Minimum meal frequency (MMF) indicator for last-born children aged 6–23 months, derived from m39, m4 and total milk feeds.",
  
  "nt_mmf", "1",
  "nt_mmf = 1 (Yes) if child aged 6–23 months meets age- and breastfeeding-status specific WHO minimum meal frequency recommendations using m39 (solid/semi-solid feeds), totmilkf and age; nt_mmf = 0 (No) otherwise.",
  
  "nt_mad", NA_character_,
  "Minimum acceptable diet (MAD) indicator combining MDD and MMF for children aged 6–23 months.",
  
  "nt_mad", "1",
  "nt_mad = 1 (Yes) if child aged 6–23 months meets WHO minimum acceptable diet criteria: breastfed children with nt_mdd = 1 AND nt_mmf = 1; non-breastfed children with sufficient solid food diversity (foodsum2 ≥ 4), nt_mmf = 1 and ≥2 milk feeds; nt_mad = 0 (No) otherwise.",
  
  # =====================================================
  # Micronutrient supplementation (iron, vitamin A)
  # =====================================================
  "iron", NA_character_,
  "Child received iron supplements in the 7 days preceding the interview, based on h42.",
  
  "iron", "1",
  "iron = 1 (Yes) if h42 = 1 (received iron in last 7 days); iron = 0 (No) otherwise.",
  
  "vit_a", NA_character_,
  "Child received vitamin A supplement in the 6 months preceding the interview, based on v008a and h34.",
  
  "vit_a", "1",
  "vit_a = 1 (Yes) if v008a < 7 (months since last vitamin A campaign ≤ 6) OR h34 = 1 (received vitamin A); vit_a = 0 (No) otherwise.",
  
  # =====================================================
  # Maternal anemia (mother_anemia, v457)
  # =====================================================
  "mother_anemia", NA_character_,
  "Four-category anemia status of the mother, recoded from v457.",
  
  "mother_anemia", "4",
  "mother_anemia = 4 (Normal) if v457 ∈ {4, 9} (no anemia or not anemic; missing recoded as normal).",
  
  "mother_anemia", "3",
  "mother_anemia = 3 (Mild anemia) if v457 = 3 (mild anemia).",
  
  "mother_anemia", "2",
  "mother_anemia = 2 (Moderate anemia) if v457 = 2 (moderate anemia).",
  
  "mother_anemia", "1",
  "mother_anemia = 1 (Severe anemia) if v457 = 1 (severe anemia).",
  
  # =====================================================
  # Maternal BMI (mother_bmi, v445)
  # =====================================================
  "mother_bmi", NA_character_,
  "Body mass index category of the mother derived from v445 (BMI × 100).",
  
  "mother_bmi", "3",
  "mother_bmi = 3 (Normal) if v445 ∈ [1850, 2499] corresponding to BMI 18.5–24.9 kg/m²; missing values recoded as normal (3).",
  
  "mother_bmi", "1",
  "mother_bmi = 1 (Very thin) if v445 ≤ 1699 (BMI < 17).",
  
  "mother_bmi", "2",
  "mother_bmi = 2 (Thin) if 1700 ≤ v445 ≤ 1849 (BMI 17.0–18.49).",
  
  "mother_bmi", "4",
  "mother_bmi = 4 (Overweight) if 2500 ≤ v445 ≤ 2999 (BMI 25.0–29.9).",
  
  "mother_bmi", "5",
  "mother_bmi = 5 (Obese) if v445 ≥ 3000 (BMI ≥ 30.0).",
  
  # =====================================================
  # Child anthropometry (height-for-age, weight-for-height, weight-for-age)
  # =====================================================
  "nt_ch_ht_age", NA_character_,
  "Child height-for-age category derived from hc70 (height-for-age z-score × 100).",
  
  "nt_ch_ht_age", "1",
  "nt_ch_ht_age = 1 (With stunting) if hv103 = 1 (child present) AND hc70 < -200 (height-for-age z-score < -2 SD).",
  
  "nt_ch_ht_age", "0",
  "nt_ch_ht_age = 0 (Normal) if hv103 = 1 AND hc70 ≥ -200; observations with hc70 ≥ 9996 are set to missing.",
  
  "nt_ch_wt_ht", NA_character_,
  "Child weight-for-height category derived from hc72 (weight-for-height z-score × 100).",
  
  "nt_ch_wt_ht", "2",
  "nt_ch_wt_ht = 2 (With wasting) if hv103 = 1 AND hc72 < -200 (weight-for-height z-score < -2 SD).",
  
  "nt_ch_wt_ht", "1",
  "nt_ch_wt_ht = 1 (Normal) if hv103 = 1 AND -200 ≤ hc72 ≤ 200 (z-score between -2 and +2 SD).",
  
  "nt_ch_wt_ht", "0",
  "nt_ch_wt_ht = 0 (Overweight) if hv103 = 1 AND hc72 > 200 & hc72 < 9990 (z-score > +2 SD). Observations with hc72 ≥ 9996 are set to missing.",
  
  "nt_ch_wt_age", NA_character_,
  "Child weight-for-age category derived from hc71 (weight-for-age z-score × 100).",
  
  "nt_ch_wt_age", "2",
  "nt_ch_wt_age = 2 (Underweight) if hv103 = 1 AND hc71 < -200 (weight-for-age z-score < -2 SD).",
  
  "nt_ch_wt_age", "1",
  "nt_ch_wt_age = 1 (Normal) if hv103 = 1 AND -200 ≤ hc71 ≤ 200 (z-score between -2 and +2 SD).",
  
  "nt_ch_wt_age", "0",
  "nt_ch_wt_age = 0 (Overweight for age) if hv103 = 1 AND hc71 > 200 & hc71 < 9990 (z-score > +2 SD).",
  
  # =====================================================
  # Currently breastfed (nt_bf_curr, m4)
  # =====================================================
  "nt_bf_curr", NA_character_,
  "Currently breastfeeding indicator derived from m4 (current breastfeeding status).",
  
  "nt_bf_curr", "1",
  "nt_bf_curr = 1 (Currently breastfed) if m4 = 95 (child is still breastfeeding); nt_bf_curr = 0 (Not currently breastfed) if m4 ∈ {93, 94, 98, 99} or child never/previously breastfed.",
  
  # =====================================================
  # Household size (hh_members, hv009)
  # =====================================================
  "hh_members", NA_character_,
  "Household size category based on hv009 (number of de jure household members).",
  
  "hh_members", "1",
  "hh_members = 1 (2–4 members) if hv009 ∈ [2, 4].",
  
  "hh_members", "2",
  "hh_members = 2 (5–7 members) if hv009 ∈ [5, 7].",
  
  "hh_members", "3",
  "hh_members = 3 (8 or more members) if hv009 ≥ 8.",
  
  # =====================================================
  # ITN-related variables
  # =====================================================
  "ml_itnhh", NA_character_,
  "Household owns at least one insecticide-treated net (ITN), based on hml10_1–hml10_7.",
  
  "ml_itnhh", "1",
  "ml_itnhh = 1 (Yes) if any of hml10_1–hml10_7 = 1 (household owns ≥1 ITN); ml_itnhh = 0 (No) otherwise.",
  
  "ml_hhaccess", NA_character_,
  "Household has more than one ITN per two household members, derived from number of ITNs and household size hv013.",
  
  "ml_hhaccess", "1",
  "ml_hhaccess = 1 (Yes) if hv013 > 0 AND 2 × (number of ITNs from hml10_1–hml10_7) ÷ hv013 ≥ 1; ml_hhaccess = 0 (No) if hv013 > 0 AND 2 × ITNs ÷ hv013 < 1.",
  
  "ml_netcat", NA_character_,
  "Child slept under an insecticide-treated net (ITN) the previous night, based on hml12.",
  
  "ml_netcat", "1",
  "ml_netcat = 1 (Yes) if hml12 ∈ {1, 2, 3} (slept under an ITN the previous night); ml_netcat = 0 (No) if hml12 = 0.",
  
  # =====================================================
  # WASH & housing: electricity, floor, rooms, cooking, land, animals
  # =====================================================
  "ph_electric", NA_character_,
  "Indicator for household electricity, based on hv206.",
  
  "ph_electric", "1",
  "ph_electric = 1 (Yes) if hv206 = 1 (household has electricity); ph_electric = 0 (No) if hv206 = 0.",
  
  "ph_floor", NA_character_,
  "Binary indicator for finished flooring material, recoded from hv213 (main floor material).",
  
  "ph_floor", "1",
  "ph_floor = 1 (Finished flooring material) if hv213 ∈ {30–35} (parquet/polished wood, vinyl, tiles, cement, carpet, etc.); ph_floor = 0 (Natural or rudimentary) if hv213 ∈ {10–22}.",
  
  "ph_rooms_sleep", NA_character_,
  "Number of rooms used for sleeping (hv216) recoded to <2 vs ≥2 rooms.",
  
  "ph_rooms_sleep", "1",
  "ph_rooms_sleep = 1 (Two or more rooms) if hv216 ≥ 2; ph_rooms_sleep = 0 (One room) if hv216 = 1.",
  
  "ph_cook_place", NA_character_,
  "Place where food is cooked, derived from hv241 and hv226, recoded as indoors/none vs outdoors/other.",
  
  "ph_cook_place", "1",
  "ph_cook_place = 1 (Inside a building / no food cooked) if hv241 ∈ {1, 2, 6} OR hv226 = 95 (no food cooked in household); corresponds to cooking inside dwelling, in a separate building, or no cooking.",
  
  "ph_cook_place", "0",
  "ph_cook_place = 0 (Outdoors) if hv241 = 3 (food cooked outdoors). Other hv241 codes are grouped as 'Other' or missing.",
  
  "ph_agriland", NA_character_,
  "Ownership of land usable for agriculture, based on hv244.",
  
  "ph_agriland", "1",
  "ph_agriland = 1 (Yes) if hv244 = 1 (household owns land usable for agriculture); ph_agriland = 0 (No) if hv244 = 0 or missing.",
  
  "ph_animals", NA_character_,
  "Household owns livestock or farm animals, based on hv246.",
  
  "ph_animals", "1",
  "ph_animals = 1 (Yes) if hv246 = 1 (household has livestock, herds or farm animals); ph_animals = 0 (No) if hv246 = 0 or missing.",
  
  "sanitation_imp", NA_character_,
  "Access to improved sanitation facility, recoded from hv205 (type of toilet facility).",
  
  "sanitation_imp", "1",
  "sanitation_imp = 1 (Yes) if hv205 ∈ {11, 12, 13, 15, 21, 22, 41}, corresponding to improved sanitation per DHS (flush/pour flush to piped sewer/septic tank/pit latrine, ventilated improved pit latrine, pit latrine with slab, etc.).",
  
  "sanitation_imp", "0",
  "sanitation_imp = 0 (No) for all other hv205 codes (unimproved facility or no facility).",
  
  "share_toilet", NA_character_,
  "Household shares toilet facilities with other households, based on v160.",
  
  "share_toilet", "1",
  "share_toilet = 1 (Shares toilet with others) if v160 ∈ {1, 7}; share_toilet = 0 (Does not share) if v160 = 0. Other codes are treated as missing.",
  
  "handwash_place", NA_character_,
  "Place where household members most often wash their hands, recoded from hv230a.",
  
  "handwash_place", "1",
  "handwash_place = 1 (At home: fixed) if hv230a = 1 (fixed place in dwelling/yard).",
  
  "handwash_place", "2",
  "handwash_place = 2 (At home: mobile) if hv230a = 2 (mobile handwashing object such as a jug or bucket).",
  
  "handwash_place", "3",
  "handwash_place = 3 (Not at home / not observed) if hv230a is any other code or handwashing place is not observed.",
  
  "soap_water", NA_character_,
  "Presence of both soap and water at the observed handwashing place, based on hv230b and hv232.",
  
  "soap_water", "1",
  "soap_water = 1 (Yes) if hv230b = 1 (water present) AND hv232 = 1 (soap or detergent present); soap_water = 0 (No) otherwise.",
  
  "handwash_dirty", NA_character_,
  "Presence of ash, mud, or sand at the handwashing place, based on hv232b.",
  
  "handwash_dirty", "1",
  "handwash_dirty = 1 (Yes) if hv232b = 1 (ash, mud, or sand present); handwash_dirty = 0 (No) otherwise.",
  
  "water_imp", NA_character_,
  "Access to improved drinking-water source, recoded from hv201 (main source of drinking water).",
  
  "water_imp", "1",
  "water_imp = 1 (Yes) if hv201 ∈ {11, 14, 21, 31, 41, 51, 61, 62, 71}, corresponding to improved water sources per DHS (piped, tube well/borehole, protected well/spring, rainwater, bottled water, etc.).",
  
  "water_imp", "0",
  "water_imp = 0 (No) for all other hv201 codes (unimproved sources).",
  
  "own_water", NA_character_,
  "Location of drinking-water source, based on hv235.",
  
  "own_water", "1",
  "own_water = 1 (Own dwelling/yard) if hv235 ∈ {1, 2} (source located in respondent’s dwelling or yard/plot); own_water = 0 (Elsewhere) otherwise.",
  
  # =====================================================
  # Mass media & internet
  # =====================================================
  "internet_fq", NA_character_,
  "Mother used the internet at least once a week, based on v171b.",
  
  "internet_fq", "1",
  "internet_fq = 1 (Yes) if v171b ∈ {2, 3} (uses internet at least once a week); internet_fq = 0 (No) otherwise.",
  
  "mass_media", NA_character_,
  "Mother exposed to radio, newspaper, or television at least once a week, based on v157, v158, and v159.",
  
  "mass_media", "1",
  "mass_media = 1 (Yes) if any of v157, v158, or v159 ∈ {2, 3} (listens/reads/watches at least once a week); mass_media = 0 (No) otherwise.",
  
  # =====================================================
  # Health-seeking & preventive care: visits, access, stool disposal, deworming
  # =====================================================
  "visit_health", NA_character_,
  "Visited a health facility in the last 12 months, based on v394.",
  
  "visit_health", "1",
  "visit_health = 1 (Yes) if v394 = 1 (visited a health facility in the last 12 months); visit_health = 0 (No) otherwise.",
  
  "access_health", NA_character_,
  "Mother reports major problems in accessing healthcare, based on v467a, v467b, v467c, v467d, and v467f.",
  
  "access_health", "1",
  "access_health = 1 (Yes) if any of v467a, v467b, v467c, v467d, or v467f = 1 (big problem in at least one access dimension); access_health = 0 (No) if all specified items ≠ 1.",
  
  "stool_disposal", NA_character_,
  "Disposal of youngest child’s stool when not using toilet, based on v465.",
  
  "stool_disposal", "1",
  "stool_disposal = 1 (Appropriate) if v465 ∈ {1, 2} (disposed into toilet/latrine or rinsed into toilet/latrine); stool_disposal = 0 (Inappropriate) for all other v465 codes (e.g., thrown into garbage, left in open, buried, rinsed into drain/ditch).",
  
  "preg_giparasites", NA_character_,
  "Mother took drugs for intestinal parasites during pregnancy, based on m60.",
  
  "preg_giparasites", "1",
  "preg_giparasites = 1 (Yes) if m60 = 1 (took drugs for intestinal parasites during pregnancy); preg_giparasites = 0 (No) otherwise.",
  
  # =====================================================
  # Place of delivery (delivery_place, m15)
  # =====================================================
  "delivery_place", NA_character_,
  "Place of delivery of the index child, recoded from m15 into three categories.",
  
  "delivery_place", "1",
  "delivery_place = 1 (Own or other's home) if m15 ∈ {10, 11, 12, 96} (respondent’s home, other home, or other home-like locations).",
  
  "delivery_place", "2",
  "delivery_place = 2 (Government facility) if m15 ∈ {20, 21, 22, 23, 26} (public/government hospital, health center, clinic or other public facility).",
  
  "delivery_place", "3",
  "delivery_place = 3 (Private facility) if m15 ∈ {30, 31, 36} (private hospital/clinic, mission/NGO facility or other private facility).",
  
  # =====================================================
  # Community-level variables (continuous 0–1 proportions)
  # =====================================================
  "comm_pov_prop", NA_character_,
  "Continuous (0–1): For each PSU (v021), survey-weighted proportion of children living in households with ph_wealth_quint ∈ {1, 2} (poor or poorer).",
  
  "comm_motheredu_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children whose mother has ≥ primary education (mother_eduyr ≥ 1).",
  
  "comm_partneredu_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children whose partner has ≥ primary education (partner_eduyr ≥ 1).",
  
  "comm_mother_occupation_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children whose mother is employed in any occupation (mother_occupation ≠ 0).",
  
  "comm_partner_occupation_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children whose partner is employed in any occupation (partner_occupation ≠ 0).",
  
  "comm_animals_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children living in households that own livestock or farm animals (ph_animals = 1).",
  
  "comm_agriland_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children living in households that own agricultural land (ph_agriland = 1).",
  
  "comm_sanitation_imp_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children living in households with improved sanitation (sanitation_imp = 1).",
  
  "comm_water_imp_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children living in households with improved drinking-water source (water_imp = 1).",
  
  "comm_floor_finished_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children living in households with finished flooring material (ph_floor = 1).",
  
  "comm_electric_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children living in households with electricity (ph_electric = 1).",
  
  "comm_massmedia_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children whose mother is exposed to mass media at least once a week (mass_media = 1).",
  
  "comm_itnhh_prop", NA_character_,
  "Continuous (0–1): For each PSU, survey-weighted proportion of children living in households owning at least one ITN (ml_itnhh = 1)."
)

## Attach to skeleton and fix 'response' text for community variables
supp_table1 <- supp_skeleton %>%
  left_join(
    recode_defs,
    by = c("recoded_variable", "response")
  ) %>%
  mutate(
    response = if_else(
      recoded_variable %in% c(
        "comm_pov_prop", "comm_motheredu_prop", "comm_partneredu_prop",
        "comm_mother_occupation_prop", "comm_partner_occupation_prop",
        "comm_animals_prop", "comm_agriland_prop", "comm_sanitation_imp_prop",
        "comm_water_imp_prop", "comm_floor_finished_prop",
        "comm_electric_prop", "comm_massmedia_prop", "comm_itnhh_prop"
      ) & is.na(response),
      "Continuous (0–1) proportion",
      response
    )
  )

colnames(supp_table1) <- c(
  "Covariate name",
  "Recoded variable",
  "Response / unit",
  "Definition and recoding rule"
)

table1s <- supp_table1

## --------------------------------------------------
## Final tables for export
## --------------------------------------------------

table1_comm <- community_block[, c("Variable", "Median proportion (25%-75%)")]
colnames(table1_comm) <- c('Covariate name', "Median proportion (25%-75%)")

table1_final <- table1[, c("Variable", "Younger children N (%)", "Older children N (%)")]
colnames(table1_final) <- c('Covariate name', "Younger children N (%)", "Older children N (%)")

table2 <- final_results_table[, c('Variable', 'Age pattern', 'PR (6–23.9 months)', 'PR (24–59 months)')]
colnames(table2) <- c('Covariate name', 'Age pattern', 'PR (6–23.9 months)', 'PR (24–59 months)')

table2s <- var_selection_table_all[, c(
  "Variable_label",
  "age_pattern",
  "step1_status",
  "step2_status",
  "backward_main_result_step3",
  "backward_int_result_step3",
  "forward_main_result_step3",
  "forward_int_result_step3"
)]

colnames(table2s) <- c(
  "Covariate name",
  "Final age pattern",
  "Step 1: LASSO & bootstrapping",
  "Step 2: GLM & bootstrapping",
  "Step 3: backward elimination - main effect",
  "Step 3: backward elimination - interaction term",
  "Step 3: forward reassessment - main effect",
  "Step 3: forward reassessment - interaction term"
)


table3s <- regression_table_with_headings[, c('variable_name', 'PR_CI', 'estimate_log','std_error', 'z_value', 'p_value_fmt' )]
names(table3s) <- c('Covariate name', 'Exponentiated coefficient (95% CI)', 'Coefficient', 'Standard error', 'Z-value', 'p-value')
table3s <- table3s[-1,]

# and if you want for Table 4s:
table4s <- model_selection_table

wb <- createWorkbook()

addWorksheet(wb, 'Table 1')
writeData(wb, 'Table 1', table1_final)

addWorksheet(wb, 'Table 1 comm')
writeData(wb, 'Table 1 comm', table1_comm)

addWorksheet(wb, 'Table 2')
writeData(wb, 'Table 2', table2)

addWorksheet(wb, 'Table 1S')
writeData(wb, 'Table 1S', table1s)

addWorksheet(wb, 'Table 2S')
writeData(wb, 'Table 2S', table2s)

addWorksheet(wb, 'Table 3S')
writeData(wb, 'Table 3S', table3s)

addWorksheet(wb, 'Table 4S')
writeData(wb, 'Table 4S', table4s)

addWorksheet(wb, "Figure 1S")
insertImage(
  wb,
  sheet      = "Figure 1S",
  file       = fig1s_path,
  startRow   = 3,
  startCol   = 2,
  width      = 10,
  height     = 5
)

addWorksheet(wb, "Figure 2S")
insertImage(
  wb,
  sheet      = "Figure 2S",
  file       = fig2s_path,
  startRow   = 3,
  startCol   = 2,
  width      = 10,
  height     = 5
)
save(
  step1_results, step2_results, step3_results,
  final_results_table, var_selection_table_all,
  table1, table1_comm, table2, table4s,
  file = file.path(RESULT_DIR, "cms_final_results.RData")
)

saveWorkbook(
  wb,
  file = file.path(OUTPUT_DIR, "cms_tables_and_figures.xlsx"),
  overwrite = TRUE
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUTPUT_DIR, "session_info.txt")
)
