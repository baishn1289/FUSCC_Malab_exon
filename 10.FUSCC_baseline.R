setwd()

rm(list = ls())
library(tableone) 

load('FUSCC_2024.Rdata')
load('FUSCC_2022.Rdata')

dup_number = na.omit(match(cl_FUSCC_2022$Tumor_Sample_Barcode,cl_FUSCC_2024$Tumor_Sample_Barcode))
dup_patients = cl_FUSCC_2024[dup_number,]
table(dup_patients$Tumor_Sample_Barcode %in% cl_FUSCC_2022$Tumor_Sample_Barcode)
cl_FUSCC_2024 = cl_FUSCC_2024[-dup_number,]

rt_FUSCC_2024 = subsetMaf(rt_FUSCC_2024,tsb=cl_FUSCC_2024$Tumor_Sample_Barcode)

rt_FUSCC_2020to2024 = merge_mafs(list(rt_FUSCC_2022,rt_FUSCC_2024))
cl_FUSCC_2020to2024 <- rbind(rt_FUSCC_2022@clinical.data,rt_FUSCC_2024@clinical.data,fill=T)
write.csv(cl_FUSCC_2020to2024,file = 'FUSCC_2020to2024.tsv',row.names = F)
save(cl_FUSCC_2020to2024,rt_FUSCC_2020to2024,file="FUSCC_2020to2024.Rdata")

#baseline
b = fread('FUSCC_2020to2024.tsv')
names(b)

data_selected <- b[,c(3:5,8,17,35)]

names(data_selected)

Vars = c('Sex',"Age","MMR_dMMR/pMMR","Tumor site", "Stage",  "Age_class" )
data_selected$Age <- as.numeric(data_selected$Age)
table1 <- CreateTableOne(vars = Vars, 
                         data = data_selected,
                         strata =   'Age_class'#, factorVars = 'Age'
                         )
print(table1)

tableCSV <- print(table1, showAllLevels = T,
      pDigits = 4, catDigits = 1,
      contDigits = 1 ,
      printToggle = FALSE,
      addOverall = TRUE, showpm = TRUE)
write.csv(tableCSV, file = "Table 1.csv")

summary(table1)
