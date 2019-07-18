% Downlaod all the mutations and copy number files as *.gz files which are small using MATLAB

curl -XPOST -i -H 'Content-Type: application/json' 'https://dcc.icgc.org/api/v1/auth/login' -d '{"username":"smsinks@gmail.com","password":"uHC-9A5-BeL-oP8"}'

curl -L -o simple_somatic_mutation.open.COAD-US.tsv.gz https://dcc.icgc.org/api/v1/download?fn=/release_16/Projects/COAD-US/simple_somatic_mutation.open.COAD-US.tsv.gz


% Use R to Run through each ICGC mutations and copy number file to convert them from a simple file to a MAF file
setwd("~/Documents/")
ldf <- list() # creates a list
listcsv <- dir(pattern = "*.csv") # creates the list of all the csv files in the directory
for (k in 1:length(listcsv)){
 ldf[[k]] <- read.csv(listcsv[k])
}
str(ldf[[1]]) 