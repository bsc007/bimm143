#' ---
#' title: "Crop Analysis Q3 2013"
#' author: "John Smith"
#' date: "May 3rd, 2014"
output: github_document
#' ---

# Class 5 R graphics

#2A. Line plot
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header=TRUE)

plot(weight$Age, weight$Weight, xlab="Age (months)", ylab="weight(kg)", 
     main="Weight of Growing Child","l", 
     pch=15, cex=1.5, 
     lwd=2, ylim=c(2,10))

#2B
features <- read.table("bimm143_05_rstats/feature_counts.txt", sep="\t",header=TRUE)
features

par(mar=c(5,11,4,2))
barplot(features$Count, horiz=TRUE, names.arg = features$Feature,las=1)
        
#3A
gender <- read.delim("bimm143_05_rstats/male_female_counts.txt", header=TRUE)
gender
par(mar=c(7,4,4,2))
barplot(gender$Count, names.arg=gender$Sample, las=2, col=c("pink","plum"))

barplot(gender$Count, names.arg=gender$Sample, las=2, col=c(1,2))

barplot(gender$Count, names.arg=gender$Sample, las=2, col=(rainbow(10)))
