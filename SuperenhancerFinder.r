library(GenomicRanges)
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)

reps = list(c('MAE2_1_F','MAE4_3_F'),c('MAS1_4_F','MAS3_2_F'),c('YE2_1_F','YE4_3_F'),c('YS1_4_F','YS3_2_F'),c('MAE2_1_M','MAE4_3_M'),c('MAS1_4_M','MAS3_2_M'),c('YE2_1_M','YE4_3_M'),c('YS1_4_M','YS3_2_M'))

cells = c('M1','M2','MSC','LSC','Dendritic','NK','Cardiomyocyte','Coronary_EC','Endocardial_EC','Lymphatic_EC','Fibroblast','Mesothelial','B_cell','T_cell','Pericyte','Smooth_muscle','Cardiac_Neuronal') 
#Cardiomyoctye misspelt in directory with ATAC data. 
# Was Pericytes, but Pericyte in directory
# Was Smooth_Muscle, but Smooth_muscle in directory

fails = c()

groups = (c(cells[1], c(unlist(reps[1]))))
for (cell in cells) {
    for (rep in reps) {
        group = c(cell, c(unlist(rep)))
        groups = rbind(groups, group)
        }
    }
groups = groups[87:nrow(groups), ]

# RefSeq genes are output of Annotation() function from Seurat
ref.gr=readRDS("/g/data/zk16/arthuss/superenhancer/RefSeq_mm10.rds")
# get TSS
tss.gr=ref.gr
# end of the + strand genes must be equalized to start pos
end(tss.gr[strand(tss.gr)=="+",])  =start(tss.gr[strand(tss.gr)=="+",])
# startof the - strand genes must be equalized to end pos
start(tss.gr[strand(tss.gr)=="-",])=end(tss.gr[strand(tss.gr)=="-",])
# remove duplicated TSSes ie alternative transcripts
# this keeps the first instance and removes duplicates
tss.gr=tss.gr[!duplicated(tss.gr),]
print("Imported TSS")

# Function for finding cutoff value for what counts as a superenhancer, slides a diagonal. Directly from ROSE source code (slightly modified). 
calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
 	inputVector <- sort(inputVector)
	inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
	slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
	xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
	
	if(drawPlot){  #if TRUE, draw the plot
		plot(1:length(inputVector), inputVector,type="l",...)
		b <- y_cutoff-(slope* xPt)
		abline(v= xPt,h= y_cutoff,lty=2,col=8)
		points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
		abline(coef=c(b,slope),col=2)
		title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",signif(y_cutoff/mean(inputVector),3),"x",sep=""))
		axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
	}
	return(y_cutoff) #return(list(absolute=y_cutoff,overMedian=y_cutoff/median(inputVector),overMean=y_cutoff/mean(inputVector), cutoff = xPt)) - modified from original ROSE script, https://bitbucket.org/young_computation/rose/src/master/ROSE_callSuper.R
}

# Also from ROSE, keeping because it didn't work without it for some reason
#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}
print("Defined slope functions")

for (i in c(2:nrow(groups))) {
    args = groups[i, ]
    print(args)

#args = commandArgs(trailingOnly=TRUE) # For testing use c("M1", "MAE2_1_F", "MAE4_3_F")

if (args[1] == "Cardiomyocyte") {
    reads = readRDS("/g/data/zk16/jack/HeartSCRNA-Jun22/analysis/2023-10-11.AllMacs/combined_only_Cardiomyoctye.rds")
    } else {
reads = readRDS(paste("/g/data/zk16/jack/HeartSCRNA-Jun22/analysis/2023-10-11.AllMacs/combined_only_", args[1], ".rds", sep = ""))
    }
print("Imported reads")

# 1.Input read count matrices for two replicates of the same treatment+sex
# 2.Exclude intervals with fewer reads than the most frequent read number for that individual (i.e., if 6 is the most common number of reads across all intervals, all intervals with <6 reads are excluded)
# 3.Take the union of intervals for the two replicate individuals (i.e., if an interval is only present in one of the replicates after the filtering above, it is kept)
individual1 = args[2]
individual2 = args[3]
if (individual1 %in% as.data.frame(table(reads$mouseID))[, 1] & individual2 %in% as.data.frame(table(reads$mouseID))[, 1]) {
    print("Individuals in reads file")
    } else {
    print("Individuals not in reads file")
    fails = append(fails, args)
    next
    }
reads.gr = reads@assays$peaks@ranges
reads.gr1 = reads@assays$peaks@ranges
reads.gr2 = reads@assays$peaks@ranges
readcounts1 = subset(x = reads, subset = mouseID == individual1)@assays$peaks@data
readcounts2 = subset(x = reads, subset = mouseID == individual2)@assays$peaks@data
ReadSums1 = as.integer(rowSums(readcounts1))
ReadSums2 = as.integer(rowSums(readcounts2))
mcols(reads.gr1)$TotalMappingReads = ReadSums1
mcols(reads.gr2)$TotalMappingReads = ReadSums2
# Two options for filtering read counts here. One is setting a hard limit for the number of reads per peak interval pseudobulked for all cells. Option 2 is keeping all read numbers above the most common read number for a given individual (commented out). 
readCutoff = 10
reads.gr1 = subset(x = reads.gr1, subset = TotalMappingReads >= readCutoff)
reads.gr2 = subset(x = reads.gr2, subset = TotalMappingReads >= readCutoff)
# reads.gr1 = subset(x = reads.gr1, subset = TotalMappingReads > 0)
# reads.gr2 = subset(x = reads.gr2, subset = TotalMappingReads > 0)
# reads.gr1 = reads.gr1[reads.gr1$TotalMappingReads >= as.numeric(levels(droplevels(as.data.frame(table(reads.gr1$TotalMappingReads))[, 1][as.data.frame(table(reads.gr1$TotalMappingReads))[, 2] == max(as.data.frame(table(reads.gr1$TotalMappingReads))[, 2])]))), ]
# reads.gr2 = reads.gr2[reads.gr2$TotalMappingReads >= as.numeric(levels(droplevels(as.data.frame(table(reads.gr2$TotalMappingReads))[, 1][as.data.frame(table(reads.gr2$TotalMappingReads))[, 2] == max(as.data.frame(table(reads.gr2$TotalMappingReads))[, 2])]))), ]
reads.gr = reads.gr[reads.gr %in% reads.gr1[!(reads.gr1 %in% reads.gr2)] | reads.gr %in% reads.gr2[!(reads.gr2 %in% reads.gr1)] | reads.gr %in% reads.gr1[reads.gr1 %in% reads.gr2]]
print("Union of replicate intervals done")
saveRDS(reads.gr, paste("/g/data/zk16/arthuss/superenhancer/SecondRun/", as.character(args[1]), "_", as.character(individual1), "_", as.character(individual2), "_readCutoff_", as.character(readCutoff), "_", "IntervalUnion.gr.rds", sep = ""))

# 4.Exclude intervals overlapping any TSS
reads.gr = reads.gr[!reads.gr %in% subsetByOverlaps(reads.gr, tss.gr)]
if (length(reads.gr) == 0) {
    print("No peaks left after excluding TSS")
    fails = append(fails, args)
    next
    }
non_tss_peak_ranges = c()
for (i in c(1:length(reads.gr))) {
    range = paste(as.character(seqnames(reads.gr[i, ])), as.character(start(reads.gr[i, ])), as.character(end(reads.gr[i, ])), sep = "-")
    non_tss_peak_ranges = append(non_tss_peak_ranges, range)
    }
print("TSS excluded")
if (length(reads.gr) == 1) {
    print("Only one interval left after excluding TSS, skipping")
    next
    }

saveRDS(reads.gr, paste("/g/data/zk16/arthuss/superenhancer/SecondRun/", as.character(args[1]), "_", as.character(individual1), "_", as.character(individual2), "_readCutoff_", as.character(readCutoff), "_", "IntervalUnionNoTSS.gr.rds", sep = ""))

# 5.Stitch intervals closer than 12500 bp together for both individuals (so there will be cases where a peak with a good number of reads in one of the replicates but not in the other is used for stitching in the replicate with the lower read number)
# The loop looks for equal intervals in the forward direction. Break for the last value means it has either been stitched with the penultimate interval, or is >StitchingWindow away, and is corectly left unstitched. 
StitchingWindow = 12500
merged.gr = reads.gr
IDs = as.numeric(c(1:length(reads.gr)))
scores = c()
for (range in IDs) {
    if (range == length(merged.gr)) {
        break
        }
    if ((start(merged.gr[range + 1, ]) - end(merged.gr[range, ])) <= StitchingWindow) {
        if (as.character(seqnames(reads.gr[range, ])) == as.character(seqnames(reads.gr[range + 1, ]))) {
            end(merged.gr[range, ]) = end(merged.gr[range + 1, ])
            } else {
            next
            }        
        }
    }
merged.gr = reduce(merged.gr)
print("Intervals stitched")

# 6.Sum read numbers for the resulting stitched clusters of peaks (same clusters for both replicates, but with different read numbers), plus record constituent peak numbers, indices, individual read numbers, standard deviation, coefficient of variation, width, etc
for (individual in c(args[2], args[3])) {
readcounts = subset(x = reads, subset = mouseID == individual)@assays$peaks@data
readcounts = readcounts[non_tss_peak_ranges, ]
if (typeof(readcounts) == 'double' & !(FALSE %in% (names(readcounts) %in% non_tss_peak_ranges))) {
ReadSums = as.integer(readcounts)
    } else if (typeof(readcounts) == 'double' & FALSE %in% (names(readcounts) %in% non_tss_peak_ranges)) {
    ReadSums = sum(readcounts)
    next
    }
    else {
    ReadSums = as.integer(rowSums(readcounts))
    }
mcols(reads.gr)$TotalMappingReads = ReadSums

overlaps = as.data.frame(findOverlaps(merged.gr, reads.gr))
PeakIndices = overlaps[1, 2]
ReadSums = c()
ReadSDs = c()
MappingReadNumbers = c()
ConstituentPeakCoords = c()
NumberOfConstituentPeaks = c()
ConstituentPeakIndices = c()
print("Individuals' readcounts entering loop")
    
for (ID in c(2:nrow(overlaps))) {
    if (overlaps[ID, 1] == overlaps[ID - 1, 1]) {
        PeakIndices = append(PeakIndices, overlaps[ID, 2])
        appended = TRUE
        } else {
        appended = FALSE
        MappingReads = c()
        PeakCoords = c()
        MappingReadsMerged = ""
        PeakCoordsMerged = ""
        PeakIndicesMerged = ""
        for (i in PeakIndices) {
            MappingReads = append(MappingReads, mcols(reads.gr)[i, 1])
            PeakCoords = append(PeakCoords, paste(as.character(seqnames(reads.gr[i, ])), as.character(start(reads.gr[i, ])), as.character(end(reads.gr[i, ])), sep = "-"))
            PeakIndicesMerged = paste(PeakIndicesMerged, as.character(i), sep = ",")
            }
        ReadSums = append(ReadSums, sum(as.numeric(MappingReads)))
        if (length(MappingReads) >= 3) {
        ReadSDs = append(ReadSDs, sd(as.numeric(MappingReads)))
            } else {
            ReadSDs = append(ReadSDs, NA)
            }
        for (number in MappingReads) {
            MappingReadsMerged = paste(MappingReadsMerged, as.character(number))
            }
        MappingReadNumbers = append(MappingReadNumbers, MappingReadsMerged)
        for (coord in PeakCoords) {
            PeakCoordsMerged = paste(PeakCoordsMerged, coord, sep = ",")
            }
        ConstituentPeakCoords = append(ConstituentPeakCoords, PeakCoordsMerged)
        NumberOfConstituentPeaks = append(NumberOfConstituentPeaks, length(PeakIndices))
        ConstituentPeakIndices = append(ConstituentPeakIndices, PeakIndicesMerged)
        PeakIndices = overlaps[ID, 2]
        }
    if (ID == nrow(overlaps)) { # Sorry I couldn't handle the edge case any other way...
        appended = FALSE
        MappingReads = c()
        PeakCoords = c()
        MappingReadsMerged = ""
        PeakCoordsMerged = ""
        PeakIndicesMerged = ""
        for (i in PeakIndices) {
            MappingReads = append(MappingReads, mcols(reads.gr)[i, 1])
            PeakCoords = append(PeakCoords, paste(as.character(seqnames(reads.gr[i, ])), as.character(start(reads.gr[i, ])), as.character(end(reads.gr[i, ])), sep = "-"))
            PeakIndicesMerged = paste(PeakIndicesMerged, as.character(i), sep = ",")
            }
        ReadSums = append(ReadSums, sum(as.numeric(MappingReads)))
        if (length(MappingReads) >= 3) {
        ReadSDs = append(ReadSDs, sd(as.numeric(MappingReads)))
            } else {
            ReadSDs = append(ReadSDs, NA)
            }
        for (number in MappingReads) {
            MappingReadsMerged = paste(MappingReadsMerged, as.character(number), sep = ",")
            }
        MappingReadNumbers = append(MappingReadNumbers, MappingReadsMerged)
        for (coord in PeakCoords) {
            PeakCoordsMerged = paste(PeakCoordsMerged, coord, sep = ",")
            }
        ConstituentPeakCoords = append(ConstituentPeakCoords, PeakCoordsMerged)
        NumberOfConstituentPeaks = append(NumberOfConstituentPeaks, length(PeakIndices))
        ConstituentPeakIndices = append(ConstituentPeakIndices, PeakIndicesMerged)
        PeakIndices = overlaps[ID, 2]
        }
    }
merged.gr$ReadSums = ReadSums
merged.gr$ReadSDs = ReadSDs
merged.gr$MappingReadNumbers = MappingReadNumbers
merged.gr$ConstituentPeakCoords = ConstituentPeakCoords
merged.gr$NumberOfConstituentPeaks = NumberOfConstituentPeaks
merged.gr$ConstituentPeakIndices = ConstituentPeakIndices

merged.gr$CoeffVar = (merged.gr$ReadSDs / (merged.gr$ReadSums / merged.gr$NumberOfConstituentPeaks))*100
merged.gr$StitchedLength = width(merged.gr)
merged.gr$ReadsPerbp = merged.gr$ReadSums / merged.gr$StitchedLength
print("Peaks stitched")

saveRDS(merged.gr, paste("/g/data/zk16/arthuss/superenhancer/SecondRun/", as.character(args[1]), "_", as.character(individual), "_stitchingWindow_", as.character(StitchingWindow), "_readCutoff_", as.character(readCutoff), "_merged.gr.rds", sep = ""))

# 7.Get the top-scoring stitched intervals using the ROSE slope function for each individual
score_cutoff = calculate_cutoff(merged.gr$ReadSums, drawPlot = FALSE)
mcols_temp = filter(as.data.frame(mcols(merged.gr)), ReadSums >= score_cutoff)
SE_candidates = makeGRangesFromDataFrame(filter(as.data.frame(merged.gr), ReadSums >= score_cutoff))
mcols(SE_candidates) = mcols_temp
saveRDS(SE_candidates, paste("/g/data/zk16/arthuss/superenhancer/SecondRun/SE_candidates_", as.character(args[1]), "_", as.character(individual), "_stitchingWindow_", as.character(StitchingWindow), "_readCutoff_", as.character(readCutoff), ".rds", sep = ""))
    }

if ((typeof(readcounts) == 'double' & FALSE %in% (names(readcounts) %in% non_tss_peak_ranges))) {
    print("One peak left after excluding TSS")
    fails = append(fails, args)
    next
    }

print("Individual SE candidates saved")

# 8.Find overlapping (i.e., identical) intervals between the two individuals and only keep them (along with metadata columns from both replicates); these are superenhancer candidates
SE_cands_rep1 = readRDS(paste("/g/data/zk16/arthuss/superenhancer/SecondRun/SE_candidates_", as.character(args[1]), "_", as.character(individual1), "_stitchingWindow_", as.character(StitchingWindow), "_readCutoff_", as.character(readCutoff), ".rds", sep = ""))
SE_cands_rep2 = readRDS(paste("/g/data/zk16/arthuss/superenhancer/SecondRun/SE_candidates_", as.character(args[1]), "_", as.character(individual2), "_stitchingWindow_", as.character(StitchingWindow), "_readCutoff_", as.character(readCutoff), ".rds", sep = ""))

SE_overlaps = as.data.frame(findOverlaps(SE_cands_rep1, SE_cands_rep2))
SE_combined = SE_cands_rep1[SE_overlaps[, 1], ]
colnames(mcols(SE_combined)) = paste(colnames(mcols(SE_combined)), individual1, sep = "_")
SE_cands_rep2 = SE_cands_rep2[SE_overlaps[, 2], ]
colnames(mcols(SE_cands_rep2)) = paste(colnames(mcols(SE_cands_rep2)), individual2, sep = "_")
mcols(SE_combined) = cbind(mcols(SE_combined), mcols(SE_cands_rep2))
SE_combined = SE_combined[!SE_combined %in% subsetByOverlaps(SE_combined, tss.gr)]
saveRDS(SE_combined, paste("/g/data/zk16/arthuss/superenhancer/SecondRun/CombinedSuperenhancerCandidates/", as.character(args[1]), "_", as.character(individual1), "_", as.character(individual2), "_stitchingWindow_", as.character(StitchingWindow), "_readCutoff_", as.character(readCutoff), "_CombinedSuperenhancerCandidates.gr.rds", sep = ""))
}

print(paste("Fails: ", as.character(fails), sep = ""))
