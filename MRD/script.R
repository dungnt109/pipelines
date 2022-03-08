### configurations

# channels
args = commandArgs(trailingOnly=TRUE)

separator <- "/"

#
library(rmarkdown)

silence <- FALSE

runmode <- "interactive"

if (length(args) > 0 && args[1] == "Silence") {
	silence <- TRUE
	runmode <- "silence"
}

cat(paste(args[1], "\n"))
cat(paste(silence, "\n"))

### to get less sensitive outlier detection, either:
# 1. increase niqr
# 2. decrease quantile.range
# suggested default: 
#niqr <- 3
#quantile.range <- 0.02

# the parameters below may not sensitive
niqr <- 5
quantile.range <- 0.005


show_outliers <- TRUE



marker.channel <- 1
alb.channel <- 2





mark.outlier <- function(intensities, range = 10) {

	#quantiles <- quantile(intensities)
	lower.quantile.value <- quantile(intensities, quantile.range)
	higher.quantile.value <- quantile(intensities, 1- quantile.range)
	iqr <- higher.quantile.value - lower.quantile.value
	upper.bound <- higher.quantile.value + range * iqr
	lower.bound <- higher.quantile.value - range * iqr
	mask <- (intensities < upper.bound) & (intensities > lower.bound)
	return(list(upper.bound = upper.bound, lower.bound = lower.bound, mask = mask))
}






## using k-means to divide positive control into 2 clusters and determine the threshold
## using c1 to represent for the positive cluster and c2 for the negative one
clustering_and_threshold <- function(intensities) {
	
	
	## initialize from max and min
	km <- kmeans(as.matrix(intensities), c(max(intensities), min(intensities)), nstart=1000)
	c1 <- intensities[km$clust==1]
	c2 <- intensities[km$clust==2]
	c1.median <- median(c1)
	c2.median <- median(c2)	
	if (c1.median < c2.median) {
		c1 <- intensities[km$clust==2]
		c2 <- intensities[km$clust==1]
		c1.median <- median(c1)
		c2.median <- median(c2)
	}
	
	solution <- (min(c1) + max(c2)) / 2
	

	## identify outlier in high cluster
	intensities.highcluster <- intensities[intensities > solution]
	upper.bound <- mark.outlier(intensities.highcluster, niqr)$upper.bound
	
	return(list(c1.median=c1.median, c2.median=c2.median, threshold=solution, upper.bound = upper.bound, mask = intensities < upper.bound))
}




##### load sample sheet
sample_sheet_file <- "/home/dungnt/Documents/Repository/pipeline/ddPCR/Exported_Run/KL1347b.txt"

folder <- dirname(sample_sheet_file)


cat(paste("Analyzing sample sheet", sample_sheet_file, "\n", sep=""))


cat("Please key in the MNC dilution factor or hit enter to use the default value (200x).\n")

mnc.alb.dilutionX <- trimws(readLines("stdin",n=1))
mnc.alb.dilutionX <- as.numeric(mnc.alb.dilutionX)

if (is.na(mnc.alb.dilutionX)) {
	mnc.alb.dilutionX <- 200
	cat("Using 200x.\n")
} else {
	cat(paste("Using user-specified value ", mnc.alb.dilutionX, "\n", sep=""))
}



cat("Please key in the albumin concentration (copies/uL, after dilution) of current MNC batch:\n");
mnc.alb.concentration <- trimws(readLines("stdin",n=1))
mnc.alb.concentration <- as.numeric(mnc.alb.concentration)




#folder <- choose.dir()

sample_sheet <- read.table(sample_sheet_file, header=TRUE, row.names=1)


row.col.combinations <- expand.grid(rownames(sample_sheet), colnames(sample_sheet))

samples <- rep("", nrow(row.col.combinations))

for (i in 1:nrow(row.col.combinations)) {
	samples[i] <- as.character(sample_sheet[row.col.combinations[i, 1], row.col.combinations[i, 2]])
}

row.col.combinations[,2] <- sprintf("%02d", as.numeric(gsub("X", "", row.col.combinations[,2])))

names(samples) <- apply(row.col.combinations, 1, paste, collapse="")



#### get the files
files <- list.files(folder, pattern = "*.csv")

well.position <- sapply(strsplit(files, "_", fixed=TRUE), "[[", 2)

names(files) <- well.position



#### get the control
## get ALB, MNC
mnc.alb <- samples[grepl("^MNC_ALB_", samples)]
mnc.alb.file <- paste(folder, files[names(mnc.alb)], sep=separator) 

mnc.alb.int <- read.csv(mnc.alb.file, header=TRUE)
mnc.alb.int <- mnc.alb.int[sample(nrow(mnc.alb.int)), ]

mnc.alb.clust <- clustering_and_threshold(mnc.alb.int[, alb.channel])

##if (!silence) {
##	plot(mnc.alb.int[,alb.channel],pch=16, cex=0.5, xlab="MNC, ALB", ylab="Amplitude")
##	lines(c(-1e6, 1e6), c(mnc.alb.clust$threshold, mnc.alb.clust$threshold), lty=2, col=2)
##}

# NTC, ALB
h2o.alb <- samples[grepl("^H2O_ALB_", samples)]
h2o.alb.file <-paste(folder, files[names(h2o.alb)], sep=separator) 
h2o.alb.int <- read.csv(h2o.alb.file, header=TRUE)
h2o.alb.int <- h2o.alb.int[sample(nrow(h2o.alb.int)), ]


### get Dx list
dx.marker.samples <- samples[grepl("^Dx_Mk_", samples)]
dx.marker.files <- paste(folder, files[names(dx.marker.samples)], sep=separator)
dx.sample.sid <- sapply(strsplit(dx.marker.samples, "_", fixed=TRUE), "[[", 3)
dx.sample.pid <- sapply(dx.sample.sid, function(x) {
	reg <- regexpr("[A-Z]+[0-9]+", x)
	substring(x, reg[1], reg[1] + attr(reg, "match.length")[1] - 1)
})
dx.sample.mid <- sapply(strsplit(dx.marker.samples, "_", fixed=TRUE), "[[", 4)
dx.sample.ng <- as.numeric(gsub("ng", "", sapply(strsplit(dx.marker.samples, "_", fixed=TRUE), "[[", 5)))


### get fu list
fu.marker.samples <- samples[grepl("^FU_Mk_", samples)]
fu.marker.files <- paste(folder, files[names(fu.marker.samples)], sep=separator)
fu.sample.sid <- sapply(strsplit(fu.marker.samples, "_", fixed=TRUE), "[[", 3)
fu.sample.pid <- sapply(fu.sample.sid, function(x) {
	reg <- regexpr("[A-Z]+[0-9]+", x)
	substring(x, reg[1], reg[1] + attr(reg, "match.length")[1] - 1)
})
fu.sample.mid <- sapply(strsplit(fu.marker.samples, "_", fixed=TRUE), "[[", 4)
fu.sample.ng <- as.numeric(gsub("ng", "", sapply(strsplit(fu.marker.samples, "_", fixed=TRUE), "[[", 5)))



### get MNC marker list
mnc.marker.samples <- samples[grepl("^MNC_Mk_", samples)]
mnc.marker.files <- paste(folder, files[names(mnc.marker.samples)], sep=separator)
mnc.sample.sid <- sapply(strsplit(mnc.marker.samples, "_", fixed=TRUE), "[[", 3)
mnc.sample.pid <- sapply(mnc.sample.sid, function(x) {
	reg <- regexpr("[A-Z]+[0-9]+", x)
	substring(x, reg[1], reg[1] + attr(reg, "match.length")[1] - 1)
})
mnc.sample.mid <- sapply(strsplit(mnc.marker.samples, "_", fixed=TRUE), "[[", 4)
mnc.sample.ng <- as.numeric(gsub("ng", "", sapply(strsplit(mnc.marker.samples, "_", fixed=TRUE), "[[", 5)))

### get H2O marker list
h2o.marker.samples <- samples[grepl("^H2O_Mk_", samples)]
h2o.marker.files <- paste(folder, files[names(h2o.marker.samples)], sep=separator)
h2o.sample.sid <- sapply(strsplit(h2o.marker.samples, "_", fixed=TRUE), "[[", 3)
h2o.sample.pid <- sapply(h2o.sample.sid, function(x) {
	reg <- regexpr("[A-Z]+[0-9]+", x)
	substring(x, reg[1], reg[1] + attr(reg, "match.length")[1] - 1)
})
h2o.sample.mid <- sapply(strsplit(h2o.marker.samples, "_", fixed=TRUE), "[[", 4)
h2o.sample.ng <- as.numeric(gsub("ng", "", sapply(strsplit(h2o.marker.samples, "_", fixed=TRUE), "[[", 5)))



### get ALB samples
alb.samples <- samples[grepl("ALB", samples)]
alb.files <- paste(folder, files[names(alb.samples)], sep=separator)








alb.samples.h2o <- alb.samples[grepl("^H2O_ALB", alb.samples)]
alb.files <- alb.files[!grepl("^H2O_ALB", alb.samples)]
alb.samples <- alb.samples[!grepl("^H2O_ALB", alb.samples)]
alb.samples.split <- strsplit(alb.samples, "_", fixed=TRUE)
alb.samples.info <- sapply(alb.samples.split, function(x) {x[1:3]})


alb.h2o.samples <- samples[grepl("^H2O_ALB", samples)]
alb.h2o.files <- paste(folder, files[names(alb.h2o.samples)], sep=separator)


cat("Analyzing albumin wells...\n")

alb.results.individual <- lapply(1:length(alb.samples), function(i) {
	ff <- alb.files[i]
	alb.int <- read.csv(ff, header=TRUE)
	alb.int <- alb.int[sample(nrow(alb.int)), ]

	
	outlier.in.silence.channel <- mark.outlier(alb.int[, marker.channel], niqr)
	alb.clust <- clustering_and_threshold(alb.int[outlier.in.silence.channel$mask, alb.channel])
	
	mask <- outlier.in.silence.channel$mask
	mask[mask] <- alb.clust$mask
	
	if (!silence) {
		windows()
		plot(alb.int[, marker.channel], alb.int[, alb.channel], col = (alb.int[, alb.channel] > alb.clust$threshold) + 1,
			pch=c(4, 16)[mask + 1],
		main=paste(alb.samples[i], names(alb.samples)[i], "\n Click to remove, and click the middle key when done." ), xlab = "Channel 1", ylab="Channel 2")
		
		lines(c(-1e8, 1e8), c(alb.clust$threshold, alb.clust$threshold), col=2, lwd=1.5)
		
		to.inverse <- identify(alb.int[,1], alb.int[,2], c("✓", "x")[mask + 1], col="blue")
		
		mask[to.inverse] <- !mask[to.inverse]
	}
	
	alb.clust2 <- clustering_and_threshold(alb.int[mask, alb.channel])
	#plot(alb.int[,marker.channel], alb.int[,alb.channel], col = (alb.int[,alb.channel] > alb.clust2$threshold) + 1,
	#		pch = c(4, 16)[mask + 1],
	#		main=paste(alb.samples[i], names(alb.samples)[i], "\n Final" ), xlab = "Channel 1", ylab="Channel 2")
	#
	#lines(c(-1e8, 1e8), c(alb.clust2$threshold, alb.clust2$threshold), col=2, lwd=1.5)
	#Sys.sleep(5)
	n.positive.droplets <- sum(alb.int[,alb.channel][mask] > alb.clust2$threshold)
	n.droplets <- sum(mask)
	n.outliers <- sum(!mask)
	concentration <- -log(1-(n.positive.droplets/n.droplets))/0.00085
	list(intensities=alb.int, threshold=alb.clust2$threshold, mask=mask, n.positive.droplets=n.positive.droplets, n.droplets = n.droplets, concentration=concentration, name=alb.samples[i], upper = alb.clust2$upper, n.outliers = n.outliers)
})

alb.h2o.results.individual <- lapply(1:length(alb.h2o.samples), function(i) {
	ff <- alb.h2o.files[i]
	alb.h2o.int <- read.csv(ff, header=TRUE)
	alb.h2o.int <- alb.h2o.int[sample(nrow(alb.h2o.int)), ]

	
	alb.mnc <- alb.results.individual[grepl("^MNC_ALB", alb.samples)]
	alb.mnc[[1]]$threshold
	
	alb.min <- min(c(alb.mnc[[1]]$intensities[, alb.channel], alb.h2o.int[, alb.channel]))
	alb.max <- max(c(alb.mnc[[1]]$intensities[, alb.channel], alb.h2o.int[, alb.channel]))
	
	
	outlier.in.silence.channel <- mark.outlier(alb.h2o.int[, marker.channel], niqr)
	mask <- outlier.in.silence.channel$mask & alb.h2o.int[, alb.channel] < alb.mnc[[1]]$upper
	
	if (!silence) {
		plot(alb.h2o.int[, marker.channel], alb.h2o.int[, alb.channel], col = (alb.h2o.int[, alb.channel] > alb.mnc[[1]]$threshold) + 1,
			pch=c(4, 16)[mask + 1], ylim = c(alb.min, alb.max),
		main=paste(alb.h2o.samples[i], names(alb.h2o.samples)[i], "\n Click to remove, and click the middle key when done." ), xlab = "Channel 1", ylab="Channel 2")
		
		lines(c(-1e8, 1e8), c(alb.mnc[[1]]$threshold, alb.mnc[[1]]$threshold), col=2, lwd=1.5)
		
		to.inverse <- identify(alb.h2o.int[,1], alb.h2o.int[,2], c("✓", "x")[mask + 1], col="blue")
		
		mask[to.inverse] <- !mask[to.inverse]
	}
	
	n.positive.droplets <- sum(alb.h2o.int[,alb.channel][mask] > alb.mnc[[1]]$threshold)
	n.droplets <- sum(mask)
	n.outliers <- sum(!mask)
	concentration <- -log(1-(n.positive.droplets/n.droplets))/0.00085
	list(intensities=alb.h2o.int, threshold=alb.mnc[[1]]$threshold, mask=mask, n.positive.droplets=n.positive.droplets, n.droplets = n.droplets, concentration=concentration, name=alb.h2o.samples[i], n.outliers = n.outliers)
})


cat("Analyzing marker wells...\n")

dx.sample.clust <- lapply(1:length(dx.marker.samples), function(i) {
	ff <- dx.marker.files[i]
	
	pid <- dx.sample.pid[i]
	dx.sid <- dx.sample.sid[i]
	
	mid <- dx.sample.mid[i]
	
	dx.marker.int <- read.csv(ff, header=TRUE)
	dx.marker.int <- dx.marker.int[sample(nrow(dx.marker.int)), ]
	
	#if (!silence) {
	#	winDialog(type = c("ok"), paste("Start processing", pid))
	#}
	### identify outlier in the other channel
	outlier.in.silence.channel <- mark.outlier(dx.marker.int[, alb.channel], niqr)
	
	dx.marker.clust <- clustering_and_threshold(dx.marker.int[outlier.in.silence.channel$mask, marker.channel])
	
	
	
	dx.mask <- outlier.in.silence.channel$mask
	dx.mask[dx.mask] <- dx.marker.clust$mask
	
	if (!silence) {
		plot(dx.marker.int[,marker.channel], dx.marker.int[,alb.channel], col = (dx.marker.int[,marker.channel] > dx.marker.clust$threshold) + 1,
				pch = c(4, 16)[dx.mask + 1],
				main=paste(dx.marker.samples[i], names(dx.marker.samples)[i], "\n Click to remove, and click the middle key when done." ), xlab = "Channel 1", ylab="Channel 2")
		
		lines(c(dx.marker.clust$threshold, dx.marker.clust$threshold), c(-1e8, 1e8), col=2, lwd=1.5)
		
		to.inverse <- identify(dx.marker.int[,1], dx.marker.int[,2], c("✓", "x")[dx.mask + 1], col="blue")
		
		
		
		
		dx.mask[to.inverse] <- !dx.mask[to.inverse]
	}
	
	dx.marker.clust2 <- clustering_and_threshold(dx.marker.int[dx.mask, marker.channel])
	
	
	#plot(dx.marker.int[,marker.channel], dx.marker.int[,alb.channel], col = (dx.marker.int[,marker.channel] > dx.marker.clust2$threshold) + 1,
	#		pch = c(4, 16)[dx.mask + 1],
	#		main=paste(dx.marker.samples[i], names(dx.marker.samples)[i], "\n Final" ), xlab = "Channel 1", ylab="Channel 2")
	#
	#lines(c(dx.marker.clust2$threshold, dx.marker.clust2$threshold), c(-1e8, 1e8), col=2, lwd=1.5)
	#
	
	n.positive.droplets <- sum(dx.marker.int[,marker.channel][dx.mask] > dx.marker.clust2$threshold)
	n.droplets <- sum(dx.mask)
	n.outliers <- sum(!dx.mask)
	concentration <- -log(1-(n.positive.droplets/n.droplets))/0.00085
	
	dx.results <- list(list(intensities=dx.marker.int, 
						threshold=dx.marker.clust2$threshold, 
						mask=dx.mask, 
						n.positive.droplets=n.positive.droplets, 
						n.droplets = n.droplets, 
						concentration=concentration, 
						name=dx.marker.samples[i], 
						n.outliers = n.outliers))
	## MNC marker files
	mnc.samples <- mnc.marker.samples[mnc.sample.pid == pid & mnc.sample.mid == mid]
	mnc.files <- mnc.marker.files[mnc.sample.pid == pid & mnc.sample.mid == mid]
	
	mnc.results <- lapply(1:length(mnc.samples), function(j) {
		mnc.file <- paste(folder, files[names(mnc.samples)[j]], sep=separator)
		mnc.marker.int <- read.csv(mnc.file, header=TRUE)
		mnc.marker.int <- mnc.marker.int[sample(nrow(mnc.marker.int)), ]


		mnc.mask <- (mnc.marker.int[, alb.channel] < outlier.in.silence.channel$upper.bound) & (mnc.marker.int[, alb.channel] > outlier.in.silence.channel$lower.bound) & (mnc.marker.int[, marker.channel] < dx.marker.clust2$upper.bound)
		
		
		if (!silence) {
		
			## get plot range
			x.min <- min(c(mnc.marker.int[, marker.channel], dx.marker.int[, marker.channel]))
			x.max <- max(c(mnc.marker.int[, marker.channel], dx.marker.int[, marker.channel]))
			y.min <- min(c(mnc.marker.int[, alb.channel], dx.marker.int[, alb.channel]))
			y.max <- max(c(mnc.marker.int[, alb.channel], dx.marker.int[, alb.channel]))
			
			
			plot(mnc.marker.int[,marker.channel], mnc.marker.int[,alb.channel], col = (mnc.marker.int[,marker.channel] > dx.marker.clust2$threshold) + 1,
					pch = c(4, 16)[mnc.mask + 1],
					main=paste(mnc.samples[j], names(mnc.samples)[j], "\n Click to remove, and click the middle key when done." ), xlab = "Channel 1", ylab="Channel 2", xlim=c(x.min, x.max), ylim=c(y.min, y.max))
			lines(c(dx.marker.clust2$threshold, dx.marker.clust2$threshold), c(-1e8, 1e8), col=2, lwd=1.5)
			mnc.to.inverse <- identify(mnc.marker.int[,1], mnc.marker.int[,2], c("✓", "x")[mnc.mask + 1], col="blue")
			mnc.mask[mnc.to.inverse] <- !mnc.mask[mnc.to.inverse]
		}
		n.positive.droplets <- sum(mnc.marker.int[,marker.channel][mnc.mask] > dx.marker.clust2$threshold)
		n.droplets <- sum(mnc.mask)
		n.outliers <- sum(!mnc.mask)
		concentration <- -log(1-(n.positive.droplets/n.droplets))/0.00085
	
		list(intensities=mnc.marker.int, 
						threshold=dx.marker.clust2$threshold, 
						mask=mnc.mask, 
						n.positive.droplets=n.positive.droplets, 
						n.droplets = n.droplets, 
						concentration=concentration, 
						name=mnc.samples[j],   #### updated from name=mnc.marker.samples[j],
						n.outliers = n.outliers)
	})
	
	
	h2o.samples <- h2o.marker.samples[h2o.sample.pid == pid & h2o.sample.mid == mid]
	h2o.files <- h2o.marker.files[h2o.sample.pid == pid & h2o.sample.mid == mid]
	
	h2o.results <- lapply(1:length(h2o.samples), function(j) {
		h2o.file <- paste(folder, files[names(h2o.samples)[j]], sep=separator)
		h2o.marker.int <- read.csv(h2o.file, header=TRUE)
		h2o.marker.int <- h2o.marker.int[sample(nrow(h2o.marker.int)), ]

		
		
		h2o.mask <- (h2o.marker.int[, alb.channel] < outlier.in.silence.channel$upper.bound) & (h2o.marker.int[, alb.channel] > outlier.in.silence.channel$lower.bound) & (h2o.marker.int[, marker.channel] < dx.marker.clust2$upper.bound)
		
		if (!silence) {
			## get plot range
			x.min <- min(c(h2o.marker.int[, marker.channel], dx.marker.int[, marker.channel]))
			x.max <- max(c(h2o.marker.int[, marker.channel], dx.marker.int[, marker.channel]))
			y.min <- min(c(h2o.marker.int[, alb.channel], dx.marker.int[, alb.channel]))
			y.max <- max(c(h2o.marker.int[, alb.channel], dx.marker.int[, alb.channel]))
			
			
			plot(h2o.marker.int[,marker.channel], h2o.marker.int[,alb.channel], col = (h2o.marker.int[,marker.channel] > dx.marker.clust2$threshold) + 1,
					pch = c(4, 16)[h2o.mask + 1],
					main=paste(h2o.samples[j], names(h2o.samples)[j], "\n Click to remove, and click the middle key when done." ), xlab = "Channel 1", ylab="Channel 2",
					xlim=c(x.min, x.max), ylim=c(y.min, y.max))
			lines(c(dx.marker.clust2$threshold, dx.marker.clust2$threshold), c(-1e8, 1e8), col=2, lwd=1.5)
			h2o.to.inverse <- identify(h2o.marker.int[,1], h2o.marker.int[,2], c("✓", "x")[h2o.mask + 1], col="blue")
			h2o.mask[h2o.to.inverse] <- !h2o.mask[h2o.to.inverse]
			#pos.droplets <- sum(h2o.marker.int[h2o.mask, marker.channel] > dx.marker.clust2$threshold)
			#n.droplets <- sum(h2o.mask)
			#list(int = h2o.marker.int, h2o.mask = h2o.mask, pos.droplets = pos.droplets)
		}
		
		n.positive.droplets <- sum(h2o.marker.int[,marker.channel][h2o.mask] > dx.marker.clust2$threshold)
		n.droplets <- sum(h2o.mask)
		n.outliers <- sum(!h2o.mask)
		concentration <- -log(1-(n.positive.droplets/n.droplets))/0.00085
	
		list(intensities=h2o.marker.int, 
						threshold=dx.marker.clust2$threshold, 
						mask=h2o.mask, 
						n.positive.droplets=n.positive.droplets, 
						n.droplets = n.droplets, 
						concentration=concentration, 
						name=h2o.samples[j],   ########### updated from name=h2o.marker.samples[j],
						n.outliers = n.outliers)
						
		
	})
	
	
	
	### search for follow up channels
	
	
	
	fu.sids <- unique(fu.sample.sid[fu.sample.pid == pid & fu.sample.mid == mid])
	
	
	for (fu.sid in fu.sids) {
	
		#dx.baseline <- readline(prompt = paste("Processing ", fu.sid, "_", mid, ". Please key in the Dx Baseline value or hit enter to use in-plate Dx baseline.\n", sep=""))
		#dx.baseline <- as.numeric(dx.baseline)
		
		cat(paste("Processing ", fu.sid, "_", mid, ". Please key in the tumour load at Dx or hit enter to use the value calculated from the current test.\n", sep=""));
		dx.baseline <- trimws(readLines("stdin",n=1))
		dx.baseline <- as.numeric(dx.baseline)
		if(is.na(dx.baseline)) {
			cat("Using in-plate Dx baseline.\n")
		} else {
			cat("Using user specified Dx baseline.\n")
		
		}
		
		fu.samples <- fu.marker.samples[fu.sample.pid == pid & fu.sample.mid == mid & fu.sample.sid == fu.sid]
		fu.files <- fu.marker.files[fu.sample.pid == pid & fu.sample.mid == mid & fu.sample.sid == fu.sid]
		fu.results <- lapply (1:length(fu.samples), function(j) {
			fu.file <- paste(folder, files[names(fu.samples)[j]], sep=separator)
			fu.marker.int <- read.csv(fu.file, header=TRUE)
			fu.marker.int <- fu.marker.int[sample(nrow(fu.marker.int)), ]

			fu.mask <- (fu.marker.int[, alb.channel] < outlier.in.silence.channel$upper.bound) & (fu.marker.int[, alb.channel] > outlier.in.silence.channel$lower.bound) & (fu.marker.int[, marker.channel] < dx.marker.clust2$upper.bound)
			
			if (!silence) {
				## get plot range
				x.min <- min(c(fu.marker.int[, marker.channel], dx.marker.int[, marker.channel]))
				x.max <- max(c(fu.marker.int[, marker.channel], dx.marker.int[, marker.channel]))
				y.min <- min(c(fu.marker.int[, alb.channel], dx.marker.int[, alb.channel]))
				y.max <- max(c(fu.marker.int[, alb.channel], dx.marker.int[, alb.channel]))
			
			
				
				plot(fu.marker.int[,marker.channel], fu.marker.int[,alb.channel], col = (fu.marker.int[,marker.channel] > dx.marker.clust2$threshold) + 1,
						pch = c(4, 16)[fu.mask + 1],
						main=paste(fu.samples[j], names(fu.samples)[j], "\n Click to remove, and click the middle key when done." ), xlab = "Channel 1", ylab="Channel 2")
				lines(c(dx.marker.clust2$threshold, dx.marker.clust2$threshold), c(-1e8, 1e8), col=2, lwd=1.5)
				fu.to.inverse <- identify(fu.marker.int[,1], fu.marker.int[,2], c("✓", "x")[fu.mask + 1], col="blue")
				fu.mask[fu.to.inverse] <- !fu.mask[fu.to.inverse]
			}
			
			pos.droplets <- sum(fu.marker.int[fu.mask, marker.channel] > dx.marker.clust2$threshold)
			
			list(int = fu.marker.int, fu.mask = fu.mask, pos.droplets = pos.droplets)
			
			
			n.positive.droplets <- sum(fu.marker.int[,marker.channel][fu.mask] > dx.marker.clust2$threshold)
			n.droplets <- sum(fu.mask)
			n.outliers <- sum(!fu.mask)
			concentration <- -log(1-(n.positive.droplets/n.droplets))/0.00085
		
			list(intensities=fu.marker.int, 
							threshold=dx.marker.clust2$threshold, 
							mask=fu.mask, 
							n.positive.droplets=n.positive.droplets, 
							n.droplets = n.droplets, 
							concentration=concentration, 
							name=fu.samples[j],           ##### updated from name=fu.marker.samples[j], 
							n.outliers = n.outliers)
				
		})
		
		
		dx.alb <- alb.results.individual[alb.samples.info[1,] == "Dx" &  alb.samples.info[3,] == dx.sid]
		fu.alb <- alb.results.individual[alb.samples.info[1,] == "FU" &  alb.samples.info[3,] == fu.sid]
		mnc.alb <- alb.results.individual[alb.samples.info[1,] == "MNC"]
		h2o.alb <- alb.h2o.results.individual
		
		dx.marker <- dx.results
		fu.marker <- fu.results
		mnc.marker <- mnc.results
		h2o.marker <- h2o.results
		
		## draw figure 
		render("/home/dungnt/Documents/Repository/pipeline/ddPCR/ddPCR_rTemplate_report.Rmd", params = list(sid = fu.sid, mid = mid, dx.alb = dx.alb, fu.alb = fu.alb, mnc.alb = mnc.alb, h2o.alb = h2o.alb, dx.marker = dx.marker, fu.marker=fu.marker, mnc.marker=mnc.marker, h2o.marker=h2o.marker, show_outliers=show_outliers, dx.baseline=dx.baseline, mnc.alb.dilutionX=mnc.alb.dilutionX, mnc.alb.concentration=mnc.alb.concentration, date=Sys.time()), 
				output_file = paste(folder, separator, fu.sid, "_", mid, "_", runmode, "_report.pdf", sep="")
		)
		## draw figure 
		render("/home/dungnt/Documents/Repository/pipeline/ddPCR/ddPCR_rTemplate_2dFigs.Rmd", params = list(sid = fu.sid, mid = mid, dx.alb = dx.alb, fu.alb = fu.alb, mnc.alb = mnc.alb, h2o.alb = h2o.alb, dx.marker = dx.marker, fu.marker=fu.marker, mnc.marker=mnc.marker, h2o.marker=h2o.marker, show_outliers=show_outliers, dx.baseline=dx.baseline, mnc.alb.dilutionX=mnc.alb.dilutionX, mnc.alb.concentration=mnc.alb.concentration, date=Sys.time()), 
				output_file = paste(folder, separator, fu.sid, "_", mid, "_", runmode, "_2dFigs.pdf", sep="")
		)
	
	}
	
	#msk <- rep(FALSE, nrow(dx.marker.int))
	
	#msk[to.exclude] <- TRUE
	#dx.marker.clust <- clustering_and_threshold(dx.marker.int[!msk, marker.channel])
	#list(dx.marker.int, dx.marker.clust, to.exclude)
	
})







