library(seqinr)
library(stringr)

gaps_cigar = function(fasta) {
	if(file.exists(fasta)) {
		seqs = read.fasta(fasta, set.attributes = TRUE, forceDNAtolower = FALSE)
	} else {
		stop(paste("file:", fasta, "not found!"))
	}
	s1 = getSequence(seqs)[[1]]
	s2 = getSequence(seqs)[[2]]

	cigar = "" 	# cigar string
	current = ""
	nxt = ""
	count = 0

	stopifnot(length(s1) == length(s2))

	for(i in 1:length(s1)) {
		 if(s1[i] == '-'){
			if(current == "I") {
				count = count + 1
				next
			} else {
			  nxt = "I"
			  }
		} else if(s2[i] == '-') {
			if(current == "D") {
				count = count + 1
				next
			} else {
			  nxt = "D"
			  }
		} else if((s1[i] != '-') & (s2[i] != '-')) {
			if(current == "M") {
				count = count + 1
				next
			} else {
			  nxt = "M"
			  }
		}

	  	if(count > 0) {
			cigar = paste(cigar,paste(count,current, sep = ""), sep = "")
	  	}
		current = nxt
		count = 1
	  
	}
	
	cigar = paste(cigar,paste(count,current,sep=""),sep="")
	print(paste0(basename(fasta), ",", cigar, ",", basename(dirname(fasta))))
}

if(!interactive()) {
	ARGS = commandArgs(trailingOnly = TRUE)
	gaps_cigar(ARGS[1])
}
