# calculate distance metric between pairwise alignments

#library(stringr)
library(seqinr)

pa_distance_main = function(filename, aligners, metric) {
  refseq = read.fasta(paste0("data/ref_alignments/", filename),
                      set.attributes = FALSE,
                      forceDNAtolower = FALSE)
  aln = c()
  d = c()
  for(aln in aligners) {
    aln_file = paste0("aln/ref/", aln, "/", filename)
    alseq = read.fasta(aln_file,set.attributes = FALSE, forceDNAtolower = FALSE)

    h1 = vector("list",2)
    h2 = vector("list",2)
    h1[[1]] = homo_set(refseq,1,metric)
    h1[[2]] = homo_set(refseq,2,metric)
    h2[[1]] = homo_set(alseq,1,metric)
    h2[[2]] = homo_set(alseq,2,metric)

    # for all metrics except dSSP:
    stopifnot(all.equal(length(h1),length(h2)))
    
    c = length(c(which(refseq[[1]] != "-"),which(refseq[[2]] != "-")))
    distance = c()
    
    for(i in 1:2) {
      distance[i] = suppressWarnings(sym_distance(h1[[i]],h2[[i]]))
    }
    
    d[length(d)+1] = 1/c * sum(distance)
  }
  
  return(d)
}

homo_set = function(al, set, option) {
  # option: 0-dSSP, 1-dSEQ, 2-dPOS, 3-dEVOL

  if(option == 1) {
	  if(set == 1) {
		return(al[[2]][which(al[[1]] != "-")])
	  } else{
		return(al[[1]][which(al[[2]] != "-")])
	  }
  } else if(option == 2) {
    if(set == 1) {
      s = 2
    } else {
      s = 1
    }
    
    cont = 0
    for(i in 1:length(al[[s]])){
      if(al[[s]][i] == "-"){
        if(i==1) {
          al[[s]][i]=0
        } else if(!is.na(suppressWarnings(as.numeric(al[[s]][i-1])))) {
          al[[s]][i] = as.numeric(al[[s]][i-1])
        } else {
          al[[s]][i] = as.numeric(cont)
        }
      } else {
        cont = cont + 1
      }
    }
     
    h_set = al[[s]][which(al[[set]] != "-")]
    return(h_set)
    
  } else {
	  cat(option," is not a valid option.")
	  break
  }

} 

sym_distance = function(h1,h2) {
  return(sum(length(which(h1 != h2)))/2)
}
