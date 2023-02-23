
num_alns = function(species) {
    dseq = read.csv(paste0("data/",species,"/dseq.csv"), header = TRUE, stringsAsFactors = FALSE)
#dseq = read.csv("../data/mouse/dseq.csv", header = TRUE, stringsAsFactors = FALSE)
#dpos = read.csv("../data/dpos.csv", header = TRUE)

    col_s = ncol(dseq)
    row_s = nrow(dseq)
#col_p = ncol(dpos)

# perfect alignments
    print("Perfect alignments")
    if(col_s == 2) {
        print(length(which(dseq[,2] == 0)))
    } else{
        print(apply(dseq[,2:col_s],2,function(x){length(which(x == 0))}))
        #apply(dpos[,2:7],2,function(x){length(which(x == 0))})
    }

# best alignments (lower d_seq)
    print("Best alignments")
    best = rep(0,col_s-1)
    for(r in 1:(row_s)) {
      for(i in which(min(dseq[r,2:col_s])==dseq[r,2:col_s])){#1:length(i)){
        best[i] = best[i] + 1
      }
    }
    print(best)

    print("Imperfect alignments")
    imp = rep(0,col_s-1)
    for(r in 1:(row_s)) {
      if(min(dseq[r,2:col_s]) == 0){
        for(i in which(dseq[r,2:col_s]>0)) {
          imp[i] = imp[i] + 1
        }
      }
    }
print(imp)
}

# check if best alignments by toy and marg are the same
#best_toy = which(dseq[,3] == 0)
#best_marg = which(dseq[,4] == 0)
#print(length(which(dseq[,3] == 0 & dseq[,4] == 0 & dseq[,2] == 0)))

if(!interactive()) {
    ARGS = commandArgs(trailing = TRUE)
    num_alns(ARGS[1])
}


