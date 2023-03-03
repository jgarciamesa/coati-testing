library(seqinr)
library(stringr)
library(TeachingDemos)


simulate_two_main = function(input, output) {
    seqs = read.fasta(input, set.attributes=FALSE, forceDNAtolower=FALSE)

	char2seed(basename(gsub(".fasta", "", output)))
	cigar_id = sample.int(n = length(cigar), size = 1)
	cigars = cigar[cigar_id]
	current_cigar = 1

    repeat {
        GG = G
        sims = seqs
        A = sims[names(sims)[[1]]][[1]]
        B = sims[names(sims)[[2]]][[1]]#$B
        stopifnot(length(A) == length(B), length(A) %% 3 == 0)

        # sample from cigar strings
        # g = sample(cigar, 1)[[1]]
		g = cigars[current_cigar][[1]]
		current_cigar = current_cigar + 1
        # Identify any flexible width matches
        bDM = names(g) %in% c("M","D")
        repeat {
            bF = names(g) == "M" & g >= GG
            if(sum(bF) != 0) {
                break
            }
            GG = (GG %/% 6)*3
        }
        # save phase
        p = g %% 3
        # total committed length
        w = sum(c(p[bF]+GG,g[bDM & !bF]))
        if(w > length(A)) {
            # this cigar is too big for our sequence, so try again
            next
        }
        # number of glue sections
        n = sum(bF)
        # sample glue widths from the remaining length
		stopifnot((length(A)-w) %% 3 == 0) #if((length(A)-w) %% 3 != 0) {next}
		r = (length(A)-w) %/% 3
        s = round(r*sort(runif(n-1)))
        x = c(s,r)-c(0,s)
        # update widths
        g[bF] = GG+p[bF]+3*x
        stopifnot(g %% 3 == p)
		stopifnot(sum(g[bDM]) == length(A)) #if(sum(g[bDM]) != length(A)) {next}

        # update alignments
        pos = cumsum(g)-g+1
        for(i in seq_along(g)) {
            if(names(g)[i] == "D") {
                o = pos[i]+seq.int(0,length.out=g[i])
                B[o] = "-"
            } else if(names(g)[i] == "I") {
                A = append(A,rep("-", g[i]), pos[i]-1)
                B = append(B,sample(nucs,g[i],prob=p_nucs,replace=TRUE),pos[i]-1)
            }
        }
        am = A != "-"
        bm = B != "-"
        h = ifelse(am, ifelse(bm, "M", "D"), ifelse(bm, "I", "X"))
        h = rle(h)
		stopifnot(h$values == names(g), h$lengths == g)# if(h$values != names(g) || h$lengths != g) {next}
        sims[names(sims)[[1]]][[1]] = A
        sims[names(sims)[[2]]][[1]] = B
        break;
    }
    write.fasta(sequences = sims, file.out = output)
    print(paste0(
        gaps_file$raw_name[cigar_id], ",",
        gaps_file$cigar[cigar_id], ",",
        gaps_file$origin[cigar_id], ",",
        basename(output)
    ))
}

if(!interactive()) {
    ARGS = commandArgs(trailing=TRUE)

	# fix issues with seqinr
	source("./scripts/write_fasta.R")

	# maximum length of a M before glue starts
	G = 99

	gaps_file = read.csv("./data/gaps_cigar.csv", header = TRUE, stringsAsFactors = FALSE)
	gaps = gaps_file$cigar
	cigar = gaps %>% str_extract_all(pattern="\\d+\\w") %>% lapply(., function(x) {
    	lens = as.integer(x %>% str_extract_all(pattern="^\\d+"))
    	chars = x %>% str_extract_all(pattern="\\w$")
    	a = setNames(lens,chars)
    	a
	})

	p_nucs = c("A"=0.308,"C"=0.185,"G"=0.199,"T"=0.308)
	nucs = names(p_nucs)

    simulate_two_main(input=ARGS[1], output=ARGS[2])
}
