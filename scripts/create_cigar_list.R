library(seqinr)
library(stringr)

if(!interactive()) {
    gaps_file = read.csv("./data/gaps_cigar.csv", header = TRUE, stringsAsFactors = FALSE)
    gaps = gaps_file$cigar
    cigar = gaps %>% str_extract_all(pattern="\\d+\\w") %>% lapply(., function(x) {
        lens = as.integer(x %>% str_extract_all(pattern="^\\d+"))
        chars = x %>% str_extract_all(pattern="\\w$")
        a = setNames(lens,chars)
        a
    })
    save(cigar, file = "data/cigar.rda")
}