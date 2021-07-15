
args <- commandArgs(T)

saf_cases <- read.table(args[1])
saf_controls <- read.table(args[2])

p0 <- apply(FUN=which.max, X=saf_controls, MAR=1) # in practice should be -1, these are indexs

for (j in 1:length(p0)) {

	H0 <- saf_cases[j,p0[j]]
        HA <- max(saf_cases[j])
        LRT <- -2*(H0-HA)

        cat(j, "\t", LRT,"\n")

}




