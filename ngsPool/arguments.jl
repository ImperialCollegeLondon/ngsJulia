
# Specific functions for ngsPool

using GZip
using ArgParse

# http://argparsejl.readthedocs.io/en/latest/argparse.html
function parse_commandline_pool()
	s = ArgParseSettings()

	@add_arg_table! s begin

		"--fin"
			help = "input file gzipped mpileup"
			arg_type = AbstractString
			required = true

		"--fout"
			help = "output file gzipped text"
			arg_type = AbstractString
			required = true #not printed in the --help in linux

		"--fsaf"
			help = "output gzipped saf file"
			arg_type = AbstractString
			default = "/dev/null"

		"--nChroms"
			help = "total number of chromosomes pooled (ploidy * number of individuals) [>0 ensables saf likelihoods]"
			arg_type = Int64
			default = 0

		"--lrtSnp"
			help = "LRT for SNP calling" #likelihood ratio test statistic
			arg_type = Float64
			default = -Inf #e.g. if 7.82 the LRT value in chisquare value in one degree of freedom p value is 0.01

		"--lrtBia"
			help = "LRT for biallelic calling"
			arg_type = Float64
			default = -Inf

		"--lrtTria"
			help = "LRT for triallelic (non) calling"
			arg_type = Float64
			default = Inf

		"--minQ"
			help = "minimum base quality in phredscore"
			arg_type = Int64
			default = 5

		"--minDepth"
			help = "minimum global depth"
			arg_type = Int64
			default = 1

		"--maxDepth"
			help = "maximum global depth"
			arg_type = Int64
			default = 100000

		"--nGrids"
			help = "grid density for grid-search estimation of allele frequencies"
			arg_type = Int64
			default = 0

		"--tol"
			help = "tolerance for GSS estimation of allele frequencies" #golden-section-search(GSS)
			arg_type = Float64
			default = 1e-5

		"--phredscale"
			help = "phredscale"
			arg_type = Int64
			default = 33

		"--verbose"
			help = "verbosity level"
			arg_type = Int64
			default = 0

		"--printSites"
			help =  "print on stdout every --printSites sites"
			arg_type = Int64
			default = 10000

	end

	return parse_args(s)

end

# -----------------------------------------
