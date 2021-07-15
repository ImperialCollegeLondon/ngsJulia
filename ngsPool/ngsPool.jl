# SNP calling, allele frequency estimation, sample allele frequency likelihoods from pooled sequencing data with gzipped mpileup files as input

# run 'julia ngsPool --help' for documentation
include("../templates.jl")
include("../functions.jl")
include("arguments.jl")

# alleles = ['A','C','G','T']

# parsing parameters
parsed_args = parse_commandline_pool()
if parsed_args["verbose"]>0
	println("Parsed args:")
	for (arg,val) in parsed_args
		println("  $arg  =>  $val")
	end
end

nSamp = parsed_args["nChroms"]
if nSamp < 0
	nSamp = 0
end

# open output file for writing
f2 = GZip.open(parsed_args["fout"], "w")#write
# write header
write(f2, join( ("##chrom", "position", "reference", "nonreference", "major",
	"minor", "lrtSNP", "lrtBia", "lrtTria", "maf", "freqMax", "freqE"), "\t"),
	"\n")

if parsed_args["fsaf"]!="/dev/null"
	f3 = GZip.open(parsed_args["fsaf"], "w")
end

# read mpileup

# Base.convert(::Type{AbstractString}, a::SubString) = absstring(a)#commented out as nowhere else called (Amend)

GZip.open(parsed_args["fin"]) do file
	# println(parsed_args["fin"], "  ", typeof(parsed_args["fin"]))
	allSites = infSites = 0

	for line in eachline(file)
		l = (split(line, "\t"))
		#1-6:
		 # chromosome, 1-based coordinate, reference base,
		 # the number of reads covering the site, read bases and base qualities

		# println("Char(l[3][1]) ", Char(l[3][1])) #extract one-element
		mySite = Site(l[1], parse(Int64, l[2]), Char(l[3][1]))
			# chrom::AbstractString
			# position::Int64
			# reference::Char
		myReads = Reads(l[5], chomp(l[6])) #chomp("Hello\n") "Hello"
			# base::AbstractString
	        # baseQuality::AbstractString

		(bases, indexDelN) = convertSyms(myReads, mySite)#translate of reads bases
		myReads = Reads(bases, myReads.baseQuality);
		# if there is a deleted base or a N, the base is marked as X
		 # and will be filtered out in the next stage

		error = 0
		# check if conversion was successful, if not print to stdout
		if length(myReads.base)!=length(myReads.baseQuality)
			println("Error! conversion:", myReads, mySite, "\n", length(myReads.base),"/",length(myReads.baseQuality))
			error = 1
		end

		# filter by base quality
		if error == 0
			newReads = filterReads(myReads; phredScale=parsed_args["phredscale"], minBaseQuality=parsed_args["minQ"])
			myReads = newReads
			newReads = 0
		end
		depth = length(myReads.base) #the num of times a given nucleotide in the genome has been read in an experiement

		# check if filtering was successful, if not print to stdout
		if length(myReads.base)!=length(myReads.baseQuality)
			if error == 0
				println("Error! filtering:", myReads, mySite, "\n", length(myReads.base),"/",length(myReads.baseQuality))
			end
			error = 1
		end

		# filter by depth
		if (error == 0 && depth >= parsed_args["minDepth"] && depth <= parsed_args["maxDepth"])

			# valid sites
			allSites += 1

			# calculate genotype likelihoods for all haploid cases
			haploid=calcGenoLogLike1(myReads, mySite)
			# We use a maximum likelihood approach to choose the major and minor alleles.
			(major, minor, minor2, minor3) = sortperm(haploid, rev=true)

			# is this biallelic?
			biaLike = calcGenoLogLike1_Bia(myReads, mySite, major, minor, phredScale=parsed_args["phredscale"])
			lrtBia = 2*(biaLike-haploid[major])

			# is this triallelic?
			triaLike = calcGenoLogLike1_Tria(myReads, mySite, major, minor, minor2, phredScale=parsed_args["phredscale"])
			lrtTria = 2*(triaLike-biaLike)

			#without SNP calling (all sites)
			if parsed_args["nGrids"] > 0 #grid-search
				freqsMLE = optimFreq_MajorMinor(myReads, major, minor, parsed_args["nGrids"]+1) #as maf=1/nGrid
			else #golden section search
				freqsMLE = optimFreq_MajorMinor_GSS(myReads, major, minor, parsed_args["tol"])
			end

			lrtSNP = snpPval_MajorMinor(myReads, freqsMLE[1], major, minor)
			#       snpPval_MajorMinor(read::Reads, maxlike::Float64, major::Int64, minor::Int64)

			#with SNP calling (only polymorphic sites)
			if (lrtSNP>parsed_args["lrtSnp"] && lrtBia>parsed_args["lrtBia"] && lrtTria<parsed_args["lrtTria"])
			 #checking whether they could be SNP, biallelic and triallelic

				# used sites
				infSites += 1

				# calculate saf (sample allele frequencies) likelihoods for
				# the non-ref allele
				freqMax = -1.0
				freqE = -1.0
				safLikes = zeros(nSamp+1)
				ref = nonref = 0
				nonref_allele ='N'

				ref = ((1:4)[mySite.reference .== alleles])[1] #read from mpileup file

				#just keep major minor allel as one of ref or nonref allele
				#and follow setting from mpile file (ref) instead of highest allele(haploid genotype) likelihood
				if ref == major
				 #major: allele with highest frequency defined in the genotype likelihood in view of haploid case
					nonref = minor
					nonref_allele = alleles[nonref]
				elseif ref == minor
					nonref = major
					nonref_allele = alleles[nonref]
				end

				# to get freqMax(maximum likelihood) and freqE(expected=p*frequency)
				 # SNP calling
				if nonref != 0 && nSamp>0
					# help = "number of samples [>0 ensables saf likelihoods]"

					for i = 0:nSamp
						safLikes[i+1] = calcFreqLogLike1_MajorMinor(myReads, ref, nonref, i/nSamp, parsed_args["phredscale"])
						# function calcFreqLogLike1_MajorMinor(read::Reads, major::Int64, minor::Int64, maf::Float64, phredScale::Int64=33)
					end

					safLikes = safLikes .- maximum(safLikes)
					# calculate maximum (and expectation)
					# freqMax = [((1:length(safLikes))[safLikes .== 0])] .- 1
					freqMax = findall(iszero, safLikes) .- 1
					freqMax = freqMax[1] #temperarily(2020/7/5)
					 #find the minor allele number (among samples) with maximum likelihood

					if freqMax > nSamp/2
						freqMax = nSamp - freqMax #ensure it's minor
					end

					# safLikes_norm = e.^safLikes
					safLikes_norm = pi.^safLikes
					safLikes_norm = safLikes_norm / sum(safLikes_norm)
					freqE = sum( (0:nSamp) .* safLikes_norm )
					if freqE > nSamp/2
						freqE = nSamp - freqE #ensure it's minor
					end

					safLikes_norm = 0
				end

				# write to .saf files
				if parsed_args["fsaf"]!="/dev/null"
					write(f3, join(safLikes,"\t"), "\n")
				end

				#notation:
				# write(f2, join( ("##chrom", "position", "reference", "nonreference", "major",
				# 	"minor", "lrtSNP", "lrtBia", "lrtTria", "maf", "naf_max", "naf_E"), "\t"),
				# 	"\n")
				# maf = freqsMLE[2] #minor allele frequency?

			#likelihood
				# naf_max = freqMax/nSamp #minor allele frequency with maximum likelihood
				# naf_E = freqE/nSamp #expected minor allele frequency
				 #naf output in individual instead of proportion (2020/5/4, Matteo suggested)

				#write to .out file
				write(f2, join( (mySite.chrom, mySite.position,
				 mySite.reference, nonref_allele, alleles[major],
				 alleles[minor], lrtSNP, lrtBia, lrtTria, freqsMLE[2],
				 freqMax[1]/nSamp, freqE/nSamp),
				 # freqMax, freqE),
				 "\t"), "\n")

				if parsed_args["verbose"]>0 && (allSites % parsed_args["printSites"])==0
					println("allSites:", allSites, "; ", "infSites:", infSites)
				end

			end
		end
	end

	GZip.close(f2)

	if parsed_args["fsaf"]!="/dev/null"
		GZip.close(f3)
	end

	if parsed_args["verbose"] > 0
		println("Done!\nallSites:", allSites, "; ", "infSites:", infSites)
	end

end #do file
