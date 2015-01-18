from __future__ import division
import argparse

def macawInterpret(args):
	# Interprets Macaw file.  NOTE: NUMBER OF MARKERS MUST BE UPDATED WHEN MARKER LIST CHANGES!
	from scipy import stats

	# set up all dictionariess, etc.
	i = 0
	linCountsPresent = {} 			# number present in the indicated group
	linCountsPresentInverse = {} 	# number present in all except the indicated group
	linCountsAbsent = {}			# number absent in the indicated group
	linCountsAbsentInverse = {} 	# number absent in all except the indicated group
	totalMarkersPerLin = {}			# total number of markers in each lineage
	markerCounts = {}				# read depths for markers in each lineage
	avgMarkerPerLin = {}			# average read depth for markers in each lineage
	linPercentTotal = {}
	linPercentEnd = {}
	fisherPValues = {} 				# p values for each lineage
	predLin = []					# predicted lineages

	# function for fisher's exact test
	def fisher(tbl1, tbl2, tbl3, tbl4):
		oddsratio, pvalue = stats.fisher_exact([[tbl1, tbl2], [tbl3, tbl4]], alternative='greater')
		return pvalue

	# open files
	f = open(args.macawFile, "r")
	savf = open(args.outputFile + ".interpreted.macaw", "w")

	# go through each line of macaw file and count absent/present/read depths
	for line in f:
		line2 = line.rstrip()
		if not line2.startswith("#"):
			i = i + 1
			if not i > 3722:
				li = line2.split("\t")
				markerCount = int(li[1])
				mutInfo = li[0].split("_")
				mutInfoLIN = mutInfo[2]
				# add marker counts to total
				try:
					totalMarkersPerLin[mutInfoLIN] += 1
				except KeyError:
					totalMarkersPerLin[mutInfoLIN] = 1
				# add marker read depth counts to total
				try:
					markerCounts[mutInfoLIN] += markerCount
				except KeyError:
					markerCounts[mutInfoLIN] = markerCount
				# add count of present/absent markers to total
				if li[3] == "present":
					try:
						linCountsPresent[mutInfoLIN] += 1
					except KeyError:
						linCountsPresent[mutInfoLIN] = 1
				elif li[3] == "absent":
					try:
						linCountsAbsent[mutInfoLIN] += 1
					except KeyError:
						linCountsAbsent[mutInfoLIN] = 1
	
	# sum up everything for analysis
	sumAllPresent = sum(linCountsPresent.values())
	sumAllAbsent = sum(linCountsAbsent.values())
	
	# make sure that each key is in each other's dictionary (for fisher's)
	for key in linCountsPresent:
		if not key in linCountsAbsent:
			linCountsAbsent[key] = 0
	for key in linCountsAbsent:
		if not key in linCountsPresent:
			linCountsPresent[key] = 0

	# calculate p values using fisher's exact test for each lineage
	for key in linCountsPresent:
		linCountsPresentInverse[key] = sumAllPresent - linCountsPresent[key]
		linCountsAbsentInverse[key] = sumAllAbsent - linCountsAbsent[key]
		roughPValue = fisher(linCountsPresent[key], linCountsPresentInverse[key], linCountsAbsent[key], linCountsAbsentInverse[key]) * len(linCountsPresent.keys())
		if roughPValue > 1:
			thePValue = 1
		else:
			thePValue = roughPValue
		fisherPValues[key] = thePValue

	# average marker read depths for each lineage
	for key in markerCounts:
		avgMarkerPerLin[key] = markerCounts[key] / totalMarkersPerLin[key]
	
	# percent of total depth
	allLinDepth = 0
	for lineage in avgMarkerPerLin:
		allLinDepth += avgMarkerPerLin[lineage]
	for lineage in avgMarkerPerLin:
		linPercentTotal[lineage] = avgMarkerPerLin[lineage] / allLinDepth
			
	for key in fisherPValues:
		if fisherPValues[key] < 0.05:
			predLin.append(key)		
	
	savf.write("## RESULTS ##" + "\n")
	
	if len(predLin) > 1:
		totalLinDepth = 0
		for lineage in predLin:
			totalLinDepth += avgMarkerPerLin[lineage]
		savf.write("Predicted group(s):\t")
		for lineage in predLin:
			linPercentEnd[lineage] = avgMarkerPerLin[lineage] / totalLinDepth
		thisI = 0
		for key in linPercentEnd:
			thisI += 1
			savf.write(str(key) + " (" + str(linPercentEnd[key]) + ")")
			if thisI < len(linPercentEnd):
				savf.write(", ")
		savf.write("\nMixed sample?:\tTRUE")
	elif len(predLin) == 0:
		savf.write("Predicted group(s):\tNONE")
		savf.write("\nMixed sample?:\tFALSE")
	else:
		savf.write("Predicted group(s):\t" + ', '.join(predLin))
		savf.write("\nMixed sample?:\tFALSE")
		
	savf.write("\n\n## DETAILS ##" + "\n")
	
	for key in fisherPValues:
		savf.write(str(key) + " Present markers:\t" + str(linCountsPresent[key]) + "\n")
		savf.write(str(key) + " Absent markers:\t" + str(linCountsAbsent[key]) + "\n")
		savf.write(str(key) + " Raw marker total:\t" + str(markerCounts[key]) + " (" + str(linPercentTotal[key]) + ")\n")
		savf.write(str(key) + " Corrected P Value:\t" + str(fisherPValues[key]) + "\n\n")

def macawCondense(args):
	import glob
	import os

	outf = open(args.outputFile, "w")
	mixedSamples = []
	
	for filename in glob.glob(os.path.join(args.macawDir, '*.interpreted.macaw')):
		f = open(filename, "r")
		nextli = 0
		for line in f:
			li = line.rstrip()
			fspli = filename.split("/")
			endf = fspli[-1]
			if nextli == 1:
				li2 = li.split("\t")
				outf.write(endf[:-18] + "\t" + li2[-1] + "\n")
				nextli = 2
			elif nextli == 2:
				li2 = li.split("\t")
				if li2[-1] == "True":
					mixedSamples.append(endf[:-18])
				nextli = 0
			elif li == "## RESULTS ##":
				nextli = 1
			else:
				pass
	
	outf.write("Mixed samples: " + ', '.join(mixedSamples))

# parser prep
parser = argparse.ArgumentParser(description="Suite of utilities for analyzing macaw files.")
subparsers = parser.add_subparsers(title='title', description='description')

# parser for condense command
parser_c = subparsers.add_parser('condense', help='condense help')
parser_c.add_argument("-d", "--directory", action="store", metavar="<input file directory>", dest="macawDir", help="Directory with all .interpreted.macaw files in it")
parser_c.add_argument("-o", "--output", action="store", metavar="<output file>", dest="outputFile", help="Output file")
parser_c.set_defaults(func=macawCondense)

# parser for interpret command
parser_i = subparsers.add_parser('interpret', help='interpret help')
parser_i.add_argument("-m", "--macaw", action="store", metavar="<marker detection file>", dest="macawFile", help="Macaw marker detection file")
parser_i.add_argument("-o", "--output", action="store", metavar="<output file>", dest="outputFile", help="Output file prefix")
parser_i.set_defaults(func=macawInterpret)

# wrap up and call the function
args = parser.parse_args()
args.func(args)
