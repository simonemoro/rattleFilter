import sys
import pandas as pd
from collections import OrderedDict
from collections import Counter
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--RATTLE_GTF", type=str, help="RATTLE transcripts mapped to reference genome in GTF format")
parser.add_argument("-f", "--RATTLE_FASTQ", type=str, help="RATTLE output 'transcriptome.fq' in FASTQ format")
parser.add_argument("-o", "--output", type=str, help="specify an output directory (default: current folder)")
parser.add_argument("-n", "--minNucleotides", type=int, help="cutoff of transcript length (nt) (default: 100)")
parser.add_argument("-i", "--intersectLength", type=float, help="minimum percentage of transcript overlap to be accepted in a group of transcripts (default: 0.5)")
parser.add_argument("-p", "--percentageReads", type=float, help="minimum percentage of associated reads that a transcript (in a group of overlapping transcripts) must have if it is not the longest one (default: 0.5)")
parser.add_argument("-m", "--minReadLength", type=float, help="cutoff of MINUMUM length that transcript selected by minumum percentage of associated reads (-p) should have. It is a length percentage of the transcript with the highest number of associated reads (default: 0.5)")
parser.add_argument("-M", "--maxReadLength", type=float, help="cutoff of MAXIMUM length that transcript selected by minumum percentage (-p) should have. It is a length percentage of the transcript with the highest number of associated reads (default: 2)")
args = parser.parse_args()

if args.RATTLE_GTF == None or args.RATTLE_FASTQ==None:
		print('filtRATTLEd.py  unable to start:    The fields -g and -f must be specified...')
		sys.exit()

RATTLE_GTF=args.RATTLE_GTF
RATTLE_FASTQ=args.RATTLE_FASTQ

if (args.output):
    outputDir=args.output
    if (outputDir[-1] != "/"):
        outputDir= outputDir+"/"
else:
    outputDir="./"

if(not args.minNucleotides):
    minNucleotides=100
if(not args.intersectLength):
    intersectLength=0.5
if(not args.percentageReads):
    percentageReads=0.5
if(not args.minReadLength):
    minReadLength=0.5
if(not args.maxReadLength):
    maxReadLength=2



RATTLE_GTF = open(RATTLE_GTF,'r') ; RATTLE_GTF = RATTLE_GTF.readlines() ; RATTLE_GTF = list(OrderedDict.fromkeys(RATTLE_GTF))
RATTLE_FASTQ = open(RATTLE_FASTQ,'r') ; RATTLE_FASTQ = RATTLE_FASTQ.readlines()

RATTLE_FASTQHeaders = []
for line in RATTLE_FASTQ:
	if line.startswith('@'):
		RATTLE_FASTQHeaders += [line]

###########################################################################################################################################

# Filtering out smaller transcripts sharing one of the two coordinates 

allLengths = []
allStarts  = []
allEnds    = []
for line in RATTLE_GTF:
	allStarts  += [line.split('\t')[3]]
	allEnds    += [line.split('\t')[4]]
	allLengths += [str(int(line.split('\t')[4]) - int(line.split('\t')[3])+1)]


sameEnds = [k for k,v in Counter(allEnds).items() if v>1] # list of duplicates (clusters with same end)

toRemove = []
for dup in sameEnds:
	indicesDup  = [h for h, w in enumerate(RATTLE_GTF) if dup+'\t' in w]  
	chromosomes = [] ; clusters = []
	for idx in indicesDup:
		chromosomes += [RATTLE_GTF[idx].split('\t')[0]]
		clusters    += [RATTLE_GTF[idx].split('"')[-2]]
	if any(element != clusters[0] for element in clusters):
		keepCluster = []
		C = Counter(clusters)   
		groupSameClusters = [[k,]*v for k,v in C.items()]
		if all(len(group) == 1 for group in groupSameClusters):
			continue
		else:
			for group in groupSameClusters:
				if len(group) == 1:
					keepCluster += [group[0]]
			for j in keepCluster:
				toPopPos = clusters.index(j) ; chromosomes.pop(toPopPos) ; clusters.pop(toPopPos)
	if all(element == clusters[0] for element in clusters) and all(element == chromosomes[0] for element in chromosomes): # if same cluster and same chromosome 
		lengths = []
		for i in range(len(indicesDup)):
			lengths += [int(allLengths[indicesDup[i]])]
		indicesLengths = set(range(len(lengths))) ; indicesLengths.remove(lengths.index(max(lengths)))
		for idxLength in indicesLengths:
			toRemove += [indicesDup[idxLength]]


sameStarts = [k for k,v in Counter(allStarts).items() if v>1]

for dup in sameStarts:
	indicesDup  = [h for h, w in enumerate(RATTLE_GTF) if dup+'\t' in w]  
	chromosomes = [] ; clusters = []
	for idx in indicesDup:
		chromosomes += [RATTLE_GTF[idx].split('\t')[0]] ; clusters += [RATTLE_GTF[idx].split('"')[-2]]
	if any(element != clusters[0] for element in clusters):
		keepCluster = []
		C = Counter(clusters)  
		groupSameClusters = [[k,]*v for k,v in C.items()] 
		if all(len(group) == 1 for group in groupSameClusters):
			continue
		else:
			for group in groupSameClusters:
				if len(group) == 1:
					keepCluster += [group[0]]
			for j in keepCluster:
				toPopPos = clusters.index(j) ; chromosomes.pop(toPopPos) ; clusters.pop(toPopPos)
	if all(element == clusters[0] for element in clusters) and all(element == chromosomes[0] for element in chromosomes):
		lengths = []
		for i in range(len(indicesDup)):
			lengths += [int(allLengths[indicesDup[i]])]
		indicesLengths = set(range(len(lengths))) ; indicesLengths.remove(lengths.index(max(lengths)))
		for idxLength in indicesLengths:
			toRemove += [indicesDup[idxLength]]

counter = -1
for line in RATTLE_GTF:
	counter += 1
	if (int(line.split('\t')[4]) - int(line.split('\t')[3]) + 1) <= minNucleotides :
		toRemove += [counter]  

toRemove = list(set(toRemove))
for toRemove in sorted(toRemove, reverse=True):
	del RATTLE_GTF[toRemove] 


###########################################################################################################################################

#  Create 3 dictionaries:
# 			clusters_features : features of overlapping transcripts ordered in clusters (name, # of associated reads, transcript length)  
#			clusters_union : number coordinates of the union of all the overlapping transcripts
#			clusters_chromosome : chromosome number of transcripts in cluster

counter = 1
clusters_features   = {}
clusters_union  = {}
clusters_chromosome = {}
for line in RATTLE_GTF:
	coords      = set(range(int(line.split('\t')[3]),int(line.split('\t')[4]) + 1))
	length      = max(coords)-min(coords)+1
	cluster_x   = line.split('"')[-2]
	chromosome  = line.split('\t')[0]
	position    = [h for h, w in enumerate(RATTLE_FASTQHeaders) if cluster_x+' ' in w]
	for j in position:
		readsNumber = int(RATTLE_FASTQHeaders[j].split('=')[-1].split('\n')[0])
		if counter == 1:
			clusters_union[counter] = coords ; clusters_features[counter] = [[cluster_x,readsNumber,length]] ; clusters_chromosome[counter] = chromosome
			counter += 1
		else:
			noOverlap = True
			keysChromosome = [k for k,v in clusters_chromosome.items() if v == chromosome]
			clusters_union_1chrom = { x: clusters_union[x] for x in keysChromosome } 
			for group_union in clusters_union_1chrom.values():  
				if len(group_union & coords) > 0 and len(group_union & coords) >= length*intersectLength: # if trascripts overlap and does it with at least 50% of its length (default)
					noOverlap = False
					for key, value in clusters_union_1chrom.items():  
						if value == group_union:
							keyUnion = key
					clusters_union[keyUnion] = group_union.union(coords) 
					clusters_features[keyUnion].append([cluster_x,readsNumber,length])
			if noOverlap == True:
				counter += 1
				clusters_union[counter] = coords ; clusters_features[counter]  = [[cluster_x,readsNumber,length]] ; clusters_chromosome[counter] = chromosome

deleteItem = []
for key,val in clusters_features.items():
	if len(val) == 1:
		deleteItem.append(key)

for i in deleteItem:
	del clusters_features[i] ; del clusters_chromosome[i]
###########################################################################################################################################
				
# Keeping the best transcript 

totalReads = {}
for cluster in clusters_features.keys():
	unique = set(tuple(x) for x in clusters_features[cluster])
	clusters_features[cluster] = [ list(x) for x in unique ] 
	associatedReads = [number[1] for number in clusters_features[cluster]]; total_clusterReads = sum(associatedReads) ; length = [number[2] for number in clusters_features[cluster]]
	pos_max_associatedReads = associatedReads.index(max(associatedReads)) ; pos_max_length = length.index(max(length))
	if pos_max_associatedReads == pos_max_length:
		totalReads[cluster] = clusters_features[cluster][pos_max_length] + [total_clusterReads] + [clusters_chromosome[cluster]]
		clusters_features[cluster].remove(clusters_features[cluster][pos_max_length]) ; clusters_features[cluster].append(clusters_chromosome[cluster])
	else:
		associatedReads_p = [i for i in sorted(associatedReads,reverse=True) if i >= max(associatedReads)*percentageReads]
		if len(associatedReads_p) == 1: # if there are no other transcripts with at least the half of the maximum number of reads associated
			totalReads[cluster] = clusters_features[cluster][pos_max_associatedReads] + [total_clusterReads] + [clusters_chromosome[cluster]]
			clusters_features[cluster].remove(clusters_features[cluster][pos_max_associatedReads]) ; clusters_features[cluster].append(clusters_chromosome[cluster])
		if len(associatedReads_p) > 1:
			indices = [associatedReads.index(i) for i in associatedReads_p[1:len(associatedReads_p)]]  
			lengthIndices = [length[i] for i in indices]
			electedTranscript = None
			for j in range(len(lengthIndices)):   # if element length is in interval: {half length of transcript with max associated reads : double of his length} (default)
				if lengthIndices[j] > clusters_features[cluster][pos_max_associatedReads][2]*minReadLength and lengthIndices[j] < clusters_features[cluster][pos_max_associatedReads][2]*maxReadLength:
					electedTranscript = lengthIndices[j]
					break
			if electedTranscript != None:
				totalReads[cluster] = clusters_features[cluster][indices[lengthIndices.index(electedTranscript)]] + [total_clusterReads] + [clusters_chromosome[cluster]]
				clusters_features[cluster].remove(clusters_features[cluster][indices[lengthIndices.index(electedTranscript)]]) ; clusters_features[cluster].append(clusters_chromosome[cluster]) # remove element in position of maximum length transcript of group of above 50% reads of maximum reads transcript
				length = [number[2] for number in clusters_features[cluster][:-1]] ; pos_max_length = length.index(max(length))
				if totalReads[cluster][:-2] != clusters_features[cluster][pos_max_length]:
					totalReads[cluster][3] = totalReads[cluster][3] - clusters_features[cluster][pos_max_length][1] # substract reads of cluster exiting the group
					totalReads[cluster+0.1] = clusters_features[cluster][pos_max_length] 
					clusters_features[cluster].remove(clusters_features[cluster][pos_max_length])
			if electedTranscript == None:
				totalReads[cluster] = clusters_features[cluster][pos_max_associatedReads] + [total_clusterReads] + [clusters_chromosome[cluster]]
				clusters_features[cluster].remove(clusters_features[cluster][pos_max_associatedReads]) ; clusters_features[cluster].append(clusters_chromosome[cluster])
				length = [number[2] for number in clusters_features[cluster][:-1]] ; pos_max_length = length.index(max(length))
				if totalReads[cluster][:-2] != clusters_features[cluster][pos_max_length]:
					totalReads[cluster][3] = totalReads[cluster][3] - clusters_features[cluster][pos_max_length][1]
					totalReads[cluster+0.1] = clusters_features[cluster][pos_max_length]
					clusters_features[cluster].remove(clusters_features[cluster][pos_max_length])

deleteItem = []
for key,val in clusters_features.items():
	if len(val) <= 0:
		deleteItem.append(key)

for i in deleteItem:
	del clusters_features[i] ; del clusters_chromosome[i]


GTFlinesRemove = []
for cluster in clusters_features.keys():
	chromosome = clusters_features[cluster][-1]
	reducedRATTLEgtf = []
	for line in RATTLE_GTF:
		if line.startswith(chromosome):
			reducedRATTLEgtf += [line]
	for transcript in clusters_features[cluster]:
		if type(transcript) == list:
			clusterName   = transcript[0]
			clusterLength = transcript[2]
			for line in reducedRATTLEgtf:
				GTFlineLength = int(line.split('\t')[4]) - int(line.split('\t')[3]) + 1
				if clusterName in line and clusterLength == GTFlineLength:
					GTFlinesRemove += [line]

GTFlinesRemove = list(set(GTFlinesRemove))
for deleteLine in GTFlinesRemove:
	RATTLE_GTF.remove(deleteLine)




final_GTF = open(outputDir+'filtered.gtf', 'w') ; final_GTF.writelines(RATTLE_GTF) ; final_GTF.close()

exit(0)
