#! /usr/bin/python

######################################################################################
## Python scripts to simulate 3D tumor growth and SMALT lineage tracing system	  ##
## using an agent-based model. Deme subdivision is assumed in order to model cell   ##
## mixing and spatial contraint.													##
##																				  ## 
## *Spatial model: pripheral growth												 ##
## Author: Zheng Hu																 ##
## Date: 11/06/2022
## Modified by Duo Xie
## duo.xie@siat.ac.cn
######################################################################################
########conda activate py2
########conda deactivate
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO, AlignIO
import sys,os,math,random
import numpy as np
from collections import Counter
#import sets
import matplotlib.pyplot as plt
#from scipy import stats

class deme():
	def __init__(self):
		self.present= 0		 ## whether the deme is empty or occupied: 0-empty;1-occupied
		self.neutral = []		## the neutral founder lineage after tumor tranformation
		self.advant = []		## the advantageous cells having one driver mutation
		self.empty = 27		 ## number of empty sites in the neighbourhood

def createLattice(d):
	"""
	Create a 3D cubic lattice with side length of 2d+1 where each site contains a empty deme.
	"""
	lattice = {}
	for x in range(0,2*d+1):
		for y in range(0,2*d+1):
			for z in range(0,2*d+1):
				lattice[(x,y,z)] = deme()
	return lattice


def neighbor26(xxx_todo_changeme):
	"""
	Moore neighbourhood: 26 neighbour sites of (a,b,c)
	"""
	(a,b,c) = xxx_todo_changeme
	neighbor = [(a+i,b+j,c+k)
				for i in [-1,0,1]
				for j in [-1,0,1]
				for k in [-1,0,1]
				if not (i==0 and j==0 and k==0)]
	
	return neighbor


#def neighbor6((a,b,c)):
#	"""
#	von Neumann neighbourhood: 6 neighbour sites of (a,b,c).
#	"""
#	neighbor = [(a-1, b, c),(a+1, b, c),(a, b-1, c),(a, b+1, c),(a, b, c-1),(a, b, c+1)]
#	return neighbor


def localNeighbor(xxx_todo_changeme1,r):
	"""
	A function to search the local neighbour sites of (a,b,c) within an area of radius r in the 3D cubic lattice.
	"""
	(a,b,c) = xxx_todo_changeme1
	neighbor = []
	for x in range(-r,r+1):
		for y in range(-r,r+1):
			for z in range(-r,r+1):
				if pow(x,2)+pow(y,2)+pow(z,2) < pow(r+1,2):
					neighbor += [(a+x,b+y,c+z)]
	return neighbor


def traceLineage(mlineage,mutid):
	"""
	A function to obtain the mutational lineage of a cell from the mutation id of the most recently occurred mutation in the cell. 
	For example, the input ID (most recently occurred mutation) of target cell is "100" and the output is "1-12-35-56-100", which is the mutation lineage of the cell
	
	mlineage - the list that could be used to recover the mutational lineage given the most recent mutation id of a lineage
	mutid - the mutation ID of the most recently occurred mutation in the cell
	"""
	recent_muts = mutid.split(',')  # it is possible that multiple mutations occur during in a cell division. For instance, the mutation id of most recently occurred mutations is "100,101"
	recent_muts = [int(t) for t in recent_muts]
	first_mut = recent_muts[0]	  # the first mutation id in a multi-mutation event
	trace = []
	while first_mut > 0:
		trace += recent_muts
		recent_muts = mlineage[first_mut].split(',')
		recent_muts = [int(t) for t in recent_muts]
		first_mut = recent_muts[0]
	return trace

	
def lowerORupper(value):
	"""
	A function to choose the upper or lower integral value given a non-integral number
	"""
	lower_int = int(value)
	upper_int = lower_int+1
	if random.random() < value-lower_int:
		return upper_int
	else:
		return lower_int


def initiateFirstDeme_t1(maxsize,lineage,current_id,current_driver_id,sfit,barcode_muts):
	"""
	The growth of the initial deme from a single transformed tumor cell via a random discrete-time birth-death process
	t1 - one-tier driver model
	maxsize - size limit of a deme
	lineage - a list that stores the lineage information of mutations
	current_id - the starting mutation ID
	sfit - selection fitness of advantageous mutations
	"""
	#print(barcode_muts)
#	neu_list = ["0" for i in range(0, rand_int[0])]
 #   neu_list = ",".join(neu_list)
  #  neu_list = [neu_list]
	neu_list = [str(current_id)]
	adv_list = []
	current_deme_size = 1
	#total_mut_rate2 = 0
	total_mut_rate2 = total_mut_rate
	while current_deme_size < maxsize:
		n1,n2 = len(neu_list),len(adv_list)
		n1_double = 0
		n2_double = 0
		neu_div_list_double = []
		adv_div_list_double = []
		neu_qui_list = []
		adv_qui_list = []
		if n1 > 0:
			neu_qui_number = int(n1*quies_rate)
			neu_div_number = int(n1*birth_rate+1)							#number of dividing cells in this generation
			neu_pass_number = neu_qui_number+neu_div_number
			random.shuffle(neu_list)
			neu_qui_list = neu_list[0:neu_qui_number]
			neu_div_list = neu_list[neu_qui_number:neu_pass_number]
			neu_div_list_double = neu_div_list*2
			n1_double = len(neu_div_list_double)
		
		if n2 > 0:
			adv_qui_number = lowerORupper(n2*(quies_rate-sfit))		#number of dividing cells in this generation		
			adv_div_number = lowerORupper(n2*(birth_rate+sfit))
			adv_pass_number = adv_qui_number+adv_div_number
			if adv_pass_number > n2:
				adv_qui_number = int(n2*(quies_rate-sfit))
				adv_pass_number = adv_qui_number+adv_div_number
			random.shuffle(adv_list)
			adv_qui_list = adv_list[0:adv_qui_number]
			adv_div_list = adv_list[adv_qui_number:adv_pass_number]
			adv_div_list_double = adv_div_list*2
			n2_double = len(adv_div_list_double)
		if n1_double > 0:
			new_mut1 = np.random.poisson(total_mut_rate2*n1_double)
			mut_assig1 = Counter(np.random.choice(n1_double,new_mut1))
			for x1 in list(mut_assig1.keys()):
				nmut = mut_assig1[x1]
				new_mut1 = list(range(current_id+1,current_id+1+nmut)) 
				mut_str = ",".join(map(str,new_mut1))
				#if nmut > 1:
				#	for t in new_mut1:
				#		multi_events[str(t)] = mut_str
				for xn in range(0,nmut):
					current_id += 1
					lineage += [neu_div_list_double[x1]]
					if random.random() < hotspot_chance:
						barcode_mut_pos = random.choice(list(range(0,hotspot_number)))+1
						barcode_muts += [barcode_mut_pos]
					else:
						barcode_mut_pos = random.choice(list(range(hotspot_number,non_hotspot_number+hotspot_number)))+1
						barcode_muts +=  [barcode_mut_pos]
				neu_div_list_double[x1] = mut_str
				
		if n2_double > 0:
			new_mut2 = np.random.poisson(total_mut_rate2*n2_double)
			mut_assig2 = Counter(np.random.choice(n2_double,new_mut2))
			for x2 in list(mut_assig2.keys()):
				nmut = mut_assig2[x2]
				new_mut2 = list(range(current_id+1,current_id+1+nmut))
				mut_str = ",".join(map(str,new_mut2))
				#if nmut > 1:
				#	for t in new_mut2:
				#		multi_events[str(t)] = mut_str
				for xn in range(0,nmut):
					current_id += 1
					lineage += [adv_div_list_double[x2]]
					if random.random() < hotspot_chance:
						barcode_mut_pos = random.choice(list(range(0,hotspot_number)))+1
						barcode_muts += [barcode_mut_pos]
					else:
						barcode_mut_pos = random.choice(list(range(hotspot_number,non_hotspot_number+hotspot_number)))+1
						barcode_muts += [barcode_mut_pos]
					#cut_type = random.choice(range(0,cut_type_number))+1
				adv_div_list_double[x2] = mut_str
		
		if random.random() < adv_rate*n1_double:
			current_driver_id += 1
			#current_n1 = len(neu_list)
			#lineage += [str(neu_list[current_n1-1])]
			adv_div_list_double += [neu_div_list_double[n1_double-1]]
			neu_div_list_double = neu_div_list_double[0:n1_double-1]
		
		neu_list=neu_qui_list+neu_div_list_double
		adv_list=adv_qui_list+adv_div_list_double
		current_deme_size=len(neu_list)+len(adv_list)
	return neu_list,adv_list,current_id,current_driver_id,lineage,barcode_muts


def demeGrowthFission_t1(neu_list,adv_list,lineage,current_id,current_driver_id,current_deme_number,sfit,barcode_muts):
	"""
	A function to simulate deme growth and fission and keep track of the mutational lineages
	t1 - one-tier driver model
	"""
	current_deme_size = len(neu_list)+len(adv_list)
	while current_deme_size < 2*deme_size:
		n1,n2 = len(neu_list),len(adv_list)
		n1_double = 0
		n2_double = 0
		neu_div_list_double = []
		adv_div_list_double = []
		neu_qui_list = []
		adv_qui_list = []
		if n1 > 0:
			neu_qui_number = lowerORupper(n1*quies_rate)
			neu_div_number = lowerORupper(n1*birth_rate)							#number of dividing cells in this generation
			neu_pass_number = neu_qui_number+neu_div_number
			if neu_pass_number > n1:
				neu_qui_number = int(n1*quies_rate)
				neu_pass_number = neu_qui_number+neu_div_number
			random.shuffle(neu_list)
			neu_qui_list = neu_list[0:neu_qui_number]
			neu_div_list = neu_list[neu_qui_number:neu_pass_number]
			neu_div_list_double = neu_div_list*2
			n1_double = len(neu_div_list_double)

		if n2 > 0:
			adv_qui_number = lowerORupper(n2*(quies_rate-sfit))		#number of dividing cells in this generation		
			adv_div_number = lowerORupper(n2*(birth_rate+sfit))
			adv_pass_number = adv_qui_number+adv_div_number
			if adv_pass_number > n2:
				adv_qui_number = int(n2*(quies_rate-sfit))
				adv_pass_number = adv_qui_number+adv_div_number
			random.shuffle(adv_list)
			adv_qui_list = adv_list[0:adv_qui_number]
			adv_div_list = adv_list[adv_qui_number:adv_pass_number]
			adv_div_list_double = adv_div_list*2
			n2_double = len(adv_div_list_double)
		
		if n1_double > 0:
			new_mut1 = np.random.poisson(total_mut_rate*n1_double)
			mut_assig1 = Counter(np.random.choice(n1_double,new_mut1))
			for x1 in list(mut_assig1.keys()):
				nmut = mut_assig1[x1]
				new_mut1 = list(range(current_id+1,current_id+1+nmut))
				mut_str = ",".join(map(str,new_mut1))
				#if nmut > 1:
				#	for t in new_mut1:
				#		multi_events[str(t)] = mut_str
				for xn in range(0,nmut):
					current_id += 1
					lineage += [neu_div_list_double[x1]]
					if random.random() < hotspot_chance:
						barcode_mut_pos = random.choice(list(range(0,hotspot_number)))+1
						barcode_muts += [barcode_mut_pos]
					else:
						barcode_mut_pos = random.choice(list(range(hotspot_number,non_hotspot_number+hotspot_number)))+1
						barcode_muts +=  [barcode_mut_pos]
				neu_div_list_double[x1] = mut_str
		
		if n2_double > 0:
			new_mut2 = np.random.poisson(total_mut_rate*n2_double)
			mut_assig2 = Counter(np.random.choice(n2_double,new_mut2))
			for x2 in list(mut_assig2.keys()):
				nmut = mut_assig2[x2]
				new_mut2 = list(range(current_id+1,current_id+1+nmut))
				mut_str = ",".join(map(str,new_mut2))
				#if nmut > 1:
				#	for t in new_mut2:
				#		multi_events[str(t)] = mut_str
				for xn in range(0,nmut):
					current_id += 1
					lineage += [adv_div_list_double[x2]]
					if random.random() < hotspot_chance:
						barcode_mut_pos = random.choice(list(range(0,hotspot_number)))+1
						barcode_muts += [barcode_mut_pos]
					else:
						barcode_mut_pos = random.choice(list(range(hotspot_number,non_hotspot_number+hotspot_number)))+1
						barcode_muts +=  [barcode_mut_pos]
				adv_div_list_double[x2] = mut_str
		
		if random.random() < adv_rate*n1_double:
			current_driver_id += 1
			#current_n1 = len(neu_list)
			#lineage += [str(neu_list[current_n1-1])]
			adv_div_list_double += [neu_div_list_double[n1_double-1]]
			neu_div_list_double = neu_div_list_double[0:n1_double-1]
		
		neu_list=neu_qui_list+neu_div_list_double
		adv_list=adv_qui_list+adv_div_list_double
		current_deme_size=len(neu_list)+len(adv_list)
	
	random.shuffle(neu_list)
	if len(neu_list) > 0:
		offspring_neu = np.random.binomial(len(neu_list),0.5)
	else:
		offspring_neu = 0
	neu_list1=neu_list[0:offspring_neu]
	neu_list2=neu_list[offspring_neu:len(neu_list)]
	
	random.shuffle(adv_list)
	if len(adv_list) > 0:
		offspring_adv = np.random.binomial(len(adv_list),0.5)
	else:
		offspring_adv = 0
	adv_list1=adv_list[0:offspring_adv]
	adv_list2=adv_list[offspring_adv:len(adv_list)]
	
	return neu_list1,neu_list2,adv_list1,adv_list2,current_id,current_driver_id,lineage,barcode_muts


def deme2cells(sp,sample_keys):
	"""
	Get the cell from a pool/bulk of demes
	"""
	all_cur_id = []
	for key in sample_keys:
		all_cells = list(sp[key].neutral + sp[key].advant)
		all_cur_id += all_cells

	all_cur_id2 = [x[0] for x in all_cur_id]
	
	return all_cur_id2


def seqProcessing(sp,sample_keys,mlineage,size_par,mean_depth,purity):
	"""
	Model the random sampling process in NGS and report the sequencing allele frequencies in a sample of cells
	
	sp- the lattice space
	sample_keys- the locations for the demes in a bulk sample
	size_par- variance parameter for negative-binomial distribution
	mean_depth- the mean depth of the sequencing
	purity- tumor purity
	"""
	all_cur_id = []									 # all most recently occurred mutations
	all_mut_id = []									 # all mutations in the sampled cells
	for key in sample_keys:
		smuts = list(sp[key].neutral + sp[key].advant)
		all_cur_id += smuts
	sample_size = 10000								 # the number of cells for sequencing analysis
	sample_id = random.sample(all_cur_id,sample_size)
	id_count = Counter(sample_id)
	for x in list(id_count.keys()):
		xlineage = traceLineage(mlineage,x)
		all_mut_id += xlineage*id_count[x]
	mut_count = Counter(all_mut_id)
	prob_par=size_par*1.0/(size_par+mean_depth)
	sampleAF = {}									   # a dictionary storing the mutation IDs and corresponding depth and allele frequency the seq data
	for x in list(mut_count.keys()):
		true_af = mut_count[x]*0.5*purity/sample_size   # the true allele frequency in the sample
		if true_af > 0.001:							 # filter mutations with very low frequency that is not detectable by ~100X sequencing depth
			site_depth = np.random.negative_binomial(size_par,prob_par)
			if site_depth >= 15:						# seq depth cutoff for "calling" a mutation
				var_reads = np.random.binomial(site_depth,true_af)
				seq_af = var_reads*1.0/site_depth
				if var_reads >= 4:					  # variant reads cutof for "calling" a mutation
					sampleAF[str(x)] = (site_depth,seq_af)
	return sampleAF


def highMuts(sp,position,mlineage,cutoff):
	"""
	Obtain the high-frequency mutations (vaf>cutoff) in a particular deme
	
	sp - the lattice space
	position - the location of the deme
	mlineage - mutation lineage dictionary
	cutoff - the VAF cutoff for a "high-frequency" mutation, e.g. 0.4
	"""
	all_cur_id = sp[position].neutral + sp[position].advant
	id_count = Counter(all_cur_id)
	
	all_mut_id = []
	sample_size = len(all_cur_id)
	#sample_id = random.sample(all_cur_id,sample_size)
	#id_count = Counter(sample_id)
	for y in list(id_count.keys()):
		xlineage = traceLineage(mlineage,y)
		all_mut_id += xlineage*id_count[y]
	mut_count = Counter(all_mut_id)
	highAF_muts = []
	for x in list(mut_count.keys()):
		allele_freq = mut_count[x]*1.0/sample_size
		if allele_freq > cutoff:
			highAF_muts += [(int(x),round(allele_freq,3))]
			#highAF_muts[str(x)] = allele_freq
	
	return highAF_muts


def highMuts_bc(sp,position,mlineage,barcode_muts,cutoff):
	"""
	Obtain the high-frequency mutations (vaf>cutoff) in a particular deme
	
	sp - the lattice space
	position - the location of the deme
	mlineage - mutation lineage dictionary
	cutoff - the VAF cutoff for a "high-frequency" mutation, e.g. 0.4
	"""
	all_cur_id = sp[position].neutral + sp[position].advant
	id_count = Counter(all_cur_id)
	
	all_mut_id = []
	sample_size = len(all_cur_id)
	#sample_id = random.sample(all_cur_id,sample_size)
	#id_count = Counter(sample_id)
	for y in list(id_count.keys()):
		xlineage = traceLineage(mlineage,y)
		all_mut_id += xlineage*id_count[y]
	all_bc_muts = [barcode_muts[i] for i in all_mut_id]
	mut_count = Counter(all_mut_id)
	bc_mut_count = Counter(all_bc_muts)
	highAF_bc_muts = []
	for x in list(bc_mut_count.keys()):
		allele_freq = bc_mut_count[x]*1.0/sample_size
		if allele_freq > cutoff:
			highAF_bc_muts += [(int(x),round(allele_freq,3))]
	
	return highAF_bc_muts


def highDriverMuts(sp,position,cutoff):
	"""
	Obtain the high-frequency driver mutations (vaf>cutoff) in a particular deme
	
	sp - the lattice space
	position - the location of the deme
	mlineage - mutation lineage dictionary
	cutoff - the VAF cutoff for a "high-frequency" driver mutation, e.g. 0.1
	"""
	all_cur_id = sp[position].neutral + sp[position].advant
	#sample_size = 100
	#sample_id = random.sample(all_cur_id,sample_size)
	all_mut_id1 = [x[1] for x in all_cur_id]
	all_mut_id2 = [] 
	for k in all_mut_id1:
		all_mut_id2 += list(k.split(","))
	
	mut_count = Counter(all_mut_id2)
	highAF_muts = []
	for x in list(mut_count.keys()):
		allele_freq = mut_count[x]*1.0/len(all_cur_id)
		if allele_freq > cutoff:
			highAF_muts += [int(x)]
	
	highAF_muts = sorted(highAF_muts)
	if highAF_muts[0]==0 and len(highAF_muts)>1:
		highAF_muts = highAF_muts[1:len(highAF_muts)]
	
	return highAF_muts


def pubMutGenerator(n,size_par,mean_depth,purity):
	"""
	A function to generate the public clonal mutations occured during the multi-step tumorigenesis before transformation.
	
	n- number of clonal mutations
	size_par- variation parameter in the negative binomial distribution
	mean_death- mean seq depth
	"""
	prob_par=size_par*1.0/(size_par+mean_depth)
	mean_af = 0.5*purity
	depth_pub = []
	maf_pub = []
	for k in range(0,n):
		correct = 0
		while correct == 0:
			site_depth = np.random.negative_binomial(size_par,prob_par)
			if site_depth >= 15:
				correct =1
		var_reads = np.random.binomial(site_depth,mean_af)
		site_maf = var_reads*1.0/site_depth
		depth_pub += [site_depth]
		maf_pub += [site_maf]
	return depth_pub,maf_pub


def localSampling(region,sample_number,cutoff):
	"""
	A function to sampling the locations of multiple bulk samples in a local region.
	"""
	success = 0
	while success == 0:
		locations = random.sample(region,sample_number)
		repeat = sample_number*(sample_number-1)
		minall = 999
		for x in range(0,repeat):
			rs = random.sample(locations,2)
			min_distance = min([abs(rs[0][0]-rs[1][0]),abs(rs[0][1]-rs[1][1]),abs(rs[0][2]-rs[1][2])])
			if min_distance < minall:
				minall = min_distance
		if min_distance > 2*cutoff:
			success = 1
	return locations


def bulkTissueSampling(sp,location,radius):
	"""
	A function to sampling a bulk sample in a local region.
	"""
	local_region = localNeighbor(location,radius)
	bulk_tissue = []
	for x in local_region:
		if sp[x].present == 1:
			bulk_tissue += [x]
	return bulk_tissue


def lineageDashLink(mlist):
	"""
	Transform the mutation lineage from list (e.g [1,3,10,20]) to dash-linked string (e.g. 1-3-10-20)
	"""
	if len(mlist) > 0:
		dstring = str(mlist[0])
		for x in mlist[1:len(mlist)]:
			dstring += "-"
			dstring += str(x)
		return dstring
	else:
		return "0"


def StringLink(mlist):
	"""
	[1,2,3,4] -> "1,2,3,4"
	"""
	if len(mlist) > 0:
		dstring = str(mlist[0])
		for x in mlist[1:len(mlist)]:
			dstring += ":"
			dstring += str(x)
		return dstring
	else:
		return "0"


def missingDepth(vafdata,absent_muts,mean_depth):
	"""
	Randomly generate the sequencing depth for the mutation-absent sites across samples
	"""
	for x in absent_muts:
		done = 0
		while done == 0:
			missing_depth = np.random.negative_binomial(2,2.0/(2+mean_depth))
			if missing_depth >= 15:
				done = 1
		vafdata[str(x)] = (missing_depth,0)
	
	return vafdata


def cut2allele(cutted,number_sites):
	allele = [0 for i in range(0,number_sites)]
	for ct in cutted:
		cposition = int(ct[0])
		ctype = int(ct[1])
		allele[cposition-1] = ctype

	return allele


############# main script to simulate a tumor and SLOTH lineage tracing data #########
###parameter intiation###
processID = os.getpid()
repl = int(sys.argv[1])			 # replication of simulation
hotspot_number = int(sys.argv[2])				 # the number of AID hotspot editing sites in each cell
s_coef = float(sys.argv[3])						# selection coefficient e.g. 0 (neutral), 0.01, 0.02, ...
if sys.argv[3] == '0':
	s_coef = 0
else:
	s_coef = float(sys.argv[3])
ncells = 1000
rd = 30							 # the radius of the pre-created 3D space
deme_size = 1000					
final_tumor_size = pow(10,7)		# the number of cells in the final tumor
final_deme_number = int(final_tumor_size/deme_size)	# the final number of demes in the tumor
birth_rate = 0.4					# the birth probability at each cell generation during tumor growth
quies_rate = 0.5					# the quiescent probability at each cell generation during tumor growth
death_rate = round(1-birth_rate-quies_rate, 1)
print("birth rate, quiescent rate and death rate=",birth_rate,quies_rate,death_rate)

#mut_rate = 0.0001				   # the mutation rate by AID at each site
mut_rate = pow(10, -3.5)				   # the mutation rate by AID at each site
non_hotspot_number = 1000		   # the number of AID editing sites in each cell 
hotspot_increase = 20
total_mut_rate = mut_rate*(non_hotspot_number+hotspot_number*hotspot_increase)
hotspot_chance = hotspot_number*hotspot_increase*1.0/(hotspot_number*hotspot_increase+non_hotspot_number)

adv_rate = pow(10,-4) if s_coef > 0 else 0

#s_birth_rate = birth_rate+s_coef
#s_quies_rate = quies_rate-s_coef
#s_coef_adj = round(math.log(1+s_coef)/math.log(2*birth_rate),3)
#print s_coef,s_coef_adj
#percentage = int(s_coef*100)	   # the percentage form of the selection

mut_id = 0						  # the ids of each mutation event, 1,2,3,4,5,...
driver_mut_id = 0				   # the ids of advantageous driver mutations (modeling clonal selection)
mutlineage = ['0']				  # the lineage tracer
barcode_mut_site = [0]			  # the mutation positions on general SMALT barcode sequence
filename = f'mu_{mut_rate}_birth_rate_{birth_rate}_death_rate_{death_rate}_s_coef_{s_coef}_{processID}.fasta'
low = 2
high = 20
size = 1
rand_int = np.random.randint(low=low, high=high, size=size)
#	print(type(rand_int))
clonal_mut_site = []
for i in range(0,rand_int[0]):
	if random.random() < hotspot_chance:
		barcode_mut_pos = random.choice(list(range(0,hotspot_number)))+1
		clonal_mut_site += [barcode_mut_pos]
	else:
		barcode_mut_pos = random.choice(list(range(hotspot_number,non_hotspot_number+hotspot_number)))+1
		clonal_mut_site +=  [barcode_mut_pos]
print(clonal_mut_site)
######################################################################################
first_neu,first_adv,mut_id,driver_mut_id,mutlineage,barcode_mut_site = initiateFirstDeme_t1(deme_size,mutlineage,mut_id,driver_mut_id,s_coef,barcode_mut_site)  #the growth of the fisrt deme from single transformed cell

print("# of neutral and advantageous mutations in the first deme:",len(first_neu),len(first_adv))

space = createLattice(rd)
space[(rd,rd,rd)].present = 1
space[(rd,rd,rd)].empty = 26
space[(rd,rd,rd)].neutral = list(first_neu)
space[(rd,rd,rd)].advant = list(first_adv)
current_keys = [(rd,rd,rd)]
current_deme_number = 1								 #current deme number
surface_keys = [(rd,rd,rd)]
surface_times = [1.0]
surface_deme_number = 1
current_time = 0

while current_deme_number < final_deme_number:
	deme_index = surface_times.index(min(surface_times))
	ckey = surface_keys[deme_index]
	
	ctime = surface_times[deme_index]
	current_time += ctime
	surface_times = [xt-ctime for xt in surface_times]
	#if current_deme_number%100 == 0:
#		print(current_time,current_deme_number,surface_deme_number,mut_id,len(mutlineage)-1,len(barcode_mut_site)-1,driver_mut_id)
	
	nei_sites = neighbor26(ckey)  # neighbor sites of (rx,ry,rz)
	empty_sites = [key for key in nei_sites if space[key].present == 0]					# the empty neighbor sites
	
	if len(empty_sites) > 0:
		num_neu = len(space[ckey].neutral)
		num_advant = len(space[ckey].advant)
		faction_advant = num_advant*1.0/(num_neu+num_advant)
		next_time = np.random.exponential(1.0/(1+faction_advant))
		#next_time = np.random.exponential(1.0)
		surface_times[deme_index] = next_time
		
		rand_prob = random.random()
		if rand_prob < 1-math.exp(-len(empty_sites)*0.25): # the probability for a deme to grow and divide is proportional to the # of empty neighbor sites
			pre_neu = list(space[ckey].neutral)
			pre_adv = list(space[ckey].advant)
			post_neu_l1,post_neu_l2,post_adv_l1,post_adv_l2,mut_id,driver_mut_id,mutlineage,barcode_mut_site = demeGrowthFission_t1(pre_neu,pre_adv,mutlineage,mut_id,driver_mut_id,current_deme_number,s_coef,barcode_mut_site)
			space[ckey].neutral = list(post_neu_l1)
			space[ckey].advant = list(post_adv_l1)
			
			nkey = random.choice(empty_sites)
			space[nkey].neutral = list(post_neu_l2)
			space[nkey].advant = list(post_adv_l2)
			space[nkey].present = 1
			current_keys += [nkey]
			current_deme_number += 1
			
			next_nei_sites = neighbor26(nkey)
			next_empty_sites = []
			for key in next_nei_sites:
				if space[key].present == 1:
					space[key].empty -= 1
					if space[key].empty == 0:
						remove_index = surface_keys.index(key)
						del surface_keys[remove_index]
						del surface_times[remove_index]
						surface_deme_number -= 1
				else:
					next_empty_sites += [key]
			if len(next_empty_sites) > 0:
				num_neu = len(space[nkey].neutral)
				num_advant = len(space[nkey].advant)
				faction_advant = num_advant*1.0/(num_neu+num_advant)
				next_time = np.random.exponential(1.0/(1+faction_advant))
				#next_time = np.random.exponential(1.0)
				
				surface_keys += [nkey]
				surface_deme_number += 1
				surface_times += [next_time]

			space[nkey].empty = len(next_empty_sites)

			#if current_deme_number == met_timing:
			#	met_ancestor_list = list(post_neu_l2+post_adv_l2+post_adv2_l2)
			#	met_founder = random.choice(met_ancestor_list)
		
	else:
		print("something is wrong!")
		break

####visulization of spatial clonal structure in the central slice###
sample_demes = random.sample(current_keys,10)

sample_deme_cells = []
for key in sample_demes:
#	print(key)
	cells0 = space[key].neutral+space[key].advant
#	print(cells0)
	sample_deme_cells += cells0

sample_cells = random.sample(sample_deme_cells,ncells)
mut_array_pos = []
barcodeseq = np.zeros([ncells+1, 3004])
count = 0
repl = str(repl)
print("barcode mutation number, mutation lineage, mutation positions")
file_name= f'./mu_{mut_rate}_birth_rate_{birth_rate}_death_rate_{death_rate}_s_coef_{s_coef}_hotspot_{hotspot_number}_{processID}.lineage'
with open(file_name, 'w') as file:
	for k in sample_cells:
		mut_lin = traceLineage(mutlineage,k)
		mut_position = [barcode_mut_site[x] for x in mut_lin]
		mut_position = mut_position + clonal_mut_site
		mut_position = np.array(mut_position)
		mut_position = np.unique(mut_position)
		mut_position = mut_position.astype(int)
		barcodeseq[count][mut_position] = 1
		count += 1
		#print(len(mut_lin),mut_lin,mut_position,file=sys.stderr)
		mut_lin.sort()
		file.write(f"P_{repl}_{count}\t{mut_lin}\n")
barcodeseq = barcodeseq.astype(int)
numbers = np.arange(1, ncells+1)
names = np.char.add(np.char.add(np.char.add("P", str(repl)), "_"), numbers.astype(str))
names = np.append(names, 'ref')
sequences = ["".join([str(x) for x in row]) for row in barcodeseq]
records = [SeqRecord(Seq(seq), id=name, description="") for seq, name in zip(sequences, names)]


alignment = MultipleSeqAlignment(records)

filename= f'./mu_{mut_rate}_birth_rate_{birth_rate}_death_rate_{death_rate}_s_coef_{s_coef}_hotspot_{hotspot_number}_{processID}.fasta'
with open(filename, "w") as handle:
	AlignIO.write(alignment, handle, "fasta")

'''
#########visulization of clones########
sample_slice = [key for key in current_keys if key[2] == rd]

#####visualization for general (non-barcode) mutations
slice_muts = {} #high frequency mutations in each deme
mutposVAF = {}
uniq_high_muts = []
for key in sample_slice:
	highAFmuts = highMuts(space, key, mutlineage, 0.05)
	uniq_high_muts += [mut[0] for mut in highAFmuts]
	for mut in highAFmuts:
		if mut[0] not in list(mutposVAF.keys()):
			mutposVAF[mut[0]] = [(key[0],key[1],mut[1])]
		else:
			mutposVAF[mut[0]] += [(key[0],key[1],mut[1])]
	#VAF = mutposVAF(mutpos, hi
	slice_muts[key] = [x[0] for x in highAFmuts]
	#print key,len(highAFmuts)

uniq_high_muts = list(set(uniq_high_muts))
#print("number of uniq high mutations=",len(uniq_high_muts),len(mutposVAF))

wholeVAF = {}
whole_size = 0
for key in sample_slice:
	deme_cell_number = len(space[key].neutral+space[key].advant)
	whole_size += deme_cell_number

for mut in list(mutposVAF.keys()):
	clone_size = 0
	for k in mutposVAF[mut]:
		key1 = (int(k[0]),int(k[1]),rd)
		deme_cell_number = len(space[key1].neutral+space[key1].advant)
		clone_size = deme_cell_number*k[2]

	wholeVAF[mut] = clone_size*1.0/whole_size
	print(mut,clone_size*1.0/whole_size)

wholeVAF = sorted(list(wholeVAF.items()), key = lambda kv: kv[1], reverse=True)
print((len(wholeVAF)))

fig = plt.figure(figsize=(6, 6))
axes = fig.subplots(3, 3, sharex=True, sharey=True)
#fig.subplots_adjust(hspace=0.05, wspace=0.05)
top9_muts = []
for i in range(9):
	mut = wholeVAF[i][0]
	top9_muts += [mut]
	x = []
	y = []
	freq = []
	mutfreq = {}
	for key in sample_slice:
		x.append(key[0])
		y.append(key[1])
		if mut in slice_muts[key]:
			cmut_info = mutposVAF[mut]
			for z in cmut_info:
				if z[0] == key[0] and z[1] == key[1]:
					freq.append(z[2])
			#for mut2 in highAFmuts1:
			#	mutfreq[mut2[0]] = mut2[1]
		else:
			freq.append(0)
	print(max(freq),mut,barcode_mut_site[mut])
	subfig = axes[i//3, i%3].scatter(x, y, c=freq, vmin=0, vmax=1, cmap='Reds', s=10)
	axes[i//3, i%3].set_title(mut)
	axes[i//3, i%3].axis('off')

cbar = fig.add_axes([0.95, 0.35, 0.01, 0.3]) 
fig.colorbar(subfig, cax=cbar)
#plt.show()
fig.savefig("top9_general_muts.pdf") 


#####visualization for barcode mutations
slice_muts = {} #high frequency mutations in each deme
mutposVAF = {}
uniq_high_muts = []
for key in sample_slice:
	highAFmuts = highMuts_bc(space, key, mutlineage, barcode_mut_site, 0.05)
	uniq_high_muts += [mut[0] for mut in highAFmuts]
	for mut in highAFmuts:
		if mut[0] not in list(mutposVAF.keys()):
			mutposVAF[mut[0]] = [(key[0],key[1],mut[1])]
		else:
			mutposVAF[mut[0]] += [(key[0],key[1],mut[1])]
	#VAF = mutposVAF(mutpos, hi
	slice_muts[key] = [x[0] for x in highAFmuts]
	#print key,len(highAFmuts)

uniq_high_muts = list(set(uniq_high_muts))
print("number of uniq high barcode mutations=",len(uniq_high_muts),len(mutposVAF))

wholeVAF = {}
whole_size = 0
for key in sample_slice:
	deme_cell_number = len(space[key].neutral+space[key].advant)
	whole_size += deme_cell_number

for mut in list(mutposVAF.keys()):
	clone_size = 0
	for k in mutposVAF[mut]:
		key1 = (int(k[0]),int(k[1]),rd)
		deme_cell_number = len(space[key1].neutral+space[key1].advant)
		clone_size = deme_cell_number*k[2]

	wholeVAF[mut] = clone_size*1.0/whole_size
	print(mut,clone_size*1.0/whole_size)

wholeVAF = sorted(list(wholeVAF.items()), key = lambda kv: kv[1], reverse=True)
print((len(wholeVAF)))

fig = plt.figure(figsize=(6, 6))
axes = fig.subplots(3, 3, sharex=True, sharey=True)
#fig.subplots_adjust(hspace=0.05, wspace=0.05)
for i in range(9):
	#mut = wholeVAF[i][0]
	mut0 = top9_muts[i]
	mut = barcode_mut_site[mut0]
	x = []
	y = []
	freq = []
	mutfreq = {}
	for key in sample_slice:
		x.append(key[0])
		y.append(key[1])
		if mut in slice_muts[key]:
			cmut_info = mutposVAF[mut]
			for z in cmut_info:
				if z[0] == key[0] and z[1] == key[1]:
					freq.append(z[2])
			#for mut2 in highAFmuts1:
			#	mutfreq[mut2[0]] = mut2[1]
		else:
			freq.append(0)
	print(max(freq))
	subfig = axes[i//3, i%3].scatter(x, y, c=freq, vmin=0, vmax=1, cmap='Reds', s=10)
	axes[i//3, i%3].set_title(mut)
	axes[i//3, i%3].axis('off')

cbar = fig.add_axes([0.95, 0.35, 0.01, 0.3]) 
fig.colorbar(subfig, cax=cbar)
#plt.show()
fig.savefig("top9_barcode_muts.pdf") 
'''
