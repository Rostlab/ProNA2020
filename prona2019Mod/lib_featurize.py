import math
import sys	#For sys.maxint

"""
"""

class InputError(Exception): pass


#Biochemical properties of aas
#Mass and volume taken from http://prowl.rockefeller.edu/aainfo/contents.htm,
#hyd(dropathy) according to Kyte-Doolittle (e.g. http://en.wikipedia.org/wiki/Hydropathy_index),
#cbeta(branching) according to http://www.russell.embl-heidelberg.de/aas/cbb.html,
#a helix breaker is currently only proline, perhaps incorporate something like chou-fassman for helix and strand formers/breakers 
#charge according to side chain charge
#Note: All features are normalized linearly to fit [0,1]
#mass: clear
#vol: clear
#hyd: hydrophobicity
#cbeta: c beta branching
#hbreaker: hbond breaker
#charge: 0->-, 0.5->0, 1->+
aa_properties_normalized = {	'A' : {'mass':0.109,	'vol':0.170,	'hyd':0.700,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 				'R' : {'mass':0.767,	'vol':0.676,	'hyd':0.000,	'cbeta':0,	'hbreaker':0,	'charge':1},	
				 				'N' : {'mass':0.442,	'vol':0.322,	'hyd':0.111,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 			 	'D' : {'mass':0.450,	'vol':0.304,	'hyd':0.111,	'cbeta':0,	'hbreaker':0,	'charge':0},
				 				'C' : {'mass':0.357,	'vol':0.289,	'hyd':0.778,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 				'Q' : {'mass':0.550,	'vol':0.499,	'hyd':0.111,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 				'E' : {'mass':0.558,	'vol':0.467,	'hyd':0.111,	'cbeta':0,	'hbreaker':0,	'charge':0},
				 				'G' : {'mass':0.000,	'vol':0.000,	'hyd':0.456,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 				'H' : {'mass':0.620,	'vol':0.555,	'hyd':0.144,	'cbeta':0,	'hbreaker':0,	'charge':1},
				 				'I' : {'mass':0.434,	'vol':0.636,	'hyd':1.000,	'cbeta':1,	'hbreaker':0,	'charge':0.5},
				 				'L' : {'mass':0.434,	'vol':0.636,	'hyd':0.922,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 				'K' : {'mass':0.550,	'vol':0.647,	'hyd':0.067,	'cbeta':0,	'hbreaker':0, 	'charge':1},
				 				'M' : {'mass':0.574,	'vol':0.613,	'hyd':0.711,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 				'F' : {'mass':0.698,	'vol':0.774,	'hyd':0.811,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 				'P' : {'mass':0.310,	'vol':0.314,	'hyd':0.322,	'cbeta':0,	'hbreaker':1,	'charge':0.5},
				 				'S' : {'mass':0.233,	'vol':0.172,	'hyd':0.411,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 				'T' : {'mass':0.341,	'vol':0.334,	'hyd':0.422,	'cbeta':1,	'hbreaker':0,	'charge':0.5},
				 				'W' : {'mass':1.000,	'vol':1.000,	'hyd':0.400,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 				'Y' : {'mass':0.822,	'vol':0.796,	'hyd':0.356,	'cbeta':0,	'hbreaker':0,	'charge':0.5},
				 				'V' : {'mass':0.326,	'vol':0.476,	'hyd':0.967,	'cbeta':1,	'hbreaker':0,	'charge':0.5},
				 				'X' : {'mass':0,		'vol':0,		'hyd':0,		'cbeta':0,	'hbreaker':0,	'charge':0}
				 			}

aas_standard = 'ARNDCQEGHILKMFPSTWYV'

"""
Global protein features
"""
def feature_aa_composition(seq,fea_win):
	"""
	This calculates a 20 residue feature which holds the global information about the protein's
	amino acid composition.
	"""
	features = []
	featnames = []
	
	aas= {}
	for aa in seq:
		if aa not in aas_standard:
			continue
		try:
			aas[aa] += 1
		except KeyError:
			aas[aa] = 1

	if 'composition' not in fea_win:
		win=0
	else:
		win=fea_win['composition']

	for n in  range(win):
		for aa in aas_standard:
			try:
				val = float(aas[aa])/len(seq)
			except KeyError:
				val = 0
			features.append(val)
			featnames.append('%s_composition'%aa)
	
	return {'composition':featnames},{'composition':features}


def feature_protein_length(seq,fea_win):
	"""
	This essentially gives four features, corresponding to the sequence length:
	1-60
	61-120
	121-180
	181-240
	The specific feature is set to 0.5 if the protein length falls within that category and
	everything below that to 1. 
	"""

	tmp=[]	
	for i in range(1,5):
		tmp.append('length_category%s'%i)


	N = len(seq)
	
	if 'length_category1' in fea_win:
		win1 = fea_win['length_category1']
	else:
		win1 = 0
	if 'length_category2' in fea_win:
		win2 = fea_win['length_category2']
	else:
		win2 = 0
	if 'length_category3' in fea_win:
                win3 = fea_win['length_category3']
	else:
		win3 = 0
	if 'length_category4' in fea_win:
		win4 = fea_win['length_category4']
	else:
		win4 = 0

	
	if N >=1 and N <= 60:
		features1 = [0.5]*win1
		features2 = [0]*win2
		features3 = [0]*win3
		features4 = [0]*win4
	elif N >= 61 and N <= 120:
		features1 = [1]*win1
		features2 = [0.5]*win2
		features3 = [0]*win3
		features4 = [0]*win4
	elif N >= 121 and N <= 180:
		features1 = [1]*win1
		features2 = [1]*win2
		features3 = [0.5]*win3
		features4 = [0]*win4
	elif N >= 181 and N <= 240:
		features1 = [1]*win1
		features2 = [1]*win2
		features3 = [1]*win3
		features4 = [0.5]*win4
	elif N >= 241:
		features1 = [1]*win1
		features2 = [1]*win2
		features3 = [1]*win3
		features4 = [1]*win4
	
	return {'length_category1':[tmp[0]]*win1,'length_category2':[tmp[1]]*win2,'length_category3':[tmp[2]]*win3,'length_category4':[tmp[3]]*win4},{'length_category1':features1,'length_category2':features2,'length_category3':features3,'length_category4':features4}

def feature_sec_strct_composition(sec_strct,fea_win):
	"""
	This accepts a string of secondary structure, composed of letters H,E,L
	and returns a 3x3 feature, each tripel representing a category of 
	the percentage amount of that particular sec strct.
	We do that in a cumulative way, like in feature_protein_length()
	"""
	
	sec = {'H':0,'E':0,'L':0}
	#Count the occurrences
	for s in sec_strct:
		if s not in 'HEL':
			raise InputError('The secondary structure string looks improper to me. Only allowed letters are H,E,L!')
		sec[s] += 1
		
	h = float(sec['H'])/len(sec_strct)
	e = float(sec['E'])/len(sec_strct)
	l = float(sec['L'])/len(sec_strct)
	#Handle helix first

	if 'helix_composition1' in fea_win:
		win1 = fea_win['helix_composition1']
	else:
		win1 = 0
	if 'helix_composition2' in fea_win:
		win2 = fea_win['helix_composition2']
	else:
		win2 = 0
	if 'helix_composition3' in fea_win:
		win3 = fea_win['helix_composition3']
	else:
		win3 = 0

	if h >= 0 and h <= 0.25:
		featuresh1 = [0.5]*win1
		featuresh2 = [0]*win2
		featuresh3 = [0]*win3
	elif h > 0.25 and h <= 0.5:
		featuresh1 = [1]*win1
		featuresh2 = [0.5]*win2
		featuresh3 = [0]*win3
	elif h > 0.5 and h <= 0.75:
		featuresh1 = [1]*win1
		featuresh2 = [1]*win2
		featuresh3 = [0.5]*win3
	elif h > 0.75 and h <= 1:
		featuresh1 = [1]*win1
		featuresh2 = [1]*win2
		featuresh3 = [1]*win3

	name_helix_composition1=['helix_composition1']*win1
	name_helix_composition2=['helix_composition2']*win2
	name_helix_composition3=['helix_composition3']*win3


	#Strand

	if 'strand_composition1' in fea_win:
		win1 = fea_win['strand_composition1']
	else:
		win1 = 0
	if 'strand_composition2' in fea_win:
		win2 = fea_win['strand_composition2']
	else:
		win2 = 0
	if 'strand_composition3' in fea_win:
		win3 = fea_win['strand_composition3']
	else:
		win3 = 0

	if e >= 0 and e <= 0.25:
		featuress1 = [0.5]*win1
		featuress2 = [0]*win2
		featuress3 = [0]*win3
	elif e > 0.25 and e <= 0.5:
		featuress1 = [1]*win1
		featuress2 = [0.5]*win2
		featuress3 = [0]*win3
	elif e > 0.5 and e <= 0.75:
		featuress1 = [1]*win1
		featuress2 = [1]*win2
		featuress3 = [0.5]*win3
	elif e > 0.75 and e <= 1:
		featuress1 = [1]*win1
		featuress2 = [1]*win2
		featuress3 = [1]*win3

	name_strand_composition1=['strand_composition1']*win1
	name_strand_composition2=['strand_composition2']*win2
	name_strand_composition3=['strand_composition3']*win3

	#Loop

	if 'loop_composition1' in fea_win:
		win1 = fea_win['loop_composition1']
	else:
		win1 = 0
	if 'loop_composition2' in fea_win:
		win2 = fea_win['loop_composition2']
	else:
		win2 = 0
	if 'loop_composition3' in fea_win:
		win3 = fea_win['loop_composition3']
	else:
		win3 = 0

	if l >= 0 and l <= 0.25:
		featuresl1 = [0.5]*win1
		featuresl2 = [0]*win2
		featuresl3 = [0]*win3
	elif l > 0.25 and l <= 0.5:
		featuresl1 = [1]*win1
		featuresl2 = [0.5]*win2
		featuresl3 = [0]*win3
	elif l > 0.5 and l <= 0.75:
		featuresl1 = [1]*win1
		featuresl2 = [1]*win2
		featuresl3 = [0.5]*win3
	elif l > 0.75 and l <= 1:
		featuresl1 = [1]*win1
		featuresl2 = [1]*win2
		featuresl3 = [1]*win3
	
	name_loop_composition1=['loop_composition1']*win1
	name_loop_composition2=['loop_composition2']*win2
	name_loop_composition3=['loop_composition3']*win3

	
	return {'helix_composition1':name_helix_composition1,'helix_composition2':name_helix_composition2,'helix_composition3':name_helix_composition3,
		'strand_composition1':name_strand_composition1,'strand_composition2':name_strand_composition2,'strand_composition3':name_strand_composition3,
		'loop_composition1':name_loop_composition1,'loop_composition2':name_loop_composition2,'loop_composition3':name_loop_composition3},\
		{'helix_composition1':featuresh1,'helix_composition2':featuresh2,'helix_composition3':featuresh3,
		'strand_composition1':featuress1,'strand_composition2':featuress2,'strand_composition3':featuress3,
		'loop_composition1':featuresl1,'loop_composition2':featuresl2,'loop_composition3':featuresl3}


def feature_solv_acc_composition(solv_acc,fea_win):
	"""
	This accepts a string of solvent accesibility, composed of letters e,i,b
	and returns a 3x3 feature, each tripel representing a category of 
	the percentage amount of that particular solv acc.
	We do that in a cumulative way, like in feature_protein_length() and feature_sec_strct_composition()
	"""
	
	sec = {'e':0,'i':0,'b':0}
	#Count the occurrences
	for s in solv_acc:
		if s not in 'eib':
			raise InputError('The solvent accessibility string looks improper to me. Only allowed letters are e,i,b!')
		sec[s] += 1
		
	e = float(sec['e'])/len(solv_acc)
	i = float(sec['i'])/len(solv_acc)
	b = float(sec['b'])/len(solv_acc)
	'''
	#Handle exposed first
	if e >= 0 and e <= 0.25:
		features.extend([0.5,0,0])
	elif e > 0.25 and e <= 0.5:
		features.extend([1,0.5,0])
	elif e > 0.5 and e <= 0.75:
		features.extend([1,1,0.5])
	elif e > 0.75 and e <= 1:
		features.extend([1,1,1])
	#intermediate
	if i >= 0 and i <= 0.25:
		features.extend([0.5,0,0])
	elif i > 0.25 and i <= 0.5:
		features.extend([1,0.5,0])
	elif i > 0.5 and i <= 0.75:
		features.extend([1,1,0.5])
	elif i > 0.75 and i <= 1:
		features.extend([1,1,1])
	#buried
	if b >= 0 and b <= 0.25:
		features.extend([0.5,0,0])
	elif b > 0.25 and b <= 0.5:
		features.extend([1,0.5,0])
	elif b > 0.5 and b <= 0.75:
		features.extend([1,1,0.5])
	elif b > 0.75 and b <= 1:
		features.extend([1,1,1])
	'''

	#Handle exposed first

	if 'exposed_composition1' in fea_win:
		win1 = fea_win['exposed_composition1']
	else:
		win1 = 0
	if 'exposed_composition2' in fea_win:
		win2 = fea_win['exposed_composition2']
	else:
		win2 = 0
	if 'exposed_composition3' in fea_win:
		win3 = fea_win['exposed_composition3']
	else:
		win3 = 0

	if e >= 0 and e <= 0.25:
		featurese1 = [0.5]*win1
		featurese2 = [0]*win2
		featurese3 = [0]*win3
	elif e > 0.25 and e <= 0.5:
		featurese1 = [1]*win1
		featurese2 = [0.5]*win2
		featurese3 = [0]*win3
	elif e > 0.5 and e <= 0.75:
		featurese1 = [1]*win1
		featurese2 = [1]*win2
		featurese3 = [0.5]*win3
	elif e > 0.75 and e <= 1:
		featurese1 = [1]*win1
		featurese2 = [1]*win2
		featurese3 = [1]*win3

	name_exposed_composition1=['exposed_composition1']*win1
	name_exposed_composition2=['exposed_composition2']*win2
	name_exposed_composition3=['exposed_composition3']*win3


	#intermediate

	if 'intermediate_composition1' in fea_win:
		win1 = fea_win['intermediate_composition1']
	else:
		win1 = 0
	if 'intermediate_composition2' in fea_win:
		win2 = fea_win['intermediate_composition2']
	else:
		win2 = 0
	if 'intermediate_composition3' in fea_win:
		win3 = fea_win['intermediate_composition3']
	else:
		win3 = 0

	if i >= 0 and i <= 0.25:
		featuresi1 = [0.5]*win1
		featuresi2 = [0]*win2
		featuresi3 = [0]*win3
	elif i > 0.25 and i <= 0.5:
		featuresi1 = [1]*win1
		featuresi2 = [0.5]*win2
		featuresi3 = [0]*win3
	elif i > 0.5 and i <= 0.75:
		featuresi1 = [1]*win1
		featuresi2 = [1]*win2
		featuresi3 = [0.5]*win3
	elif i > 0.75 and i <= 1:
		featuresi1 = [1]*win1
		featuresi2 = [1]*win2
		featuresi3 = [1]*win3

	name_intermediate_composition1=['intermediate_composition1']*win1
	name_intermediate_composition2=['intermediate_composition2']*win2
	name_intermediate_composition3=['intermediate_composition3']*win3

	#buried

	if 'buried_composition1' in fea_win:
		win1 = fea_win['buried_composition1']
	else:
		win1 = 0
	if 'buried_composition2' in fea_win:
		win2 = fea_win['buried_composition2']
	else:
		win2 = 0
	if 'buried_composition3' in fea_win:
		win3 = fea_win['buried_composition3']
	else:
		win3 = 0

	if b >= 0 and b <= 0.25:
		featuresb1 = [0.5]*win1
		featuresb2 = [0]*win2
		featuresb3 = [0]*win3
	elif b > 0.25 and b <= 0.5:
		featuresb1 = [1]*win1
		featuresb2 = [0.5]*win2
		featuresb3 = [0]*win3
	elif b > 0.5 and b <= 0.75:
		featuresb1 = [1]*win1
		featuresb2 = [1]*win2
		featuresb3 = [0.5]*win3
	elif b > 0.75 and b <= 1:
		featuresb1 = [1]*win1
		featuresb2 = [1]*win2
		featuresb3 = [1]*win3
	
	name_buried_composition1=['buried_composition1']*win1
	name_buried_composition2=['buried_composition2']*win2
	name_buried_composition3=['buried_composition3']*win3


	
	return {'exposed_composition1':name_exposed_composition1,'exposed_composition2':name_exposed_composition2,'exposed_composition3':name_exposed_composition3,
                'intermediate_composition1':name_intermediate_composition1,'intermediate_composition2':name_intermediate_composition2,'intermediate_composition3':name_intermediate_composition3,
                'buried_composition1':name_buried_composition1,'buried_composition2':name_buried_composition2,'buried_composition3':name_buried_composition3},\
                {'exposed_composition1':featurese1,'exposed_composition2':featurese2,'exposed_composition3':featurese3,
                'intermediate_composition1':featuresi1,'intermediate_composition2':featuresi2,'intermediate_composition3':featuresi3,
                'buried_composition1':featuresb1,'buried_composition2':featuresb2,'buried_composition3':featuresb3}

"""
Local, i.e. window-based protein features
"""
def feature_chem_props(seq, central_pos, fea_win):
	"""
	This is a feature that is not derived from pp, but might be useful for structure-oriented predictions.
	"""
	N = len(seq)
	#The lists that collect the features per position

	window_={'chemprop_mass':[],'chemprop_vol':[],'chemprop_hyd':[],'chemprop_cbeta':[],'chemprop_hbreaker':[],'chemprop_charge':[],'chemprop_mass':[],'position':[]}

	#The feature names

	featnames_={'chemprop_mass':[],'chemprop_vol':[],'chemprop_hyd':[],'chemprop_cbeta':[],'chemprop_hbreaker':[],'chemprop_charge':[],'chemprop_mass':[],'position':[]}

	


	for f in window_:
		if f in fea_win:
			window_length=fea_win[f]
		else:
			continue 

		if window_length % 2 == 0:
                	raise InputError("To me, it makes more sense to use uneven window lengths.")

		start_pos = central_pos - (window_length-1)/2
		stop_pos = central_pos + (window_length-1)/2
		cnt = int(-window_length/2)	#For designating the position prefix in the feature names
		for i in range(int(start_pos), int(stop_pos)+1):
			#We hit the left or right boundary of the sequence with our window...
			if i < 0 or i >= len(seq):
				window_[f].append(0)
			else:
				props = aa_properties_normalized[seq[i]]
				if f=='position':
					window_[f].append((float(i)+1)/N)
				else:
					window_[f].append(props[f.split('_')[1]])

			featnames_[f].append("%s_"%cnt+f)

		
			cnt += 1


	return featnames_,window_

'''		
def feature_sequence(seq, central_pos, window_length):
	"""
	"""
	if window_length % 2 == 0:
		raise InputError("To me, it makes more sense to use uneven window lengths.")
	
	raw_featnames = 'ARNDCQEGHILKMFPSTWYV'
	#The list that collects the features per position
	window_seq = []
	#The feature names
	featnames_seq = []
	#Domains
	domains_seq = []
	
	start_pos = central_pos - (window_length-1)/2
	stop_pos = central_pos + (window_length-1)/2
	cnt = int(-window_length/2)	#For designating the position prefix in the feature names
	for i in range(int(start_pos), int(stop_pos)+1):
		#We hit the left or right boundary of the sequence with our window...
		if i < 0 or i >= len(seq):
			aa_seq = 20*[0]
		else:
			aa = seq[i]
			if aa not in raw_featnames:	#Like it's the case with an 'X' 
				aa_seq = 20*[0]	
			else:
				aa_seq = 20*[0]
				pos = raw_featnames.index(aa)
				aa_seq[pos] = 1
		
		featnames_seq.append( ["%s_%s_SEQ"%(cnt,column) for column in 'ARNDCQEGHILKMFPSTWYV' ] )
		domains_seq.append( ["{0,1}" for column in 'ARNDCQEGHILKMFPSTWYV' ] )
		
		window_seq.append(aa_seq)
		cnt += 1
	
	return {'names_seq':featnames_seq}, {'win_seq':window_seq}, {'domains_seq':domains_seq}
'''

def feature_blast(dict_blast, central_pos, fea_win, norm_raw=True):
	"""
	"""
	#Get the input parts
	pssm = dict_blast['pssm']
	perc = dict_blast['perc']
	inf_per_pos = dict_blast['inf_per_pos']
	rel_weight = dict_blast['rel_weight']

	#The lists that collect the features per position
	window_={'pssm':[],'perc':[],'infPP':[],'relW':[]}

	#The feature names
	featnames_={'pssm':[],'perc':[],'infPP':[],'relW':[]}



	for f in window_:
		if f in fea_win:
			window_length=fea_win[f]
		else:
			continue

		if window_length % 2 == 0:
			raise InputError("To me, it makes more sense to use uneven window lengths.")

		start_pos = central_pos - (window_length-1)/2
		stop_pos = central_pos + (window_length-1)/2
		cnt = int(-window_length/2)     #For designating the position prefix in the feature names

	
	#cnt = int(-window_length/2)	#For designating the position prefix in the feature names
		for i in range(int(start_pos), int(stop_pos)+1):
			#We hit the left or right boundary of the sequence with our window...
			if i < 0 or i >= len(pssm):
				if f=='pssm':
					aa_pssm = 20*[0]		#The 20 standard aas
				elif f=='perc':
					aa_perc = 20*[0]		#The 20 standard aas
				elif f=='infPP':
					aa_inf_per_pos = 0		#As in SNAP
				elif f=='relW':
					aa_rel_weight = 0		#As in SNAP
			#...and completely within
			else:
				if f=='pssm':
					aa_pssm = pssm[i]		#The 20 standard aas
				elif f=='perc':
					aa_perc = perc[i]
				elif f=='infPP':
					aa_inf_per_pos = inf_per_pos[i]
				elif f=='relW':
					aa_rel_weight = rel_weight[i]
		
			#Normalize 
			if norm_raw:
				if f=='pssm': 
					data = [ float(1)/(1+math.exp(-int(t))) for t in aa_pssm ]	#Logistically, according to Jones' PSI-Pred paper
				elif f=='perc':
					data = [ float(t)/100 for t in aa_perc ]	#Linearly
				elif f=='infPP':
					data = float(1)/(1+math.exp(-float(aa_inf_per_pos)))	#Logistically, since we do not know the max value needed for linearly normalization
				elif f=='relW':
					data = float(1)/(1+math.exp(-float(aa_rel_weight)))	#Logistically, since we do not know the max value needed for linearly normalization
				
			#Append the current residue's feautures to the feature spaces
			window_[f].append(data)
		
			#Finally handle the feature names too
			if f=='pssm':
				featnames_['pssm'].append( ["%s_%s_pssm"%(cnt,column) for column in 'ARNDCQEGHILKMFPSTWYV' ] )
			elif f=='perc':
				featnames_['perc'].append( ["%s_%s_perc"%(cnt,column) for column in 'ARNDCQEGHILKMFPSTWYV' ] )
			elif f=='infPP':
				featnames_['infPP'].append( "%s_infPP"%cnt )
			elif f=='relW':
				featnames_['relW'].append( "%s_relW"%cnt )
			
		
			cnt += 1
	
	#We return a 2-tuple, each containing a dictionary:
	#1. The feature names
	#2. The actual feature values, window-position dependent
	#3. The domains of each feature, i.e. numeric, nominal, string
	return featnames_,window_

'''
def feature_disis(dict_disis, central_pos, window_length, norm_raw=True):
	"""
	"""
	if window_length % 2 == 0:
		raise InputError("To me, it makes more sense to use uneven window lengths.")
	#Get the input parts
	bin = dict_disis['prd_bin']
	raw = dict_disis['prd_raw']
	#The lists that collect the features per position
	#We introduce two binary features here, one for a positive and one for negative prdct. A 1 in the features list tells which one it is
	window_bin_plus = []	
	window_bin_minus = []	
	window_raw = []
	#The feature names
	featnames_bin_plus = []
	featnames_bin_minus = []
	featnames_raw = []
	#The domains
	domains_bin_plus = []
	domains_bin_minus = []
	domains_raw = []
	
	start_pos = central_pos - (window_length-1)/2
	stop_pos = central_pos + (window_length-1)/2
	cnt = int(-window_length/2)	#For designating the position prefix in the feature names
	for i in range(int(start_pos), int(stop_pos)+1):
		#We hit the left or right boundary of the sequence with our window...
		if i < 0 or i >= len(bin):
			window_bin_plus.append(0)
			window_bin_minus.append(0)
			window_raw.append(0)
		else:
			if bin[i] == '+':
				window_bin_plus.append(1)
				window_bin_minus.append(0)
			elif bin[i] == '-':
				window_bin_plus.append(0)
				window_bin_minus.append(1)
			window_raw.append(raw[i])
		
		featnames_bin_plus.append( "%s_disis_plus"%cnt )
		featnames_bin_minus.append( "%s_disis_minus"%cnt )
		featnames_raw.append( "%s_disis_raw"%cnt )
		
		domains_bin_minus.append( "{0,1}" )
		domains_bin_plus.append( "{0,1}" )
		domains_raw.append( "NUMERIC" )
		
		cnt += 1
	
	#Normalize if wanted
	if norm_raw: window_raw = [ (float(t)+100)/200 for t in window_raw]
	
	return {'names_bin_plus':featnames_bin_plus, 'names_bin_minus':featnames_bin_minus, 'names_raw':featnames_raw}, {'win_bin_plus':window_bin_plus, 'win_bin_minus':window_bin_minus, 'win_raw':window_raw}, {'domains_bin_plus':domains_bin_plus, 'domains_bin_minus':domains_bin_minus, 'domains_raw':domains_raw}
'''	

def feature_md(dict_md, central_pos, fea_win, norm_raw=True):
	"""
	"""
	#Get the input parts
	#For now, we skip everything else that is in the parsed md output (i.e. the other helper-predictions).
	#Those shouldn't be of great interest anyway.
	bin = dict_md['prd_bin']
	raw = dict_md['prd_raw']
	ri = dict_md['prd_ri']
	#The lists that collect the features per position
	#We introduce two binary features here, one for a positive and one for negative prdct. A 1 in the features list tells which one it is
	window_ = {'md_plus':[],'md_minus':[],'md_raw':[],'md_ri':[]}
	#The feature names
	featnames_ = {'md_plus':[],'md_minus':[],'md_raw':[],'md_ri':[]}
	#Domains


	for f in window_:
		if f in fea_win:
			window_length=fea_win[f]
		else:
			continue

		if window_length % 2 == 0:
			raise InputError("To me, it makes more sense to use uneven window lengths.")


		start_pos = central_pos - (window_length-1)/2
		stop_pos = central_pos + (window_length-1)/2
		cnt = int(-window_length/2)	#For designating the position prefix in the feature names
		for i in range(int(start_pos), int(stop_pos)+1):
			#We hit the left or right boundary of the sequence with our window...
			if i < 0 or i >= len(bin):
				window_[f].append(0)
			else:
				if bin[i] == '+':
					if f=='md_plus':
						window_[f].append(1)
					elif f=='md_minus':
						window_[f].append(0)
				elif bin[i] == '-':
					if f=='md_plus':
						window_[f].append(0)
					elif f=='md_minus':
						window_[f].append(1)
				if f=='md_raw':
					window_[f].append(raw[i])
				elif f=='md_ri':
					window_[f].append(ri[i])
	
			if f=='md_plus':	
				featnames_[f].append( "%s_md_plus"%cnt )
			elif f=='md_minus':
				featnames_[f].append( "%s_md_minus"%cnt )
			elif f=='md_raw':
				featnames_[f].append( "%s_md_raw"%cnt )
			elif f=='md_ri':
				featnames_[f].append( "%s_md_ri"%cnt  )
		
			
			cnt += 1
	
	#Normalize if wanted
	#The raw output is already between 0 and 1 hence we norm only the ri here
	if norm_raw: window_['md_ri'] = [ float(t)/9 for t in window_['md_ri']]
	
	return featnames_,window_	

'''	
def feature_isis(dict_isis, central_pos, window_length, norm_raw=True):
	"""
	"""
	if window_length % 2 == 0:
		raise InputError("To me, it makes more sense to use uneven window lengths.")
	#Get the input parts
	bin = dict_isis['prd_bin']
	raw = dict_isis['prd_raw']
	#The lists that collect the features per position
	#We introduce two binary features here, one for a positive and one for negative prdct. A 1 in the features list tells which one it is
	window_bin_plus = []	
	window_bin_minus = []	
	window_raw = []
	#The feature names
	featnames_bin_plus = []
	featnames_bin_minus = []
	featnames_raw = []
	#Domains
	domains_bin_plus = []
	domains_bin_minus = []
	domains_raw = []
	
	start_pos = central_pos - (window_length-1)/2
	stop_pos = central_pos + (window_length-1)/2
	cnt = int(-window_length/2)	#For designating the position prefix in the feature names
	for i in range(int(start_pos), int(stop_pos)+1):
		#We hit the left or right boundary of the sequence with our window...
		if i < 0 or i >= len(bin):
			window_bin_plus.append(0)
			window_bin_minus.append(0)
			window_raw.append(0)
		else:
			if bin[i] == '+':
				window_bin_plus.append(1)
				window_bin_minus.append(0)
			elif bin[i] == '-':
				window_bin_plus.append(0)
				window_bin_minus.append(1)
			window_raw.append(raw[i])
		
		featnames_bin_plus.append( "%s_isis_plus"%cnt )
		featnames_bin_minus.append( "%s_isis_minus"%cnt )
		featnames_raw.append( "%s_isis_raw"%cnt )
		
		domains_bin_plus.append("{0,1}")
		domains_bin_minus.append("{0,1}")
		domains_raw.append("NUMERIC")
		
		cnt += 1
	
	#Normalize if wanted
	if norm_raw: window_raw = [ (float(t)+100)/200 for t in window_raw]
	
	return {'names_bin_plus':featnames_bin_plus, 'names_bin_minus':featnames_bin_minus, 'names_raw':featnames_raw}, {'win_bin_plus':window_bin_plus, 'win_bin_minus':window_bin_minus, 'win_raw':window_raw}, {'domains_bin_plus':domains_bin_plus, 'domains_bin_minus':domains_bin_minus, 'domains_raw':domains_raw}
'''	

def feature_profbval(dict_bval, central_pos, fea_win, norm_raw=True):
	"""
	"""
	#Get the input parts
	raw1 = dict_bval['prd_raw1']
	raw2 = dict_bval['prd_raw2']
	data = {'profbval_raw1':raw1,'profbval_raw2':raw2}
	#The lists that collect the features per position
	window_ = {'profbval_raw1':[],'profbval_raw2':[]}
	#The feature names
	featnames_ = {'profbval_raw1':[],'profbval_raw2':[]}	


	for f in window_:
		if f in fea_win:
			window_length=fea_win[f]
		else:
			continue

		if window_length % 2 == 0:
			raise InputError("To me, it makes more sense to use uneven window lengths.")



		start_pos = central_pos - (window_length-1)/2
		stop_pos = central_pos + (window_length-1)/2
		cnt = int(-window_length/2)	#For designating the position prefix in the feature names
		for i in range(int(start_pos), int(stop_pos)+1):
			#We hit the left or right boundary of the sequence with our window...
			if i < 0 or i >= len(raw1):
				window_[f].append(0)
			else:
				window_[f].append(data[f][i])
		
			featnames_[f].append( "%s_"%cnt+f )
		
		
			cnt += 1
	
		#Normalize if wanted
		if norm_raw: 
			window_[f] = [ float(t)/100 for t in window_[f]]
		
	return featnames_,window_

def feature_profsecacc(dict_profsecacc, central_pos, fea_win, norm_raw=True):
	"""
	"""

	#Get the input parts
	#We don't use everything here, only what makes sense and was used in the past in the group, e.g. in SNAP
	sec = dict_profsecacc['PHEL']		#The 3-states secondary structure, H,E,L
	sec_OtH = dict_profsecacc['OtH']	#Raw NN output for Helix
	sec_OtE = dict_profsecacc['OtE']	#Raw NN output for Strand
	sec_OtL = dict_profsecacc['OtL']	#Raw NN output for Loop
	sec_ri = dict_profsecacc['RI_S']	#Reliability index for secondary structure prediction, 0-9
	
	acc = dict_profsecacc['Pbie']		#The 3-states solvent accessibility, e,b,i
	acc_A = dict_profsecacc['PREL']		#Relative solvent accessibility, 0-100
	acc_ri = dict_profsecacc['RI_A']	#Reliability index for solvent accessibility prediction, 0-9
	
	#The lists that collect the features per position

	window_ = {'e':[],'i':[],'b':[],'helix':[],'strand':[],'loop':[],'OtH':[],'OtE':[],'OtL':[],'ri_sec':[],'ri_acc':[],'rel_acc':[]}
	
	#The feature names
	
	featnames_ = {'e':[],'i':[],'b':[],'helix':[],'strand':[],'loop':[],'OtH':[],'OtE':[],'OtL':[],'ri_sec':[],'ri_acc':[],'rel_acc':[]}



	for f in window_:
		if f in fea_win:
			window_length=fea_win[f]
		else:
			continue

		if window_length % 2 == 0:
			raise InputError("To me, it makes more sense to use uneven window lengths.")


	
		start_pos = central_pos - (window_length-1)/2
		stop_pos = central_pos + (window_length-1)/2
		cnt = int(-window_length/2)	#For designating the position prefix in the feature names
	
		for i in range(int(start_pos), int(stop_pos)+1):
			#We hit the left or right boundary of the sequence with our window...
			if i < 0 or i >= len(sec):
				window_[f].append(0)
			else:
				if sec[i] == 'H':
					if f=='helix':
						window_['helix'].append(1)
					if f=='strand':
						window_['strand'].append(0)
					if f=='loop':
						window_['loop'].append(0)
				elif sec[i] == 'E':
					if f=='helix':
						window_['helix'].append(0)
					if f=='strand':
						window_['strand'].append(1)
					if f=='loop':
						window_['loop'].append(0)
				elif sec[i] == 'L':
					if f=='helix':
						window_['helix'].append(0)
					if f=='strand':
						window_['strand'].append(0)
					if f=='loop':
						window_['loop'].append(1)
			
				if f=='OtH':
					window_['OtH'].append( sec_OtH[i] )
				elif f=='OtE':
					window_['OtE'].append( sec_OtE[i] )
				elif f=='OtL':
					window_['OtL'].append( sec_OtL[i] )
				elif f=='ri_sec':
					window_['ri_sec'].append( sec_ri[i] )
			
				if acc[i] == 'e':
					if f=='e':
						window_['e'].append(1)
					if f=='i':
						window_['i'].append(0)
					if f=='b':
						window_['b'].append(0)
				elif acc[i] == 'i':
					if f=='e':
						window_['e'].append(0)
					if f=='i':
						window_['i'].append(1)
					if f=='b':
						window_['b'].append(0)
				elif acc[i] == 'b':
					if f=='e':
						window_['e'].append(0)
					if f=='i':
						window_['i'].append(0)
					if f=='b':
						window_['b'].append(1)
				if f=='rel_acc':
					window_['rel_acc'].append(acc_A[i])
				elif f=='ri_acc':
					window_['ri_acc'].append(acc_ri[i])

			featnames_[f].append("%s_"%cnt+f)		
		
			cnt += 1

		#Normalize if wanted
		if norm_raw:
			if f=='OtH': 
				window_['OtH'] = [ float(t)/100 for t in window_['OtH']]
			if f=='OtE':
				window_['OtE'] = [ float(t)/100 for t in window_['OtE']]
			if f=='OtL':
				window_['OtL'] = [ float(t)/100 for t in window_['OtL']]
			if f=='ri_sec':
				window_['ri_sec'] = [ float(t)/9 for t in window_['ri_sec'] ]
			if f=='rel_acc':
				window_['rel_acc'] = [ float(t)/100 for t in window_['rel_acc'] ]
			if f=='ri_acc':
				window_['ri_acc'] = [ float(t)/9 for t in window_['ri_acc'] ]
	
	return featnames_,window_

def feature_psic(dict_psic, central_pos, window_length, norm_raw=True):
	"""
	"""
	if window_length % 2 == 0:
		raise InputError("To me, it makes more sense to use uneven window lengths.")
	
	#Get the input parts
	psic = dict_psic['psic']
	numSeq = dict_psic['NumSeq']
	#The lists that collect the features per position
	window_psic = []
	window_numSeq = []
	#The feature names
	featnames_psic = []
	featnames_numSeq = []
	#Domains
	domains_psic = []
	domains_numSeq = []
	
	start_pos = central_pos - (window_length-1)/2
	stop_pos = central_pos + (window_length-1)/2
	cnt = int(-window_length/2)	#For designating the position prefix in the feature names
	for i in range(int(start_pos), int(stop_pos)+1):
		#We hit the left or right boundary of the sequence with our window
		if i < 0 or i >= len(psic):
			aa_psic = 20*[0]
			window_numSeq.append(0)
		else :
			aa_psic = psic[i]		#The 20 standard aas, not represented here
			window_numSeq.append(numSeq[i])

		#Normalize according to Jones' PSI-Pred paper
		if norm_raw: 
			aa_psic = [ float(1)/(1+math.exp(-float(t))) for t in aa_psic ]
		
		window_psic.append(aa_psic)
		featnames_psic.append(["%s_%s_psic"%(cnt,column) for column in 'ARNDCQEGHILKMFPSTWYV' ])
		featnames_numSeq.append("%s_psic_numSeq"%cnt)
		
		domains_psic.append( ["NUMERIC" for column in 'ARNDCQEGHILKMFPSTWYV' ] )
		domains_numSeq.append( "NUMERIC" )
		
		cnt += 1
		
	if norm_raw: 
		window_numSeq = [ float(t)/300 for t in window_numSeq ] #The max value is 300 due to the restriction in the blast call that collects the sequences (s. SNAP paper)
		
	return {'names_psic':featnames_psic, 'names_numSeq':featnames_numSeq}, {'win_psic':window_psic, 'win_numSeq':window_numSeq}, {'domains_psic':domains_psic, 'domains_numSeq':domains_numSeq}


def feature_pfam(dict_pfam, central_pos, window_length):
	"""
	"""
	if window_length % 2 == 0:
		raise InputError("To me, it makes more sense to use uneven window lengths.")
	
	try:
		residues_annotated = set(dict_pfam.keys())
	#In that case, dict_pfam is None, hence _no_ residue is annotated.
	except AttributeError:
		residues_annotated = set()
		
	#The lists that collect the features per position
	window_within_domain = []		#Is the specific residue within a domain?
	window_domain_conservation = []	#How conserved is the domain?
	window_domain_fit = []			#How well does the current residue fit to the model's position?
	window_pp = []					#The posterior probability for the match
	#The lists that collect the features per position
	featnames_within_domain = []
	featnames_domain_conservation = []
	featnames_domain_fit = []
	featnames_pp = []
	#Feature domains
	domains_within_domain = []
	domains_domain_conservation = []
	domains_domain_fit = []
	domains_pp = []
	
	start_pos = central_pos - (window_length-1)/2
	stop_pos = central_pos + (window_length-1)/2
	cnt = int(-window_length/2)	#For designating the position prefix in the feature names
	for i in range(int(start_pos), int(stop_pos)+1):
		#Is the current residue within a domain?
		if i in residues_annotated:
			window_within_domain.append(1)
		else:
			window_within_domain.append(0)
		#What about everything else?
		try:
			res_domains = dict_pfam[i].values()
			min_eval = sys.maxint	#We capture only the highest scoring domain, i.e. with the lowest evalue
			for mdl_doms in res_domains:
				for dom in mdl_doms:
					eval = float(dom[1])	#The domain's evalue
					conservation = dom[2]	#The conservation of the current position, indicated by upper,lower letter or '.' (insertion)
					residue_fit = dom[3]	#How does the current residue fit to the current domain position?
					pp = dom[4]				#The posterior probability for that match
					if eval < min_eval:
						min_eval = eval
						#Now handle the other features
						#1. Domain conservation
						if conservation == '.': cons = 0		#Residue is not aligned to a domain, i.e. insertion
						elif conservation.isupper(): cons = 1	#Domain is highly conserved at that position
						elif conservation.islower(): cons = 0.5	#Not so well conserved 
						else: raise ValueError('There is another possibility')
						#2. Residue fit to the current domain
						if residue_fit == ' ': fit = 0.3	#No fit at all
						elif residue_fit == '+': fit = 0.6	#some fit
						else: fit = 1						#Best fit, could be lower or upper letter in that line, we do not distinguish between both
						#3. Finally the posterior probability of the fit
						if pp == '*': prob = 1
						else: prob = int(pp)/11.0
		
		#The current residue isn't aligned to a domain at all	
		except KeyError:
			cons = 0
			fit = 0
			prob = 0
		#No domains at all were found for the current sequence.
		#In that case dict_pfam is None
		except TypeError:
			cons = 0
			fit = 0
			prob = 0
		
		window_domain_conservation.append(cons)
		window_domain_fit.append(fit)
		window_pp.append(prob)
		
		#Handle the other lists as well
		featnames_within_domain.append("%s_pfam_within_domain"%cnt)
		featnames_domain_conservation.append("%s_pfam_domain_conservation"%cnt)
		featnames_domain_fit.append("%s_pfam_residue_fit"%cnt)
		featnames_pp.append("%s_pfam_pp"%cnt)
		
		domains_within_domain.append("{0,1}")
		domains_domain_conservation.append("{0,0.5,1}")
		domains_domain_fit.append("{0,0.3,0.6,1}")
		domains_pp.append("NUMERIC")
		
		#Continue with next residue within the window
		cnt += 1

	return {
			'names_within_domain':featnames_within_domain,
			'names_domain_conservation':featnames_domain_conservation,
			'names_residue_fit':featnames_domain_fit,
			'names_pp':featnames_pp
			}, {
			'win_within_domain':window_within_domain,
			'win_domain_conservation':window_domain_conservation,
			'win_residue_fit':window_domain_fit,
			'win_pp':window_pp
			}, {
			'domains_within_domain':domains_within_domain,
			'domains_domain_conservation':domains_domain_conservation,
			'domains_residue_fit':domains_domain_fit,
			'domains_pp':domains_pp
			}
	

def feature_prosite(list_prosite, central_pos, window_length):
	"""
	"""
	if window_length % 2 == 0:
		raise InputError("To me, it makes more sense to use uneven window lengths.")
	
	#We here simply note whether the current amino acid is part of a prosite pattern or not
	window_prosite = []		
	featnames_prosite = []
	domains_prosite = []
	
	start_pos = central_pos - (window_length-1)/2
	stop_pos = central_pos + (window_length-1)/2
	cnt = int(-window_length/2)	#For designating the position prefix in the feature names
	for i in range(int(start_pos), int(stop_pos)+1):
		part_of_prosite = False
		for prosite in list_prosite:
			#Note that we _implicitly_handle the left and right boundaries of the sequence here:
			#In that case the following clause never evaluates to True, hence those off-sequence-positions
			#are annotated as 0!
			if i >= prosite[0] and i <= prosite[1]:
				part_of_prosite = True
		
		if part_of_prosite:
			window_prosite.append(1)
		else:
			window_prosite.append(0)
		
		featnames_prosite.append("%s_prosite"%cnt)
		domains_prosite.append("{0,1}")
		
		cnt += 1
	
	return {'names_prosite':featnames_prosite},{'win_prosite':window_prosite},{'domains_prosite':domains_prosite}
	

if __name__ == '__main__':
	from lib_parser import *
	import os
	pp_path = '/mnt/home/schaefer/SNAPv2/pp/'
	chains = os.listdir(pp_path)
	for chain in chains:
		#d_hmmer = open(pp_path+chain+"/query.hmm3pfam").read()
		#print 200*'='
		#print d_hmmer
		d_in = open(pp_path+chain+"/query.in").read()
		d_fasta = open(pp_path+chain+"/query.fasta").read()
		d_prosite = open(pp_path+chain+"/query.prosite").read()
		
		seq = parse_sequence(d_in, d_fasta)['seq']
		l_prosite = parse_prosite(d_prosite)
		#dict_pfam = parse_pfam_annotations(d_hmmer)
		#print feature_pfam(dict_pfam, 100, 21)[1]
		#print feature_chem_props(seq, 100, 21)[1]
		print (feature_prosite(l_prosite, 100,11)[1])
		print (200*'+')
