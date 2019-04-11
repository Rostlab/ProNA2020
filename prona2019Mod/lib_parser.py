import math
import sys
import re


class ParseError(Exception): pass		#To indicate a parsing error
class EmptyError(Exception): pass		#To indicate empty data passed to the parser
class NoResultError(Exception): pass	#To indicate that a method did not feel like producing a result, used in parse_psic()


def parse_sequence(d_in, d_fasta):
	"""
	pp returns two sequence files: query.in and query.fasta. No idea why.
	Here we check that both are the same and return the sequence.
	"""
	seq_in = ''
	seq_fasta = ''
	for line in d_in.split('\n')[1:]:
		if not line: continue
		seq_in += line
	for line in d_fasta.split('\n')[1:]:
		if not line: continue
		seq_fasta += ''.join(line.split())	#Get rid of those strange whitespaces within the sequence!
	
	if seq_in != seq_fasta:
		sys.exit("Error!!!ProNA2019 can not be done for protein %s.\nProtein sequence of *in and * fasta are not identical.\npp seems to work with different sequences.\n" % d_fasta.split('\n')[0][1:])
	
	return {'seq':seq_in},d_fasta.split('\n')[0][1:]
		

def parse_blast(d_blast):
	""" 
	Note that we do not parse out the weighted observed percentages part.
	Meaning of pssm columns:
	A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
	
	Returns a dictionary with keys as follows:
	'seq':	The sequence as blast sees it
	'pssm':	pssm matrix as a list of lists. Each sublist represents a row in the PSSM matrix.
	'perc': perc matrix
	'inf_per_pos': The second last column in the blast output
	'rel_weight': The last column
	"""
	if d_blast == '':
		raise EmptyError('Empty pssm file!')
	
	pssm_mat = []
	perc_mat = []
	inf_per_pos = []
	rel_weight = []
	pssm_seq = ''
	#First turn pssm into a matrix we can handle
	for line in d_blast.split('\n'):
		tokens = line.split()
		if len(tokens) == 40 and line.strip() != 'A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V':
			raise ParseError("It seems that we have an issue now. Blast produces columns with altering meanings!")
		
		if len(tokens) != 44: continue
		
		pssm_seq += tokens[1] 
		inf_per_pos.append( float(tokens[42]) )	#The second last column in the blast output
		rel_weight.append( float(tokens[43]) )	#The very last column
		
		#The first matrix i.e. pssm
		pssm_mat_row = []
		for t in tokens[2:22]:
			pssm_mat_row.append(int(t))
		#The second one, i.e. the percentages
		perc_mat_row = []
		for t in tokens[22:42]:
			perc_mat_row.append(int(t))

		#Check if we are really dealing with 20 values here!
		if len(pssm_mat_row) != 20 or len(perc_mat_row) != 20:
			raise ParseError("It seems that we have a situation now. The expected amount of columns is 20, found: %s!" % len(pssm_mat_row))
		
		pssm_mat.append(pssm_mat_row)
		perc_mat.append(perc_mat_row)
	
	#Further consistency check...
	if len(pssm_mat) != len(pssm_seq) != len(perc_mat) != len(inf_per_pos) != len(rel_weight):
		raise ParseError("It seems that we have an issue now. Something went wrong during parsing the pssm matrix!")

	return {'seq':pssm_seq, 'pssm':pssm_mat, 'perc':perc_mat, 'inf_per_pos':inf_per_pos, 'rel_weight':rel_weight}


def parse_psic(d_psic):
	"""
	Unfortunately, psic returns no sequence. 
	Meaning of psic's columns:
	A      R      N      D      C      Q      E      G      H      I      L      K      M      F      P      S      T      W      Y      V NumSeq
	This is exactly what could be found in the sublist of each residue.
	
	Returns a dictionary with keys as follows:
	'psic':	psic matrix as a list of lists. Each sublist represents a row in the psic matrix.
	'NumSeq': the very last column, denoting NumSeq i.e. number of aligned sequences at that pos 
	"""
	if d_psic == '':
		raise EmptyError('Empty psic file!')
	elif d_psic.startswith('sequence too short'):
		raise NoResultError('Sequence seems to be too short for psic. No psic output found.')
		
	psic_mat = []
	numseq = []
	for line in d_psic.split('\n'):
		if line.startswith('Pos') or line == '': continue
		tokens = line.split()
		if len(tokens) != 22:
			raise ParseError('"It seems that we have a situation now. The expected amount of columns is 22, found: %s!" % len(tokens)')
		psic_mat_row = [ float(t) for t in tokens[1:21] ]
		numseq.append( int(tokens[21]) )	#The last column is an integer denoting the amount of aligned seqs at that pos.		
		
		psic_mat.append(psic_mat_row)		#Now glue the current column to the matrix
	
	#Check!
	if len(psic_mat) != len(numseq):
		raise ParseError("It seems that we have an issue now. Something went wrong during parsing the psic matrix!")
	
	return {'psic':psic_mat, 'NumSeq':numseq}
	

def parse_disis(d_disis):
	"""
	Returns a dictionary with keys as follows:
	'seq':		The sequence as disis sees it
	'prd_bin':	binary prdct
	'prd_raw':	raw prdct
	"""
	if d_disis == '':
		raise EmptyError('Empty disis file!')
	
	disis_seq_binprd = []	#Sequence parsed out of the binary prediction part
	disis_seq_rawprd = []	#...parsed out of the raw (numeric) predictions
	disis_prd_bin = []		#Binary predictions
	disis_prd_raw = []		#Raw numeric predictions
	
	cnt = 0
	for line in d_disis.split('\n'):
		if line == '': continue
		tokens = line.split()
		if len(tokens) == 1:	#We are in the upper part of disis' output, i.e. the binary predictions
			if cnt % 2 == 0:
				disis_seq_binprd.extend( list(line) )
			else:
				disis_prd_bin.extend( list(line.replace('P','+')) )
		elif len(tokens) == 2:	#Now we are in the lower part, i.e. the numeric outputs of disis
			disis_seq_rawprd.append( tokens[0] )
			disis_prd_raw.append( int(tokens[1]) )
			
		cnt += 1
	
	#Now do some consistency checks
	if disis_seq_binprd != disis_seq_rawprd:
		raise ParseError("It seems that we have an issue now. Disis returns different sequences in the upper and lower part!")
	if len(disis_seq_binprd) != len(disis_prd_bin) != len(disis_prd_raw):
		raise ParseError("It seems that we have an issue now. Parsed datastructures have different lengths!")
		
	return {'seq':''.join(disis_seq_binprd), 'prd_bin':disis_prd_bin, 'prd_raw':disis_prd_raw}


def parse_isis(d_isis):
	"""
	Returns a dictionary with keys as follows:
	'seq':		The sequence as isis sees it
	'prd_bin':	binary prdct
	'prd_raw':	raw prdct
	"""
	if d_isis == '':
		raise EmptyError('Empty isis file!')
	
	isis_seq_binprd = []	#Sequence parsed out of the binary prediction part
	isis_seq_rawprd = []	#...parsed out of the raw (numeric) predictions
	isis_prd_bin = []		#Binary predictions
	isis_prd_raw = []		#Raw numeric predictions
	
	cnt = 0
	for line in d_isis.split('\n'):
		if line == '' or line.startswith('>'): continue
		tokens = line.split()
		if len(tokens) == 1:	#We are in the upper part of disis' output, i.e. the binary predictions
			if cnt % 2 == 0:
				isis_seq_binprd.extend( list(line) )
			else:
				isis_prd_bin.extend( list(line.replace('P','+')) )
		elif len(tokens) == 3:	#Now we are in the lower part, i.e. the numeric outputs of disis
			isis_seq_rawprd.append( tokens[1] )
			isis_prd_raw.append( int(tokens[2]) )
			
		cnt += 1
	
	#Now do some consistency checks
	if isis_seq_binprd != isis_seq_rawprd:
		raise ParseError("It seems that we have an issue now. Isis returns different sequences in the upper and lower part!")
	if len(isis_seq_binprd) != len(isis_prd_bin) != len(isis_prd_raw):
		raise ParseError("It seems that we have an issue now. Parsed datastructures have different lengths!")
	
	return {'seq':''.join(isis_seq_binprd), 'prd_bin':isis_prd_bin, 'prd_raw':isis_prd_raw}	


def parse_md(d_md):
	"""
	Returns a dictionary with keys as follows:
	'seq': 			sequence as MD sees it
	'norsnet_raw':	raw norsnet prdct
	'norsnet_bin':	binary norsnet prdct
	'bval_raw':		raw bval prdct
	'bval_bin':		binary bval prdct 
	'ucon_raw':		raw ucon prdct
	'ucon_bin':		binary ucon prdct
	'prd_raw':		MD's raw prdct
	'prd_ri':		MD's reliability index 
	'prd_bin':		MD's binary prdct
	"""
	if d_md == '':
		raise EmptyError('Empty md file!')
	md_seq = []
	md_norsnet_raw = []
	md_norsnet_bin = []
	md_bval_raw = []
	md_bval_bin = []
	md_ucon_raw = []
	md_ucon_bin = []
	md_raw = []
	md_ri = []
	md_bin = []
	
	for line in d_md.split('\n'):
		if line.startswith('Number'): continue	#The header
		if line == '': break	#We reached the end of the output block
		tokens = line.split()
		if len(tokens) != 11:
			raise ParseError("It seems that we have an issue now. MD returned an unexpected number of columns!")
		md_seq.append( tokens[1] )
		md_norsnet_raw.append( float(tokens[2]) )
		md_norsnet_bin.append( tokens[3].replace('D','+') )
		md_bval_raw.append( float(tokens[4]) )
		md_bval_bin.append( tokens[5].replace('D','+') )
		md_ucon_raw.append( float(tokens[6]) )
		md_ucon_bin.append( tokens[7].replace('D','+') )
		md_raw.append( float(tokens[8]) )
		md_ri.append( int(tokens[9]) )
		md_bin.append( tokens[10].replace('D','+') )
	
	#Check it!
	if len(md_seq) != len(md_norsnet_raw) != len(md_norsnet_bin) != len(md_bval_raw) != len(md_bval_bin) != len(md_ucon_raw) != len(md_ucon_bin) != len(md_raw) != len(md_ri) != len(md_bin): 
		raise ParseError("It seems that we have an issue now. MD returned unequal column lengths!")
		
	return {'seq':''.join(md_seq), 'norsnet_raw':md_norsnet_raw, 'norsnet_bin':md_norsnet_bin, 'bval_raw':md_bval_raw, 'bval_bin':md_bval_bin, 'ucon_raw':md_ucon_raw, 'ucon_bin':md_ucon_bin, 'prd_raw':md_raw, 'prd_ri':md_ri, 'prd_bin':md_bin}  


def parse_profsecacc(d_prof):
	"""
	Returns a dictionary where keys have the same designation as the column names in prof's tabular output.
	Values hold lists of per-residue predictions.
	
	AA
	OHEL    
	PHEL    
	RI_S    
	OACC    
	PACC    
	OREL    
	PREL    
	RI_A    
	pH      
	pE      
	pL      
	Obe     
	Pbe     
	Obie    
	Pbie    
	OtH     
	OtE     
	OtL     
	Ot0     
	Ot1     
	Ot2     
	Ot3     
	Ot4     
	Ot5     
	Ot6     
	Ot7    
	Ot8
	Ot9
	
	Their meaning (taken from prof's output):
	# NOTATION BODY      : PROFsec
	# NOTATION OHEL      : observed secondary structure: H=helix, E=extended (sheet), blank=other (loop)
	# NOTATION PHEL      : PROF predicted secondary structure: H=helix, E=extended (sheet), blank=other (loop) PROF = PROF: Profile network prediction HeiDelberg
	# NOTATION RI_S      : reliability index for PROFsec prediction (0=lo 9=high) Note: for the brief presentation strong predictions marked by '*'
	# NOTATION pH        : 'probability' for assigning helix (1=high, 0=low)
	# NOTATION pE        : 'probability' for assigning strand (1=high, 0=low)
	# NOTATION pL        : 'probability' for assigning neither helix, nor strand (1=high, 0=low)
	# NOTATION OtH       : actual neural network output from PROFsec for helix unit
	# NOTATION OtE       : actual neural network output from PROFsec for strand unit
	# NOTATION OtL       : actual neural network output from PROFsec for 'no-regular' unit
	# 
	# ------------------------------------------------------------------------
	# NOTATION BODY      : PROFacc
	# NOTATION OACC      : observed solvent accessibility (acc) in square Angstroem (taken from DSSP: W Kabsch and C Sander, Biopolymers, 22, 2577-2637, 1983)
	# NOTATION PACC      : PROF predicted solvent accessibility (acc) in square Angstroem
	# NOTATION OREL      : observed relative solvent accessibility (acc) in 10 states: a value of n (=0-9) corresponds to a relative acc. of between n*n % and (n+1)*(n+1) % (e.g. for n=5: 16-25%).
	# NOTATION PREL      : PROF predicted relative solvent accessibility (acc) in 10 states: a value of n (=0-9) corresponds to a relative acc. of between n*n % and (n+1)*(n+1) % (e.g. for n=5: 16-25%).
	# NOTATION RI_A      : reliability index for PROFacc prediction (0=low to 9=high) Note: for the brief presentation strong predictions marked by '*'
	# NOTATION Obe       : observerd relative solvent accessibility (acc) in 2 states: b = 0-16%, e = 16-100%.
	# NOTATION Pbe       : PROF predicted  relative solvent accessibility (acc) in 2 states: b = 0-16%, e = 16-100%.
	# NOTATION Obie      : observerd relative solvent accessibility (acc) in 3 states: b = 0-9%, i = 9-36%, e = 36-100%.
	# NOTATION Pbie      : PROF predicted relative solvent accessibility (acc) in 3 states: b = 0-9%, i = 9-36%, e = 36-100%.
	# NOTATION Ot4       : actual neural network output from PROFsec for unit 0 coding for a relative solvent accessibility of 4*4 - 5*5 percent (16-25%). Note: OtN, with N=0-9 give the same information for the other output units!
	# 
	"""
	if d_prof == '':
		raise EmptyError('Empty prof file!')
	
	ret = {}
	for line in d_prof.split('\n'):
		if not line.startswith('#') and line != '':
			#First parse the column header
			if line.startswith('No'):
				column_names = re.split('\s+', line)
			#Now the predicted values per line
			else:
				value_tokens = re.split('\s+', line)
				for i in range(len(value_tokens)):
					#Get its specific column name
					col = column_names[i]
					#Try to convert the current value into an integer.
					#If that fails we are dealing with a string
					try: 
						val = int(value_tokens[i])
					except ValueError: 
						val = value_tokens[i]
					#Now append it
					try: 
						ret[col].append(val)
					except KeyError:
						ret[col] = [val]
	
	#Do some final consistency checks: Has everything the same length?
	l = len(list(ret.values())[0])
	for listt in ret.values():
		if len(listt) != l:
			raise ParseError("Something happened! profsecacc returns different column lengths!")
	
	#Add an additional entry containing the concatenated aa sequence
	seq = ''.join(ret['AA'])
	ret['seq'] = seq
	
	return ret


def parse_profbval(d_bval):
	"""
	Returns a dictionary with keys and values as follows:
	
	'prd_raw1': list of integers corresponding to first output node
	'prd_raw2': list of integers corresponding to second output node
	
	Unfortunately, there is neither sequence information nor a binary prediction to be found in the output.
	"""
	if d_bval == '':
		raise EmptyError('Empty bval file!')
	
	out1 = []
	out2 = []
	region_of_interest = False
	for line in d_bval.split('\n'):
		if region_of_interest:
			tokens = line.split()
			if len(tokens) == 0: continue
			out1.append( int(tokens[1]) )
			out2.append( int(tokens[2]) )
			
		if line.startswith('* out vec:'):
			region_of_interest = True
	
	#Consistency check
	if len(out1) != len(out2):
		raise ParseError("Something happened! profbval returns different column lengths!")
	
	return {'prd_raw1':out1, 'prd_raw2':out2}
	
	
def parse_pfam_annotations(d_hmmer):
	"""
	This method performs residue-wise domain annotations according to aligned pfam domains. Search against the PfamA database should be performed
	by means of the new hmmer3 suite, so using the old hmmer2 is strongly discouraged, due to its different output style!
	The parsing depends on hmmer-3.0rc1 output (as of February 2010), so check that before running a newer hmmer3!!!
	
	Each sequence position of the query seq is annotated in a dictionary in the following way:
	
	{	...
		i: {mdl#:[(query_i,domain_i-eval,consensus_i,match_i,pp_i),(...)],  mdl#:[(...),(...),...] } , 
		i+j: {mdl#:[(...),(...),...], ...},
		... 
	}, 
	
	where 
		i: the i-th position in the query seq (starting at 0!!),
		mdl#: the number of the model
		query_i: the residue of the query sequence, intended for checking purposes for the function caller
		domain_i-eval: the domain's i-evalue 
		consensus_i: the aligned residue in the consensus pfam domain
		match_i: the information of the conservation grade ,
		pp_i: the posterior probability of that specific aligned residue (new to hmmer3)
	
	Note, the hierarchy of hmmer output:
	A query sequence could match to different Pfam models, each consisting of several domains. Furthermore, a residue could be aligned to 
	more than one domain _within_ a model, hence the assigned list to each model number in the nested dictionary:
	Each entry essentially refers to one domain where that specific residue i is aligned to.
	
	
	A sample hmmer output against PfamA could look like this:
	-------------------------------------------------------------------------------------------------------------------------------------
	Query:       query  [L=386]
	Scores for complete sequence (score includes all domains):
	   --- full sequence ---   --- best 1 domain ---    -#dom-
	    E-value  score  bias    E-value  score  bias    exp  N  Model      Description
	    ------- ------ -----    ------- ------ -----   ---- --  --------   -----------
	    3.1e-94  315.1   5.1      1e-78  264.1   0.4    2.7  2  PF00224.14 Pyruvate kinase, barrel domain
	    5.2e-26   90.1   6.7      6e-26   89.9   3.7    1.8  1  PF02887.9  Pyruvate kinase, alpha/beta domain
	
	
	Domain annotation for each model (and alignments):
	>> PF00224.14  Pyruvate kinase, barrel domain
	   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
	 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
	   1 !   53.2   0.0   2.2e-18   1.3e-14       2      72 ..      24      89 ..      23      90 .. 0.93
	   2 !  264.1   0.4   1.7e-82     1e-78     173     344 ..      89     259 ..      88     263 .. 0.96
	
	  Alignments for each domain:
	  == domain 1    score: 53.2 bits;  conditional E-value: 2.2e-18
	                --SEEEEEE--TTTSHHHHHHHHH----EEEEETT---HHHHHHHHHHHHHHHHCTTTS-EEEEE------ CS
	  PF00224.14  2 rrtkivctlGPasesvekleklieaGlnvvRlnfshGsheehkeridnvreaeeklgkkvaillDtkGpei 72
	                ++t+ivctlGPa +sve+l kli+aG+++ R+n    she+hke  +nv +a+ +l   +++llDtkGp i
	       query 24 KKTHIVCTLGPACKSVETLVKLIDAGMDICRFN----SHEDHKEMFNNVLKAQ-ELRCLLGMLLDTKGPPI 89
	                89******************************9....789*********9986.56788**********76 PP
	
	  == domain 2    score: 264.1 bits;  conditional E-value: 1.7e-82
	                 SS-HHHHHHHH---TT.-SEEEETTE-SHHHHHHHHHHHHHTTTTSEEEEEE-S----TTHHHHHHH----EEE-------S-GGGHHHHHHHHHHHCCC-----EEESSTTGGGGTSSS--HHHHHHHHHHHH----EEEE---------HHHHHHHHHHHHHHHHCTS-H CS
	  PF00224.14 173 alsekDkadlkfgvkqgvdliaasfvRkaedvkevRevleekgkeikiiakienqegvenldeileasdgimvaRGDlGieipaekvvlaqkllikkcnlagkpvitatqmlesmiknPrptRaevsDvanavldGaDavmLsgetakGkyPveavkamaevaleaekalke 344
	                  +sekDk+d+   +    ++iaasf+ +a+dv+ +R++l+++g++ikii kien eg+ ++d+il +sdgim+aRGDlG+ei  ekv+laqkl+i+kcnl gkp+itatqmlesm+knPrptRaev+DvanavldG+D+vmLsgeta Gk+Pveav++m++++leae+ +++
	       query  89 IISEKDKNDILNFAIPMCNFIAASFIQSADDVRLIRNLLGPRGRHIKIIPKIENIEGIIHFDKILAESDGIMIARGDLGMEISPEKVFLAQKLMISKCNLQGKPIITATQMLESMTKNPRPTRAEVTDVANAVLDGTDCVMLSGETA-GKFPVEAVTIMSKICLEAEACIDY 259
	                 69******9765555579********************************************************************************************************************************8.*******************99986 PP
	
	>> PF02887.9  Pyruvate kinase, alpha/beta domain
	   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
	 ---   ------ ----- --------- --------- ------- -------    ------- -------    ------- -------    ----
	   1 !   89.9   3.7     1e-29     6e-26       2     116 ..     278     383 ..     277     384 .. 0.94
	
	  Alignments for each domain:
	  == domain 1    score: 89.9 bits;  conditional E-value: 1e-29
	                HHHHHHHHHHHHH----EEEEE-----HHHHHHCC---..EEEEE----HHH---EEE---TT---HHHHCHHHHHHHHHCCHHH-----SSS-EEEE--....-------EEEE CS
	  PF02887.9   2 eaiaeaaveaAkelgakaIvvltesGstarlvskyrpgvpIlavtpseetarqlalvwGvhplvgkeraistdeviaealraalkkglikkgdevvvtaglpfgtaggtntikvv 116
	                ea+a++ave+A++++a+ I++lte+G+tarl++ky+p++ Ila++ s++t + l++++Gv+++ + +    td vi++a+++a++++++k gd v++++g       +tn++kvv
	      query 278 EAVARSAVETAESIQASLIIALTETGYTARLIAKYKPSCTILALSASDSTVKCLNVHRGVTCIKVGSF---TDIVIRNAIEIAKQRNMAKVGDSVIAIHG------IKTNLMKVV 383
	                99************************************************************544444...59***************************......589999998 PP
	-------------------------------------------------------------------------------------------------------------------------------------
	
	Each model is introduced by an '>>', each model could have several domains, introduced by an '=='.
	Mind e.g. query residue i=88 in the first model (89 in the output above): It is annotated in both domains. Hence its annotation in the return dictionary would
	look like:
	
	88:{0:[('I','1.3e-14', 'i', 'i', '6'), ('I','1e-78', 'a', ' ', '6')]}
	
	If it would align in a domain of the second model, that annotation would accur as another entry in the sub-dictionary, introduced by a 1.
	Here, you can also see what is actually used as annotation: first the i-evalue of the domain (1.3e-14 or 1e-78) followed by the subject (consensus) residue, the
	conservation letter (line between query and subject) and the posterior probability (beneath the query line).
	There could be other information to be extracted (like bit score, start stop positions...). Perhaps in the future.
	"""
	if 'No hits detected that satisfy reporting thresholds' in d_hmmer:
		#raise NoResultError('No significant hit found')
		raise NoResultError('hmmer3 did not detect any hits.')
	#First we split up into models
	#Look for the '>>' at the beginning of the line.
	rgx = re.compile('^>>', re.M)
	models_tmp = rgx.split(d_hmmer)
	models = []
	for model in models_tmp[1:-1]:	#The first one is the hmmscan header and general information, we aren't interested in that one; the last one
		models.append(model)		#needs to be purged off the footer
	#Get rid of the last model's footer and append it to models
	rgx = re.compile('^Internal pipeline statistics summary', re.M)
	models.append(rgx.split(models_tmp[-1])[0])
	#Now handle each single domain within the models  and save models and their domains into model_list
	rgx = re.compile('^  ==', re.M)	#Each domain starts with this string
	mdl_cnt = 0	#How many models are we dealing with? Remember, each model is made up of at least one domain
	#model_list = {}	#Here we store each model along with their domains
	residues = {}
	for model in models:
		if 'No individual domains that satisfy reporting thresholds' in model:
			continue
		domains = rgx.split(model)
		domains_header = domains[0]
		domains_aligns = domains[1:]
		#Parse the header first for domain-specific information
		#   #    score  bias  c-Evalue  i-Evalue hmmfrom  hmm to    alifrom  ali to    envfrom  env to     acc
		dom = []
		for line in domains_header.split('\n'):	
			if re.match('\s+\d+\s+', line):
				tokens = line.strip().split()
				assert tokens[1] == '!'		#If this fails, we didn't fully understand what's going on during parsing!
				score = tokens[2]			#See the docu about the significance of both evalues. i-eval (independent) is the one
				c_eval = tokens[4]			#we're probably interested in later, as well as the bitscore.
				i_eval = tokens[5]			#Yes we are, since pfam website only reports i_eval as _the_ evalue
				hmm_endpoints = tokens[8]	#Is either '..', '.]', '[.' or '[]', indicating if alignment ended internally a domain ('.') or end exactly with domain boundaries ([/])
				start_model = tokens[6]
				end_model = tokens[7]
				start_query = tokens[9]
				end_query = tokens[10]
				
				#Not all of this is currently needed.
				info = {'score':score, 'i_eval':i_eval, 'start_model':start_model, 'end_model':end_model, 'start_query':start_query, 'end_query':end_query, 'hmm_endpoints':hmm_endpoints}
				dom.append(info)
		
		#Now handle the alignments in each domain
		i = 0
		feature_string = ''
		for algn in domains_aligns:
			lines = algn.strip().split('\n')
			#There could be up to two additional annotation lines present above the actual alignment (s. hmmer docu p.18). Get rid of those!
			if len(lines) == 5:
				lines = lines[1:]
			elif len(lines) == 6:
				lines = lines[2:]
			elif len(lines) == 7:
				lines = lines[3:]
			else:
				raise ParseError('Well, that is indeed interesting. Something went terribly wrong during assuming the amount of possible lines per alignment! I think I will be dying now!')
			
			line_model = lines[0]		#The line containing the consensus of the domain sequence
			line_match = lines[1]		
			line_target = lines[2]		#Our target sequence for which we just found a homologous domain
			line_pp = lines[3]
			
			name_model = line_model.split()[0]
			start_query = int(dom[i]['start_query'])
			end_query = int(dom[i]['end_query'])
			
			seq_model = line_model.split()[2]		#The domain consensus sequence
			seq_query = line_target.split()[2]		#The query sequence
			#We need the start index of the match sequence which is the same as from all the others, e.g. seq_model
			m_start = line_model.index(seq_model)
			seq_match = line_match[m_start:]		#The match sequence between both
			seq_pp = line_pp.lstrip().split()[0]	#The posterior probability sequence
			#Some semantics checks: each string length has to be the same. Otherwise, something went wrong during parsing!
			assert len(seq_model) == len(seq_match) == len(seq_query) == len(seq_pp)
			#Now do the mapping
			actual_pos = 0
			for pos in range(len(seq_query)):
				if seq_query[pos] == '-':
					continue
				
				try:
					residues[actual_pos+start_query-1][mdl_cnt].append( (seq_query[pos],dom[i]['i_eval'],seq_model[pos], seq_match[pos], seq_pp[pos]) )
				except KeyError:
					try:
						residues[actual_pos+start_query-1][mdl_cnt] = [ (seq_query[pos],dom[i]['i_eval'],seq_model[pos], seq_match[pos], seq_pp[pos]) ]
					except KeyError:
						residues[actual_pos+start_query-1] = {mdl_cnt: [ (seq_query[pos],dom[i]['i_eval'],seq_model[pos], seq_match[pos], seq_pp[pos]) ] }
				
				actual_pos += 1
			
			#A further consistency check!
			assert end_query == actual_pos+start_query-1
			#Proceed with the next residue within the alignment
			i += 1
		
		#Proceed with the next model
		mdl_cnt += 1	
		
	return residues


def parse_prosite(d_prosite):
	"""
	"""
	if d_prosite == '':
		raise EmptyError('Empty prosite file!')
	
	stretches = []	#Here we store the prosite matches: A list of 3-lists, i.e. one per prosite match: start,stop,stretch
	within = False
	for line in d_prosite.split('\n'):
		if line.startswith('Pattern:'):
			within = True
			continue
		if line.startswith('Pattern-ID:') or line.strip() == '':
			within = False
			continue
		if within:
			tokens = line.strip().split()
			start = int(tokens[0])	#Prosite starts counting at 1!
			stretch = tokens[1]
			stop = start + len(stretch) - 1
			#We return sequence positions 0-based!
			stretches.append( [start-1,stop-1,stretch] )	
	
	return stretches


if __name__ == '__main__':
	import os
	from lib_parser import *
	import sys
	
	pp_path = '/mnt/home/schaefer/SNAPv2/pp/'
	
	
	chains = os.listdir(pp_path)
	i = 0
	N = len(chains)
	mn = 0
	mx = 0
	for chain in chains:
		#print chain
		i += 1
		#print chain, i, N
		try:
			#d_blast = open(pp_path+chain+"/query.blastPsiMat").read()
			#d_disis = open(pp_path+chain+"/query.disis").read()
			#d_isis = open(pp_path+chain+"/query.isis").read()
			#d_md = open(pp_path+chain+"/query.mdisorder").read()
			#d_prof = open(pp_path+chain+"/query.profRdb").read()
			#d_bval = open(pp_path+chain+"/query.profbval").read()
			#d_psic = open(pp_path+chain+"/query.psic").read()
			d_in = open(pp_path+chain+"/query.in").read()
			d_fasta = open(pp_path+chain+"/query.fasta").read()
			#d_hmmer = open(pp_path+chain+"/query.hmm3pfam").read()
			d_prosite = open(pp_path+chain+"/query.prosite").read()
		except NoResultError:
			#print 'too short for psic'
			continue
		except IOError:
			print('file not found')
			continue
		
		#print d_hmmer
		seq = parse_sequence(d_in, d_fasta)['seq']
		for stretch in parse_prosite(d_prosite):
			print (stretch[2])
			print (seq[stretch[0]:stretch[1]+1])
			assert stretch[2] == seq[stretch[0]:stretch[1]+1]
