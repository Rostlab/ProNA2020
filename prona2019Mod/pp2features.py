#!/usr/bin/python3
import os
import os.path
import sys
from prona2019Mod.lib_parser import *
from prona2019Mod.lib_featurize import *
import random
import numpy as np

np.random.seed(0)
random.seed( 0 )

def expand_list(nested_list):
    for item in nested_list:
        if isinstance(item, (list, tuple)):
            for sub_item in expand_list(item):
                yield sub_item
        else:
            yield item



def extract_features(pp_path ,feature_file):
	if pp_path[-1] != '/':
		pp_path += '/'

	feature_space = []
	fea_win = {}

	for l in open(feature_file):
		l = l.rstrip()
		t = l.split('\t')
		fea_win[t[0]] = int(t[1])	

	'''
	#First look for the wanted sequence positions plus additional features like class label,...
	try:
		#seq_positions, additional_feature_header = parse_classfile( open("%s/%s/%s" % (pp_path,dir,classfile)).read() )
		add_names, pos_and_features, add_domains = parse_obligatory_arff( open("%s%s/%s" % (pp_path,dir,classfile)).read() )
	except IOError:
		 f_e.write("%s: %s not found, skipping...\n" % (dir,classfile) )
		 f_e.flush()
		 ignored += 1
		 continue
	'''
	#Look which pp features are defined in the config and parse those
	#0. pp's sequences?
	try:
		d_in = open(pp_path+"query.in").read()
		d_fasta = open(pp_path+"query.fasta").read()
		pp_seq,pro = parse_sequence(d_in, d_fasta)
		seq = pp_seq['seq']
	except IOError:
		sys.exit("Error!!!ProNA2019 can not be done for protein %s.\nFile *.in and/or *.fasta not found...\n" % pro )
		
	'''	
	#1. Raw sequence?
	if 'raw_seq' in wanted:
		parsed_seq = seq
	else:
		parsed_seq = None 
	'''
	
	#2. blast?
	try: 
		d_blast = open(pp_path+"query.blastPsiMat").read()
	except IOError: 
		sys.exit("Error!!!ProNA2019 can not be done for protein %s.\nFile *.blastPsiMat not found...\n" % pro)
	parsed_blast = parse_blast(d_blast)
	if seq != parsed_blast['seq']:
		sys.exit("Error!!!ProNA2019 can not be done for protein %s.\nDisagreement between pp sequence and blast sequence...\n" % pro)
	blast_norm = True
	
	#5. PROFbval?
	try:
		d_bval = open(pp_path+"query.profbval").read()
	except IOError:
		sys.exit("Error!!!ProNA2019 can not be done for protein %s.\nFile *.profbval not found...\n" % pro)
	parsed_bval = parse_profbval(d_bval)
	#We can only check here for identical lengths since bval doesn't return a sequence
	if len(seq) != len(parsed_bval['prd_raw1']) != len(parsed_bval['prd_raw2']):
		sys.exit("Error!!!ProNA2019 can not be done for protein %s.\nDisagreement between pp sequence and PROFbval sequence...\n" % pro)
	
	bval_norm = True
	
	
	#6. MD?
	try:
		d_md = open(pp_path+"query.mdisorder").read()
	except IOError:
		sys.exit("Error!!!ProNA2019 can not be done for protein %s.\nFile *.mdisorder not found...\n" % pro)
	parsed_md = parse_md(d_md)
	if seq != parsed_md['seq']:
		sys.exit("Error!!!ProNA2019 can not be done for protein %s.\nDisagreement between pp sequence and md sequence...\n" % pro)
		
	md_norm = True
	
	
	
	#7. profphd? sec strct comp? solv acc comp?
	try:
		d_profsecacc = open(pp_path+"query.profRdb").read()
	except IOError:
		sys.exit("Error!!!ProNA2019 can not be done for protein %s.\nFile *.profRdb not found...\n" % pro)
	parsed_profsecacc = parse_profsecacc(d_profsecacc)
	if seq != parsed_profsecacc['seq']:
		sys.exit("Error!!!ProNA2019 can not be done for protein %s.\nDisagreement between pp sequence and profseccacc sequence...\n" % pro)
	profsecacc_norm = True
	parsed_sec = ''.join(parsed_profsecacc['PHEL'])
	parsed_acc = ''.join(parsed_profsecacc['Pbie'])
	'''
	#8. psic?
	if not wanted.isdisjoint( set(['psic','num_seq']) ):
		flag_done = False
		try:
			d_psic = open(pp_path+dir+"/query.psic").read()
		#In that case, there is no query.psic due to whatever reason. Instead we create a position-independent psic 
		#filled with default values grabbed from SNAP's extractAllSingle.pl script.
		#No idea where those value came from originally.
		#That's the only case currently where we use default values when a query file is not there.
		except IOError:
			psic = len(seq)*[[0.85786,0.45676,0.47306,0.58022,0.18036,0.37722,0.59724,0.81155,0.21639,0.52944,0.81156,0.58717,0.21109,0.39946,0.48178,0.63047,0.60835,0.14256,0.36310,0.68436]]
			NumSeq = len(seq)*[0]
			parsed_psic = {'psic':psic,'NumSeq':NumSeq}
			f_e.write("%s: query.psic not found, using default values...\n" % dir)
			f_e.flush()
			flag_done = True
			
		try:	
			if not flag_done: parsed_psic = parse_psic(d_psic)
		#In that case, the sequence was too short for psic. Instead we create a position-independent psic 
		#filled with default values grabbed from SNAP's extractAllSingle.pl script.
		#No idea where those value came from originally.
		except NoResultError:
			psic = len(seq)*[[0.85786,0.45676,0.47306,0.58022,0.18036,0.37722,0.59724,0.81155,0.21639,0.52944,0.81156,0.58717,0.21109,0.39946,0.48178,0.63047,0.60835,0.14256,0.36310,0.68436]]
			NumSeq = len(seq)*[0]
			parsed_psic = {'psic':psic,'NumSeq':NumSeq}
			f_e.write("%s: sequence too short for psic, using default values...\n" % dir)
			f_e.flush()
			
		if len(seq) != len(parsed_psic['psic']):
			f_e.write("%s: Disagreement between pp sequence and psic sequence length, skipping...\n" % dir)
			f_e.flush()
			ignored += 1
			continue
		psic_norm = 'psic_normalize' in wanted
	else:
		parsed_psic = None
		
	
	#9. pfam?
	if not wanted.isdisjoint( set(['pfam_within_domain','pfam_dom_cons','pfam_residue_fit','pfam_pp']) ):
		try:
			d_pfam = open(pp_path+dir+"/query.hmm3pfam").read()
		except IOError:
			f_e.write("%s: query.hmm3pfam not found, skipping...\n" % dir)
			f_e.flush()
			ignored += 1
			continue
		#In that case hmmer3 did not find anything of significance.
		#We simply create an empty dictionary which means that no residues are annotated at all.
		#The featurizer handles that correctly by designating default values.
		try:
			parsed_pfam = parse_pfam_annotations(d_pfam)
		except NoResultError:
			f_e.write("%s: no significant pfam hit, using default values...\n" % dir)
			f_e.flush()
			parsed_pfam = {}
	else:
		parsed_pfam = None
	
	'''
	#10. Chemical properties?
	chemprop = True
		
	'''
	#11. Prosite?
	if not wanted.isdisjoint( set(['prosite_part']) ):
		try:
			d_prosite = open(pp_path+dir+"/query.prosite").read()
		except IOError:
			f_e.write("%s: query.prosite not found, skipping...\n" % dir)
			f_e.flush()
			ignored += 1
			continue
		parsed_prosite = parse_prosite(d_prosite)
	else:
		parsed_prosite = None
	'''
	#Global features
	#12. aa comp?
	glbl_aa_comp = True
	
	#13. protein length?
	glbl_length = True
	'''	
	#14. sec strct comp?
	try:
		d_profsecacc = open(pp_path+"query.profRdb").read()
	except IOError:
		sys.exit("%s: *.profRdb not found, skipping...\n" % pro)
	parsed_profsec = parse_profsecacc(d_profsecacc)
	if seq != parsed_profsec['seq']:
		sys.exit("%s: Disagreement between pp sequence and profseccacc sequence, skipping...\n" % pro)
	parsed_sec = ''.join(parsed_profsec['PHEL'])
	
	#15. solv acc comp?
	try:
		d_profsecacc = open(pp_path+"query.profRdb").read()
	except IOError:
		sys.exit("%s: *.profRdb not found, skipping...\n" % pro)
	parsed_profacc = parse_profsecacc(d_profsecacc)
	if seq != parsed_profacc['seq']:
		sys.exit("%s: Disagreement between pp sequence and profseccacc sequence, skipping...\n" % pro)
	parsed_acc = ''.join(parsed_profacc['Pbie'])
	'''	

	
	#Now, handle each sequence position
	for pos in range(len(seq)):
		#for instance in add_features:
		instance_feature_names = {}		#We would need that only once, but this way it's easier
		instance_features = {}


		'''
		#0. Protein identifier
		if 'prot_id' in wanted:
			instance_feature_names.append('ID_pos')
			instance_features.append(pro+'_'+str(pos))
			instance_feature_domains.append("STRING")
		'''
		'''	
		#1. Raw sequence
		if parsed_seq:
			featnames_seq, window_seq = feature_sequence(parsed_seq, pos, fea_win)
			#First flatten the list of lists
			names = [j for i in featnames_seq['names_seq'] for j in i]
			feats = [j for i in window_seq['win_seq'] for j in i]
			doms = [j for i in domains_seq['domains_seq'] for j in i]
			#Extend the names, features and domains
			instance_feature_names.extend(names)
			instance_features.extend(feats)
			instance_feature_domains.extend(doms)
		'''

		#2. Blast features
		if parsed_blast:
			featnames_blast, window_blast = feature_blast(parsed_blast, pos, fea_win, norm_raw=blast_norm)
			instance_feature_names.update(featnames_blast)
			instance_features.update(window_blast)	
		'''	
		#3. disis features
		if parsed_disis:
			win_disis = windows['disis_window']
			featnames_disis, window_disis, domains_disis = feature_disis(parsed_disis, pos, win_disis, norm_raw=disis_norm)
			if 'disis_bin' in wanted:
				#Extend the names and features
				instance_feature_names.extend(featnames_disis['names_bin_plus'])
				instance_feature_names.extend(featnames_disis['names_bin_minus'])
				instance_features.extend(window_disis['win_bin_plus'])
				instance_features.extend(window_disis['win_bin_minus'])
				instance_feature_domains.extend(domains_disis['domains_bin_plus'])
				instance_feature_domains.extend(domains_disis['domains_bin_minus'])
			if 'disis_raw' in wanted:
				#Extend the names and features
				instance_feature_names.extend(featnames_disis['names_raw'])
				instance_features.extend(window_disis['win_raw'])
				instance_feature_domains.extend(domains_disis['domains_raw'])
		
		
		#4. isis features
		if parsed_isis:
			win_isis = windows['isis_window']
			featnames_isis, window_isis, domains_isis = feature_isis(parsed_isis, pos, win_isis, norm_raw=isis_norm)
			if 'isis_bin' in wanted:
				#Extend the names and features
				instance_feature_names.extend(featnames_isis['names_bin_plus'])
				instance_feature_names.extend(featnames_isis['names_bin_minus'])
				instance_features.extend(window_isis['win_bin_plus'])
				instance_features.extend(window_isis['win_bin_minus'])
				instance_feature_domains.extend(domains_isis['domains_bin_plus'])
				instance_feature_domains.extend(domains_isis['domains_bin_minus'])
			if 'isis_raw' in wanted:
				#Extend the names and features
				instance_feature_names.extend(featnames_isis['names_raw'])
				instance_features.extend(window_isis['win_raw'])
				instance_feature_domains.extend(domains_isis['domains_raw'])
		'''
		#5. PROFbval features
		if parsed_bval:
			featnames_bval, window_bval = feature_profbval(parsed_bval, pos, fea_win, norm_raw=bval_norm)
			instance_feature_names.update(featnames_bval)
			instance_features.update(window_bval)	
	
		
		#6. MD features
		if parsed_md:
			featnames_md, window_md = feature_md(parsed_md, pos, fea_win, norm_raw=md_norm)
			instance_feature_names.update(featnames_md)
			instance_features.update(window_md)	
	
		
		#7. profphd features
		if parsed_profsecacc:
			featnames_profphd, window_profphd = feature_profsecacc(parsed_profsecacc, pos, fea_win, norm_raw=profsecacc_norm)
			instance_feature_names.update(featnames_profphd)
			instance_features.update(window_profphd)  			
	
		'''
		#8. psic features
		if parsed_psic:
			win_psic = windows['psic_window']
			featnames_psic, window_psic, domains_psic = feature_psic(parsed_psic, pos, win_psic, norm_raw=psic_norm)
			if 'psic' in wanted:
				#First flatten the list of lists
				names = [j for i in featnames_psic['names_psic'] for j in i]
				feats = [j for i in window_psic['win_psic'] for j in i]
				doms = [j for i in domains_psic['domains_psic'] for j in i]
				instance_feature_names.extend(names)
				instance_features.extend(feats)
				instance_feature_domains.extend(doms)
			if 'num_seq' in wanted:
				instance_feature_names.extend(featnames_psic['names_numSeq'])
				instance_features.extend(window_psic['win_numSeq'])
				instance_feature_domains.extend(domains_psic['domains_numSeq'])
				
		
		#9. pfam features
		if parsed_pfam != None:
			win_pfam = windows['pfam_window']
			featnames_pfam, window_pfam, domains_pfam = feature_pfam(parsed_pfam, pos, win_pfam)
			if 'pfam_within_domain' in wanted:
				instance_feature_names.extend(featnames_pfam['names_within_domain'])
				instance_features.extend(window_pfam['win_within_domain'])
				instance_feature_domains.extend(domains_pfam['domains_within_domain'])
			if 'pfam_dom_cons' in wanted:
				instance_feature_names.extend(featnames_pfam['names_domain_conservation'])
				instance_features.extend(window_pfam['win_domain_conservation'])
				instance_feature_domains.extend(domains_pfam['domains_domain_conservation'])
			if 'pfam_residue_fit' in wanted:
				instance_feature_names.extend(featnames_pfam['names_residue_fit'])
				instance_features.extend(window_pfam['win_residue_fit'])
				instance_feature_domains.extend(domains_pfam['domains_residue_fit'])
			if 'pfam_pp' in wanted:
				instance_feature_names.extend(featnames_pfam['names_pp'])
				instance_features.extend(window_pfam['win_pp'])
				instance_feature_domains.extend(domains_pfam['domains_pp'])
		
		'''
		#10. chem props
		if chemprop:
			featnames_chemprop, window_chemprop = feature_chem_props(seq, pos, fea_win)
			instance_feature_names.update(featnames_chemprop)
			instance_features.update(window_chemprop)	
		'''
		#11. Prosite
		if parsed_prosite != None:
			win_prosite = windows['prosite_window']
			featnames_prosite, window_prosite, domains_prosite = feature_prosite(parsed_prosite, pos, win_prosite)
			if 'prosite_part' in wanted:
				instance_feature_names.extend(featnames_prosite['names_prosite'])
				instance_features.extend(window_prosite['win_prosite'])
				instance_feature_domains.extend(domains_prosite['domains_prosite'])
		'''
		#The global fatures
		#12. aa composition
		if glbl_aa_comp:
			featnames_aa_comp, features_aa_comp = feature_aa_composition(seq, fea_win)
			instance_feature_names.update(featnames_aa_comp)
			instance_features.update(features_aa_comp)	

		#13. Protein length
		if glbl_length:
			featnames_length, features_length = feature_protein_length(seq, fea_win)
			instance_feature_names.update(featnames_length)
			instance_features.update(features_length) 
	
		#14 sec strct comp
		if parsed_sec != None:
			featnames_sec_comp, features_sec_comp = feature_sec_strct_composition(parsed_sec, fea_win)
			instance_feature_names.update(featnames_sec_comp)
			instance_features.update(features_sec_comp) 
	
		#15. solv acc comp
		if parsed_acc != None:
			featnames_acc_comp, features_acc_comp = feature_solv_acc_composition(parsed_acc, fea_win)
			instance_feature_names.update(featnames_acc_comp)
			instance_features.update(features_acc_comp)	
		
		
	#Finally, append the current instance to the great picture
		temp=[]
		for l in open(feature_file):
			l = l.rstrip()
			x = l.split('\t')[0]
			temp.append(instance_features[x])
		tempp = list(expand_list(temp))
		feature_space.append(tempp)
		
	return feature_space,pro,seq





	
