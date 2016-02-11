"""
Do motif cassette discovery on sequences - a pipelined version.

This pipeline trys to be as agnostic as possible about the sequence data that
are used within it. To lay out in plain English the expectations for this analysis:

* A 'data' directory, containing:
	* A 'sequences' directory, containing:
		* A 'control' directory, containing:
			* A single sequence file
		* An 'experimental' dircetory, containing:
			* A single sequence file
		* Optionally, a 'comparative' directory, containing:
			* Any number of sequence files

'build' (scratch) and 'results' directories will be created by the program.

Various notes:

* The names of files and directories are important in as much as they are used
	to name datafiles and columns later in the procedure. Thus, name your data
	in a useful way.
* The control data is always called - regardless of the what the input file is
	named - 'control'. Likewise, the experimental data is always called
	'experimental'.
* All sequence files are FASTA-formatted.
* The initial step (MEME analysis) uses just 100 sequences, regardless of how
	many are supplied. These will be taken as the first 100 of the control and
	experimental sequence files. They should therefore be ordered appropriately.

"""


### IMPORTS

from os import path
from os import environ
from glob import glob
import json
import csv
import bisect

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from statsmodels.sandbox.stats.multicomp import multipletests

# these are the core packages we require
import mcda


### CONSTANTS & DEFINES

# load snakemake job configfile which should have bulk of customizable steps
configfile: "find-cassettes.config.yaml"

## Directory structure
# various directories seperating our work
DATA_DIR = 'data'
BUILD_DIR = 'build'
RESULTS_DIR = 'results'

## Location of sequence data
# where to find the original sequence data
SEQ_DATA_DIR = path.join (DATA_DIR, 'seqs')
EXP_SEQ_DATA_DIR = path.join (SEQ_DATA_DIR, 'experimental')
EXP_SEQ_DATA = glob (path.join (EXP_SEQ_DATA_DIR, '*.fasta'))[0]
CNTRL_SEQ_DATA_DIR = path.join (SEQ_DATA_DIR, 'control')
CNTRL_SEQ_DATA = glob (path.join (CNTRL_SEQ_DATA_DIR, '*.fasta'))[0]
COMP_SEQ_DATA_DIR = path.join (SEQ_DATA_DIR, 'comparative')

# where to find the working sequences
SEQ_WORK_DIR = path.join (BUILD_DIR, 'seqs')

COMP_SEQ_WORK_DIR = path.join (SEQ_WORK_DIR, 'comparative')
ALL_COMP_SEQ_DATA = glob (path.join (COMP_SEQ_WORK_DIR, '*.fasta'))

ALL_CONTROL_SEQS = path.join (SEQ_WORK_DIR, 'all-control.fasta')
MTF_DISC_CONTROL_SEQS = path.join (SEQ_WORK_DIR, 'first100-control.fasta')
NONDISC_CONTROL_SEQS = path.join (COMP_SEQ_WORK_DIR, 'control-nd.fasta')

ALL_EXP_SEQS = path.join (SEQ_WORK_DIR, 'all-experimental.fasta')
MTF_DISC_EXP_SEQS = path.join (SEQ_WORK_DIR, 'first100-experimental.fasta')
NONDISC_EXP_SEQS = path.join (COMP_SEQ_WORK_DIR, 'experimental-nd.fasta')


SEQ_CNT_PTH = path.join (BUILD_DIR, 'seq_cnts.csv')

## External executables
PSPGEN_EXE = config['exe']['pspgen']
MEME_EXE = config['exe']['meme']
MAST_EXE = config['exe']['mast']

## Meme variables
# how many seqs to use in meme analysis
MEME_SEQ_CNT = config['meme']['num_seqs']

# what number of and size motif to look for
NUM_MOTIFS = config['meme']['num_motifs']
MIN_MOTIF_WIDTH = config['meme']['min_motif_width']
MAX_MOTIF_WIDTH = config['meme']['max_motif_width']

# temp work directories
PSPGEN_WORK_DIR = path.join (BUILD_DIR, 'pspgen')
MEME_PRIORS = path.join (PSPGEN_WORK_DIR, 'meme.priors')

MEME_WORK_DIR = path.join (BUILD_DIR, 'meme')
MEME_RESULTS = path.join (MEME_WORK_DIR, 'meme.xml')

## Cassette mining
CASS_WORK_DIR = path.join (BUILD_DIR, 'cassettes')
RAW_CASSETTES = path.join (CASS_WORK_DIR, 'raw-cassettes.csv')
FILTERED_CASSETTES = path.join (CASS_WORK_DIR, 'filtered-cassettes.csv')
CASSETTE_GRAPH = path.join (CASS_WORK_DIR, 'cassettes.graphml')
CASSETTE_GRAPH_PIC = path.join (RESULTS_DIR, 'cassettes.png')

CASSETTE_DATA_FIELDS = ('datafile', 'cutoff', 'max_pval', 'pattern',
	'complement', 'freq', 'elements')

# max pval of motif to accept when making cassettes
MAX_MOTIF_PVAL = config['cassettes']['max_pval']

# summary of motifs found
MOTIF_SUMMARY = path.join (RESULTS_DIR, 'motif-summary.csv')


## Cassette support

CASS_SUPPORT_WORK_DIR = path.join (BUILD_DIR, 'cassette-support')
CASS_SUPPORT_MAST_XML = path.join (CASS_SUPPORT_WORK_DIR, 'mast.xml')
CASS_SUPPORT_MAST_JSON = CASS_SUPPORT_MAST_XML.replace ('.xml', '.json')
CASSETTES_WITH_SUPPORT = path.join (RESULTS_DIR, 'discovery-cassettes.csv')

RANDOMIZATION_CNT = config['cassettes']['randomizations']

CASSETTE_SUPPORT_FIELDS = ('pattern', 'complement', 'frequency',
	'randomizations', 'ramdomization_rank', 'randomization_support')

## Mast search for cassettes
MAST_WORK_DIR = path.join (BUILD_DIR, 'mast')

EXP_MAST_RES_DIR = path.join (MAST_WORK_DIR, 'experimental')
CONTROL_MAST_RES_DIR = path.join (MAST_WORK_DIR, 'control')
COMPARATIVE_MAST_RES_DIR = path.join (MAST_WORK_DIR, 'comparative')

EXP_MAST_RES = path.join (EXP_MAST_RES_DIR, 'mast.xml')
CONTROL_MAST_RES = path.join (CONTROL_MAST_RES_DIR, 'mast.xml')

EXP_MAST_JSON = EXP_MAST_RES.replace ('.xml', '.json')
CONTROL_MAST_JSON = CONTROL_MAST_RES.replace ('.xml', '.json')

EXP_CASS_ALL = path.join (CASS_WORK_DIR, 'experimental.all.csv')
CONTROL_CASS_ALL = path.join (CASS_WORK_DIR, 'control.all.csv')
COMPARATIVE_CASS_DIR = path.join (CASS_WORK_DIR, 'comparative')

CASS_HDRS = ['pattern', 'dir', 'seq_name', 'start', 'stop', 'len', 'max_gap',
	'pvals', 'combined_pval']

## Summary tables
# list of experimental seqs and which cass they contain
EXP_CASS_SEQS = path.join (RESULTS_DIR, 'exp-seqs-and-cassettes.csv')

MAX_CASS_GAP = config['cassettes']['max_gap']
MAX_CASS_PVAL = config['cassettes']['max_pval']

EXP_CASS_SUMMARY = EXP_CASS_ALL.replace ('.all.', '.summary.')
CONTROL_CASS_SUMMARY = CONTROL_CASS_ALL.replace ('.all.', '.summary.')

CASS_SUMMARY_HDRS = ['pattern', 'freq', 'freq_seq']

OVERALL_CASS_TABLE = path.join (CASS_WORK_DIR, 'all-cassettes-table.csv')
OVERALL_CASS_TABLE_WITH_STATS = path.join (RESULTS_DIR, 'all-cassettes.with-stats.csv')
CASS_ENRICHMENT_TABLE = OVERALL_CASS_TABLE_WITH_STATS.replace ('.with-stats.', '.enrichment.')
CASS_COUNT_TABLE = OVERALL_CASS_TABLE_WITH_STATS.replace ('.with-stats.', '.counts.')

## Exemplar sequences

# how many examples of each to gather
NUM_EXEMPLARS = config['exemplars']['num_seqs']

EXEMPLAR_CASS= path.join (RESULTS_DIR, 'exemplars.csv')
EXEMPLAR_SEQS = path.join (RESULTS_DIR, 'exemplars.fasta')

## Reporting

REPORT_PTH = path.join (RESULTS_DIR, 'report.html')

REPORT_DATA_DIR = path.join (path.dirname (mcda.__file__), 'data')



### CODE ###

## Utils
# Mainly just assertions for checking values and arguments

def count_seqs_in_file (pth):
	"""
	Largely just for checking we're processing data right.
	"""
	all_seqs = mcda.read_seqs (pth)
	return len (all_seqs)


def assert_num_seqs_in_file (pth, cnt):
	"""
	Checking we have the expected number of sequences in each file.
	"""
	assert count_seqs_in_file (pth) == cnt, \
		"incorrect number of sequences in file '%s', should be %s" % (
			pth, cnt)


def get_names_of_seqs_in_file (pth):
	all_seqs = mcda.read_seqs (pth)
	return [s.name for s in all_seqs]



## Snakemake rules

rule all:
	input:
		REPORT_PTH

rule clean:
	message: "Delete all built / produced files, leave original data"
	shell:
		"""
		rm -rf {BUILD_DIR}
		rm -rf {RESULTS_DIR}
		"""


rule prep_data:
	message: "Check & prep all the data we need for this pipeline."
	input:
		all_exp_seqs=EXP_SEQ_DATA,
		all_cntrl_seqs=CNTRL_SEQ_DATA,
		comp_seq_data_dir=COMP_SEQ_DATA_DIR,
	output:
		all_exp_seqs=ALL_EXP_SEQS,
		mtf_disc_exp_seqs=MTF_DISC_EXP_SEQS,
		nondisc_exp_seqs=NONDISC_EXP_SEQS,
		all_cntrl_seqs=ALL_CONTROL_SEQS,
		mtf_disc_cntrl_seqs=MTF_DISC_CONTROL_SEQS,
		nondisc_control_seqs=NONDISC_CONTROL_SEQS,
		comp_seq_work_dir=COMP_SEQ_WORK_DIR,
	run:
		## Preconditions:
		# have they provide dinput sequences?
		assert path.exists (input.all_exp_seqs), \
			"no experimental sequences provided"
		assert path.exists (input.all_cntrl_seqs), \
			"no control sequences provided"

		# have they provided only 1 experimental & control sequence
		assert len (glob (path.join (EXP_SEQ_DATA_DIR, '*.fasta'))) == 1, \
			"need 1 and only 1 experimental sequence file in fasta format"
		assert len (glob (path.join (CNTRL_SEQ_DATA_DIR, '*.fasta'))) == 1, \
			"need 1 and only 1 control sequence file in fasta format"

		# have they provided enough of each type of sequence?
		assert MEME_SEQ_CNT < count_seqs_in_file (input.all_exp_seqs), \
			"need at least %s experimental sequences" % MEME_SEQ_CNT
		assert MEME_SEQ_CNT < count_seqs_in_file (input.all_cntrl_seqs), \
			"need at least %s control sequences" % MEME_SEQ_CNT

		# check seq names
		raw_exp_seq_names = get_names_of_seqs_in_file (input.all_exp_seqs)
		exp_seq_names = frozenset (raw_exp_seq_names)
		assert len (raw_exp_seq_names) == len (exp_seq_names), \
			"there are duplicate sequence names in the experimental set"

		raw_cntrl_seq_names = get_names_of_seqs_in_file (input.all_cntrl_seqs)
		cntrl_seq_names = frozenset (raw_cntrl_seq_names)
		assert len (raw_cntrl_seq_names) == len (cntrl_seq_names), \
			"there are duplicate sequence names in the control set"

		overlap = exp_seq_names.intersection (cntrl_seq_names)
		assert len (overlap) == 0, \
			"some sequences appear in experimental and control sets: %s" % \
			', '.join (overlap)

		## Main:
		# create required directories
		snakemake.utils.makedirs (BUILD_DIR)
		snakemake.utils.makedirs (RESULTS_DIR)
		snakemake.utils.makedirs (SEQ_WORK_DIR)
		snakemake.utils.makedirs (COMP_SEQ_WORK_DIR)

		# copy & rename controls, sample first 100 & make non-discovery set
		shell ("cp {input.all_cntrl_seqs} {output.all_cntrl_seqs}")
		all_control_seqs = mcda.read_seqs (input.all_cntrl_seqs)
		mcda.write_seqs (all_control_seqs[:MEME_SEQ_CNT],
			output.mtf_disc_cntrl_seqs)
		mcda.write_seqs (all_control_seqs[MEME_SEQ_CNT:],
			output.nondisc_control_seqs)

		# copy & rename experimental, take top 100 & make non-discovery set
		# NOTE: assumes exp seqs are in order of effect
		shell ("cp {input.all_exp_seqs} {output.all_exp_seqs}")
		all_exp_seqs = mcda.read_seqs (input.all_exp_seqs)
		mcda.write_seqs (all_exp_seqs[:MEME_SEQ_CNT],
			output.mtf_disc_exp_seqs)
		mcda.write_seqs (all_exp_seqs[MEME_SEQ_CNT:],
			output.nondisc_exp_seqs)

		# copy all comparative seqs across
		shell ("cp {input.comp_seq_data_dir}/*.fasta {output.comp_seq_work_dir}")


rule count_seqs:
	message: "Count the number of sequences in every file"
	input:
		all_exp_seqs=ALL_EXP_SEQS,
		all_cntrl_seqs=ALL_CONTROL_SEQS,
	output:
		seq_cnt_pth=SEQ_CNT_PTH
	run:
		# XXX: this is where a lot of naming stuff can come acrop
		ALL_SEQ_FILES = ALL_COMP_SEQ_DATA + [
			input.all_exp_seqs,
			input.all_cntrl_seqs,
		]

		cnts = []
		for f in ALL_SEQ_FILES:
			all_seqs = mcda.read_seqs (f)
			cnts.append ({
				'path': f,
				'name': mcda.get_file_basename (f),
				'count': len (all_seqs)
			})

		mcda.write_csv_as_dicts (cnts, output.seq_cnt_pth,
			hdr_flds=cnts[0].keys())


rule prep_discrim_meme_search:
	message: "Making priors for discrimination meme search"
	input:
		top100_exp_seqs=MTF_DISC_EXP_SEQS,
		random100_cntrl_seqs=MTF_DISC_CONTROL_SEQS,
	output:
		meme_priors=MEME_PRIORS,
	log: 'foo.log'
	run:
		## Main:
		mcda.run_pspgen (input.top100_exp_seqs, input.random100_cntrl_seqs,
			output.meme_priors, minw=MIN_MOTIF_WIDTH, maxw=MAX_MOTIF_WIDTH,
			exe=PSPGEN_EXE)


rule meme_search:
	message: "Looking for motifs with MEME using previous priors"
	input:
		top100_exp_seqs=MTF_DISC_EXP_SEQS,
		meme_priors=MEME_PRIORS,
	output:
		meme_results=MEME_RESULTS,
	run:
		## Preconditions:
		assert_num_seqs_in_file (input.top100_exp_seqs, MEME_SEQ_CNT)

		## Main:
		mcda.run_meme (input.top100_exp_seqs, MEME_WORK_DIR,
			input.meme_priors,
			nmotifs=NUM_MOTIFS, minw=MIN_MOTIF_WIDTH, maxw=MAX_MOTIF_WIDTH,
				exe=MEME_EXE)


rule mine_for_cassettes:
	message: "Look through MEME results for motif cassettes"
	input:
		meme_results=MEME_RESULTS,
	output:
		raw_cassettes=RAW_CASSETTES,
		motif_summary=MOTIF_SUMMARY,
	run:
		## Constants & defines:
		# XXX: this allows us to step up thresholds in a granular way, by
		# trying everything from 10 to 50. Should maybe parametrize?
		CUTOFFS = [float(x)/100.0 for x in range (10, 50, 1)]

		## Main:
		# read motifs in
		motif_data = mcda.parse_meme_results (input.meme_results)

		# init file for cassette results
		with open (output.raw_cassettes, 'w') as hndl:
			csv_wrtr = csv.DictWriter (hndl, fieldnames=CASSETTE_DATA_FIELDS)
			csv_wrtr.writeheader()

			gsp_search = mcda.motifsearch.MotifSearch (motif_data, MAX_MOTIF_PVAL)
			for c in CUTOFFS:
				gsp_res = gsp_search.search (c)
				for r in gsp_res:
					csv_wrtr.writerow ({
						'datafile': input.meme_results,
						'cutoff': c,
						# XXX: shitty column naming
						'max_pval': MAX_MOTIF_PVAL,
						'pattern': ''.join (r[0]),
						'freq': r[1],
					})

		# record some summary data for later reporting
		mtf_logos = list (snakemake.utils.listfiles (path.join (MEME_WORK_DIR,
			'logo_rc{num}.png')))

		mtf_summary = []
		for i, m in enumerate (motif_data.motifs):
			new_mtf = {
				'id': m.id,
				'name': m.name,
				'width': m.width,
				'sites': m.sites,
				'evalue': m.e_value,
				'logo_path': mtf_logos[i][0],
			}
			mtf_summary.append (new_mtf)
		mcda.write_csv_as_dicts (mtf_summary, output.motif_summary,
			hdr_flds=['id', 'name', 'width', 'sites', 'evalue', 'logo_path']
		)



rule distill_cassette_results:
	message: "Filter down cassettes to most common occurrences & calculate parts"
	input:
		raw_cassettes=RAW_CASSETTES,
	output:
		filtered_cassettes=FILTERED_CASSETTES,
	run:
		# TODO: do we need to set a minimal frequency for cassettes

		# read in cassettes from previous step
		# NOTE: this whole scheme relies on cassettes being in rising frequency
		all_cassettes = mcda.read_csv_as_dicts (input.raw_cassettes)

		# filter to the uniques
		all_cassettes.reverse()
		uniq_cassettes = []
		for r in all_cassettes:
			if 2 < len (r['pattern']):
				curr_is_uniq = True
				for u in uniq_cassettes:
					if mcda.is_eq_or_complement (r['pattern'], u['pattern']):
						curr_is_uniq = False
						break
				if curr_is_uniq:
					uniq_cassettes.append (r)

		# postprocess for some useful fields
		for c in uniq_cassettes:
			c['elements'] = mcda.extract_motifs_from_cassette (c['pattern'])
			c['complement'] = mcda.complement (c['pattern'])

		# add sort them by length
		uniq_cassettes.sort (key=lambda x: mcda.sort_key (x['pattern']))

		# try to get them all in the same orientation:
		# we check every cassette against every previous and rotate them
		# so the patterns match
		for j in range (len (uniq_cassettes)):
			curr_cass = uniq_cassettes[j]
			for prev_cass in uniq_cassettes[:j]:
				if prev_cass['pattern'] in curr_cass['pattern']:
					# if orientations match, do nothing
					break
				if prev_cass['pattern'] in curr_cass['complement']:
					# if matches with complement, swap
					curr_cass['pattern'], curr_cass['complement'] = \
						curr_cass['complement'], curr_cass['pattern']
					break

		# now record them all
		mcda.write_csv_as_dicts (uniq_cassettes, output.filtered_cassettes,
			hdr_flds=CASSETTE_DATA_FIELDS)



rule calc_cassette_support:
	message: "Randomize motif sequences to calculate cassette support"
	input:
		filtered_cassettes=FILTERED_CASSETTES,
		meme_results=MEME_RESULTS,
		mtf_disc_exp_seqs=MTF_DISC_EXP_SEQS,
	output:
		cassettes_with_support=CASSETTES_WITH_SUPPORT,
		work_dir=CASS_SUPPORT_WORK_DIR,
		mast_xml=CASS_SUPPORT_MAST_XML,
		mast_json=CASS_SUPPORT_MAST_JSON,
	run:
		## Preparation:
		snakemake.utils.makedirs (CASS_SUPPORT_WORK_DIR)

		## Main:
		# read in & organise cassettes
		all_cassettes = mcda.read_csv_as_dicts (input.filtered_cassettes)
		cass_dict = {c['pattern']:c for c in all_cassettes}

		# run mast over all meme data
		mcda.run_mast (input.meme_results, input.mtf_disc_exp_seqs,
			output.work_dir, exe=MAST_EXE)

		# convert to json
		mast_data = mcda.parse_mast_results (output.mast_xml)
		mast_data_str = json.dumps (mast_data, indent=3)
		with open (output.mast_json, 'w') as hndl:
			hndl.write (mast_data_str)

		# randomize and look for cassettes
		cass_support_dict = {c['pattern']:[] for c in all_cassettes}
		uniq_cassettes = list (cass_support_dict.keys())

		for i in range (RANDOMIZATION_CNT):
			# randomize cassette
			shuff_mast_json = mcda.randomize_mast_json (mast_data)
			# search for every cassette in randomize data
			for u in uniq_cassettes:
				hits = mcda.find_cassettes_in_mast_json (u, shuff_mast_json,
					max_gap=MAX_CASS_GAP, max_pval=MAX_CASS_PVAL)
				# record the number of sequences with a hit
				seq_names = [h['seq_name'] for h in hits]
				cass_support_dict[u].append (len (set (seq_names)))

		# sort support list
		for k,v in cass_support_dict.items():
			cass_support_dict[k] = sorted (v)

		# calculate support
		for c in all_cassettes:
			pattern = c['pattern']
			c['randomizations'] = RANDOMIZATION_CNT
			rnd_list = [x / MEME_SEQ_CNT for x in cass_support_dict[pattern]]
			cass_freq = c['cutoff']
			# this will give the number of randomizations lower than it
			randomization_rank = bisect.bisect_left (rnd_list, float (cass_freq))
			randomization_support = randomization_rank / RANDOMIZATION_CNT
			c['ramdomization_rank'] = randomization_rank
			c['randomization_support'] = randomization_support

		mcda.write_csv_as_dicts (all_cassettes, output.cassettes_with_support,
			hdr_flds=CASSETTE_SUPPORT_FIELDS)



rule graph_cassettes:
	message: "Draw a diagram of the cassette relationship"
	input:
		filtered_cassettes=FILTERED_CASSETTES,
	output:
		cass_graph=CASSETTE_GRAPH,
	run:
		## Main:
		import networkx as nx
		import matplotlib.pyplot as plt

		# read in cassettes from previous step
		all_cassettes = mcda.read_csv_as_dicts (input.filtered_cassettes)

		# build graph
		patts = [r['pattern'] for r in all_cassettes]
		grf = mcda.graph_cassettes (patts)

		# save in various forms
		nx.write_graphml (grf, output.cass_graph)
		nx.write_gml (grf, output.cass_graph.replace ('.graphml', '.gml'))

		from networkx.readwrite import json_graph
		data = json_graph.node_link_data (grf)
		import json
		json.dumps (data, output.cass_graph.replace ('.graphml', '.json'))

		# build a list of frequencies
		freqs = {}
		for c in all_cassettes:
			freqs[c['pattern']] = c['freq']

		# draw & save
		nx.draw_circular (grf,
			arrows=True,
			labels={k:k for k in grf.nodes()},
			node_size=[float (freqs[n]) * 15 for n in grf],
		)
		plt.show (block=False)
		plt.savefig (output.cass_graph.replace ('.graphml', '.png'), format="PNG")

		mcda.save_graph_as_image (grf, freqs, CASSETTE_GRAPH_PIC)



rule mast_search_for_elements:
	input:
		filtered_cassettes=FILTERED_CASSETTES,
		all_exp_seqs=ALL_EXP_SEQS,
		all_control_seqs=ALL_CONTROL_SEQS,
		meme_results=MEME_RESULTS,
	output:
		EXP_MAST_RES_DIR,
		CONTROL_MAST_RES_DIR,
		EXP_MAST_RES,
		CONTROL_MAST_RES,
	run:
		## Constants:
		# build list of data to analyse & where results go
		ANALYSES = [
			(input.all_exp_seqs, EXP_MAST_RES_DIR),
			(input.all_control_seqs, CONTROL_MAST_RES_DIR),
		]

		for pth in glob (path.join (COMP_SEQ_WORK_DIR, '*.fasta')):
			seq_dir, seq_name, seq_ext = mcda.path_to_dir_name_ext (pth)
			ANALYSES.append ((pth, path.join (COMPARATIVE_MAST_RES_DIR, seq_name)))

		## Main:
		# calculate elements needed
		all_cassettes = mcda.read_csv_as_dicts (input.filtered_cassettes)
		all_motifs_seen = ''.join ([r['elements'] for r in all_cassettes])
		uniq_motifs_seen = set (all_motifs_seen)

		for a in ANALYSES:
			seq_file, output_dir = a
			mcda.run_mast (input.meme_results, seq_file, output_dir, motif_ids=uniq_motifs_seen, exe=MAST_EXE)



rule postprocess_mast_results:
	message: "Manipulate the MAST results into a more readable form"
	input:
		exp_mast_results=EXP_MAST_RES,
		control_mast_results=CONTROL_MAST_RES,
	output:
		exp_mast_json=EXP_MAST_JSON,
		control_mast_json=CONTROL_MAST_JSON,
	run:
		## Constants:
		# build list of data to analyse & where results go
		ANALYSES = [
			(input.exp_mast_results, output.exp_mast_json),
			(input.control_mast_results, output.control_mast_json),
		]

		for pth in glob (path.join (COMP_SEQ_WORK_DIR, '*.fasta')):
			seq_dir, seq_name, seq_ext = mcda.path_to_dir_name_ext (pth)
			res_path = path.join (COMPARATIVE_MAST_RES_DIR, seq_name,
				'mast.xml')
			json_pth = res_path.replace ('.xml', '.json')
			ANALYSES.append ((res_path, json_pth))

		## Main:
		# do actual analysis
		for f, j in ANALYSES:
			mcda.progress_msg ("Mapping %s to %s" % (f, j))
			mast_data = mcda.parse_mast_results (f)
			mast_json = json.dumps (mast_data, indent=3)
			with open (j, 'w') as hndl:
				hndl.write (mast_json)



rule list_cassettes:
	message: "Process mast results to list cassettes in total sequence sets"
	input:
		filtered_cassettes=FILTERED_CASSETTES,
		exp_mast_json=EXP_MAST_JSON,
		control_mast_json=CONTROL_MAST_JSON,
	output:
		exp_cass_all=EXP_CASS_ALL,
		control_cass_all=CONTROL_CASS_ALL,
		exp_cass_seqs=EXP_CASS_SEQS,
	run:
		#
		## Constants:
		# build list of data to analyse & where results go
		ANALYSES = [
			(input.exp_mast_json, output.exp_cass_all),
			(input.control_mast_json, output.control_cass_all),
		]

		for pth in glob (path.join (COMPARATIVE_MAST_RES_DIR, '*/mast.json')):
			mast_dir, mast_file, mast_ext = mcda.path_to_dir_name_ext (pth)
			par_dir = mast_dir.split ('/')[-1]
			summary_pth = path.join (COMPARATIVE_CASS_DIR, '%s.all.csv' % par_dir)
			ANALYSES.append ((pth, summary_pth))

		## Main:
		# do actual analysis
		uniq_cassettes = mcda.read_csv_as_dicts (input.filtered_cassettes)
		uniq_patterns = [r['pattern'] for r in uniq_cassettes]

		mcda.ensure_dir_exists (COMPARATIVE_CASS_DIR)

		hdrs = ['seq'] + uniq_patterns

		all_seqs = {}
		for j, c in ANALYSES:
			with open (j, 'r') as hndl:
				json_data = json.load (hndl)
			all_hits = []
			for u in uniq_patterns:
				# find the hits
				hits = mcda.find_cassettes_in_mast_json (u, json_data,
					max_gap=MAX_CASS_GAP, max_pval=MAX_CASS_PVAL)
				# add to the global list
				all_hits.extend (hits)
				# add to the per-sequence list
				for h in hits:
					seq_name = h['seq_name']
					r = all_seqs.get (seq_name, {'seq': seq_name})
					r[u] = h['combined_pval']
					all_seqs[seq_name] = r

			# combine pvalues for storage
			for h in all_hits:
				pvals = [str(x) for x in h['pvals']]
				h['pvals'] = ' '.join (pvals)

			# save the global hits
			mcda.write_csv_as_dicts (all_hits, c, hdr_flds=CASS_HDRS)
			# save the local (per-sequence) hits
			sorted_recs = sorted (list (all_seqs.values()), key=lambda x: x['seq'])
			mcda.write_csv_as_dicts (sorted_recs, output.exp_cass_seqs,
				hdr_flds=hdrs, rest_val=0)



rule count_cassettes:
	message: "Process cassette hit listings to a summary form"
	input:
		filtered_cassettes=FILTERED_CASSETTES,
		exp_cass_all=EXP_CASS_ALL,
		control_cass_all=CONTROL_CASS_ALL,
	output:
		exp_cass_summary=EXP_CASS_SUMMARY,
		control_cass_summary=CONTROL_CASS_SUMMARY,
	run:
		## Constants:
		# build list of data to analyse & where results go
		ANALYSES = [
			(input.exp_cass_all, output.exp_cass_summary),
			(input.control_cass_all, output.control_cass_summary),
		]

		for pth in glob (path.join (COMPARATIVE_CASS_DIR, '*.all.csv')):
			summary_pth = pth.replace ('.all.', '.summary.')
			ANALYSES.append ((pth, summary_pth))

		## Main:
		# do actual analysis
		uniq_cassettes = mcda.read_csv_as_dicts (input.filtered_cassettes)
		uniq_patterns = [r['pattern'] for r in uniq_cassettes]

		for a, s in ANALYSES:
			all_cass = mcda.read_csv_as_dicts (a)
			summary_recs = []
			for u in uniq_patterns:
				hits = [x['seq_name'] for x in all_cass if x['pattern'] == u]
				# NOTE: very important! we count the frequency (all hits) and the
				# 'freq_seq' or the frequency of sequences containing at least 1
				# example cassette
				summary_recs.append ({
					'pattern': u,
					'freq': len (hits),
					'freq_seq': len (set (hits)),
				})
			mcda.write_csv_as_dicts (summary_recs, s, CASS_SUMMARY_HDRS)



rule tabulate_cassette_counts:
	message: "Tabulate cassette hits into a single table"
	input:
		filtered_cassettes=FILTERED_CASSETTES,
		exp_cass_summary=EXP_CASS_SUMMARY,
		control_cass_summary=CONTROL_CASS_SUMMARY,
	output:
		overall_cass_table=OVERALL_CASS_TABLE,
	run:
		## Constants:
		# build list of data to analyse
		comp_summaries = [x for x in glob (path.join (COMPARATIVE_CASS_DIR,
			'*.summary.csv'))]
		all_summaries = [input.control_cass_summary, input.exp_cass_summary] + \
			comp_summaries

		## Main:
		# get list of all patterns
		uniq_cassettes = mcda.read_csv_as_dicts (input.filtered_cassettes)
		uniq_patterns = [r['pattern'] for r in uniq_cassettes]
		uniq_patterns.sort (key=mcda.sort_key)

		# get counts for all of these from summaries
		sources = []
		all_freqs = {}
		for f in all_summaries:
			src = mcda.get_file_name (f).replace ('.summary', '')
			sources.append (src)
			recs = mcda.read_csv_as_dicts (f)
			# NOTE: this uses the frequency of sequences containing a cassette
			# A sequence may contain more than one example of a cassette
			freqs = {r['pattern']:r['freq_seq'] for r in recs}
			all_freqs[src] = freqs

		# tabulate
		tab_recs = []
		for p in uniq_patterns:
			curr_rec = {'pattern': p}
			for s in sources:
				curr_rec[s] = all_freqs[s][p]
			tab_recs.append (curr_rec)

		# save
		hdrs = ['pattern'] + sources
		mcda.write_csv_as_dicts (tab_recs, output.overall_cass_table,
			hdrs)



rule make_summary_table:
	message: "Calculate enrichment and pvals for cassettes"
	input:
		overall_cass_table=OVERALL_CASS_TABLE,
		seq_cnts=SEQ_CNT_PTH,
	output:
		table_with_stats=OVERALL_CASS_TABLE_WITH_STATS,
		enrichment_table=CASS_ENRICHMENT_TABLE,
		count_table=CASS_COUNT_TABLE,
	run:
		## Constants:
		# define the form of headers for convenience
		CNT_FLD_TMPL = '%s_seq_cnt'
		CASS_FLD_TMPL = '%s_cass_cnt'
		FRAC_FLD_TMPL = '%s_frac'
		ENR_FLD_TMPL = '%s_enr'
		PVAL_FLD_TMPL = '%s_pval'
		QVAL_FLD_TMPL = '%s_qval'

		CNTRL_CNT_FLD = CNT_FLD_TMPL % 'control'
		CNTRL_CASS_FLD = CASS_FLD_TMPL % 'control'
		CNTRL_FRAC_FLD = FRAC_FLD_TMPL % 'control'

		EXP_CNT_FLD = CNT_FLD_TMPL % 'experimental'
		EXP_CASS_FLD = CASS_FLD_TMPL % 'experimental'
		EXP_FRAC_FLD = FRAC_FLD_TMPL % 'experimental'
		EXP_ENR_FLD = ENR_FLD_TMPL % 'experimental'
		EXP_PVAL_FLD = PVAL_FLD_TMPL % 'experimental'
		EXP_QVAL_FLD = QVAL_FLD_TMPL % 'experimental'

		## Main:
		# get counts of all types of sequences
		seq_cnt_recs = mcda.read_csv_as_dicts (input.seq_cnts)
		seq_cnts = {}
		for r in seq_cnt_recs:
			name = r['name'].replace ('all-', '')
			seq_cnts[name] = r['count']

		# build table of comparative results
		# this should be the the # seqs, # hits, calculated fraction for all
		# plus the enrichment & pval for experimental seqs

		# read the number of cassettes found in each seq set
		cass_summaries = mcda.read_csv_as_dicts (input.overall_cass_table)

		# build list of fields for output (order is important)
		fld_list = ['cassette',
			CNTRL_CNT_FLD, CNTRL_CASS_FLD, CNTRL_FRAC_FLD,
			EXP_CNT_FLD, EXP_CASS_FLD, EXP_FRAC_FLD, EXP_ENR_FLD, EXP_PVAL_FLD,
			EXP_QVAL_FLD
		]
		for k in list (seq_cnts.keys()):
			if k not in ['control', 'experimental']:
				fld_list.extend ([
					CNT_FLD_TMPL % k,
					CASS_FLD_TMPL % k,
					FRAC_FLD_TMPL % k,
					ENR_FLD_TMPL % k,
					PVAL_FLD_TMPL % k,
				])
		#mcda.prettyprint (fld_list)

		# init records to output with seq & cassette counts & patterns
		cass_details = []
		for r in cass_summaries:
			# give it all the fields it needs
			new_rec = {}
			for f in fld_list:
				new_rec[f] = None

			# give it the pattern/cassette & cassette counts
			for k in list (r.keys()):
				if k == 'pattern':
					new_rec['cassette'] = r['pattern']
				else:
					new_rec[CASS_FLD_TMPL % k] = int (r[k])

			for k in list (seq_cnts.keys()):
				new_rec[CNT_FLD_TMPL % k] = int (seq_cnts[k])

			cass_details.append (new_rec)
		# mcda.prettyprint (cass_details[:2])

		# annotate with seq counts & generate fraction
		for src, cnt in seq_cnts.items():
			cnt_fld = CNT_FLD_TMPL % src
			cass_fld = CASS_FLD_TMPL % src
			frac_fld = FRAC_FLD_TMPL % src

			for r in cass_details:
				try:
					r[cnt_fld] = int (cnt)
					r[frac_fld] = r[cass_fld] / r[cnt_fld]
				except:
					print (cnt_fld, r[cnt_fld])
					print (cass_fld, r[cass_fld])
					print (r)
					raise

		# calculate enrichment & pval
		for src in seq_cnts.keys():
			if src != 'control':
				cnt_fld = CNT_FLD_TMPL % src
				cass_fld = CASS_FLD_TMPL % src
				frac_fld = FRAC_FLD_TMPL % src
				enr_fld = ENR_FLD_TMPL % src
				pval_fld = PVAL_FLD_TMPL % src

				for r in cass_details:
					r[enr_fld] = r[frac_fld] / r[CNTRL_FRAC_FLD]
					r[pval_fld] = mcda.enrichment_pval_via_betabinomial (r[cass_fld],
						r[cnt_fld], r[CNTRL_CASS_FLD], r[CNTRL_CNT_FLD])

		# calculate qval?
		for c in ('experimental', 'experimental-nd'):
			in_col = "%s_pval" % c
			out_col = "%s_qval" % c
			all_pvals = [r[in_col] for r in cass_details]
			adj_pvals = multipletests (all_pvals, alpha=0.05, method='fdr_by')
			qvals = adj_pvals[1]
			for i in range (len (cass_details)):
				cass_details[i][out_col] = qvals[i]

		# save this all
		mcda.write_csv_as_dicts (cass_details, output.table_with_stats,
			fld_list)
		#print (fld_list)

		# make sub-tables for presentation purposes
		# process data to make it better formatted
		for f in fld_list:
			if f.endswith ('_frac') or f.endswith ('_enr'):
				for r in cass_details:
					r[f] = '%0.3f' % r[f]
			if f.endswith ('_pval') or f.endswith ('_qval'):
				for r in cass_details:
					r[f] = '%0.3e' % r[f]

		# first one with the seqs & counts & frac
		cnt_fld_list = [f for f in fld_list if f.split('_')[-1] not in ('pval',
			'enr')]
		mcda.write_csv_as_dicts (cass_details, output.count_table,
			cnt_fld_list)

		# then one with the enr and pvals
		enr_fld_list = [f for f in fld_list if f.split('_')[-1] not in ('cnt')]
		mcda.write_csv_as_dicts (cass_details, output.enrichment_table,
			enr_fld_list)



rule extract_exemplars:
	message: "Extract best examples of cassette sequences from experimental data"
	input:
		all_exp_seqs=ALL_EXP_SEQS,
		all_exp_cass=EXP_CASS_ALL,
	output:
		exemplar_cass=EXEMPLAR_CASS,
		exemplar_seqs=EXEMPLAR_SEQS,
	run:
		# read in & index all seqs
		with open (input.all_exp_seqs, 'r') as hndl:
			exp_seqs_dict = SeqIO.to_dict (SeqIO.parse (hndl, 'fasta'))

		# read in all cassette hits & gather list of all cassettes
		patt_list = []
		all_exp_cass = mcda.read_csv_as_dicts (input.all_exp_cass)
		for c in all_exp_cass:
			c['start'] = int (c['start'])
			c['stop'] = int (c['stop'])
			c['len'] = int (c['len'])
			c['max_gap'] = int (c['max_gap'])
			c['combined_pval'] = float (c['combined_pval'])
			if c['pattern'] not in patt_list:
				patt_list.append (c['pattern'])

		# get n best examples:
		exemplar_data = []
		for p in patt_list:
			all_cass = [c for c in all_exp_cass if c['pattern'] == p]
			sorted_cass = sorted (all_cass, key=lambda x: x['combined_pval'])
			exemplar_data.extend (sorted_cass[:NUM_EXEMPLARS])
		mcda.write_csv_as_dicts (exemplar_data, output.exemplar_cass,
			CASS_HDRS)

		# now extract seqs
		exemplar_seqs = []
		for r in exemplar_data:
			patt = r['pattern']
			patt_str = patt.replace ('+', 'p').replace ('-', 'm')
			seq = exp_seqs_dict[r['seq_name']]
			sub_seq = seq[r['start']-1:r['stop']]
			if r['dir'] == '+':
				dir_str = 'p'
			else:
				dir_str = 'm'
				sub_seq = sub_seq.reverse_complement (id=True, name=True)
			seq.name = seq.id = "%s_%s_%s" % (r['seq_name'], patt_str, dir_str)
			sub_seq.description = \
				"cassette %s, strand %s, positions %s-%s, length %s, combined p-value %s" % (
					r['pattern'], r['dir'], r['start'], r['stop'], r['len'],
					r['combined_pval']
				)
			exemplar_seqs.append (sub_seq)

		# record it
		mcda.write_seqs (exemplar_seqs, EXEMPLAR_SEQS, 'fasta')



rule report_results:
	message: "Produce a nicely formatted report, summarising what has been found"
	input:
		EXP_SEQ_DATA,
		CNTRL_SEQ_DATA,
		SEQ_CNT_PTH,
		FILTERED_CASSETTES,
		MOTIF_SUMMARY,
		OVERALL_CASS_TABLE,
		OVERALL_CASS_TABLE_WITH_STATS,
		EXEMPLAR_SEQS,
		EXEMPLAR_CASS,
		EXP_CASS_SEQS,
		CASSETTE_GRAPH,
		CASSETTES_WITH_SUPPORT,
	output:
		html=REPORT_PTH
	run:
		# make some data for inclusion
		import datetime
		curr_time = datetime.datetime.now().isoformat()

		# need to draw all the logos
		mtf_info = mcda.read_csv_as_dicts (MOTIF_SUMMARY)
		logo_rst_tmpl = """
.. figure:: ../%(logo_path)s

	Motif %(name)s (found in %(num_seq)s sequences, width %(width)s bases, e-value %(eval)s)

		"""
		logo_rst_incls = [logo_rst_tmpl % {
				'logo_path': m['logo_path'],
				'name': m['name'],
				'num_seq': m['sites'],
				'width': m['width'],
				'eval': m['evalue'],
			} for m in mtf_info
		]
		logo_rst_incl_txt = '\n\n'.join (logo_rst_incls)
		# NOTE: the path to the logo file is a complete hack. Has to be a more
		# more robust way of doing this

		snakemake.utils.report ("""
=======================
Motif cassette analysis
=======================

**Date:** {curr_time}

Input data
----------

* The experimental data was ``{EXP_SEQ_DATA}``.

* The control data was ``{CNTRL_SEQ_DATA}``.

* The remainder of the control and experimental sequence sets not used for motif discovery is saved as "-nd".

* The number of sequences in all files (including those to be used for comparative purposes) were:

	.. csv-table:: All input sequences & counts
		:file: {SEQ_CNT_PTH}


Analysis
--------

Motif discovery
~~~~~~~~~~~~~~~

Motifs were discovered in a discriminative analysis (experimental vs. control sequence data) with the following parameters:

* The executables used were::

	{MEME_EXE}
	{MAST_EXE}

* The number of motifs searched was {NUM_MOTIFS}.

* The size of motif searched for was {MIN_MOTIF_WIDTH} to {MAX_MOTIF_WIDTH}.

* The first {MEME_SEQ_CNT} sequences of experimental and control sequences was used in the search.

* The following motifs were found:

	.. csv-table:: All motifs found
		:file: {MOTIF_SUMMARY}

%(logo_text)s

Cassette mining
~~~~~~~~~~~~~~~

* The discovery (MEME) results were examined for repeating motif patterns, using a motif p-value threshold of {MAX_MOTIF_PVAL}.

* The following cassettes were found:

	.. csv-table:: All cassettes found
		:file: {FILTERED_CASSETTES}

* In a graph, where node height is proportional to frequency and cassettes are connected to the next largest cassette that contains them:

	.. image:: ../{CASSETTE_GRAPH_PIC}

* By randomizing the sequence of motifs and searching for the same cassettes, a support score was calculated:

	.. csv-table:: Cassette support scores
		:file: {CASSETTES_WITH_SUPPORT}



Cassette detection in all sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* MAST (executable {MAST_EXE}) was used to search for motifs in all sequence datasets.

* MAST results were mined for cassettes:

	* The maximum gap allows between motifs was ('None' means unused): {MAX_CASS_GAP}

	* The maximum p-value threshold for cassettes (based on combined p-value of constituent motifs) was: {MAX_CASS_PVAL}

* The raw number of sequences with specific cassettes each dataset was:

	.. csv-table:: Cassettes in all sequence sets
		:file: {OVERALL_CASS_TABLE}


Cassette enrichment & exemplar sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* For each cassette and every sequence dataset, these are the number of sequences in the dataset, the raw number and fraction of sequences carrying at least one cassette:

	.. csv-table:: Cassette counts
		:file: {CASS_COUNT_TABLE}

* From this, we can calculate an enrichment value (over the control data) and an associated p-value:

	.. csv-table:: Cassette enrichment
		:file: {CASS_ENRICHMENT_TABLE}

* {NUM_EXEMPLARS} exemplar sequences (the sequence for a given cassette, taken from with sequence with the best p-value) were extracted. These are stored in ``{EXEMPLAR_SEQS}``. A summary of each cassette and the corresponding best sequences can be found in ``{EXEMPLAR_CASS}``:

	.. csv-table:: Best sequence for each cassette
		:file: {EXEMPLAR_CASS}

		""" % {'logo_text': logo_rst_incl_txt},
			output.html,
			defaultenc='utf8',
			stylesheet=path.join(REPORT_DATA_DIR, 'report.css'),
			files=input
		)



### END ###
