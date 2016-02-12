"""
Perform custom data validation that is not covered by the root script.

This script should be adapted (or commented out) for custom workflows.
"""


def get_seq_names (pth):
	"""
	Return a set of the names of sequenecs found in a file.
	"""
	return [s.name for s in SeqIO.parse (open (pth, 'r'), 'fasta')]

# base file checks that experimental and control seqs don't overlap
# the downreg and upreg sequences respectively
exp_names = get_seq_names (snakemake.input.exp_seq_data)
cntrl_names = get_seq_names (snakemake.input.cntrl_seq_data)

# these should all be a subset of the experimental seqs
phago_names = get_seq_names ('data/seqs/comparative/phago.fasta')
transmemb_names = get_seq_names ('data/seqs/comparative/transmemb.fasta')
immune_names = get_seq_names ('data/seqs/comparative/immune.fasta')

assert exp_names.issuperset (phago_names), "some phago not in down"
assert exp_names.issuperset (transmemb_names), "some transmemb not in down"
assert exp_names.issuperset (immune_names), "some immune not in down"

# these should not overlap with the up- or down-reg
unchanged_names = get_seq_names ('data/seqs/comparative/unchanged.fasta')

assert exp_names.isdisjoint (unchanged_names), "some unchanged in down"
assert cntrl_names.isdisjoint (unchanged_names), "some unchanged in up"

### END ###
