# don't use the version 33 because it lacks many important domains such as msgs
# 32 also is off for msg at least the main domain

rule search:
	input:
		pcanA="/Users/cisseoh/Documents/TMP/SP_DIV/temp_docs_191225/Orphans_wgs/data/raw/Proteomes/PcanA.faa",
		pcanCk1="/Users/cisseoh/Documents/TMP/SP_DIV/temp_docs_191225/Orphans_wgs/data/raw/Proteomes/PcanCk1.faa",
		pcanCk2="/Users/cisseoh/Documents/TMP/SP_DIV/temp_docs_191225/Orphans_wgs/data/raw/Proteomes/PcanCk2.faa",
		pc="/Users/cisseoh/Documents/TMP/SP_DIV/temp_docs_191225/Orphans_wgs/data/raw/Proteomes/Pcar.faa",
		pj="/Users/cisseoh/Documents/TMP/SP_DIV/temp_docs_191225/Orphans_wgs/data/raw/Proteomes/Pjir.faa",
		pmac="/Users/cisseoh/Documents/TMP/SP_DIV/temp_docs_191225/Orphans_wgs/data/raw/Proteomes/Pmac.faa",
		pm="/Users/cisseoh/Documents/TMP/SP_DIV/temp_docs_191225/Orphans_wgs/data/raw/Proteomes/Pmur.faa",
		po="/Users/cisseoh/Documents/TMP/SP_DIV/temp_docs_191225/Orphans_wgs/data/raw/Proteomes/Pory.faa",
		pw="/Users/cisseoh/Documents/TMP/SP_DIV/temp_docs_191225/Orphans_wgs/data/raw/Proteomes/Pwk.faa",
		hmm="../data/external/Pfam/31.0/Pfam-A.hmm"
	output:
		pcanA="../data/processed/pfam/31.0/PcanA.out",
		pcanCk1="../data/processed/pfam/31.0/PcanCk1.out",
		pcanCk2="../data/processed/pfam/31.0/PcanCk2.out",
		pc="../data/processed/pfam/31.0/Pcar.out",
		pj="../data/processed/pfam/31.0/Pjir.out",
		pmac="../data/processed/pfam/31.0/Pmac.out",
		pm="../data/processed/pfam/31.0/Pmur.out",
		po="../data/processed/pfam/31.0/Pory.out",
		pw="../data/processed/pfam/31.0/Pwk.out"

	threads: 2

	params: e="1e-5",

	run:
		shell("hmmsearch -E {params.e} --noali --domtblout {output.pcanA} --cpu {threads} {input.hmm} {input.pcanA}")
		shell("hmmsearch -E {params.e} --noali --domtblout {output.pcanCk1} --cpu {threads} {input.hmm} {input.pcanCk1}")
		shell("hmmsearch -E {params.e} --noali --domtblout {output.pcanCk2} --cpu {threads} {input.hmm} {input.pcanCk2}")
		shell("hmmsearch -E {params.e} --noali --domtblout {output.pc} --cpu {threads} {input.hmm} {input.pc} ")
		shell("hmmsearch -E {params.e} --noali --domtblout {output.pj} --cpu {threads} {input.hmm} {input.pj} ")
		shell("hmmsearch -E {params.e} --noali --domtblout {output.pmac} --cpu {threads} {input.hmm} {input.pmac} ")
		shell("hmmsearch -E {params.e} --noali --domtblout {output.pm} --cpu {threads} {input.hmm} {input.pm} ")
		shell("hmmsearch -E {params.e} --noali --domtblout {output.po} --cpu {threads} {input.hmm} {input.po} ")
		shell("hmmsearch -E {params.e} --noali --domtblout {output.pw} --cpu {threads} {input.hmm} {input.pw} ")

# I going to count one domain per protein 
rule parse_results:
	input:
		pcanA="../data/processed/pfam/31.0/PcanA.out",
		pcanCk1="../data/processed/pfam/31.0/PcanCk1.out",
		pcanCk2="../data/processed/pfam/31.0/PcanCk2.out",
		pc="../data/processed/pfam/31.0/Pcar.out",
		pj="../data/processed/pfam/31.0/Pjir.out",
		pmac="../data/processed/pfam/31.0/Pmac.out",
		pm="../data/processed/pfam/31.0/Pmur.out",
		po="../data/processed/pfam/31.0/Pory.out",
		pw="../data/processed/pfam/31.0/Pwk.out"

	output:
		pcanA=temp("../data/processed/pfam/31.0/PcanA.p"),
		pcanCk1=temp("../data/processed/pfam/31.0/PcanCk1.p"),
		pcanCk2=temp("../data/processed/pfam/31.0/PcanCk2.p"),
		pc=temp("../data/processed/pfam/31.0/Pcar.p"),
		pj=temp("../data/processed/pfam/31.0/Pjir.p"),
		pmac=temp("../data/processed/pfam/31.0/Pmac.p"),
		pm=temp("../data/processed/pfam/31.0/Pmur.p"),
		po=temp("../data/processed/pfam/31.0/Pory.p"),
		pw=temp("../data/processed/pfam/31.0/Pwk.p"),
		merge="../data/processed/pfam/31.0/all_pneumo_pfam_merged.csv"
	params:
		c="1e-5"
	run:
		shell("perl scripts/parse_hmmsearch.pl --cutoff {params.c} {input.pj} > {output.pj}")
		shell("perl scripts/parse_hmmsearch.pl --cutoff {params.c} {input.pmac} > {output.pmac}")
		shell("perl scripts/parse_hmmsearch.pl --cutoff {params.c} {input.po} > {output.po}")
		shell("perl scripts/parse_hmmsearch.pl --cutoff {params.c} {input.pcanCk1} > {output.pcanCk1}")
		shell("perl scripts/parse_hmmsearch.pl --cutoff {params.c} {input.pcanCk2} > {output.pcanCk2}")
		shell("perl scripts/parse_hmmsearch.pl --cutoff {params.c} {input.pcanA} > {output.pcanA}")
		shell("perl scripts/parse_hmmsearch.pl --cutoff {params.c} {input.pm} > {output.pm}")
		shell("perl scripts/parse_hmmsearch.pl --cutoff {params.c} {input.pc} > {output.pc}")
		shell("perl scripts/parse_hmmsearch.pl --cutoff {params.c} {input.pw} > {output.pw}")
		shell("perl scripts/merge_hmmsearch_parsed.pl {output.pj} {output.pmac} {output.po} {output.pcanCk1} {output.pcanCk2} {output.pcanA} {output.pc} {output.pm} {output.pw} > {output.merge}")

# cleaned version as all_pneumo_pfam_merged_clean.xslx
# version pfam33.0 don't have importat domains revert to 32.0 (DON'T USE 33.0)

# msg are not significant anymore but the counts are more conform with Ma et al. 2016
# focus on domains that could have link to host specificity, no need to shwo domains that are housekeeping

