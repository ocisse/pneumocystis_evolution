###
# makeblastdb -in tcdb.fasta -input_type fasta -dbtype prot -out tcdb.db
rule all:
	input:
		"../data/processed/tcdb/all_species_tcdb.csv"

rule blastp:
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
	threads: 2
	output:
		pcanA="../data/processed/tcdb/PcanA.out",
		pcanCk1="../data/processed/tcdb/PcanCk1.out",
		pcanCk2="../data/processed/tcdb/PcanCk2.out",
		pc="../data/processed/tcdb/Pcar.out",
		pj="../data/processed/tcdb/Pjir.out",
		pmac="../data/processed/tcdb/Pmac.out",
		pm="../data/processed/tcdb/Pmur.out",
		po="../data/processed/tcdb/Pory.out",
		pw="../data/processed/tcdb/Pwk.out"
	run:
		shell("blastp -db ../data/external/tcdb.db -out {output.pcanA} -query {input.pcanA} -outfmt 7 -evalue 1e-5 -num_threads {threads}")
		shell("blastp -db ../data/external/tcdb.db -out {output.pcanCk1} -query {input.pcanCk1} -outfmt 7 -evalue 1e-5 -num_threads {threads}")
		shell("blastp -db ../data/external/tcdb.db -out {output.pcanCk2} -query {input.pcanCk2} -outfmt 7 -evalue 1e-5 -num_threads {threads}")
		shell("blastp -db ../data/external/tcdb.db -out {output.pc} -query {input.pc} -outfmt 7 -evalue 1e-5 -num_threads {threads}")
		shell("blastp -db ../data/external/tcdb.db -out {output.pj} -query {input.pj} -outfmt 7 -evalue 1e-5 -num_threads {threads}")
		shell("blastp -db ../data/external/tcdb.db -out {output.pmac} -query {input.pmac} -outfmt 7 -evalue 1e-5 -num_threads {threads}")
		shell("blastp -db ../data/external/tcdb.db -out {output.pm} -query {input.pm} -outfmt 7 -evalue 1e-5 -num_threads {threads}")
		shell("blastp -db ../data/external/tcdb.db -out {output.po} -query {input.po} -outfmt 7 -evalue 1e-5 -num_threads {threads}")
		shell("blastp -db ../data/external/tcdb.db -out {output.pw} -query {input.pw} -outfmt 7 -evalue 1e-5 -num_threads {threads}")

rule parse_blastp:
	input:
		pj=rules.blastp.output.pj,pmac=rules.blastp.output.pmac,po=rules.blastp.output.po,
		pcanCk1=rules.blastp.output.pcanCk1,pcanCk2=rules.blastp.output.pcanCk2,pcanA=rules.blastp.output.pcanA,
		pc=rules.blastp.output.pc,pm=rules.blastp.output.pm,pw=rules.blastp.output.pw
	output:
		pcanA="../data/processed/tcdb/PcanA.tp",
		pcanCk1="../data/processed/tcdb/PcanCk1.tp",
		pcanCk2="../data/processed/tcdb/PcanCk2.tp",
		pc="../data/processed/tcdb/Pcar.tp",
		pj="../data/processed/tcdb/Pjir.tp",
		pmac="../data/processed/tcdb/Pmac.tp",
		pm="../data/processed/tcdb/Pmur.tp",
		po="../data/processed/tcdb/Pory.tp",
		pw="../data/processed/tcdb/Pwk.tp"
	run:
		shell("perl scripts/parse_blastp.pl {input.pj} > {output.pj}")
		shell("perl scripts/parse_blastp.pl {input.pmac} > {output.pmac}")
		shell("perl scripts/parse_blastp.pl {input.po} > {output.po}")
		shell("perl scripts/parse_blastp.pl {input.pcanCk1} > {output.pcanCk1}")
		shell("perl scripts/parse_blastp.pl {input.pcanA} > {output.pcanA}")
		shell("perl scripts/parse_blastp.pl {input.pcanCk2} > {output.pcanCk2}")
		shell("perl scripts/parse_blastp.pl {input.pc} > {output.pc}")
		shell("perl scripts/parse_blastp.pl {input.pm} > {output.pm}")
		shell("perl scripts/parse_blastp.pl {input.pw} > {output.pw}")

rule merge_in_on_table:
	input:
		t="../data/external/tcdb.mapping.txt",
		pj=rules.parse_blastp.output.pj,pmac=rules.parse_blastp.output.pmac,po=rules.parse_blastp.output.po,
		pcanCk1=rules.parse_blastp.output.pcanCk1,pcanCk2=rules.parse_blastp.output.pcanCk2,pcanA=rules.parse_blastp.output.pcanA,
		pc=rules.parse_blastp.output.pc,pm=rules.parse_blastp.output.pm,pw=rules.parse_blastp.output.pw
	output:
		"../data/processed/tcdb/all_species_tcdb.csv"
	shell:
		"perl scripts/merge_in_on_table.pl {input.pj} {input.pmac} {input.po} {input.pcanCk1} {input.pcanCk2} {input.pcanA} {input.pc} "
		"{input.pm} {input.pw} {input.t} > {output}"









