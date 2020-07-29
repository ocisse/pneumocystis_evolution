### rule all

rule build_bt2_idx:
	input:
		rat="../../data/external/db/rat/GCF_000001895.5_Rnor_6.0_genomic.fna", 
		pc="../../data/external/db/pcarinii/GCA_001477545.1_Pneu_cari_B80_V3_genomic.fasta",

	output:
		rat="../../data/external/db/rat/rat",
		pc="../../data/external/db/pcarinii/pcar"

	shell:
		"echo > {output} && bowtie2-build {input} {output}"

## HYBRID assembly using all Pwk datasets  - regular Spades (individual assemblies not shown here)

rule remove_rat_from_raw_reads:
	input:
		ref=rules.build_bt2_idx.output.rat,
		pwrnf="../../data/raw/Pwk/reads/Reads/Pwk/RN/Pw_RNA_1.fastq", # NCBI accession: SRR11794282
		pwrnr="../../data/raw/Pwk/reads/Reads/Pwk/RN/Pw_RNA_2.fastq", # NCBI accession: SRR11794282
		pw1af="../../data/raw/Pwk/reads/Reads/Pwk/Pw1A/Pw1A_USD16081096_HGC3CALXX_L1_1.fq", # SRR11794285
		pw1ar="../../data/raw/Pwk/reads/Reads/Pwk/Pw1A/Pw1A_USD16081096_HGC3CALXX_L1_1.fq",	# SRR11794285
		pw2af="../../data/raw/Pwk/reads/Reads/Pwk/Pw2A/Pw2A_DSW31654_HF5KLALXX_L7_1.fq", # SRR11794283 
		pw2ar="../../data/raw/Pwk/reads/Reads/Pwk/Pw2A/Pw2A_DSW31654_HF5KLALXX_L7_2.fq", # SRR11794283
		pw3af="../../data/raw/Pwk/reads/Reads/Pwk/Pw3A/Pw3A_USD16081097_HGC3CALXX_L2_1.fq", # SRR11794284
		pw3ar="../../data/raw/Pwk/reads/Reads/Pwk/Pw3A/Pw3A_USD16081097_HGC3CALXX_L2_2.fq", # SRR11794284
		pwC1f="../../data/raw/Pwk/reads/Reads/Pwk/PwC1/PwC1_DSW31653_HF5KLALXX_L7_1.fq", # SRR11794286
		pwC1r="../../data/raw/Pwk/reads/Reads/Pwk/PwC1/PwC1_DSW31653_HF5KLALXX_L7_2.fq", # SRR11794286

	output:
		m="../../data/processed/cleaning/Pwk/allreads_mapped2rat.sam",
		u="../../data/processed/cleaning/Pwk/allreads_unmapped2rat"
		
	conda:
		"envs/Pwk_config.yaml"

	threads: 12

	run:
		shell("echo > {output.u}")
		shell("bowtie2 -x {input.ref} -1 {input.pw1af} {input.pw2af} {input.pw3af} {input.pwC1f} "
			  "-2 {input.pw1ar} {input.pw2ar} {input.pw3ar} {input.pwC1r} -S {output.m}  --threads {threads} "
			  "--very-sensitive-local --un {output.u}")

### rename and gzip

rule remove_Pcarinii_from_raw_reads:
	input:
			ref=rules.build_bt2_idx.pc,
			r1="../../data/processed/cleaning/allreads_unmapped2rat.1.fq.gz",
			r2="../../data/processed/cleaning/allreads_unmapped2rat.2.fq.gz",
	output:
			m="../../data/processed/cleaning/Pwk/allreads_mapped2rat_pc.sam",
			u="../../data/processed/cleaning/Pwk/allreads_unmapped2rat_unmap2pc"
	conda:
		"envs/Pwk_config.yaml"

	threads: 12
	params: a="-D 5 -R 1 -N 0 -L 32 -i S,0,2.50 --end-to-end"

	run:
		shell("echo > {output.u}")
		shell("bowtie2 -x {input.ref} -1 {input.r1} -2 {input.r2} -S {output.m} --threads {threads} {params.a} --un {output.un}")


run assembly:
	input:
		r1="../../data/processed/cleaning/Pwk/allreads_unmapped2rat_unmap2pc.1.fq.gz",
		r2="../../data/processed/cleaning/Pwk/allreads_unmapped2rat_unmap2pc.2.fq.gz",

	threads: 24
	
	resources:
        mem_mb=250000

	output:
			directory("../../data/processed/assembly/Pwk/ass_v1")
	shell:
		"spades.py -t {threads} -1 {input.r1} -2 {input.r2} -o {output}"


run assem_cleanup1:
	input:
		rat=rules.build_bt2_idx.input.rat,
		s="../../data/processed/assembly/Pwk/ass_v1/scaffolds.fasta"

	output:
		rat="../../data/processed/assembly/Pwk/ass_v1/scaffolds_2rat.psl",
		rathits=temp("../../data/processed/assembly/Pwk/ass_v1/tmp1"),
		s="../../data/processed/assembly/Pwk/ass_v1/scaffolds.f1.fas"

	params: l="500"
	run:
		shell("blat {input.rat} {input.s} {output.rat}")
		# remove everything that map - becareful, might contains ambigous sequences + seq < 500
		shell("perl filter_psl {output.rat} | cut -f 1 | sort -u > {output.rathits}")
		shell("seqtk subseq {input.s} {output.rathits} -l 500 > ") 

# Quick check: 
run augustus_1:
	input:
		rules.assem_cleanup1.output.s
	output:
		"../../data/processed/assembly/Pwk/ass_v1/scaffolds.f1.gff"

	run:
		shell("augustus --gff3=on --species=pneumocystis {input} > {output}")
		shell("getAnnoFasta.pl {output} --seqfile={input}")

#~/utils/augustus-3.2.1/bin/augustus --gff3=on --species=pneumocystis scaffolds.f2.fas > scaffolds.f2.annot.gff
#perl ~/utils/augustus-3.2.1/scripts/getAnnoFasta.pl scaffolds.f2.annot.gff --seqfile=scaffolds.f2.fas

# manual examination
# I know P.wk is closer to P. murina than Pc
#cp NCBI/GCF_001477545.1_Pneu_cari_B80_V3_protein.faa Pc.faa
#cp NCBI/GCF_000349005.1_Pneumo_murina_B123_V2_protein.faa Pm.faa 

# change the here
#perl -pi.old -E 's/>/>Pc_/' Pc.faa
#perl -pi.old -E 's/>/>Pm_/' Pm.faa
#cat Pc.faa Pm.faa > db.faa
#~/utils/oldBlastall/blast-2.2.11/bin/blastall -p blastp -d db.faa -i scaffolds.f2.annot.aa -v 1 -b 1 -e 1e-20 -m 8 -o comp.bl &
# grep genes that hits first Pm 
#grep 'Pm_XP' comp.bl | cut -f 1 | sort -u > geneBhPm.txt

# extract scafs that contains those genes 
#cat geneBhPm.txt | while read f; do grep $f scaffolds.f2.annot.gff | grep 'transcript' >> tmp.2.txt;done
#cat tmp.2.txt | cut -f 1 | sort -u > tmp.3
#perl  ~/utils/retrieveFasta.pl tmp.3 scaffolds.f2.fas > scaffolds.f3.fas
# removed NODE_25196_length_513_cov_3440.89 by hand # looks like Pc

# quick megablast if - too similar to Pc
#~/utils/oldBlastall/blast-2.2.11/bin/megablast -i scaffolds.f3.fas \
#-d NCBI/GCA_001477545.1_Pneu_cari_B80_V3_genomic.fasta -W 1000 -v 1 -b 1 -e 1e-150 -m 8 -o scaffolds.f3.fas.2Pc.mgbl

#cat scaffolds.f3.fas.2Pc.mgbl | cut -f 1 | sort -u > x1
#grep '^>' scaffolds.f3.fas | cut -f 2 -d '>' | sort -u > x2
#diff -y x2 x1 | grep '<' | cut -f 1 > x3
#perl ~/utils/retrieveFasta.pl x3 scaffolds.f3.fas > scaffolds.f4.fas # 7.2 Mb genome

# verification via phylogeny 
#~/utils/augustus-3.2.1/bin/augustus --gff3=on --species=pneumocystis scaffolds.f4.fas > scaffolds.f4.annot.gff
#perl ~/utils/augustus-3.2.1/scripts/getAnnoFasta.pl scaffolds.f4.annot.gff --seqfile=scaffolds.f4.fas

#perl ~/utils/proteinortho_v2.3.0.pl -a=5 \
#scaffolds.f4.annot.aa \
#GCF_001477545.1_Pneu_cari_B80_V3_protein.faa \
#GCF_001477535.1_Pneu_jiro_RU7_V2_protein.faa \
#GCF_000349005.1_Pneumo_murina_B123_V2_protein.faa \
#> out

# generate csv for R
#./generate_csv.pl out

# 674 orthologs are missing in the current Pwk genome assembly
# try to retrieve scafs in which these genes are found 
#grep -v '#' out | grep '^\*' | cut -f 2,3,4 | grep -v '\*' | cut -f 3 > Pm_Conserved_missing_in_f4.txt
#perl  ~/utils/retrieveFasta.pl Pm_Conserved_missing_in_f4.txt \
#/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_000349005.1_Pneumo_murina_B123_V2_protein.faa \
#> Pm_Conserved_missing_in_f4.faa

#~/utils/oldBlastall/blast-2.2.11/bin/blastall -p tblastn \
#-d scaffolds.f3.fas \
#-i Pm_Conserved_missing_in_f4.faa -e 1e-100 -v 1 -b 1 -m 8 -o Pm_Conserved_missing_in_f4.tbln

#cat Pm_Conserved_missing_in_f4.tbln | cut -f 2 | sort -u > tmp.1
#perl ~/utils/retrieveFasta.pl tmp.1 scaffolds.f3.fas > tmp.1.fa
# find to
#~/utils/oldBlastall/blast-2.2.11/bin/megablast -i tmp.1.fa \
#-d NCBI/Pneumocystis_Pj_Pc_Pm_concat_WGS.fasta -W 1000 -v 1 -b 1 -e 1e-150 -m 8 -o tmp.1.fa.mgbl

# all those map to Pc, so it;s sure thpse are carinii = > remove them form tmp.1.fa
#cat tmp.1.fa.mgbl | cut -f 1 | sort -u > x1
#grep '^>' tmp.1.fa | cut -f 2 -d '>' | sort -u > x2
#diff x2 x1 | grep '<' | cut -f 2 -d '<' | sort -u > x3
#perl ~/utils/retrieveFasta.pl x3 tmp.1.fa > tmp.2.fa

#cat scaffolds.f4.fas tmp.2.fa > scaffolds.f5.fas # 13 Mb genome size

# annot
#~/utils/augustus-3.2.1/bin/augustus --gff3=on --species=pneumocystis scaffolds.f5.fas > scaffolds.f5.annot.gff
#perl ~/utils/augustus-3.2.1/scripts/getAnnoFasta.pl scaffolds.f5.annot.gff --seqfile=scaffolds.f5.fas

#perl ~/utils/proteinortho_v2.3.0.pl -a=5 -e=1e-4 \
#scaffolds.f5.annot.aa \
#/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_001477545.1_Pneu_cari_B80_V3_protein.faa \
#/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_001477535.1_Pneu_jiro_RU7_V2_protein.faa \
#/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_000349005.1_Pneumo_murina_B123_V2_protein.faa \
#> out2

#./generate_csv.pl out2

# this script ~/utils/proteinortho_v2.3.0.pl is leaving most of the genes behind 
# do it my self - parse the bla files to look for hits for genes missing in Pwk
# plus Pwk not found in other pneumo

#./update_bbh.pl # will generate out2.updated
#./generate_csv.pl out2.updated

# seems a bunch no hits in other pnuemocystis - but really a lot of misclassification or the best hit is 
# already taken by another protein

# simple blast - how many are really not found in other pneumocystis
#cat \
#/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_001477545.1_Pneu_cari_B80_V3_protein.faa \
#/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_001477535.1_Pneu_jiro_RU7_V2_protein.faa \
#/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_000349005.1_Pneumo_murina_B123_V2_protein.faa \
#> ref.pneumo.faa
#~/utils/oldBlastall/blast-2.2.11/bin/blastall -p blastp -i scaffolds.f5.annot.aa \
#-d ref.pneumo.faa -v 1 -b 1 -e 1e-20 -m 8 -o scaffolds.f5.annot.bl
#cat scaffolds.f5.annot.bl | cut -f 1 | sort -u | wc -l # 5491 
