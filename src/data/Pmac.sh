# modules
module load minimap2/2.10
module load SAMtools/1.7-goolf-1.7.20
module load BEDTools/2.23.0-goolf-1.7.20
module load pilon/1.22
module load bwa/0.7.17-goolf-1.7.20

#
cat \
    /Ousmane_Data/DATA/Pneumocystis_data/Reads/Pmacacae/NAID_2018/NANOPORE_DATA/run3/LM_PCP_D2A*/fastq_*/*.fastq > FAK53620_allReads.fastq
cat \
    /Ousmane_Data/DATA/Pneumocystis_data/Reads/Pmacacae/NAID_2018/NANOPORE_DATA/run3/LM_PCP_D2_35K/fastq_*/*.fastq > FAK53629_allReads.fastq

# Mapping to macaque genome
minimap2 -ax map-ont \
    /Ousmane_Data/DATA/genomesDB/monkey/GCF_000772875.2_Mmul_8.0.1_genomic.fna \
    FAK53620_allReads.fastq.gz > FAK53620_allReads2mac.sam
minimap2 -ax map-ont \
    /Ousmane_Data/DATA/genomesDB/monkey/GCF_000772875.2_Mmul_8.0.1_genomic.fna \
    FAK53629_allReads.fastq.gz > FAK53629_allReads2mac.sam

# macaque
samtools view -h -S -F 4 FAK53620_allReads2mac.sam > FAK53620_mapped2macaque.sam
samtools view -h -S -F 4 FAK53629_allReads2mac.sam > FAK53629_mapped2macaque.sam

grep -v '^@' FAK53620_mapped2macaque.sam | cut -f 1 | sort -u > tmp.1 # 471533 (82.8% of 569427)
grep -v '^@' FAK53629_mapped2macaque.sam | cut -f 1 | sort -u > tmp.2 # 361436 (77.6% of 465214)

# retrieve reads that don't map to macaque genome
samtools view -h -S -f 4 FAK53620_allReads2mac.sam > FAK53620_notMap2Macaque.sam
grep -v '^@' FAK53620_notMap2Macaque.sam | cut -f 1 | sort -u > tmp.3 # 97894
~/utils/seqtk-master/seqtk subseq FAK53620_allReads.fastq.gz tmp.3 > FAK53620_notMap2Macaque.fq

samtools view -h -S -f 4 FAK53629_allReads2mac.sam > FAK53629_notMap2Macaque.sam
grep -v '^@' FAK53629_notMap2Macaque.sam | cut -f 1 | sort -u > tmp.4 # 103778
~/utils/seqtk-master/seqtk subseq FAK53629_allReads.fastq.gz tmp.4 > FAK53629_notMap2Macaque.fq

# map 2 Pmk
/Ousmane_Data/DATA/softs/ngmlr-0.2.7/ngmlr -t 4 \
    -r /Ousmane_Data/projects/Pmacaca_genome/NAID_2018/Assembly_v2_nanopore/Pmac_genome_merge_v2b.fasta \
    -q FAK53620_notMap2Macaque.fq.gz -o FAK53620_notMap2Macaque_2PmacV2b_ngmlr.sam -x ont 

/Ousmane_Data/DATA/softs/ngmlr-0.2.7/ngmlr -t 4 \
    -r /Ousmane_Data/DATA/Pneumocystis_data/Mitochondria/genomes/PMA_H835_mitogenome_0306.fasta \
    -q FAK53620_notMap2Macaque.fq.gz -o FAK53620_notMap2Macaque_2PmacMt.sam -x ont # Done (407 reads mapped (0.42%), 97487 reads not mapped)

/Ousmane_Data/DATA/softs/ngmlr-0.2.7/ngmlr -t 4 \
    -r /Ousmane_Data/projects/Pmacaca_genome/NAID_2018/Assembly_v2_nanopore/Pmac_genome_merge_v2b.fasta \
    -q FAK53629_notMap2Macaque.fq.gz -o FAK53629_notMap2Macaque_2PmacV2b_ngmlr.sam -x ont

/Ousmane_Data/DATA/softs/ngmlr-0.2.7/ngmlr -t 4 \
    -r /Ousmane_Data/DATA/Pneumocystis_data/Mitochondria/genomes/PMA_H835_mitogenome_0306.fasta \
    -q FAK53629_notMap2Macaque.fq.gz -o FAK53629_notMap2Macaque_2PmacMt.sam -x ont # Done (291 reads mapped (0.28%), 103487 reads not mapped)

samtools view -h -S -F 4 FAK53620_notMap2Macaque_2PmacV2b_ngmlr.sam > FAK53620_Pmac.sam
samtools view -h -S -F 4 FAK53629_notMap2Macaque_2PmacV2b_ngmlr.sam > FAK53629_Pmac.sam
cat FAK53620_Pmac.sam | grep -v '^@' | cut -f 1 | sort -u > tmp.5
cat FAK53629_Pmac.sam | grep -v '^@' | cut -f 1 | sort -u > tmp.6
~/utils/seqtk-master/seqtk subseq FAK53620_allReads.fastq.gz tmp.5 > FAK53620_Pmac.fq
~/utils/seqtk-master/seqtk subseq FAK53629_allReads.fastq.gz tmp.6 > FAK53629_Pmac.fq

# get coverage
samtools view -Sb FAK53620_Pmac.sam > FAK53620_Pmac.bam && samtools sort FAK53620_Pmac.bam -O bam -o FAK53620_Pmac_sorted.bam 
samtools view -Sb FAK53629_Pmac.sam > FAK53629_Pmac.bam && samtools sort FAK53629_Pmac.bam -O bam -o FAK53629_Pmac_sorted.bam
bedtools genomecov -ibam FAK53620_Pmac_sorted.bam -g Pmac_genome_merge_v2b.genome -bg > FAK53620_Pmac_sorted.cov
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' FAK53620_Pmac_sorted.cov # 10.9386
bedtools genomecov -ibam FAK53629_Pmac_sorted.bam -g Pmac_genome_merge_v2b.genome -bg > FAK53629_Pmac_sorted.cov
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' FAK53629_Pmac_sorted.cov # 17.2797

# merge the two bam to get the overal coverage
samtools merge FAK53620_AND_FAK53629_merged.bam FAK53620_Pmac_sorted.bam FAK53629_Pmac_sorted.bam 
bedtools genomecov -ibam FAK53620_AND_FAK53629_merged.bam -g Pmac_genome_merge_v2b.genome -bg > FAK53620_AND_FAK53629_merged.cov
awk '{ sum += $4; n++ } END { if (n > 0) print sum / n; }' FAK53620_AND_FAK53629_merged.cov # 26.2486
# find reads that map neither to macqque or Pneumocysts

### Assembly  (nuclear)
# Canu
module load canu/1.8.0
module load java/1.8.0_45
canu -p FAK53620_Canu_assembly genomeSize=8.0m useGrid=false -nanopore-raw FAK53620_notMap2Macaque.fq.gz # failed -> too many short reads.  Check your reads!
canu -p FAK53629_Canu_assembly genomeSize=8.0m useGrid=false -nanopore-raw FAK53629_notMap2Macaque.fq.gz # failed -> too many short reads.  Check your reads!

zcat FAK53620_notMap2Macaque.fq.gz FAK53629_notMap2Macaque.fq.gz > FAK53620_AND_FAK53629_notMap2Macaque.fq
canu -p FAK53620_AND_FAK53629_Canu_assembly genomeSize=8.0m useGrid=false -nanopore-raw FAK53620_AND_FAK53629_notMap2Macaque.fq.gz # failed -> too many short reads.  Check your reads!

zcat FAK53620_notMap2Macaque.fq.gz FAK53629_notMap2Macaque.fq.gz \
    ../Nanopore_run2/FAK51336_allReads2mac_unmapped.fastq.gz | gzip - >  FAK53620_AND_FAK53629_FAK51336_notMap2Macaque.fq.gz
canu -p Run2_3_combined_Canu_assembly genomeSize=8.0m useGrid=false -nanopore-raw FAK53620_AND_FAK53629_FAK51336_notMap2Macaque.fq.gz # failed too many short reads.  Check your reads!

zcat FAK53620_Pmac.fq.gz FAK53629_Pmac.fq.gz | gzip - > FAK53620_AND_FAK53629_Pmac.fq.gz
zcat FAK53620_Pmac.fq.gz FAK53629_Pmac.fq.gz ../Nanopore_run2/FAK51336_PmacReads.fastq.gz | gzip - > FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz 

canu -p Canu_run2_3_combined_clean genomeSize=8.0m useGrid=false -nanopore-raw FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz stopOnLowCoverage=10 # too many short reads.  Check your reads!

# see the length distribution of these reads (seems like they contains some long reads, will dig later)
# gunzip -c file.fqz | awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)>10000{print a"\n"b"\n"c"\n"$0;}' - | gzip -c - > result.fqz
 gunzip -c FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz \
     | awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)>5000{print a"\n"b"\n"c"\n"$0;}' - | gzip -c - > FAK53620_AND_FAK53629_FAK51336_Pmac_sup5kb.fq.gz

canu -p Canu_run2_3_combined_clean2  genomeSize=8.0m useGrid=false -nanopore-raw FAK53620_AND_FAK53629_FAK51336_Pmac_sup5kb.fq.gz # failed
canu -p Canu_run2_3_combined_clean2  genomeSize=8.0m useGrid=false -nanopore-raw FAK53620_AND_FAK53629_FAK51336_Pmac_sup5kb.fq.gz stopOnLowCoverage=10 # failed
canu -assemble -p Canu_run2_3_combined_clean2  genomeSize=8.0m useGrid=false -nanopore-raw FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz # failed
canu -p Canu_run2_3_combined_clean2  genomeSize=8.0m useGrid=false -nanopore-raw FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz  minReadLength=5000 # failed
canu -nanopore_raw -p Canu_run2_3_combined_clean2 -d test_canu FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz genomeSize=8000000 gnuplotTested=true # failed
canu -p Canu_run2_3_combined_clean2 useGrid=false -nanopore-raw FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz genomeSize=8000000 minReadLength=1100 stopOnLowCoverage=10 correctedErrorRate=0.10 # run


# MiniAsm (error - ava-pb is for pacbio - not used at the end anyway)
minimap2  -x ava-pb -t8 FAK53620_notMap2Macaque.fq.gz FAK53620_notMap2Macaque.fq.gz  | gzip -1 > FAK53620_notMap2Macaque.paf.gz
~/utils/miniasm/miniasm -f FAK53620_notMap2Macaque.fq.gz FAK53620_notMap2Macaque.paf.gz > FAK53620_notMap2Macaque.gfa
awk '/^S/{print ">"$2"\n"$3}' FAK53620_notMap2Macaque.gfa | fold > FAK53620_notMap2Macaque.miniasm.fa # Number of Contigs=56, Total bp=2456822

minimap2  -x ava-pb -t8 FAK53620_Pmac.fq.gz FAK53620_Pmac.fq.gz  | gzip -1 > FAK53620_Pmac.paf.gz
~/utils/miniasm/miniasm -f FAK53620_Pmac.fq.gz FAK53620_Pmac.paf.gz > FAK53620_Pmac.gfa
awk '/^S/{print ">"$2"\n"$3}' FAK53620_Pmac.gfa | fold > FAK53620_Pmac.miniasm.fa # Number of Contigs=52, Total bp=2343966

minimap2  -x ava-pb -t8 FAK53629_notMap2Macaque.fq.gz FAK53629_notMap2Macaque.fq.gz | gzip -1 > FAK53629_notMap2Macaque.paf.gz
~/utils/miniasm/miniasm -f FAK53629_notMap2Macaque.fq.gz FAK53629_notMap2Macaque.paf.gz > FAK53629_notMap2Macaque.gfa
awk '/^S/{print ">"$2"\n"$3}' FAK53629_notMap2Macaque.gfa | fold > FAK53629_notMap2Macaque.miniasm.fa # Number of Contigs=82, Total bp=4687177

minimap2  -x ava-pb -t8 FAK53629_Pmac.fq.gz FAK53629_Pmac.fq.gz  | gzip -1 > FAK53629_Pmac.paf.gz
~/utils/miniasm/miniasm -f FAK53629_Pmac.fq.gz FAK53629_Pmac.paf.gz > FAK53629_Pmac.gfa
awk '/^S/{print ">"$2"\n"$3}' FAK53629_Pmac.gfa | fold > FAK53629_Pmac.miniasm.fa # Number of Contigs=77, Total bp=4472276

minimap2  -x ava-pb -t8 FAK53620_AND_FAK53629_notMap2Macaque.fq.gz FAK53620_AND_FAK53629_notMap2Macaque.fq.gz  | gzip -1 >  FAK53620_AND_FAK53629_notMap2Macaque.paf.gz
~/utils/miniasm/miniasm -f FAK53620_AND_FAK53629_notMap2Macaque.fq.gz FAK53620_AND_FAK53629_notMap2Macaque.paf.gz > FAK53620_AND_FAK53629_notMap2Macaque.gfa
awk '/^S/{print ">"$2"\n"$3}' FAK53620_AND_FAK53629_notMap2Macaque.gfa | fold >  FAK53620_AND_FAK53629_notMap2Macaque.miniasm.fa

minimap2  -x ava-pb -t8 FAK53620_AND_FAK53629_Pmac.fq.gz FAK53620_AND_FAK53629_Pmac.fq.gz | gzip -1 >  FAK53620_AND_FAK53629_Pmac.paf.gz
~/utils/miniasm/miniasm -f FAK53620_AND_FAK53629_Pmac.fq.gz FAK53620_AND_FAK53629_Pmac.paf.gz > FAK53620_AND_FAK53629_Pmac.gfa
awk '/^S/{print ">"$2"\n"$3}' FAK53620_AND_FAK53629_Pmac.gfa | fold >  FAK53620_AND_FAK53629_Pmac.miniasm.fa # Number of Contigs=59, Total bp=7340125

minimap2  -x ava-pb -t8 FAK53620_AND_FAK53629_FAK51336_notMap2Macaque.fq.gz FAK53620_AND_FAK53629_FAK51336_notMap2Macaque.fq.gz | gzip -1 > FAK53620_AND_FAK53629_FAK51336_notMap2Macaque.paf.gz
~/utils/miniasm/miniasm -f FAK53620_AND_FAK53629_FAK51336_notMap2Macaque.fq.gz FAK53620_AND_FAK53629_FAK51336_notMap2Macaque.paf.gz >  FAK53620_AND_FAK53629_FAK51336_notMap2Macaque.gfa
awk '/^S/{print ">"$2"\n"$3}' FAK53620_AND_FAK53629_FAK51336_notMap2Macaque.gfa  | fold >  FAK53620_AND_FAK53629_FAK51336_notMap2Macaque.miniasm.fa

minimap2  -x ava-pb -t8 FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz | gzip -1 > FAK53620_AND_FAK53629_FAK51336_Pmac.paf.gz
~/utils/miniasm/miniasm -f FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz FAK53620_AND_FAK53629_FAK51336_Pmac.paf.gz > FAK53620_AND_FAK53629_FAK51336_Pmac.gfa
awk '/^S/{print ">"$2"\n"$3}' FAK53620_AND_FAK53629_FAK51336_Pmac.gfa | fold > FAK53620_AND_FAK53629_FAK51336_Pmac.miniasm.fa # Number of Contigs=25, Total bp=7893027

## plot
module load gnuplot/5.0.3-foss-2016a
module load MUMmer/4.0.0
nucmer --maxmatch -c 100 -p V2b_vs_Nanopore_run3 Pmac_genome_merge_v2b.fasta FAK53620_AND_FAK53629_FAK51336_Pmac.miniasm.fa 
show-coords -c V2b_vs_Nanopore_run3.delta > V2b_vs_Nanopore_run3.coords
mummerplot --fat --filter --png --large -p V2b_vs_Nanopore_run3 V2b_vs_Nanopore_run3.delta -R Pmac_genome_merge_v2b.fasta -Q FAK53620_AND_FAK53629_FAK51336_Pmac.miniasm.fa 
gnuplot V2b_vs_Nanopore_run3.gp

nucmer --maxmatch -c 100 -p V2b_vs_CanuMapped Pmac_genome_merge_v2b.fasta Pmac_canu_assembly/Pmac.contigs.fasta
show-coords -c V2b_vs_CanuMapped.delta > V2b_vs_CanuMapped.coords
mummerplot --fat --filter --png --large -p V2b_vs_CanuMapped V2b_vs_CanuMapped.delta -R Pmac_genome_merge_v2b.fasta -Q Pmac_canu_assembly/Pmac.contigs.fasta
gnuplot V2b_vs_CanuMapped.gp

nucmer --maxmatch -c 100 -p V2b_vs_CanuUnmapped Pmac_genome_merge_v2b.fasta Unmap_canu_assembly/Unmap.contigs.fasta
show-coords -c V2b_vs_CanuUnmapped.delta > V2b_vs_CanuUnmapped.coords
mummerplot --fat --filter --png --large -p V2b_vs_CanuUnmapped V2b_vs_CanuUnmapped.delta -R Pmac_genome_merge_v2b.fasta -Q Unmap_canu_assembly/Unmap.contigs.fasta
gnuplot V2b_vs_CanuUnmapped.gp

nucmer --maxmatch -c 100 -p V2_vs_miniasm Pmac_genome_merge_v2b.fasta FAK51336_PmacReads.miniasm.fa
show-coords -c V2_vs_miniasm.delta > V2_vs_miniasm.coords
mummerplot --fat --filter --png --large -p V2_vs_miniasm V2_vs_miniasm.delta -R Pmac_genome_merge_v2b.fasta -Q FAK51336_PmacReads.miniasm.fa
gnuplot V2_vs_miniasm.gp

nucmer --maxmatch -c 100 -p V2_vs_Canu Pmac_genome_merge_v2b.fasta Canu_run2_3_combined_clean2.contigs.fasta
mummerplot --fat --filter --png --large -p V2_vs_Canu V2_vs_Canu.delta -R Pmac_genome_merge_v2b.fasta -Q Canu_run2_3_combined_clean2.contigs.fasta
gnuplot V2_vs_Canu.gp # looks better than Miniasm 

nucmer --maxmatch -c 100 -p Canu_vs_Miniasm Canu_run2_3_combined_clean2.contigs.fasta FAK53620_AND_FAK53629_FAK51336_Pmac.miniasm.fa
mummerplot --fat --filter --png --large -p Canu_vs_Miniasm Canu_vs_Miniasm.delta  -R Canu_run2_3_combined_clean2.contigs.fasta -Q FAK53620_AND_FAK53629_FAK51336_Pmac.miniasm.fa
gnuplot Canu_vs_Miniasm.gp # Canu is better

# best assembly : Canu_run2_3_combined_clean2.contigs.fasta
# but Canu_run2_3_combined_clean2.unassembled.fasta is pretty large and there are clearly Pneumocystis (not used for now_
# How many of these unassembled unitigs are missing the Canu assembly?
# minimap2 -x asm5 -t 4 Canu_run2_3_combined_clean2.contigs.fasta  Canu_run2_3_combined_clean2.unassembled.fasta > tmp.1.paf
# cat tmp.1.paf | cut -f 6 | sort -u | wc -l # only 22
# minimap2 -x asm20 -t 4 Canu_run2_3_combined_clean2.contigs.fasta  Canu_run2_3_combined_clean2.unassembled.fasta > tmp.2.paf
#cat tmp.2.paf | cut -f 6 | sort -u | wc -l # 22
# OK, so most of them are not in the assembly which makes sense https://github.com/marbl/canu/issues/103
# I don't think these reads are useful, but I cannot exclude that these reads contain something useful

## Polishing
# difficult to install https://github.com/nanoporetech/ont-assembly-polish/blob/master/analysis.mk
# so I reproduce manually

minimap2  -x ava-ont -t 4 Canu_run2_3_combined_clean2.contigs.fasta FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz | gzip -1 > minimap_overlap.paf.gz
/Ousmane_Data/DATA/softs/racon/build/bin/racon \
    -t 4 FAK53620_AND_FAK53629_FAK51336_Pmac.fq.gz minimap_overlap.paf.gz Canu_run2_3_combined_clean2.contigs.fasta > racon_contigs.fasta
bwa index racon_contigs.fasta
bwa mem -M -t 6 racon_contigs.fasta P2C_allSamples_merged_1.fq.gz P2C_allSamples_merged_2.fq.gz | samtools view -S -b -u - | samtools sort - -o racon_contigs.bam
samtools index racon_contigs.bam
java -Xmx16G -jar $EBROOTPILON/pilon-1.22.jar --threads 4 --genome racon_contigs.fasta --bam racon_contigs.bam --outdir pilon --output pilon.contigs.fasta

# annotation w Funannote
### Problem is that 236 genes are missing compared to v2b version
# Well v2b is not annotated but v2a (which is little bit contamined is)
cd /Ousmane_Data/projects/Pmacaca_genome/NAID_2018/Nanopore_run3/FUNANN/predict_results

module load blast/2.2.22-Linux_x86_64
~/utils/proteinortho_v2.3.0.pl Pjir.fa v2a_proteins.fa Pneumocystis_NAID_NANO_v1.proteins.fa -a=4 > rbbh.out
perl find_missing_genes.pl rbbh.out # find the missing genes, locate the scaffolds where they are located in the old assembly and write the file :'Scaffolds_to_look_in_the_newVersion.txt' and 'oldGenes_missing_in_New.txt'

perl ~/utils/retrieveFasta.pl Scaffolds_to_look_in_the_newVersion.txt v2a_genome.fa >  Scaffolds_to_look_in_the_newVersion.fa # too big
perl ~/utils/retrieveFasta.pl OldGenes_missing_in_New.txt v2a_cds.fasta > OldGenes_missing_in_New.fa

dnadiff Pneumocystis_NAID_NANO_v1.scaffolds.fa OldGenes_missing_in_New.fa -p OldGenes_missing_in_New # looks like everything is mapped
dnadiff Pneumocystis_NAID_NANO_v1.scaffolds.fa Pmac_genome_merge_v2b.fasta -p Nano_v1_vs_V2b # looks like everything is mapped (looks like collapsed)
# conclusion, not why there are differences, but it's not due to missing sequences, looks like some sequences were expanded or duplicated in the Illumina v2b
# final version is Pneumocystis_NAID_NANO_v1.scaffolds.fa
# clean up and move tmp files here Investigation

__END__

