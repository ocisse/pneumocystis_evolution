#Mon Apr 17 17:25:10 EDT 2017

# modules
module load Bowtie2/2.2.5-goolf-1.7.20
module load SAMtools/1.2-goolf-1.7.20-HTSlib-1.2.1
module load BEDTools/2.23.0-goolf-1.7.20

# Reads
export RLIB1=/Ousmane_Data/DATA/Pneumocystis_data/Reads/P_dogs/raw_data_LS03061701
export RLIB2=/Ousmane_Data/DATA/Pneumocystis_data/Reads/P_dogs/170523
export RLIB3=/Ousmane_Data/DATA/Pneumocystis_data/Reads/P_dogs/Raw_data_DOF-F2-PE150
export RLIB4=/Ousmane_Data/DATA/Pneumocystis_data/Reads/P_dogs/Raw_data_DOG-F2-PE250

export WRKDIR=/Ousmane_Data/projects/P_dogs
export DB=/Ousmane_Data/DATA/genomesDB
export PROTDB=/Ousmane_Data/DATA/proteomesDB

# building db
bowtie2-build $DB/dog/GCF_000002285.3_CanFam3.1_genomic.fna $DB/dog/dog

# clean 1
cd $WRKDIR

bowtie2 -x $DB/dog/dog \
-1 $RLIB1/DOG_F_1.fq.gz,$RLIB2/DOG_F_1.fq.gz,$RLIB3/DOG_F2_1.fq.gz,$RLIB4/DOG_F2_1.fq.gz \
-2 $RLIB1/DOG_F_2.fq.gz,$RLIB2/DOG_F_2.fq.gz,$RLIB3/DOG_F2_2.fq.gz,$RLIB4/DOG_F2_2.fq.gz \
-S $WRKDIR/mapped2dog.sam --threads 12 --sensitive-local --un-conc AllLib_unmapped

# clean 2
bowtie2 -x $DB/Arabidopsis/Atha -1 AllLib_unmapped.1.fq.gz -2 AllLib_unmapped.2.fq.gz \
-S $WRKDIR/mapped2atha.sam --threads 12 --sensitive-local --un-conc $WRKDIR/AllLib_unmappedV2

# 1.5% reads mapp to Arabidopsis genome
# clean 3  removeal puccina
bowtie2 -x $DB/Puccina/pgra -1 AllLib_unmappedV2.1.fq.gz -2 AllLib_unmappedV2.2.fq.gz \
-S $WRKDIR/mapped2puccina.sam --threads 12 --very-fast --un-conc $WRKDIR/AllLib_unmappedV3

# Assembly unmapped
spades.py --careful -t 24 -1 $WRKDIR/AllLib_unmappedV3.1.fq.gz -2 $WRKDIR/AllLib_unmappedV3.2.fq.gz \
-o $WRKDIR/Assembly_cleanV3

# clean up
module load BLAT/3.5-goolf-1.7.20
blat \
-dots=10000 \
$DB/dog/GCF_000002285.3_CanFam3.1_genomic.fna \
$WRKDIR/Assembly_cleanV3/scaffolds.fasta \
$WRKDIR/Assembly_cleanV3/scaffolds.fasta2CanFam3.psl

blat \
-dots=10000 \
$DB/Arabidopsis/GCF_000001735.3_TAIR10_genomic.fna \
$WRKDIR/Assembly_cleanV3/scaffolds.fasta \
$WRKDIR/Assembly_cleanV3/scaffolds.fasta2Atha.psl

# removing these scafs
cat $WRKDIR/Assembly_cleanV3/scaffolds.fasta2CanFam3.psl \
$WRKDIR/Assembly_cleanV3/scaffolds.fasta2CanFam3.psl \
> $WRKDIR/Assembly_cleanV3/tmp.1

perl remove_longMatches.pl tmp.1 > bad_scafs.txt 

cat bad_scafs.txt | cut -f 1 | sort -u > tmp.2
grep '^>' scaffolds.fasta | cut -f 2 -d '>' | sort -u > tmp.3
diff -y tmp.3 tmp.2 | grep '<' | cut -f 1 | sort -u > tmp.4
perl ~/utils/retrieveFasta.pl tmp.4 scaffolds.fasta > scaffolds.cl1.fasta
rm tmp.1 tmp.2 tmp.3 tmp.4

# remove bacteria
~/utils/oldBlastall/blast-2.2.11/bin/blastall -p blastx \
-d $PROTDB/bacteria/uniprot_sprot_bacteria.fa \
-i $WRKDIR/Assembly_cleanV3/scaffolds.cl1.fasta \
-e 1e-15 -v 1 -b 1 -a 6 -m 8 \
-o $WRKDIR/Assembly_cleanV3/scaffolds.cl1.2bact.blast # not a good idea, many contigs w hits are infact Pneumocystis

# manual check - looks like there some dog microsat seque
blat -dots=1000 Dog_micro_satellites.fas scaffolds.cl1.fasta scaffolds.cl1.2Dogmt.psl
perl remove_longMatches.pl scaffolds.cl1.2Dogmt.psl > dog_microsat.txt
grep '^NODE' dog_microsat.txt | cut -f 1 | sort -u > tmp.1
grep '^>' scaffolds.cl1.fasta | cut -f 2 -d '>' | sort -u > tmp.2
diff -y tmp.2 tmp.1 | grep '<' | cut -f 1 | sort -u > tmp.3
perl ~/utils/retrieveFasta.pl tmp.3 scaffolds.cl1.fasta > scaffolds.cl2.fasta
rm tmp.1 tmp.2 tmp.3

# remove contigs w GC% > 55% and/or length  500
# remove High GC contigs - contam 
~/utils/get_fasta_stats.pl -a scaffolds.cl2.fasta > tmp.1
perl remove_highGC.pl tmp.1 | sort -u > tmp.2
grep '^>' scaffolds.cl2.fasta | cut -f 2 -d '>' | sort -u > tmp.3
diff -y tmp.3 tmp.2 |  grep '<' | cut -f 1 | sort -u > tmp.4
perl ~/utils/retrieveFasta.pl tmp.4 scaffolds.cl2.fasta > scaffolds.cl3.fasta
rm tmp.1 tmp.2 tmp.3 tmp.4

# remove scafs < 500
~/utils/get_fasta_stats.pl -a scaffolds.cl3.fasta > tmp.1
perl select_longScafs.pl tmp.1 > tmp.2
perl ~/utils/retrieveFasta.pl tmp.2 scaffolds.cl3.fasta > scaffolds.cl4.fasta
rm tmp.1 tmp.2 

# how many belongs to MSGs
~/utils/oldBlastall/blast-2.2.11/bin/blastall -p blastx \
-d /Ousmane_Data/DATA/Pneumocystis_data/MSGs/merge.prot \
-i /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl4.fasta \
-e 1e-5 -v 1 -b 1 -a 6 -m 8 \
-o /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl4.putative_msgs.blast
# not much only 298 hits => the vast vajority of contigs < 5 kb are not MSGs
# MSGs are not the problem

# run redundans
module load Python/2.7.12-foss-2016b
python ~/utils/redundans/redundans.py -v -t 15 \
-i AllLib_unmappedV3.1.fq.gz AllLib_unmappedV3.2.fq.gz \
-f /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl4.fasta \
-o Redundans_ass/run1
# based on Redundans report most of the small contigs are variants of regions integrated in large contigs (i.e. different pops)
# redundans will for haplotigs analysis

ln -s Redundans_ass/run1/contigs.reduced.fa scaffolds.cl5.fasta
# not much difference w redundans scaffolding

# first remove all contigs that map 2 pneu
~/utils/oldBlastall/blast-2.2.11/bin/megablast \
-d /Ousmane_Data/DATA/Pneumocystis_data/NCBI/Pneumocystis_Pj_Pc_Pm_concat_WGS.fasta \
-i scaffolds.cl5.fasta \
-m 8 -e 1e-5 -a 10 \
-o scaffolds.cl5.2otherPneumo.mgblast

cat scaffolds.cl5.2otherPneumo.mgblast | cut -f 1 | sort -u > scaffolds.cl5.tokeep.txt # only 251 contigs
perl ~/utils/retrieveFasta.pl scaffolds.cl5.tokeep.txt scaffolds.cl5.fasta > scaffolds.cl5.tokeep.fasta

# extract the contigs without hits in Pneumo
grep '^>' scaffolds.cl5.fasta | cut -f 2 -d '>' | sort -u > tmp.1
diff -y tmp.1 scaffolds.cl5.tokeep.txt |  grep '<' | cut -f 1 | sort -u > tmp.2
perl ~/utils/retrieveFasta.pl tmp.2 scaffolds.cl5.fasta >  scaffolds.cl5.unkn.fasta
rm tmp.1 tmp.2 
# this contains mostly dog sequences, but need to make sure there no fungi in here

~/utils/oldBlastall/blast-2.2.11/bin/blastall -p blastx \
-d /Ousmane_Data/DATA/Maker2/fun_db.faa \
-i /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl5.unkn.fasta \
-e 1e-20 -v 1 -b 1 -m 8 -a 10 \
-o /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl5.unkn.2fun_db.blastx

# some pneumocystis are here, but divergent, I rescue 373 contigs
cat scaffolds.cl5.unkn.2fun_db.blastx | cut -f 1  | sort -u > tmp.1
perl ~/utils/retrieveFasta.pl tmp.1 scaffolds.cl5.unkn.fasta > scaffolds.cl5.unkn.rescued.fasta

# merge w 
cat scaffolds.cl5.tokeep.fasta scaffolds.cl5.unkn.rescued.fasta >  scaffolds.cl6.fasta # final version 8 Mb; 624 scafs


# Extract large contigs 
perl ~/utils/retrieveFasta.pl pj.ref.Large_scafs_only /Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCA_001477535.1_Pneu_jiro_RU7_V2_genomic.fasta > pj.ref.Large_scafs_only.fas
perl ~/utils/retrieveFasta.pl pc.ref.Large_scafs_only /Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCA_001477545.1_Pneu_cari_B80_V3_genomic.fasta > pc.ref.Large_scafs_only.fas


# GAP closure - primers design
# Extract the first 300 np and last 300 bp from each scaffolds
cd /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/PRIMERS
./search_primers.pl
# seems primer cannot run outside so I copy in test fold
cp /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/PRIMERS/fasta/*.txt \
~/utils/primer3-2.3.7/test
cd ~/utils/primer3-2.3.7/test
ls NODE*.txt | while read f; do \
~/utils/primer3-2.3.7/src/primer3_core < $f;done
mv NODE*.for /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/PRIMERS/fasta/
mv NODE*.rev /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/PRIMERS/fasta/

cd /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/PRIMERS/fasta/
# proocess primer3 file - pick one per file/scaf and rename to make it easy to rememe
ls *.for | while read f; do perl pick_onePerFile.pl $f >> primers_raw.txt;done
ls *.rev | while read f; do perl pick_onePerFile.pl $f >> primers_raw.txt;done
cat primers_raw.txt | cut -f 4 -d '_' > len.txt
paste primers_raw.txt len.txt > primers_w_len.txt


# quick annotation to see if something is missing
~/utils/augustus-3.2.1/bin/augustus --gff3=on --species=pneumocystis scaffolds.cl6.fasta > scaffolds.cl6.annot.gff
perl ~/utils/augustus-3.2.1/scripts/getAnnoFasta.pl scaffolds.cl6.annot.gff --seqfile=scaffolds.cl6.fasta
cd ~/utils/ORTHOMCLV1.4
ln -s /Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_001477545.1_Pneu_cari_B80_V3_protein.faa pc.faa
ln -s /Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_001477535.1_Pneu_jiro_RU7_V2_protein.faa pj.faa
ln -s /Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_000349005.1_Pneumo_murina_B123_V2_protein.faa pm.faa
ln -s /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl6.annot.aa

perl orthomcl.pl --mode 1 --fa_files pc.faa,pj.faa,pm.faa,scaffolds.cl6.annot.aa

### Haplotype analysis
cd /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/Haplotype_analysis

# how many of the filtered assembly have distinct haplotypes
grep '^>' scaffolds.cl6.fasta | cut -f 2 -d '>' \
| while read f; do grep $f \
/Ousmane_Data/projects/P_dogs/Assembly_cleanV3/Redundans_ass/run1/contigs.reduced.fa.hetero.tsv \
>> /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/Haplotype_analysis/scaffolds.cl6.fasta.hetero.tsv;done

./group_haplotigs.pl
grep 'WITH_HAPLOTYPE' Haplotype_report | cut -f 1 > x
perl ~/utils/retrieveFasta.pl x scaffolds.cl6.fasta > WITH_HAPLOTYPE.fa
rm x

grep 'WITHOUT_HAPLOTYPE' Haplotype_report | cut -f 1 > x
perl ~/utils/retrieveFasta.pl x scaffolds.cl6.fasta > WITHOUT_HAPLOTYPE.fa
rm x

~/utils/augustus-3.2.1/bin/augustus --gff3=on --species=pneumocystis WITH_HAPLOTYPE.fa  > WITH_HAPLOTYPE.annot.gff
perl ~/utils/augustus-3.2.1/scripts/getAnnoFasta.pl WITH_HAPLOTYPE.annot.gff --seqfile=WITH_HAPLOTYPE.fa

~/utils/augustus-3.2.1/bin/augustus --gff3=on --species=pneumocystis WITHOUT_HAPLOTYPE.fa > WITHOUT_HAPLOTYPE.annot.gff
perl ~/utils/augustus-3.2.1/scripts/getAnnoFasta.pl WITHOUT_HAPLOTYPE.annot.gff --seqfile=WITHOUT_HAPLOTYPE.fa

# clearly, there are multiple pops: 90% genome is heterogenous

# recovering the alternative genome assembly
cd /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/Haplotype_analysis
./export_alternative_contigs.pl 
perl ~/utils/retrieveFasta.pl Haplotype_report.contigs /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.fasta > ALTERNATIVE_UNORDER.fa
ln -s /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/Haplotype_analysis/ALTERNATIVE_REORDER/superscaffolds.fasta scaffolds.cl6.alter.fasta


# re-order using  Satsuma2

/nethome/cisseoh/utils/satsuma2/bin/Chromosemble \
-q /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/Haplotype_analysis/WITH_HAPLOTYPE.fa \
-t /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl6.fasta \
-o /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/Haplotype_analysis/ChromosembleDIR_WITH_HAPLOTYPE \
-n 12

/nethome/cisseoh/utils/satsuma2/bin/Chromosemble \
-q /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/Haplotype_analysis/ALTERNATIVE_UNORDER.fa \
-t /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl6.fasta \
-o /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/Haplotype_analysis/ALTERNATIVE_REORDER \
-n 12

### latest
/nethome/cisseoh/utils/satsuma2/bin/Chromosemble \
-t /Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCA_001477535.1_Pneu_jiro_RU7_V2_genomic.fasta \
-q /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl6.fasta \
-o /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/ChromoAssemble_scaffolds.cl6 \
-n 6 -pseudochr -s

/nethome/cisseoh/utils/satsuma2/bin/Chromosemble \
-t /Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCA_001477535.1_Pneu_jiro_RU7_V2_genomic.fasta \
-q /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl6.alter.fasta \
-o /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/ChromoAssemble_scaffolds_cl6_alt \
-n 6 -pseudochr -s

# annotation
cd /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/MAKER2_ann_RefPOP/Pdg_ref_scaffolds.f6_v2.maker.output # same
cd /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/MAKER2_ann_AltPOP/

# seems like there is some contam - Pcanis_strainCK_003359-RA - this is mice
module load BLAT
blat -t=prot -q=prot \
/Ousmane_Data/DATA/genomesDB/dog/GCF_000002285.3_CanFam3.1_protein.faa \
Pcanis_strainCK_Ref.all.maker.proteins.fasta \
Pcanis_strainCK_Ref.all.maker.proteins.vs.GCF_000002285.3_CanFam3.1_protein.psl

cat \
/Ousmane_Data/DATA/genomesDB/dog/GCF_000002285.3_CanFam3.1_protein.faa \
/Ousmane_Data/DATA/Pneumocystis_data/NCBI/GCF_001477535.1_Pneu_jiro_RU7_V2_protein.SP_ID_ADDED.faa \
> tmp.1.faa

module load hmmer
phmmer -E 1e-20 --cpu 6 \
--tblout Pcanis_strainCK_Ref.all.maker.proteins.vs.GCF_000002285.3_CanFam3.1_protein.tab \
Pcanis_strainCK_Ref.all.maker.proteins.fasta \
tmp.1.faa
# OK problem there are some serious hits against dog  proteome
perl identify_problematic.pl \
Pcanis_strainCK_Ref.all.maker.proteins.vs.GCF_000002285.3_CanFam3.1_protein.tab \
> Pcanis_strainCK_Ref.all.maker.proteins.vs.GCF_000002285.3_CanFam3.1_protein.report

# manually inspected Pcanis_strainCK_Ref.all.maker.proteins.vs.GCF_000002285.3_CanFam3.1_protein.report
# most of them are false positive create by local homology
# possible real contam
Pcanis_strainCK_003352-RA
# seems many are mitochdondrion

Pcanis_strainCK_003224-RA # removed from Pcanis_strainCK_Ref.all.maker.proteins.vs.GCF_000002285.3_CanFam3.1_protein.report
Pcanis_strainCK_003220-RA # removed from Pcanis_strainCK_Ref.all.maker.proteins.vs.GCF_000002285.3_CanFam3.1_protein.report

##### - ASSEMBLY SECOND SAMPLE DOG_A_Austria
# Thu Feb  1 13:36:14 EST 2018
# Reads
export RLIBDA=/Ousmane_Data/DATA/Pneumocystis_data/Reads/P_dogs/DOG_A_Austria/
export WRKDIR2=/Ousmane_Data/projects/P_dogs/dogA_Austria

cd $WRKDIR2
bowtie2 -x \
/Ousmane_Data/DATA/genomesDB/dog/dog \
-1 $RLIBDA/DOG_A_USD16083955L-A52_HF5MWCCXY_L1_1.fq.gz,$RLIBDA/DOG_A_USD16083955L_HF72CCCXY_L4_1.fq.gz,$RLIBDA/DOG_A_USD16083955L_HFCG2CCXY_L4_1.fq.gz \
-2 $RLIBDA/DOG_A_USD16083955L-A52_HF5MWCCXY_L1_2.fq.gz,$RLIBDA/DOG_A_USD16083955L_HF72CCCXY_L4_2.fq.gz,$RLIBDA/DOG_A_USD16083955L_HFCG2CCXY_L4_2.fq.gz \
-S $WRKDIR2/mapped2dog.sam --threads 12 --sensitive-local --un-conc AllLib_unmapped


spades.py --careful -t 24 -1 AllLib_unmapped.1.fq.gz -2 AllLib_unmapped.2.fq.gz -o Ass_V1

# now simply filter using the first assembly
~/utils/oldBlastall/blast-2.2.11/bin/megablast \
-i /Ousmane_Data/projects/P_dogs/dogA_Austria/Ass_V1/scaffolds.fasta \
-d /Ousmane_Data/projects/P_dogs/Assembly_cleanV3/scaffolds.cl6.fasta \
-m 8 -v 1 -b 1 -e 1e-20 \
-o /Ousmane_Data/projects/P_dogs/dogA_Austria/Ass_V1/scaffolds.blast

cat scaffolds.blast | cut -f 1 | sort -u > tmp.txt
perl ~/utils/retrieveFasta.pl tmp.txt scaffolds.fasta > P_dog_Austria_assembly_v1.fasta


