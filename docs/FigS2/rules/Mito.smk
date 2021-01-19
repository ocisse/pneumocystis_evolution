# first to extract genes from all
# convert mfannot tbl into genbank
#../tools/mac.tbl2asn \
#-M n -J -c w -euk \
#-t ../../data/raw/template.sbt.txt \
#-gaps-min 10 -l paired-ends \
#-j "[organism=Pneumocystis jirovecii] [isolate=SE8]" \
#-i ../../data/raw/Pjirovecii/SE8.fasta \
#-f ../../data/raw/Pjirovecii/MFannot_SE8/SE8.fasta.new.tbl \
#-o ../../data/raw/Pjirovecii/MFannot_SE8/SE8.OC.sqn -Z -V b
# not working

rule extract_cds:
	input:
		pj1="../../data/raw/Pjirovecii/Ar_25kb.gb",
		pj2="../../data/raw/Pjirovecii/CH26_full.gb",
		pj3="../../data/raw/Pjirovecii/DK07_25kb.gb",
		pj4="../../data/raw/Pjirovecii/GR11b_full.gb",
		pj5="../../data/raw/Pjirovecii/GR5_full.gb",
		pj6="../../data/raw/Pjirovecii/IT918_25kb.gb",
		pj7="../../data/raw/Pjirovecii/JX855936.1.gb",
		pj8="../../data/raw/Pjirovecii/JX855937.1.gb",
		pj9="../../data/raw/Pjirovecii/JX855938.1.gb",
		pj10="../../data/raw/Pjirovecii/S567_25kb.gb",
		pj11="../../data/raw/Pjirovecii/SW1_full.gb",
		pj12="../../data/raw/Pjirovecii/SW4_full.gb",
		pj13="../../data/raw/Pjirovecii/SW7_full.gb",
		pmac="../../data/raw/Pmacacae/MFannot_H835/PMA_H835_mitogenome_0306.gb",
		pory="../../data/raw/Poryctolagi/MFannot_MT/Poryc_MT_genome_v2.gb",
		pcan1="../../data/raw/Pcanis/MFannot/P_dog_Ck_assembly.Mt_v1.NODE_35.gb",
		pcan2="../../data/raw/Pcanis/MFannot/P_dog_Ck_assembly.Mt_v1.NODE_34.gb",
		pcan3="../../data/raw/Pcanis/MFannot/P_dog_Austria_assembly_v2.Mt.gb",
		pcar1="../../data/raw/Pcarinii/JX499145.1.gb",
		pcar2="../../data/raw/Pcarinii/NC_013660.1.gb",
		pcar3="../../data/raw/Pcarinii/NC_013660.2.gb",
		pmur="../../data/raw/Pmurina/JX499144.1.gb",
		pwk="../../data/raw/Pwakefieldiae/MFannot/Pwk_Mt_genome_v1_renamed.gb"



	run:
		shell("perl ../scripts/cds_extractor.pl -i {input.pj1} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj2} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj3} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj4} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj5} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj6} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj7} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj8} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj9} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj10} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj11} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj12} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pj13} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pmac} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pory} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pcan1} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pcan2} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pcan3} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pcar1} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pcar2} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pcar3} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pmur} -n")
		shell("perl ../scripts/cds_extractor.pl -i {input.pwk} -n")

# 
# bedtools getfasta -fi ../SE8.fasta -bed SE8.cds.bed -fo test -s -name

rule add_species_tag:
	input:
		pj1="../../data/raw/Pjirovecii/Ar_25kb.ffn",
		pj2="../../data/raw/Pjirovecii/CH26_full.ffn",
		pj3="../../data/raw/Pjirovecii/DK07_25kb.ffn",
		pj4="../../data/raw/Pjirovecii/GR11b_full.ffn",
		pj5="../../data/raw/Pjirovecii/GR5_full.ffn",
		pj6="../../data/raw/Pjirovecii/IT918_25kb.ffn",
		pj7="../../data/raw/Pjirovecii/JX855936.1.ffn",
		pj8="../../data/raw/Pjirovecii/JX855937.1.ffn",
		pj9="../../data/raw/Pjirovecii/JX855938.1.ffn",
		pj10="../../data/raw/Pjirovecii/S567_25kb.ffn",
		pj11="../../data/raw/Pjirovecii/SW1_full.ffn",
		pj12="../../data/raw/Pjirovecii/SW4_full.ffn",
		pj13="../../data/raw/Pjirovecii/SW7_full.ffn",
		pj14="../../data/raw/Pjirovecii/MFannot_SE8revcomp/SE8.revcomp.cds.fasta",
		pj15="../../data/raw/Pjirovecii/MFannot_Pj46/Pj46_mito_full.cds.fasta",
		pj16="../../data/raw/Pjirovecii/MFannot_Pj54c/Pj54c_mito_full.cds.fasta",
		pj17="../../data/raw/Pjirovecii/MF_annot_YUrevcomp/YU.mito.revcomp.cds.fasta",
		pmac1="../../data/raw/Pmacacae/MFannot_H835/PMA_H835_mitogenome_0306.ffn",
		pmac2="../../data/raw/Pmacacae/MFannot_P2C/Pmac_P2C_mitogenome.cds.fasta",
		pmac3="../../data/raw/Pmacacae/MFannot_ER17/Pmac_ER17_mitogenome.cds.fasta",
		pmac4="../../data/raw/Pmacacae/MFannot_UC86/Pmac_UC86_mitogenome.cds.fasta",
		pory1="../../data/raw/Poryctolagi/MFannot_MT/Poryc_MT_genome_v2.ffn",
		pory2="../../data/raw/Poryctolagi/MFannot_Prab1/Prab1_mitogenome.cds.fasta",
		pory3="../../data/raw/Poryctolagi/MFannot_PrabF/PrabF_mitogenome.cds.fasta",
		pory4="../../data/raw/Poryctolagi/MFannot_PrabM/PrabM_mitogenome.cds.fasta",
		pcan1="../../data/raw/Pcanis/MFannot/P_dog_Ck_assembly.Mt_v1.NODE_35.ffn",
		pcan2="../../data/raw/Pcanis/MFannot/P_dog_Ck_assembly.Mt_v1.NODE_34.ffn",
		pcan3="../../data/raw/Pcanis/MFannot/P_dog_Austria_assembly_v2.Mt.ffn",
		pcar1="../../data/raw/Pcarinii/JX499145.1.ffn",
		pcar2="../../data/raw/Pcarinii/NC_013660.1.ffn",
		pcar3="../../data/raw/Pcarinii/NC_013660.2.ffn",
		pmur="../../data/raw/Pmurina/JX499144.1.ffn",
		pwk="../../data/raw/Pwakefieldiae/MFannot/Pwk_Mt_genome_v1_renamed.ffn"

	output:
		pj1="../../data/processed/Pjirovecii/Ar_25kb.ffn",
		pj2="../../data/processed/Pjirovecii/CH26_full.ffn",
		pj3="../../data/processed/Pjirovecii/DK07_25kb.ffn",
		pj4="../../data/processed/Pjirovecii/GR11b_full.ffn",
		pj5="../../data/processed/Pjirovecii/GR5_full.ffn",
		pj6="../../data/processed/Pjirovecii/IT918_25kb.ffn",
		pj7="../../data/processed/Pjirovecii/JX855936.1.ffn",
		pj8="../../data/processed/Pjirovecii/JX855937.1.ffn",
		pj9="../../data/processed/Pjirovecii/JX855938.1.ffn",
		pj10="../../data/processed/Pjirovecii/S567_25kb.ffn",
		pj11="../../data/processed/Pjirovecii/SW1_full.ffn",
		pj12="../../data/processed/Pjirovecii/SW4_full.ffn",
		pj13="../../data/processed/Pjirovecii/SW7_full.ffn",
		pj14="../../data/processed/Pjirovecii/SE8.ffn",
		pj15="../../data/processed/Pjirovecii/Pj46.ffn",
		pj16="../../data/processed/Pjirovecii/Pj54c.ffn",
		pj17="../../data/processed/Pjirovecii/YU.ffn",
		pmac1="../../data/processed/Pmacacae/PMA_H835_mitogenome_0306.ffn",
		pmac2="../../data/processed/Pmacacae/Pmac_P2C_mitogenome.ffn",
		pmac3="../../data/processed/Pmacacae/Pmac_ER17_mitogenome.ffn",
		pmac4="../../data/processed/Pmacacae/Pmac_UC86_mitogenome.ffn",
		pory1="../../data/raw/Poryctolagi/MFannot/Poryc_MT_genome_v2.ffn",
		pory2="../../data/raw/Poryctolagi/MFannot/Prab1_mitogenome.ffn",
		pory3="../../data/raw/Poryctolagi/MFannot/PrabF_mitogenome.ffn",
		pory4="../../data/raw/Poryctolagi/MFannot/PrabM_mitogenome.ffn",
		pcan1="../../data/processed/Pcanis/P_dog_Ck_assembly.Mt_v1.NODE_35.ffn",
		pcan2="../../data/processed/Pcanis/P_dog_Ck_assembly.Mt_v1.NODE_34.ffn",
		pcan3="../../data/processed/Pcanis/P_dog_Austria_assembly_v2.Mt.ffn",
		pcar1="../../data/processed/Pcar/JX499145.1.ffn",
		pcar2="../../data/processed/Pcar/NC_013660.1.ffn",
		pcar3="../../data/processed/Pcar/NC_013660.2.ffn",
		pmur="../../data/processed/Pmur/JX499144.1.ffn",
		pwk="../../data/processed/Pwk/Pwk_Mt_genome_v1_renamed.ffn"

	run:
		shell("sed \'s/>/>Pj_Ar_/\' {input.pj1} > {output.pj1}")
		shell("sed \'s/>/>Pj_CH26_/\' {input.pj2} > {output.pj2}")
		shell("sed \'s/>/>Pj_DK07_/\' {input.pj3} > {output.pj3}")
		shell("sed \'s/>/>Pj_GR11b_/\' {input.pj4} > {output.pj4}")
		shell("sed \'s/>/>Pj_GR5_/\' {input.pj5} > {output.pj5}")
		shell("sed \'s/>/>Pj_IT918_/\' {input.pj6} > {output.pj6}")
		shell("sed \'s/>/>Pj_JX855936_/\' {input.pj7} > {output.pj7}")
		shell("sed \'s/>/>Pj_JX855937_/\' {input.pj8} > {output.pj8}")
		shell("sed \'s/>/>Pj_JX855938_/\' {input.pj9} > {output.pj9}")
		shell("sed \'s/>/>Pj_S567_/\' {input.pj10} > {output.pj10}")
		shell("sed \'s/>/>Pj_SW1_/\' {input.pj11} > {output.pj11}")
		shell("sed \'s/>/>Pj_SW4_/\' {input.pj12} > {output.pj12}")
		shell("sed \'s/>/>Pj_SW7_/\' {input.pj13} > {output.pj13}")
		shell("sed \'s/>/>Pj_SE8_/\' {input.pj14} > {output.pj14}")
		shell("sed \'s/>/>Pj_46_/\' {input.pj15} > {output.pj15}")
		shell("sed \'s/>/>Pj_54c_/\' {input.pj16} > {output.pj16}")
		shell("sed \'s/>/>Pj_U_/\' {input.pj17} > {output.pj17}")	
		shell("sed \'s/>/>Pmac_H835_/\' {input.pmac1} > {output.pmac1}")
		shell("sed \'s/>/>Pmac_P2C_/\' {input.pmac2} > {output.pmac2}")
		shell("sed \'s/>/>Pmac_ER17_/\' {input.pmac3} > {output.pmac3}")
		shell("sed \'s/>/>Pmac_UC86_/\' {input.pmac4} > {output.pmac4}")	
		shell("sed \'s/>/>Po_CS1_/\' {input.pory1} > {output.pory1}")
		shell("sed \'s/>/>Po_rab1_/\' {input.pory2} > {output.pory2}")
		shell("sed \'s/>/>Po_rabF_/\' {input.pory3} > {output.pory3}")
		shell("sed \'s/>/>Po_rabM_/\' {input.pory4} > {output.pory4}")
		shell("sed \'s/>/>Pcan_1_/\' {input.pcan1} > {output.pcan1}")
		shell("sed \'s/>/>Pcan_2_/\' {input.pcan2} > {output.pcan2}")
		shell("sed \'s/>/>Pcan_3_/\' {input.pcan3} > {output.pcan3}")
		shell("sed \'s/>/>Pcar_1_/\' {input.pcar1} > {output.pcar1}")
		shell("sed \'s/>/>Pcar_2_/\' {input.pcar2} > {output.pcar2}")
		shell("sed \'s/>/>Pcar_3_/\' {input.pcar3} > {output.pcar3}")
		shell("sed \'s/>/>Pmur_1_/\' {input.pmur} > {output.pmur}")
		shell("sed \'s/>/>Pwk_1_/\' {input.pwk} > {output.pwk}")
rule combine_cds:
	input:
		pj1=rules.add_species_tag.output.pj1,pj2=rules.add_species_tag.output.pj2,pj3=rules.add_species_tag.output.pj3,
		pj4=rules.add_species_tag.output.pj4,pj5=rules.add_species_tag.output.pj5,pj6=rules.add_species_tag.output.pj6,
		pj7=rules.add_species_tag.output.pj7,pj8=rules.add_species_tag.output.pj8,pj9=rules.add_species_tag.output.pj9,
		pj10=rules.add_species_tag.output.pj10,pj11=rules.add_species_tag.output.pj11,pj12=rules.add_species_tag.output.pj12,
		pj13=rules.add_species_tag.output.pj13,pj14=rules.add_species_tag.output.pj14,pj15=rules.add_species_tag.output.pj15,
		pj16=rules.add_species_tag.output.pj16,pj17=rules.add_species_tag.output.pj17,
		pmac1=rules.add_species_tag.output.pmac1,pmac2=rules.add_species_tag.output.pmac2,pmac3=rules.add_species_tag.output.pmac3,
		pmac4=rules.add_species_tag.output.pmac4,
		pory1=rules.add_species_tag.output.pory1,pory2=rules.add_species_tag.output.pory2,
		pory3=rules.add_species_tag.output.pory3,pory4=rules.add_species_tag.output.pory4,
		pcan1=rules.add_species_tag.output.pcan1,pcan2=rules.add_species_tag.output.pcan2,pcan3=rules.add_species_tag.output.pcan3,
		pcar1=rules.add_species_tag.output.pcar1,pcar2=rules.add_species_tag.output.pcar2,pcar3=rules.add_species_tag.output.pcar3,
		pmur=rules.add_species_tag.output.pmur,
		pwk=rules.add_species_tag.output.pwk,
	output:
		a="../../data/processed/all_genes.cds.fa",
		rnlt=temp("../../data/processed/all_genes.cds.rnl.tmp"),rnl="../../data/processed/all_genes.cds.rnl.fa",
		rnst=temp("../../data/processed/all_genes.cds.rns.tmp"),rns="../../data/processed/all_genes.cds.rns.fa",
		rnpBt=temp("../../data/processed/all_genes.cds.rnpB.tmp"),rnpB="../../data/processed/all_genes.cds.rnpB.fa",
		cox1t=temp("../../data/processed/all_genes.cds.cox1.tmp"),cox1="../../data/processed/all_genes.cds.cox1.fa",
		cox2t=temp("../../data/processed/all_genes.cds.cox2.tmp"),cox2="../../data/processed/all_genes.cds.cox2.fa",
		cox3t=temp("../../data/processed/all_genes.cds.cox3.tmp"),cox3="../../data/processed/all_genes.cds.cox3.fa",	
		atp6t=temp("../../data/processed/all_genes.cds.atp6.tmp"),atp6="../../data/processed/all_genes.cds.atp6.fa",
		atp8t=temp("../../data/processed/all_genes.cds.atp8.tmp"),atp8="../../data/processed/all_genes.cds.atp8.fa",
		atp9t=temp("../../data/processed/all_genes.cds.atp9.tmp"),atp9="../../data/processed/all_genes.cds.atp9.fa",		
		nad1t=temp("../../data/processed/all_genes.cds.nad1.tmp"),nad1="../../data/processed/all_genes.cds.nad1.fa",
		nad2t=temp("../../data/processed/all_genes.cds.nad2.tmp"),nad2="../../data/processed/all_genes.cds.nad2.fa",
		nad3t=temp("../../data/processed/all_genes.cds.nad3.tmp"),nad3="../../data/processed/all_genes.cds.nad3.fa",
		nad4Lt=temp("../../data/processed/all_genes.cds.nad4L.tmp"),nad4L="../../data/processed/all_genes.cds.nad4L.fa",
		nad4t=temp("../../data/processed/all_genes.cds.nad4.tmp"),nad4="../../data/processed/all_genes.cds.nad4.fa",
		nad5t=temp("../../data/processed/all_genes.cds.nad5.tmp"),nad5="../../data/processed/all_genes.cds.nad5.fa",
		nad6t=temp("../../data/processed/all_genes.cds.nad6.tmp"),nad6="../../data/processed/all_genes.cds.nad6.fa",
		cobt=temp("../../data/processed/all_genes.cds.cob.tmp"),cob="../../data/processed/all_genes.cds.cob.fa"

	run:
		shell("cat {input.pj1} {input.pj2} {input.pj3} {input.pj4} {input.pj5} {input.pj6} {input.pj7} {input.pj8} {input.pj9} "
			  " {input.pj10} {input.pj11} {input.pj12} {input.pj13} {input.pj14} {input.pj15} {input.pj16} {input.pj17} "
			  " {input.pmac1} {input.pmac2} {input.pmac3} {input.pmac4} "
			  " {input.pory1} {input.pory2} {input.pory3} {input.pory4} "
			  " {input.pcan1} {input.pcan2} {input.pcan3} "
			  " {input.pcar1} {input.pcar2} {input.pcar3} {input.pmur} {input.pwk} > {output.a}")

		shell("grep \'rnl\' {output.a} | cut -f 2 -d \'>\' > {output.rnlt}")
		shell("seqtk subseq {output.a} {output.rnlt} > {output.rnl}")
		shell("grep \'rns\' {output.a} | cut -f 2 -d \'>\' > {output.rnst}")
		shell("seqtk subseq {output.a} {output.rnst} > {output.rns}")
		shell("grep \'rnpB\' {output.a} | cut -f 2 -d \'>\' > {output.rnpBt}")
		shell("seqtk subseq {output.a} {output.rnpBt} > {output.rnpB}")
		shell("grep \'cox1\' {output.a} | cut -f 2 -d \'>\' > {output.cox1t}")
		shell("seqtk subseq {output.a} {output.cox1t} > {output.cox1}")
		shell("grep \'cox2\' {output.a} | cut -f 2 -d \'>\' > {output.cox2t}")
		shell("seqtk subseq {output.a} {output.cox2t} > {output.cox2}")
		shell("grep \'cox3\' {output.a} | cut -f 2 -d \'>\' > {output.cox3t}")
		shell("seqtk subseq {output.a} {output.cox3t} > {output.cox3}")
		shell("grep \'atp6\' {output.a} | cut -f 2 -d \'>\' > {output.atp6t}")
		shell("seqtk subseq {output.a} {output.atp6t} > {output.atp6}")
		shell("grep \'atp8\' {output.a} | cut -f 2 -d \'>\' > {output.atp8t}")
		shell("seqtk subseq {output.a} {output.atp8t} > {output.atp8}")
		shell("grep \'atp9\' {output.a} | cut -f 2 -d \'>\' > {output.atp9t}")
		shell("seqtk subseq {output.a} {output.atp9t} > {output.atp9}")
		shell("grep \'nad1\' {output.a} | cut -f 2 -d \'>\' > {output.nad1t}")
		shell("seqtk subseq {output.a} {output.nad1t} > {output.nad1}")
		shell("grep \'nad2\' {output.a} | cut -f 2 -d \'>\' > {output.nad2t}")
		shell("seqtk subseq {output.a} {output.nad2t} > {output.nad2}")
		shell("grep \'nad3\' {output.a} | cut -f 2 -d \'>\' > {output.nad3t}")
		shell("seqtk subseq {output.a} {output.nad3t} > {output.nad3}")
		shell("grep -E \'nad4\s\' {output.a} | cut -f 2 -d \'>\' > {output.nad4t}")
		shell("seqtk subseq {output.a} {output.nad4t} > {output.nad4}")
		shell("grep \'nad4L\' {output.a} | cut -f 2 -d \'>\' > {output.nad4Lt}")
		shell("seqtk subseq {output.a} {output.nad4Lt} > {output.nad4L}")
		shell("grep \'nad5\' {output.a} | cut -f 2 -d \'>\' > {output.nad5t}")
		shell("seqtk subseq {output.a} {output.nad5t} > {output.nad5}")
		shell("grep \'nad6\' {output.a} | cut -f 2 -d \'>\' > {output.nad6t}")
		shell("seqtk subseq {output.a} {output.nad6t} > {output.nad6}")
		shell("grep \'cob\' {output.a} | cut -f 2 -d \'>\' > {output.cobt}")
		shell("seqtk subseq {output.a} {output.cobt} > {output.cob}")

# these two are off (some genome miss some genes)
#../../data/processed/all_genes.cds.nad4.fa:23
#../../data/processed/all_genes.cds.nad4L.fa:32 - 
# Ok I checked this and they all have the proper orientation but it's too slow so I am skipping for now
rule check_orientation:
	input:
		rnl=rules.combine_cds.output.rnl,rns=rules.combine_cds.output.rns,
		rnpB=rules.combine_cds.output.rnpB,
		cox1=rules.combine_cds.output.cox1,cox2=rules.combine_cds.output.cox2,cox3=rules.combine_cds.output.cox3,
		atp6=rules.combine_cds.output.atp6,atp8=rules.combine_cds.output.atp8,atp9=rules.combine_cds.output.atp9,
		nad1=rules.combine_cds.output.nad1,nad2=rules.combine_cds.output.nad2,nad3=rules.combine_cds.output.nad3,
		nad5=rules.combine_cds.output.nad5,nad6=rules.combine_cds.output.nad6,cob=rules.combine_cds.output.cob,
	output:
		rnl="../../data/processed/all_genes.cds.rnl.clean.fa",rns="../../data/processed/all_genes.cds.rns.clean.fa",
		rnpB="../../data/processed/all_genes.cds.rnpB.clean.fa",
		cox1="../../data/processed/all_genes.cds.cox1.clean.fa",cox2="../../data/processed/all_genes.cds.cox2.clean.fa",
		cox3="../../data/processed/all_genes.cds.cox3.clean.fa",	
		atp6="../../data/processed/all_genes.cds.atp6.clean.fa",atp8="../../data/processed/all_genes.cds.atp8.clean.fa",
		atp9="../../data/processed/all_genes.cds.atp9.clean.fa",		
		nad1="../../data/processed/all_genes.cds.nad1.clean.fa",nad2="../../data/processed/all_genes.cds.nad2.clean.fa",
		nad3="../../data/processed/all_genes.cds.nad3.clean.fa",nad5="../../data/processed/all_genes.cds.nad5.clean.fa",
		nad6="../../data/processed/all_genes.cds.nad6.clean.fa",
		cob="../../data/processed/all_genes.cds.cob.clean.fa"

	run:
		shell("perl ../scripts/SeqOrient.pl {input.rnl} > {output.rnl}")
		shell("perl ../scripts/SeqOrient.pl {input.rns} > {output.rns}")
		shell("perl ../scripts/SeqOrient.pl {input.cox1} > {output.cox1}")
		shell("perl ../scripts/SeqOrient.pl {input.cox2} > {output.cox2}")
		shell("perl ../scripts/SeqOrient.pl {input.cox3} > {output.cox3}")
		shell("perl ../scripts/SeqOrient.pl {input.atp6} > {output.atp6}")
		shell("perl ../scripts/SeqOrient.pl {input.atp8} > {output.atp8}")
		shell("perl ../scripts/SeqOrient.pl {input.atp9} > {output.atp9}")
		shell("perl ../scripts/SeqOrient.pl {input.nad1} > {output.nad1}")
		shell("perl ../scripts/SeqOrient.pl {input.nad2} > {output.nad2}")
		shell("perl ../scripts/SeqOrient.pl {input.nad3} > {output.nad3}")
		shell("perl ../scripts/SeqOrient.pl {input.nad5} > {output.nad5}")
		shell("perl ../scripts/SeqOrient.pl {input.nad6} > {output.nad6}")
		shell("perl ../scripts/SeqOrient.pl {input.cob} > {output.cob}")

rule prepare_aln_for_concat:
	input:
		rnl=rules.combine_cds.output.rnl,rns=rules.combine_cds.output.rns,
		rnpB=rules.combine_cds.output.rnpB,
		cox1=rules.combine_cds.output.cox1,cox2=rules.combine_cds.output.cox2,cox3=rules.combine_cds.output.cox3,
		atp6=rules.combine_cds.output.atp6,atp8=rules.combine_cds.output.atp8,atp9=rules.combine_cds.output.atp9,
		nad1=rules.combine_cds.output.nad1,nad2=rules.combine_cds.output.nad2,nad3=rules.combine_cds.output.nad3,
		nad5=rules.combine_cds.output.nad5,nad6=rules.combine_cds.output.nad6,cob=rules.combine_cds.output.cob,
	output:
		rnl=temp("../../data/processed/all_genes.cds.rnl.clean.spID.fa"),rns=temp("../../data/processed/all_genes.cds.rns.clean.spID.fa"),
		rnpB=temp("../../data/processed/all_genes.cds.rnpB.clean.spID.fa"),
		cox1=temp("../../data/processed/all_genes.cds.cox1.clean.spID.fa"),cox2=temp("../../data/processed/all_genes.cds.cox2.clean.spID.fa"),
		cox3=temp("../../data/processed/all_genes.cds.cox3.clean.spID.fa"),	
		atp6=temp("../../data/processed/all_genes.cds.atp6.clean.spID.fa"),atp8=temp("../../data/processed/all_genes.cds.atp8.clean.spID.fa"),
		atp9=temp("../../data/processed/all_genes.cds.atp9.clean.spID.fa"),		
		nad1=temp("../../data/processed/all_genes.cds.nad1.clean.spID.fa"),nad2=temp("../../data/processed/all_genes.cds.nad2.clean.spID.fa"),
		nad3=temp("../../data/processed/all_genes.cds.nad3.clean.spID.fa"),nad5=temp("../../data/processed/all_genes.cds.nad5.clean.spID.fa"),
		nad6=temp("../../data/processed/all_genes.cds.nad6.clean.spID.fa"),
		cob=temp("../../data/processed/all_genes.cds.cob.clean.spID.fa")
	run:
		shell("perl ../scripts/chop_headers.pl {input.rnl} > {output.rnl}")
		shell("perl ../scripts/chop_headers.pl {input.rns} > {output.rns}")
		shell("perl ../scripts/chop_headers.pl {input.rnpB} > {output.rnpB}")
		shell("perl ../scripts/chop_headers.pl {input.cox1} > {output.cox1}")
		shell("perl ../scripts/chop_headers.pl {input.cox2} > {output.cox2}")
		shell("perl ../scripts/chop_headers.pl {input.cox3} > {output.cox3}")
		shell("perl ../scripts/chop_headers.pl {input.atp6} > {output.atp6}")
		shell("perl ../scripts/chop_headers.pl {input.atp8} > {output.atp8}")
		shell("perl ../scripts/chop_headers.pl {input.atp9} > {output.atp9}")
		shell("perl ../scripts/chop_headers.pl {input.nad1} > {output.nad1}")
		shell("perl ../scripts/chop_headers.pl {input.nad2} > {output.nad2}")
		shell("perl ../scripts/chop_headers.pl {input.nad3} > {output.nad3}")
		shell("perl ../scripts/chop_headers.pl {input.nad5} > {output.nad5}")
		shell("perl ../scripts/chop_headers.pl {input.nad6} > {output.nad6}")
		shell("perl ../scripts/chop_headers.pl {input.cob} > {output.cob}")

rule align_msa:
	input:
		rnl=rules.prepare_aln_for_concat.output.rnl,rns=rules.prepare_aln_for_concat.output.rns,
		rnpB=rules.prepare_aln_for_concat.output.rnpB,
		cox1=rules.prepare_aln_for_concat.output.cox1,cox2=rules.prepare_aln_for_concat.output.cox2,cox3=rules.prepare_aln_for_concat.output.cox3,
		atp6=rules.prepare_aln_for_concat.output.atp6,atp8=rules.prepare_aln_for_concat.output.atp8,atp9=rules.prepare_aln_for_concat.output.atp9,
		nad1=rules.prepare_aln_for_concat.output.nad1,nad2=rules.prepare_aln_for_concat.output.nad2,nad3=rules.prepare_aln_for_concat.output.nad3,
		nad5=rules.prepare_aln_for_concat.output.nad5,nad6=rules.prepare_aln_for_concat.output.nad6,cob=rules.prepare_aln_for_concat.output.cob,
	output:
		rnl="../../data/processed/all_genes.cds.rnl.clean.spID.aln",
		rns="../../data/processed/all_genes.cds.rns.clean.spID.aln",
		rnpB="../../data/processed/all_genes.cds.rnpB.clean.spID.aln",
		cox1="../../data/processed/all_genes.cds.cox1.clean.spID.aln",
		cox2="../../data/processed/all_genes.cds.cox2.clean.spID.aln",
		cox3="../../data/processed/all_genes.cds.cox3.clean.spID.aln",	
		atp6="../../data/processed/all_genes.cds.atp6.clean.spID.aln",
		atp8="../../data/processed/all_genes.cds.atp8.clean.spID.aln",
		atp9="../../data/processed/all_genes.cds.atp9.clean.spID.aln",		
		nad1="../../data/processed/all_genes.cds.nad1.clean.spID.aln",
		nad2="../../data/processed/all_genes.cds.nad2.clean.spID.aln",
		nad3="../../data/processed/all_genes.cds.nad3.clean.spID.aln",
		nad5="../../data/processed/all_genes.cds.nad5.clean.spID.aln",
		nad6="../../data/processed/all_genes.cds.nad6.clean.spID.aln",
		cob="../../data/processed/all_genes.cds.cob.clean.spID.aln"

	run:
		shell("clustalo -i {input.rnl} -o {output.rnl}")
		shell("clustalo -i {input.rns} -o {output.rns}")	
		shell("clustalo -i {input.rnpB} -o {output.rnpB}")	
		shell("clustalo -i {input.cox1} -o {output.cox1}")
		shell("clustalo -i {input.cox2} -o {output.cox2}")
		shell("clustalo -i {input.cox3} -o {output.cox3}")
		shell("clustalo -i {input.atp6} -o {output.atp6}")
		shell("clustalo -i {input.atp8} -o {output.atp8}")
		shell("clustalo -i {input.atp9} -o {output.atp9}")
		shell("clustalo -i {input.nad1} -o {output.nad1}")
		shell("clustalo -i {input.nad2} -o {output.nad2}")
		shell("clustalo -i {input.nad3} -o {output.nad3}")
		shell("clustalo -i {input.nad5} -o {output.nad5}")
		shell("clustalo -i {input.nad6} -o {output.nad6}")
		shell("clustalo -i {input.cob} -o {output.cob}")

rule concat:
	input:
#		rnl=rules.align_msa.output.rnl,rns=rules.align_msa.output.rns,
#		rnpB=rules.align_msa.output.rnpB,
		cox1=rules.align_msa.output.cox1,cox2=rules.align_msa.output.cox2,cox3=rules.align_msa.output.cox3,
		atp6=rules.align_msa.output.atp6,atp8=rules.align_msa.output.atp8,atp9=rules.align_msa.output.atp9,
		nad1=rules.align_msa.output.nad1,nad2=rules.align_msa.output.nad2,nad3=rules.align_msa.output.nad3,
		nad5=rules.align_msa.output.nad5,nad6=rules.align_msa.output.nad6,cob=rules.align_msa.output.cob,
	output:
		"../../data/processed/all_mitogenes_concat.fasta"
	shell:
		"perl ../scripts/catfasta2phyml/catfasta2phyml.pl -f"
#		" {input.rnl} {input.rns} {input.rnpB}"
		" {input.cox1} {input.cox2} {input.cox3}"
		" {input.atp6} {input.atp8} {input.atp8} {input.atp9}"
		" {input.nad1} {input.nad2} {input.nad3} {input.nad5} {input.nad6} {input.cob} > {output}"

		
#rule Pomo:
#	input:
#		rules.concat.output
#	output:
		#"../../data/processed/all_mitogenes_concat.cf"
#	threads: 2
#	run:
		#shell("python ../tools/cflib/scripts/FastaToCounts.py {input} {output}") # not working yet
		#shell("iqtree -nt {threads} -s {input} -m HKY+P -bb 1000") #

#iqtree -s ../../data/processed/all_mitogenes_concat.fasta -alrt 1000 -bb 1000 -nt 2
	
		
		
		
		
		
		
		
		


### 
#    /Volumes/Shares/Papers_DATA/Pneumo_Comp_G/Mitochondria/Divergence_test/IQ-Tree-Pomo/cflib-master/scripts/FastaToCounts.py \
#    all_seqs_oriented.aln \
#    all_seqs_oriented.aln.cf
#~/utils/iqtree-1.6.beta4-Linux/bin/iqtree  -nt 6 -s all_seqs_oriented.aln.cf -m HKY+P -bb 1000
