rule all:
    input:
        "../../data/processed/minimap/scores.txt"

rule build_lastdb:
    input:
        pj="../../data/raw/Pjir/GCF_001477535.1_Pneu_jiro_RU7_V2_genomic.fna",
        pmac="../../data/processed/Pmac/Funannotate/PMAC_NAID_v2/predict_results/Pneumocystis_macacae_NIAID.scaffolds.fa",
        po="../../data/processed/Poryc/Funannotate/PORYCT_v0110_Ns_removed/predict_results/Pneumocystis_oryctolagi_MERGE.scaffolds.fa",
        ck1="../../data/processed/Pcanis/Funaannotate/Pneumocystis_canis_CK.scaffolds.fa",
        ck2="../../data/processed/Pcanis/Ck2/NCBI_sub/Pdg_scaffolds.cl6.alterv2_1kremoved_cl1_trim.fasta",
        pcanA="../../data/processed/Pcanis/dogA_Austria/Ass_V1/FUNANNOTATE/PCANB_v1/predict_results/Pneumocystis_canis_AUS.scaffolds.fa",
        pc="../../data/raw/Pcar/GCF_001477545.1_Pneu_cari_B80_V3_genomic.fna",
        pm="../../data/raw/Pmur/GCF_000349005.2_Pneumo_murina_B123_V4_genomic.fna",
        pwk="../../data/processed/Pwk/Funannotate/Pneumocystis_wakefieldiae_MERGE.scaffolds.fa"
    output:
        pj="../../data/processed/last/db/pjdb",pmac="../../data/processed/last/db/pmacdb",po="../../data/processed/last/db/podb",
        ck1="../../data/processed/last/db/pcanck1db",ck2="../../data/processed/last/db/pcanck2db",
        pcanA="../../data/processed/last/db/pcanadb",
        pc="../../data/processed/last/db/pcardb",pm="../../data/processed/last/db/pmurdb",pwk="../../data/processed/last/db/pwkdb",
    run:
        shell("echo > {output.pj}")
        shell("echo > {output.pmac}")
        shell("echo > {output.po}")
        shell("echo > {output.ck1}")
        shell("echo > {output.ck2}")
        shell("echo > {output.pcanA}")
        shell("echo > {output.pc}")
        shell("echo > {output.pm}")
        shell("echo > {output.pwk}")
        shell("lastdb -cR01 {output.pj} {input.pj}")
        shell("lastdb -cR01 {output.pmac} {input.pmac}")
        shell("lastdb -cR01 {output.po} {input.po}")
        shell("lastdb -cR01 {output.ck1} {input.ck1}")
        shell("lastdb -cR01 {output.ck2} {input.ck2}")
        shell("lastdb -cR01 {output.pcanA} {input.pcanA}")
        shell("lastdb -cR01 {output.pc} {input.pc}")
        shell("lastdb -cR01 {output.pm} {input.pm}")
        shell("lastdb -cR01 {output.pwk} {input.pwk}")

rule last_maf:
    input:
        pjref=rules.build_lastdb.output.pj, 
        pmacref=rules.build_lastdb.output.pmac,
        pjru7r1="../../data/raw/Pjir/Reads/RU7/SRR1043749_1.fastq",
        pjru7r2="../../data/raw/Pjir/Reads/RU7/SRR1043749_2.fastq",
        pj46r1="../../../Pj/Pj46/Pj46_mapped2PjRefG.1.fq",
        pj46r2="../../../Pj/Pj46/Pj46_mapped2PjRefG.2.fq",
        pj54cr1="../../../Pj/Pj54c/Pj54c_mapped2PjRefG.1.fq", 
        pj54cr2="../../../Pj/Pj54c/Pj54c_mapped2PjRefG.2.fq",
        pj55="../../../Pj/Pj55/Pj55_aligned2Ref.fq",
        pjzr1="/nethome/cisseoh/DATA/Pneumocystis_data/Reads/Pjirovecii/Genome/Pj_reads_Broad/PJ_Z_HighANDLow_merged_1.fq",
        pjzr2="/nethome/cisseoh/DATA/Pneumocystis_data/Reads/Pjirovecii/Genome/Pj_reads_Broad/PJ_Z_HighANDLow_merged_2.fq",
        pjwr1="/nethome/cisseoh/DATA/Pneumocystis_data/Reads/Pjirovecii/Genome/Pj_reads_Broad/PJ_W_HighANDLow_merged_1.fq",
        pjwr2="/nethome/cisseoh/DATA/Pneumocystis_data/Reads/Pjirovecii/Genome/Pj_reads_Broad/PJ_W_HighANDLow_merged_2.fq", 
        m2r1="/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/Reads/Pmacacae/BroadSequencing/M2_r1.fq.gz",
        m2r2="/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/Reads/Pmacacae/BroadSequencing/M2_r2.fq.gz",
        mp3r1="/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/Reads/Pmacacae/BroadSequencing/MP3_r1.fq.gz",
        mp3r2="/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/Reads/Pmacacae/BroadSequencing/MP3_r2.fq.gz",
        h835r1="/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/Reads/Pmacacae/BroadSequencing/Pneumocystis_macaca_H835_r1.fastq.gz",
        h835r2="/hpcdata/dcr_cc/Ousmane_Data/DATA/Pneumocystis_data/Reads/Pmacacae/BroadSequencing/Pneumocystis_macaca_H835_r2.fastq.gz",
        p2c1r1="../../data/raw/Pmac/DNA/Illumina/NIAID_2018/run_180809/P2C_1.fq",
        p2c1r2="../../data/raw/Pmac/DNA/Illumina/NIAID_2018/run_180809/P2C_2.fq",
        s1a2r1="../../data/raw/Pmac/DNA/Illumina/NIAID_2018/run_180809/S1A2_1.fq",
        s1a2r2="../../data/raw/Pmac/DNA/Illumina/NIAID_2018/run_180809/S1A2_2.fq",
#        p2c2r1="../../data/raw/Pmac/DNA/Illumina/NIAID_2018/data_release-LS08131802/P2C_150/P2C_1.fq.gz",
#        p2c2r2="../../data/raw/Pmac/DNA/Illumina/NIAID_2018/data_release-LS08131802/P2C_150/P2C_2.fq.gz",
#        p2c3r1="../../data/raw/Pmac/DNA/Illumina/NIAID_2018/data_release-LS08131802/P2C_250/P2C_1.fq.gz",
#        p2c3r2="../../data/raw/Pmac/DNA/Illumina/NIAID_2018/data_release-LS08131802/P2C_250/P2C_2.fq.gz",
        cj36r1="../../data/raw/Pmac/DNA/Illumina/TULANE/CJ36/CJ36_1.fq",
        cj36r2="../../data/raw/Pmac/DNA/Illumina/TULANE/CJ36/CJ36_2.fq",
        er17r1="../../data/raw/Pmac/DNA/Illumina/TULANE/ER17/ER17_allRuns_1.fq",
        er17r2="../../data/raw/Pmac/DNA/Illumina/TULANE/ER17/ER17_allRuns_2.fq",
        gl92r1="../../data/raw/Pmac/DNA/Illumina/TULANE/GL92/GL92_USD16091310L_HL2H3DSXX_L1_1.fq",
        gl92r2="../../data/raw/Pmac/DNA/Illumina/TULANE/GL92/GL92_USD16091310L_HL2H3DSXX_L1_2.fq",
        uc86r1="../../data/raw/Pmac/DNA/Illumina/TULANE/UC86/UC86_allRuns_1.fq",
        uc86r2="../../data/raw/Pmac/DNA/Illumina/TULANE/UC86/UC86_allRuns_2.fq"
    output:
        pjru7m1pj=protected("../../data/processed/last/maf/vsPj/RU7_1.maf"),
        pjru7m2pj=protected("../../data/processed/last/maf/vsPj/RU7_2.maf"),
        pjru7m1pmac=protected("../../data/processed/last/maf/vsPmac/RU7_1.maf"),
        pjru7m2pmac=protected("../../data/processed/last/maf/vsPmac/RU7_2.maf"),
        pj46m1pj=protected("../../data/processed/last/maf/vsPj/Pj46_1.maf"),
        pj46m2pj=protected("../../data/processed/last/maf/vsPj/Pj46_2.maf"),
        pj46m1pmac=protected("../../data/processed/last/maf/vsPmac/Pj46_1.maf"),
        pj46m2pmac=protected("../../data/processed/last/maf/vsPmac/Pj46_2.maf"),
        pj54cm1pj=protected("../../data/processed/last/maf/vsPj/Pj54c_1.maf"),
        pj54cm2pj=protected("../../data/processed/last/maf/vsPj/Pj54c_2.maf"),
        pj54cm1pmac=protected("../../data/processed/last/maf/vsPmac/Pj54c_1.maf"),
        pj54cm2pmac=protected("../../data/processed/last/maf/vsPmac/Pj54c_2.maf"),
        pj55m1pj=protected("../../data/processed/last/maf/vsPj/Pj55_1.maf"),
        pj55m1pmac=protected("../../data/processed/last/maf/vsPmac/Pj55_1.maf"),
        pjzm1pj=protected("../../data/processed/last/maf/vsPj/PjZ_1.maf"),
        pjzm2pj=protected("../../data/processed/last/maf/vsPj/PjZ_2.maf"),
        pjzm1pmac=protected("../../data/processed/last/maf/vsPmac/PjZ_1.maf"),
        pjzm2pmac=protected("../../data/processed/last/maf/vsPmac/PjZ_2.maf"),
        pjwm1pj=protected("../../data/processed/last/maf/vsPj/PjW_1.maf"),
        pjwm2pj=protected("../../data/processed/last/maf/vsPj/PjW_2.maf"),
        pjwm1pmac=protected("../../data/processed/last/maf/vsPmac/PjW_1.maf"),
        pjwm2pmac=protected("../../data/processed/last/maf/vsPmac/PjW_2.maf"),
        p2cm1pj=protected("../../data/processed/last/maf/vsPj/PmacP2C_1.maf"),
        p2cm2pj=protected("../../data/processed/last/maf/vsPj/PmacP2C_2.maf"),
        p2cm1pmac=protected("../../data/processed/last/maf/vsPmac/PmacP2C_1.maf"),
        p2cm2pmac=protected("../../data/processed/last/maf/vsPmac/PmacP2C_2.maf"),
        cj36m1pj=protected("../../data/processed/last/maf/vsPj/PmacCj36_1.maf"),
        cj36m2pj=protected("../../data/processed/last/maf/vsPj/PmacCj36_2.maf"),
        cj36m1pmac=protected("../../data/processed/last/maf/vsPmac/PmacCj36_1.maf"),
        cj36m2pmac=protected("../../data/processed/last/maf/vsPmac/PmacCj36_2.maf"),
        er17m1pj=protected("../../data/processed/last/maf/vsPj/PmacER17_1.maf"),
        er17m2pj=protected("../../data/processed/last/maf/vsPj/PmacER17_2.maf"),
        er17m1pmac=protected("../../data/processed/last/maf/vsPmac/PmacER17_1.maf"),
        er17m2pmac=protected("../../data/processed/last/maf/vsPmac/PmacER17_2.maf"),
        gl92m1pj=protected("../../data/processed/last/maf/vsPj/Pmacgl92_1.maf"),
        gl92m2pj=protected("../../data/processed/last/maf/vsPj/Pmacgl92_2.maf"),
        gl92m1pmac=protected("../../data/processed/last/maf/vsPmac/Pmacgl92_1.maf"),
        gl92m2pmac=protected("../../data/processed/last/maf/vsPmac/Pmacgl92_2.maf"),
        uc86m1pj=protected("../../data/processed/last/maf/vsPj/Pmacuc86_1.maf"),
        uc86m2pj=protected("../../data/processed/last/maf/vsPj/Pmacuc86_2.maf"),
        uc86m1pmac=protected("../../data/processed/last/maf/vsPmac/Pmacuc86_1.maf"),
        uc86m2pmac=protected("../../data/processed/last/maf/vsPmac/Pmacuc86_2.maf"),

    params: fq="-Q0 -r5 -q5 -a35 -b5 -i8"

    threads: 8

    run:
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pjru7r1} | last-split > {output.pjru7m1pj}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pjru7r2} | last-split > {output.pjru7m2pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pjru7r1} | last-split > {output.pjru7m1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pjru7r2} | last-split > {output.pjru7m2pmac}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pj46r1} | last-split > {output.pj46m1pj}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pj46r2} | last-split > {output.pj46m2pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pj46r1} | last-split > {output.pj46m1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pj46r2} | last-split > {output.pj46m2pmac}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pj54cr1} | last-split > {output.pj54cm1pj}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pj54cr2} | last-split > {output.pj54cm2pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pj54cr1} | last-split > {output.pj54cm1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pj54cr2} | last-split > {output.pj54cm2pmac}") 
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pj55} | last-split > {output.pj55m1pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pj55} | last-split > {output.pj55m1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pjzr1} | last-split > {output.pjzm1pj}") 
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pjzr2} | last-split > {output.pjzm2pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pjzr1} | last-split > {output.pjzm1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pjzr2} | last-split > {output.pjzm2pmac}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pjwr1} | last-split > {output.pjwm1pj}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.pjwr2} | last-split > {output.pjwm2pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pjwr1} | last-split > {output.pjwm1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.pjwr2} | last-split > {output.pjwm2pmac}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.p2c1r1} | last-split > {output.p2cm1pj}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.p2c1r2} | last-split > {output.p2cm2pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.p2c1r1} | last-split > {output.p2cm1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.p2c1r2} | last-split > {output.p2cm2pmac}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.cj36r1} | last-split > {output.cj36m1pj}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.cj36r2} | last-split > {output.cj36m2pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.cj36r1} | last-split > {output.cj36m1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.cj36r2} | last-split > {output.cj36m2pmac}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.er17r1} | last-split > {output.er17m1pj}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.er17r2} | last-split > {output.er17m2pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.er17r1} | last-split > {output.er17m1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.er17r2} | last-split > {output.er17m2pmac}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.gl92r1} | last-split > {output.gl92m1pj}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.gl92r2} | last-split > {output.gl92m2pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.gl92r1} | last-split > {output.gl92m1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.gl92r2} | last-split > {output.gl92m2pmac}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.uc86r1} | last-split > {output.uc86m1pj}")
        shell("lastal {params.fq} -P {threads} {input.pjref} {input.uc86r2} | last-split > {output.uc86m2pj}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.uc86r1} | last-split > {output.uc86m1pmac}")
        shell("lastal {params.fq} -P {threads} {input.pmacref} {input.uc86r2} | last-split > {output.uc86m2pmac}")

# module load samtools 1.9 -- create the config
rule maf2bam:
    input:
        pj="../../data/raw/Pjir/GCF_001477535.1_Pneu_jiro_RU7_V2_genomic.fna.fai",
        pmac="../../data/processed/Pmac/Funannotate/PMAC_NAID_v2/predict_results/Pneumocystis_macacae_NIAID.scaffolds.fa.fai",
        po="../../data/processed/Poryc/Funannotate/PORYCT_v0110_Ns_removed/predict_results/Pneumocystis_oryctolagi_MERGE.scaffolds.fa.fai",
        ck1="../../data/processed/Pcanis/Funaannotate/Pneumocystis_canis_CK.scaffolds.fa.fai",
        ck2="../../data/processed/Pcanis/Ck2/NCBI_sub/Pdg_scaffolds.cl6.alterv2_1kremoved_cl1_trim.fasta.fai",
        pcanA="../../data/processed/Pcanis/dogA_Austria/Ass_V1/FUNANNOTATE/PCANB_v1/predict_results/Pneumocystis_canis_AUS.scaffolds.fa.fai",
        pc="../../data/raw/Pcar/GCF_001477545.1_Pneu_cari_B80_V3_genomic.fna.fai",
        pm="../../data/raw/Pmur/GCF_000349005.2_Pneumo_murina_B123_V4_genomic.fna.fai",
        pwk="../../data/processed/Pwk/Funannotate/Pneumocystis_wakefieldiae_MERGE.scaffolds.fa.fai",
        pjru7m1pj=rules.last_maf.output.pjru7m1pj,pjru7m2pj=rules.last_maf.output.pjru7m2pj,
        pjru7m1pmac=rules.last_maf.output.pjru7m1pmac,pjru7m2pmac=rules.last_maf.output.pjru7m2pmac,
        pj46m1pj=rules.last_maf.output.pj46m1pj,pj46m2pj=rules.last_maf.output.pj46m2pj,
        pj46m1pmac=rules.last_maf.output.pj46m1pmac,pj46m2pmac=rules.last_maf.output.pj46m2pmac,
        pj54cm1pj=rules.last_maf.output.pj54cm1pj,pj54cm2pj=rules.last_maf.output.pj54cm2pj,
        pj54cm1pmac=rules.last_maf.output.pj54cm1pmac,pj54cm2pmac=rules.last_maf.output.pj54cm2pmac,
        pj55m1pj=rules.last_maf.output.pj55m1pj,
        pj55m1pmac=rules.last_maf.output.pj55m1pmac,
        pjzm1pj=rules.last_maf.output.pjzm1pj,pjzm2pj=rules.last_maf.output.pjzm2pj,
        pjzm1pmac=rules.last_maf.output.pjzm1pmac,pjzm2pmac=rules.last_maf.output.pjzm2pmac,
        pjwm1pj=rules.last_maf.output.pjwm1pj,pjwm2pj=rules.last_maf.output.pjwm2pj,
        pjwm1pmac=rules.last_maf.output.pjwm1pmac,pjwm2pmac=rules.last_maf.output.pjwm2pmac,
        p2cm1pj=rules.last_maf.output.p2cm1pj,p2cm2pj=rules.last_maf.output.p2cm2pj,
        p2cm1pmac=rules.last_maf.output.p2cm1pmac,p2cm2pmac=rules.last_maf.output.p2cm2pmac,
        cj36m1pj=rules.last_maf.output.cj36m1pj,cj36m2pj=rules.last_maf.output.cj36m2pj,
        cj36m1pmac=rules.last_maf.output.cj36m1pmac,cj36m2pmac=rules.last_maf.output.cj36m2pmac,
        er17m1pj=rules.last_maf.output.er17m1pj,er17m2pj=rules.last_maf.output.er17m2pj,
        er17m1pmac=rules.last_maf.output.er17m1pmac,er17m2pmac=rules.last_maf.output.er17m2pmac,
        gl92m1pj=rules.last_maf.output.gl92m1pj,gl92m2pj=rules.last_maf.output.gl92m2pj,
        gl92m1pmac=rules.last_maf.output.gl92m1pmac,gl92m2pmac=rules.last_maf.output.gl92m2pmac,
        uc86m1pj=rules.last_maf.output.uc86m1pj,uc86m2pj=rules.last_maf.output.uc86m2pj,
        uc86m1pmac=rules.last_maf.output.uc86m1pmac,uc86m2pmac=rules.last_maf.output.uc86m2pmac,


    output:
        pjru7pj="../../data/processed/last/maf/vsPj/PjRU7.bam",
        pj46pj="../../data/processed/last/maf/vsPj/Pj46.bam",
        pj54cpj="../../data/processed/last/maf/vsPj/Pj54c.bam",
        pj55pj="../../data/processed/last/maf/vsPj/Pj55.bam",
        pjzpj="../../data/processed/last/maf/vsPj/PjZ.bam",
        pjwpj="../../data/processed/last/maf/vsPj/PjW.bam",
        p2cpj="../../data/processed/last/maf/vsPj/PmacP2C.bam",
        cj36pj="../../data/processed/last/maf/vsPj/PmacCj36.bam",
        er17pj="../../data/processed/last/maf/vsPj/PmacER17.bam",
        gl92pj="../../data/processed/last/maf/vsPj/Pmacgl92.bam",
        uc86pj="../../data/processed/last/maf/vsPj/Pmacuc86.bam",
        pjru7pmac="../../data/processed/last/maf/vsPmac/PjRU7.bam",
        pj46pmac="../../data/processed/last/maf/vsPmac/Pj46.bam",
        pj54cpmac="../../data/processed/last/maf/vsPmac/Pj54c.bam",
        pj55pmac="../../data/processed/last/maf/vsPmac/Pj55.bam",
        pjzpmac="../../data/processed/last/maf/vsPmac/PjZ.bam",
        pjwpmac="../../data/processed/last/maf/vsPmac/PjW.bam",
        p2cpmac="../../data/processed/last/maf/vsPmac/PmacP2C.bam",
        cj36pmac="../../data/processed/last/maf/vsPmac/PmacCj36.bam",
        er17pmac="../../data/processed/last/maf/vsPmac/PmacER17.bam",
        gl92pmac="../../data/processed/last/maf/vsPmac/Pmacgl92.bam",
        uc86pmac="../../data/processed/last/maf/vsPmac/Pmacuc86.bam",

    threads: 8

    conda:
        "../envs/samtools.yaml"

    run:
        shell("perl ../scripts/buddle {input.pj} {input.pjru7m1pj} {input.pjru7m2pj} {threads} {output.pjru7pj}")
        shell("perl ../scripts/buddle {input.pmac} {input.pjru7m1pmac} {input.pjru7m2pmac} {threads} {output.pjru7pmac}")
        shell("perl ../scripts/buddle {input.pj} {input.pj46m1pj} {input.pj46m2pj} {threads} {output.pj46pj}")
        shell("perl ../scripts/buddle {input.pmac} {input.pj46m1pmac} {input.pj46m2pmac} {threads} {output.pj46pmac}")
        shell("perl ../scripts/buddle {input.pj} {input.pj54cm1pj} {input.pj54cm2pj} {threads} {output.pj54cpj}") 
        shell("perl ../scripts/buddle {input.pmac} {input.pj54cm1pmac} {input.pj54cm2pmac} {threads} {output.pj54cpmac}")
        shell("perl ../scripts/buddle2 {input.pj} {input.pj55m1pj} {threads} {output.pj55pj}")
        shell("perl ../scripts/buddle2 {input.pmac} {input.pj55m1pmac} {threads} {output.pj55pmac}")
        shell("perl ../scripts/buddle {input.pj} {input.pjzm1pj} {input.pjzm2pj} {threads} {output.pjzpj}")
        shell("perl ../scripts/buddle {input.pmac} {input.pjzm1pmac} {input.pjzm2pmac} {threads} {output.pjzpmac}")
        shell("perl ../scripts/buddle {input.pj} {input.pjwm1pj} {input.pjwm2pj} {threads} {output.pjwpj}")
        shell("perl ../scripts/buddle {input.pmac} {input.pjwm1pmac} {input.pjwm2pmac} {threads} {output.pjwpmac}")
        shell("perl ../scripts/buddle {input.pj} {input.p2cm1pj} {input.p2cm2pj} {threads} {output.p2cpj}")
        shell("perl ../scripts/buddle {input.pmac} {input.p2cm1pmac} {input.p2cm2pmac} {threads} {output.p2cpmac}")
        shell("perl ../scripts/buddle {input.pj} {input.cj36m1pj} {input.cj36m2pj} {threads} {output.cj36pj}")
        shell("perl ../scripts/buddle {input.pmac} {input.cj36m1pmac} {input.cj36m2pmac} {threads} {output.cj36pmac}")
        shell("perl ../scripts/buddle {input.pj} {input.er17m1pj} {input.er17m2pj} {threads} {output.er17pj}")
        shell("perl ../scripts/buddle {input.pmac} {input.er17m1pmac} {input.er17m2pmac} {threads} {output.er17pmac}")
        shell("perl ../scripts/buddle {input.pj} {input.gl92m1pj} {input.gl92m2pj} {threads} {output.gl92pj}")
        shell("perl ../scripts/buddle {input.pmac} {input.gl92m1pmac} {input.gl92m2pmac} {threads} {output.gl92pmac}")
        shell("perl ../scripts/buddle {input.pj} {input.uc86m1pj} {input.uc86m2pj} {threads} {output.uc86pj}")
        shell("perl ../scripts/buddle {input.pmac} {input.uc86m1pmac} {input.uc86m2pmac} {threads} {output.uc86pmac}")

rule sort:
    input:
        pjru7pj=rules.maf2bam.output.pjru7pj,
        pj46pj=rules.maf2bam.output.pj46pj,
        pj54cpj=rules.maf2bam.output.pj54cpj,
        pj55pj=rules.maf2bam.output.pj55pj,
        pjzpj=rules.maf2bam.output.pjzpj,
        pjwpj=rules.maf2bam.output.pjwpj,
        p2cpj=rules.maf2bam.output.p2cpj,
        cj36pj=rules.maf2bam.output.cj36pj,
        er17pj=rules.maf2bam.output.er17pj,
        gl92pj=rules.maf2bam.output.gl92pj,
        uc86pj=rules.maf2bam.output.uc86pj,
        pjru7pmac=rules.maf2bam.output.pjru7pmac,
        pj46pmac=rules.maf2bam.output.pj46pmac,
        pj54cpmac=rules.maf2bam.output.pj54cpmac,
        pj55pmac=rules.maf2bam.output.pj55pmac,
        pjzpmac=rules.maf2bam.output.pjzpmac,
        pjwpmac=rules.maf2bam.output.pjwpmac,
        p2cpmac=rules.maf2bam.output.p2cpmac,
        cj36pmac=rules.maf2bam.output.cj36pmac,
        er17pmac=rules.maf2bam.output.er17pmac,
        gl92pmac=rules.maf2bam.output.gl92pmac,
        uc86pmac=rules.maf2bam.output.uc86pmac,
    output:
        pjru7pj="../../data/processed/last/maf/vsPj/PjRU7.sorted.bam",
        pj46pj="../../data/processed/last/maf/vsPj/Pj46.sorted.bam",
        pj54cpj="../../data/processed/last/maf/vsPj/Pj54c.sorted.bam",
        pj55pj="../../data/processed/last/maf/vsPj/Pj55.sorted.bam",
        pjzpj="../../data/processed/last/maf/vsPj/PjZ.sorted.bam",
        pjwpj="../../data/processed/last/maf/vsPj/PjW.sorted.bam",
        p2cpj="../../data/processed/last/maf/vsPj/PmacP2C.sorted.bam",
        cj36pj="../../data/processed/last/maf/vsPj/PmacCj36.sorted.bam",
        er17pj="../../data/processed/last/maf/vsPj/PmacER17.sorted.bam",
        gl92pj="../../data/processed/last/maf/vsPj/Pmacgl92.sorted.bam",
        uc86pj="../../data/processed/last/maf/vsPj/Pmacuc86.sorted.bam",
        pjru7pmac="../../data/processed/last/maf/vsPmac/PjRU7.sorted.bam",
        pj46pmac="../../data/processed/last/maf/vsPmac/Pj46.sorted.bam",
        pj54cpmac="../../data/processed/last/maf/vsPmac/Pj54c.sorted.bam",
        pj55pmac="../../data/processed/last/maf/vsPmac/Pj55.sorted.bam",
        pjzpmac="../../data/processed/last/maf/vsPmac/PjZ.sorted.bam",
        pjwpmac="../../data/processed/last/maf/vsPmac/PjW.sorted.bam",
        p2cpmac="../../data/processed/last/maf/vsPmac/PmacP2C.sorted.bam",
        cj36pmac="../../data/processed/last/maf/vsPmac/PmacCj36.sorted.bam",
        er17pmac="../../data/processed/last/maf/vsPmac/PmacER17.sorted.bam",
        gl92pmac="../../data/processed/last/maf/vsPmac/Pmacgl92.sorted.bam",
        uc86pmac="../../data/processed/last/maf/vsPmac/Pmacuc86.sorted.bam",

    threads: 8

    conda:
        "../envs/samtools.yaml"
    params: b="bam"

    run:
        shell("samtools sort -@ {threads} -O {params.b} {input.pjru7pj} -o {output.pjru7pj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pjru7pmac} -o {output.pjru7pmac}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pj46pj} -o {output.pj46pj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pj46pmac} -o {output.pj46pmac}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pj54cpj} -o {output.pj54cpj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pj54cpmac} -o {output.pj54cpmac}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pj55pj} -o {output.pj55pj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pj55pmac} -o {output.pj55pmac}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pjzpj} -o {output.pjzpj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pjzpmac} -o {output.pjzpmac}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pjwpj} -o {output.pjwpj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.pjwpmac} -o {output.pjwpmac}")
        shell("samtools sort -@ {threads} -O {params.b} {input.p2cpj} -o {output.p2cpj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.p2cpmac} -o {output.p2cpmac}")
        shell("samtools sort -@ {threads} -O {params.b} {input.cj36pj} -o {output.cj36pj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.cj36pmac} -o {output.cj36pmac}")
        shell("samtools sort -@ {threads} -O {params.b} {input.er17pj} -o {output.er17pj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.er17pmac} -o {output.er17pmac}")
        shell("samtools sort -@ {threads} -O {params.b} {input.gl92pj} -o {output.gl92pj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.gl92pmac} -o {output.gl92pmac}")
        shell("samtools sort -@ {threads} -O {params.b} {input.uc86pj} -o {output.uc86pj}")
        shell("samtools sort -@ {threads} -O {params.b} {input.uc86pmac} -o {output.uc86pmac}")
 
rule snp_call:
    input:
        pj="../../data/raw/Pjir/GCF_001477535.1_Pneu_jiro_RU7_V2_genomic.fna",
        pmac="../../data/processed/Pmac/Funannotate/PMAC_NAID_v2/predict_results/Pneumocystis_macacae_NIAID.scaffolds.fa",
        pjru7pj=rules.sort.output.pjru7pj,
        pj46pj=rules.sort.output.pj46pj,
        pj54cpj=rules.sort.output.pj54cpj,
        pj55pj=rules.sort.output.pj55pj,
        pjzpj=rules.sort.output.pjzpj,
        pjwpj=rules.sort.output.pjwpj,
        p2cpj=rules.sort.output.p2cpj,
        cj36pj=rules.sort.output.cj36pj,
        er17pj=rules.sort.output.er17pj,
        gl92pj=rules.sort.output.gl92pj,
        uc86pj=rules.sort.output.uc86pj,
        pjru7pmac=rules.sort.output.pjru7pmac,
        pj46pmac=rules.sort.output.pj46pmac,
        pj54cpmac=rules.sort.output.pj54cpmac,
        pj55pmac=rules.sort.output.pj55pmac,
        pjzpmac=rules.sort.output.pjzpmac,
        pjwpmac=rules.sort.output.pjwpmac,
        p2cpmac=rules.sort.output.p2cpmac,
        cj36pmac=rules.sort.output.cj36pmac,
        er17pmac=rules.sort.output.er17pmac,
        gl92pmac=rules.sort.output.gl92pmac,
        uc86pmac=rules.sort.output.uc86pmac,
    output:
#        pjru7pj="../../data/processed/last/vcf/vsPj/PjRU7.vcf.gz",
        pj46pj="../../data/processed/last/vcf/vsPj/Pj46.vcf.gz",
        pj54cpj="../../data/processed/last/vcf/vsPj/Pj54c.vcf.gz",
        pj55pj="../../data/processed/last/vcf/vsPj/Pj55.vcf.gz",
        pjzpj="../../data/processed/last/vcf/vsPj/PjZ.vcf.gz",
        pjwpj="../../data/processed/last/vcf/vsPj/PjW.vcf.gz",
#        p2cpj="../../data/processed/last/vcf/vsPj/PmacP2C.vcf.gz",
#        cj36pj="../../data/processed/last/vcf/vsPj/PmacCj36.vcf.gz",
#        er17pj="../../data/processed/last/vcf/vsPj/PmacER17.vcf.gz",
#        gl92pj="../../data/processed/last/vcf/vsPj/Pmacgl92.vcf.gz",
#        uc86pj="../../data/processed/last/vcf/vsPj/Pmacuc86.vcf.gz",
#        pjru7pmac="../../data/processed/last/vcf/vsPmac/PjRU7.vcf.gz",
#        pj46pmac="../../data/processed/last/vcf/vsPmac/Pj46.vcf.gz",
#        pj54cpmac="../../data/processed/last/vcf/vsPmac/Pj54c.vcf.gz",
#        pj55pmac="../../data/processed/last/vcf/vsPmac/Pj55.vcf.gz",
#        pjzpmac="../../data/processed/last/vcf/vsPmac/PjZ.vcf.gz",
#        pjwpmac="../../data/processed/last/vcf/vsPmac/PjW.vcf.gz",
#        p2cpmac="../../data/processed/last/vcf/vsPmac/PmacP2C.vcf.gz",
        cj36pmac="../../data/processed/last/vcf/vsPmac/PmacCj36.vcf.gz",
        er17pmac="../../data/processed/last/vcf/vsPmac/PmacER17.vcf.gz",
        gl92pmac="../../data/processed/last/vcf/vsPmac/Pmacgl92.vcf.gz",
        uc86pmac="../../data/processed/last/vcf/vsPmac/Pmacuc86.vcf.gz",

    threads: 8

    conda:
        "../envs/samtools.yaml"

    params:
        p="1", c="5", q="20" , g="10000", pl="1", v="indels" # indels are creating errors for cns
    run:
#        shell("bcftools mpileup -Ou -f {input.pj} {input.pjru7pj} | bcftools call -mv -Oz --ploidy {params.pl} -o {output.pjru7pj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.pjru7pj} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pjru7pj}")
        shell("bcftools mpileup -Ou -f {input.pj} {input.pj46pj} | bcftools call -mv --threads {threads} -V {params.v} -Oz --ploidy {params.pl} -o {output.pj46pj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.pj46pj} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pj46pj}")
        shell("bcftools mpileup -Ou -f {input.pj} {input.pj54cpj} | bcftools call -mv --threads {threads} -V {params.v} -Oz --ploidy {params.pl} -o {output.pj54cpj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.pj54cpj} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pj54cpj}")
        shell("bcftools mpileup -Ou -f {input.pj} {input.pj55pj} | bcftools call -mv --threads {threads} -V {params.v} -Oz --ploidy {params.pl} -o {output.pj55pj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.pj55pj}  -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pj55pj}")
        shell("bcftools mpileup -Ou -f {input.pj} {input.pjzpj} | bcftools call -mv --threads {threads} -V {params.v} -Oz --ploidy {params.pl} -o {output.pjzpj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.pjzpj} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pjzpj}")
        shell("bcftools mpileup -Ou -f {input.pj} {input.pjwpj} | bcftools call -mv --threads {threads} -V {params.v} -Oz --ploidy {params.pl} -o  {output.pjwpj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.pjwpj} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pjwpj}") 
#        shell("bcftools mpileup -Ou -f {input.pj} {input.p2cpj} | bcftools call -mv -Oz --ploidy {params.pl} -o {output.p2cpj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.p2cpj} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.p2cpj}")
#        shell("bcftools mpileup -Ou -f {input.pj} {input.cj36pj} | bcftools call -mv -Oz --ploidy {params.pl} -o {output.cj36pj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.cj36pj} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.cj36pj}")
#        shell("bcftools mpileup -Ou -f {input.pj} {input.er17pj}  | bcftools call -mv -Oz --ploidy {params.pl} -o {output.er17pj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.er17pj} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.er17pj}")
#        shell("bcftools mpileup -Ou -f {input.pj} {input.gl92pj} | bcftools call -mv -Oz --ploidy {params.pl} -o {output.gl92pj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.gl92pj} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.gl92pj}")
#        shell("bcftools mpileup -Ou -f {input.pj} input.uc86pj} | bcftools call -mv -Oz --ploidy {params.pl} -o {output.uc86pj}")
#        shell("freebayes -f {input.pj} -C {params.c} {input.uc86pj} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.uc86pj}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.pjru7pmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pjru7pmac}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.pj46pmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pj46pmac}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.pj54cpmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pj54cpmac}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.pj55pmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pj55pmac}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.pjzpmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pjzpmac}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.pjwpmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.pjwpmac}")
#        shell("bcftools mpileup -Ou -f {input.pmac} {input.p2cpmac} | bcftools call -mv -Oz --ploidy {params.pl} -o  {output.p2cpmac}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.p2cpmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.p2cpmac}")
        shell("bcftools mpileup -Ou -f {input.pmac} {input.cj36pmac} | bcftools call -mv -V {params.v} -Oz --threads {threads} --ploidy {params.pl} -o {output.cj36pmac}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.cj36pmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.cj36pmac}")
        shell("bcftools mpileup -Ou -f {input.pmac} {input.er17pmac} | bcftools call -mv -V {params.v} -Oz --threads {threads} --ploidy {params.pl} -o {output.er17pmac}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.er17pmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.er17pmac}")
        shell("bcftools mpileup -Ou -f {input.pmac} {input.gl92pmac} | bcftools call -mv -V {params.v} -Oz --threads {threads} --ploidy {params.pl} -o {output.gl92pmac}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.gl92pmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.gl92pmac}")
        shell("bcftools mpileup -Ou -f {input.pmac} {input.uc86pmac} | bcftools call -mv -V {params.v} -Oz --threads {threads} --ploidy {params.pl} -o {output.uc86pmac}")
#        shell("freebayes -f {input.pmac} -C {params.c} {input.uc86pmac} -p {params.p} --gvcf  | vcffilter -f \"QUAL > {params.q} \" > {output.uc86pmac}")

rule normalize_indels:
    input:
        pjref=rules.snp_call.input.pj,
        pmacref=rules.snp_call.input.pmac,
        pj46pj=rules.snp_call.output.pj46pj,
        pj54cpj=rules.snp_call.output.pj54cpj,
        pj55pj=rules.snp_call.output.pj55pj,
        pjzpj=rules.snp_call.output.pjzpj,
        pjwpj=rules.snp_call.output.pjwpj,
        cj36pmac=rules.snp_call.output.cj36pmac,
        er17pmac=rules.snp_call.output.er17pmac,
        gl92pmac=rules.snp_call.output.gl92pmac,
        uc86pmac=rules.snp_call.output.uc86pmac,
    output:
        pj46pjt=temp("../../data/processed/last/vcf/vsPj/Pj46.norm.bcf"),pj46pj="../../data/processed/last/vcf/vsPj/Pj46.norm.flt-indels.bcf.gz",
        pj54cpjt=temp("../../data/processed/last/vcf/vsPj/Pj54c.norm.bcf"),pj54cpj="../../data/processed/last/vcf/vsPj/Pj54c.norm.flt-indels.bcf.gz",
        pj55pjt=temp("../../data/processed/last/vcf/vsPj/Pj55.norm.bcf"),pj55pj="../../data/processed/last/vcf/vsPj/Pj55.norm.flt-indels.bcf.gz",
        pjzpjt=temp("../../data/processed/last/vcf/vsPj/PjZ.norm.bcf"),pjzpj="../../data/processed/last/vcf/vsPj/PjZ.norm.flt-indels.bcf.gz",
        pjwpjt=temp("../../data/processed/last/vcf/vsPj/PjW.norm.bcf"),pjwpj="../../data/processed/last/vcf/vsPj/PjW.norm.flt-indels.bcf.gz",
        cj36pmact=temp("../../data/processed/last/vcf/vsPmac/PmacCj36.norm.bcf"),cj36pmac="../../data/processed/last/vcf/vsPmac/PmacCj36.norm.flt-indels.bcf.gz",
        er17pmact=temp("../../data/processed/last/vcf/vsPmac/PmacER17.norm.bcf"),er17pmac="../../data/processed/last/vcf/vsPmac/PmacER17.norm.flt-indels.bcf.gz",
        gl92pmact=temp("../../data/processed/last/vcf/vsPmac/Pmacgl92.norm.bcf"),gl92pmac="../../data/processed/last/vcf/vsPmac/Pmacgl92.norm.flt-indels.bcf.gz",
        uc86pmact=temp("../../data/processed/last/vcf/vsPmac/Pmacuc86.norm.bcf"),uc86pmac="../../data/processed/last/vcf/vsPmac/Pmacuc86.norm.flt-indels.bcf.gz"

    threads: 8

    conda:
        "../envs/bcftools.yaml"

    params:
            i="5", b="+any"
    run:
        shell("bcftools index -f {input.pj46pj}")
        shell("bcftools norm -f {input.pjref} {input.pj46pj} -m {params.b} -Ob -o {output.pj46pjt}")
        shell("bcftools filter --IndelGap {params.i} -Oz {output.pj46pjt} -o {output.pj46pj}")
        shell("bcftools index -f {output.pj46pj}")
 
        shell("bcftools index -f {input.pj54cpj}")
        shell("bcftools norm -f {input.pjref} {input.pj54cpj} -m {params.b} -Ob -o {output.pj54cpjt}")
        shell("bcftools filter --IndelGap {params.i} -Oz {output.pj54cpjt} -Oz -o {output.pj54cpj}")
        shell("bcftools index -f {output.pj54cpj}")

        shell("bcftools index -f {input.pj55pj}")
        shell("bcftools norm -f {input.pjref} {input.pj55pj} -m {params.b} -Ob -o {output.pj55pjt}")
        shell("bcftools filter --IndelGap {params.i} -Oz {output.pj55pjt} -o {output.pj55pj}")
        shell("bcftools index -f {output.pj55pj}")

        shell("bcftools index -f {input.pjzpj}")
        shell("bcftools norm -f {input.pjref} {input.pjzpj} -m {params.b} -Ob -o {output.pjzpjt}")
        shell("bcftools filter --IndelGap {params.i} -Oz {output.pjzpjt} -o {output.pjzpj}")
        shell("bcftools index -f {output.pjzpj}")

        shell("bcftools index -f {input.pjwpj}")
        shell("bcftools norm -f {input.pjref} {input.pjwpj} -m {params.b} -Ob -o {output.pjwpjt}")
        shell("bcftools filter --IndelGap {params.i} -Oz {output.pjwpjt} -o {output.pjwpj}")
        shell("bcftools index -f {output.pjwpj}")

        shell("bcftools index -f {input.cj36pmac}")
        shell("bcftools norm -f {input.pmacref} {input.cj36pmac} -m {params.b} -Ob -o {output.cj36pmact}")
        shell("bcftools filter --IndelGap {params.i} -Oz {output.cj36pmact} -o {output.cj36pmac}")
        shell("bcftools index -f {output.cj36pmac}")

        shell("bcftools index -f {input.er17pmac}")
        shell("bcftools norm -f {input.pmacref} {input.er17pmac} -m {params.b} -Ob -o {output.er17pmact}")
        shell("bcftools filter --IndelGap {params.i} -Oz {output.er17pmact} -o {output.er17pmac}")
        shell("bcftools index -f {output.er17pmac}")

        shell("bcftools index -f {input.gl92pmac}")
        shell("bcftools norm -f {input.pmacref} {input.gl92pmac} -m {params.b} -Ob -o {output.gl92pmact}")
        shell("bcftools filter --IndelGap {params.i} -Oz {output.gl92pmact} -o {output.gl92pmac}")
        shell("bcftools index -f {output.gl92pmac}")

        shell("bcftools index -f {input.uc86pmac}")
        shell("bcftools norm -f {input.pmacref} {input.uc86pmac} -m {params.b} -Ob -o {output.uc86pmact}")
        shell("bcftools filter --IndelGap {params.i} -Oz {output.uc86pmact} -o {output.uc86pmac}")
        shell("bcftools index -f {output.uc86pmac}")

rule consensus:
    input:
        pjref=rules.snp_call.input.pj,
        pmacref=rules.snp_call.input.pmac,
        pj46pj=rules.normalize_indels.output.pj46pj,
        pj54cpj=rules.normalize_indels.output.pj54cpj,
        pj55pj=rules.normalize_indels.output.pj55pj,
        pjzpj=rules.normalize_indels.output.pjzpj,
        pjwpj=rules.normalize_indels.output.pjwpj,
        cj36pmac=rules.normalize_indels.output.cj36pmac,
        er17pmac=rules.normalize_indels.output.er17pmac,
        gl92pmac=rules.normalize_indels.output.gl92pmac,
        uc86pmac=rules.normalize_indels.output.uc86pmac,
    output:
        pj46pj="../../data/processed/last/vcf/vsPj/Pj46.cns.fa",
        pj54cpj="../../data/processed/last/vcf/vsPj/Pj54c.cns.fa",
        pj55pj="../../data/processed/last/vcf/vsPj/Pj55.cns.fa",
        pjzpj="../../data/processed/last/vcf/vsPj/PjZ.cns.fa",
        pjwpj="../../data/processed/last/vcf/vsPj/PjW.cns.fa",
        cj36pmac="../../data/processed/last/vcf/vsPmac/PmacCj36.cns.fa",
        er17pmac="../../data/processed/last/vcf/vsPmac/PmacER17.cns.fa",
        gl92pmac="../../data/processed/last/vcf/vsPmac/Pmacgl92.cns.fa",
        uc86pmac="../../data/processed/last/vcf/vsPmac/Pmacuc86.cns.fa"

    params: h="A"
    run:
        shell("cat {input.pjref} | bcftools consensus -H {params.h} {input.pj46pj} > {output.pj46pj}")
        shell("cat {input.pjref} | bcftools consensus -H {params.h} {input.pj54cpj} > {output.pj54cpj}")
        shell("cat {input.pjref} | bcftools consensus -H {params.h} {input.pj55pj} > {output.pj55pj}")
        shell("cat {input.pjref} | bcftools consensus -H {params.h} {input.pjzpj} > {output.pjzpj}")
        shell("cat {input.pjref} | bcftools consensus -H {params.h} {input.pjwpj} > {output.pjwpj}")
        shell("cat {input.pmacref} | bcftools consensus -H {params.h} {input.cj36pmac} > {output.cj36pmac}")
        shell("cat {input.pmacref} | bcftools consensus -H {params.h} {input.er17pmac} > {output.er17pmac}")
        shell("cat {input.pmacref} | bcftools consensus -H {params.h} {input.gl92pmac} > {output.gl92pmac}")
        shell("cat {input.pmacref} | bcftools consensus -H {params.h} {input.uc86pmac} > {output.uc86pmac}")

# finally succeed!
rule minimap:
    input:
        pjru7="../../data/raw/Pjir/GCF_001477535.1_Pneu_jiro_RU7_V2_genomic.fna",
        pjse8="../../data/processed/Pjir/SE8/Mk2/Pj_SE8_assembly.maker.output/Pj_SE8_assembly.fasta",
        pjse2178="../../data/processed/Pjir/SE2178/Mk2/Pneumocystis_jirovecii_StrE2178.fasta",
        pj46pj=rules.consensus.output.pj46pj,
        pj54cpj=rules.consensus.output.pj54cpj,
        pj55pj=rules.consensus.output.pj55pj,
        pjzpj=rules.consensus.output.pjzpj,
        pjwpj=rules.consensus.output.pjwpj,
        pmac="../../data/processed/Pmac/Funannotate/PMAC_NAID_v2/predict_results/Pneumocystis_macacae_NIAID.scaffolds.fa",
        pmaccj36=rules.consensus.output.cj36pmac,
        pmacer17=rules.consensus.output.er17pmac,
        pmacgl92=rules.consensus.output.gl92pmac,
        pmacuc86=rules.consensus.output.uc86pmac, 
        po="../../data/processed/Poryc/Funannotate/PORYCT_v0110_Ns_removed/predict_results/Pneumocystis_oryctolagi_MERGE.scaffolds.fa",
        ck1="../../data/processed/Pcanis/Funaannotate/Pneumocystis_canis_CK.scaffolds.fa",
        ck2="../../data/processed/Pcanis/Ck2/NCBI_sub/Pdg_scaffolds.cl6.alterv2_1kremoved_cl1_trim.fasta",
        pcanA="../../data/processed/Pcanis/dogA_Austria/Ass_V1/FUNANNOTATE/PCANB_v1/predict_results/Pneumocystis_canis_AUS.scaffolds.fa",
        pc="../../data/raw/Pcar/GCF_001477545.1_Pneu_cari_B80_V3_genomic.fna",
        pm="../../data/raw/Pmur/GCF_000349005.2_Pneumo_murina_B123_V4_genomic.fna",
        pwk="../../data/processed/Pwk/Funannotate/Pneumocystis_wakefieldiae_MERGE.scaffolds.fa"

    output:
        s="../../data/processed/minimap/all_files.txt",
        tmp=temp("../../data/processed/minimap/tmp"),
        r="../../data/processed/minimap/scores.txt"
    run:
        shell("echo \"{input.pjru7},Pjir\" > {output.s}")
        shell("echo \"{input.pjse8},Pjir\" >> {output.s}")
        shell("echo \"{input.pjse2178},Pjir\" >> {output.s}")
        shell("echo \"{input.pj46pj},Pjir\" >> {output.s}")
        shell("echo \"{input.pj54cpj},Pjir\" >> {output.s}")
        shell("echo \"{input.pj55pj},Pjir\" >> {output.s}")
        shell("echo \"{input.pjzpj},Pjir\" >> {output.s}")
        shell("echo \"{input.pjwpj},Pjir\" >> {output.s}")
        shell("echo \"{input.pmac},Pmac\" >> {output.s}")
        shell("echo \"{input.pmaccj36},Pmac\" >> {output.s}")
        shell("echo \"{input.pmacer17},Pmac\" >> {output.s}")
        shell("echo \"{input.pmacgl92},Pmac\" >> {output.s}")
        shell("echo \"{input.pmacuc86},Pmac\" >> {output.s}")
        shell("echo \"{input.po},Pory\" >> {output.s}")
        shell("echo \"{input.ck1},PcanCk1\" >> {output.s}")
        shell("echo \"{input.ck2},PcanCk2\" >> {output.s}")
        shell("echo \"{input.pcanA},PcanA\" >> {output.s}")
        shell("echo \"{input.pc},Pcar\" >> {output.s}")
        shell("echo \"{input.pm},Pmur\" >> {output.s}")
        shell("echo \"{input.pwk},Pwk\" >> {output.s}")
        shell("echo > {output.tmp}")
        shell("perl ../scripts/pairwise_score.pl {output.s} {output.tmp} > {output.r}") 
