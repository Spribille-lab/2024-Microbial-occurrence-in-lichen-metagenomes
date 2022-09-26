#!/usr/bin/env bash

while getopts a:c: flag
do
    case "${flag}" in
        a) AAI=${OPTARG};;
        c) ALT=${OPTARG};;
    esac
done
echo "AAI: $AAI";
echo "ALT: $ALT";

all(){
    set_variables
    create_mmseqDB
    run_linclust
    extract_reps
    #add_seqinfo
}


set_variables(){
    DB=mgy_seqs
    INPUT_FASTA=/hps/research/finn/escameron/databases/mgnify_proteins/emg_peptidedb_protein_lichenappend.fasta
    TEMPDIR=./temp
    THREADS=16
    RUN=aai-${AAI}-alt-${ALT}
}

create_mmseqDB(){
	printf "\n=====Create mmseq2 database=====\n"
	mmseqs createdb ${INPUT_FASTA} ${DB}.mmseqs
}

run_linclust(){
	printf "\n=====Cluster Seqences with linclust=====\n"
	mkdir -p ${RUN}
	mmseqs linclust ${DB}.mmseqs ${RUN}/${DB}-${RUN}.cluster ${TEMPDIR} --min-seq-id ${AAI} -c ${ALT} --cov-mode 1 --threads ${THREADS}
}

extract_reps(){
	printf "\n=====Extract representative sequences=====\n"
	mmseqs createsubdb ${RUN}/${DB}-${RUN}.cluster ${DB}.mmseqs ${RUN}/${DB}-${RUN}.cluster_rep
	mmseqs convert2fasta ${RUN}/${DB}-${RUN}.cluster_rep ${RUN}/${DB}-${RUN}.cluster_rep.fasta
	mmseqs createtsv ${DB}.mmseqs ${DB}.mmseqs ${RUN}/${DB}-${RUN}.cluster ${RUN}/${DB}-${RUN}.cluster.tsv
}

add_seqinfo(){
	printf "\n=====Add sequence info=====\n"
	mmseqs createseqfiledb ${DB}.mmseqs ${RUN}/${DB}-${RUN}.cluster ${RUN}/${DB}-${RUN}.cluster_seq
	mmseqs result2flat ${DB}.mmseqs ${DB}.mmseqs ${RUN}/${DB}-${RUN}.cluster_seq ${DB}.cluster_seq.fasta
}

all
