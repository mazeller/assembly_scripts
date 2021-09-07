#!/bin/bash
set -o errexit

#Required programs: fastqc, BBDuk, BBNorm, Kaiju, seqtq, SPADES
#TODO - modify thread counts for other machines
echo "
   __    ___  ___  ____  __  __  ____  __   _  _    ____  ____  ____  ____  __    ____  _  _  ____ 
  /__\  / __)/ __)( ___)(  \/  )(  _ \(  ) ( \/ )  (  _ \(_  _)(  _ \( ___)(  )  (_  _)( \( )( ___)
 /(__)\ \__ \\\__ \ )__)  )    (  ) _ < )(__ \  /    )___/ _)(_  )___/ )__)  )(__  _)(_  )  (  )__) 
(__)(__)(___/(___/(____)(_/\/\_)(____/(____)(__)   (__)  (____)(__)  (____)(____)(____)(_)\_)(____)
MZ 25 Aug 2021 v0.01
This script is very conservative.
"

#Help text
help() { 
	echo "Usage: $0 -a <read1> -b <read2>
	
The following are run:
fastqc - Initial quality check.
seqtk - Filtering classified sequences
BBDuk - Trim adapters and low quality reads.
fastqc - Trim quality check.
BBNorm - Normalize the read representation. Do not use for quantitative studies.
fastqc - Normalize quality check
Kaiju - Taxonomic binning of reads to isolate reads of interest.
Kraken - Second taxonomic binning step

*** pause

SPADES assembly on select set

FLAGS
-n	No quality controls. Skips FastQC steps
-q	Skip normalizationstep. Allows analyssi to remain quantitative.
	" 1>&2; 
	exit 1; 
}

#GLOBALS
SCRIPTDIR=`dirname "$0"`
FASTQCDIR="fastqc_results"
BBDUKADAPTER="~/bbmap/resources/adapters.fa"
KAIJUDBDIR="/mnt/c/Users/mazeller.NUSSTF/Desktop/Lab/kaijudb/viruses"
QCSKIP=0
NORMSKIP=0

#COLORS
GREEN='\033[1;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

#get opts
while getopts ":hnqa:b:" o; do
	case "${o}" in
		h)
			help
			;;
		a)
			input1=${OPTARG}
			FILE1=`basename "$input1"`
			;;
		b)
			input2=${OPTARG}
			FILE2=`basename "$input2"`
			;;
		n)
			QCSKIP=1
			echo -e "${RED}Skipping all FastQC steps!${NC}"
			;;
		q)
			NORMSKIP=1
			echo -e "${RED}Skipping read normalization!${NC}"
			;;
	esac
done
#shift $((OPTIND-1))

#if [ -z "${a}" ] || [ -z "${b}" ]; then
#	help
#fi

#Initial quality controls
if [ ${QCSKIP} != 1 ]; then
	echo -e "${GREEN}Initial quality checks${NC}"
	mkdir -p ${FASTQCDIR}
	mkdir -p results
	fastqc -t 8 -f fastq  ${input1} -noextract -o ${FASTQCDIR} 
	fastqc -t 8 -f fastq  ${input2} -noextract -o ${FASTQCDIR} 
fi

#Trim adapters and QC
echo -e "${GREEN}Trimming adapters and quality${NC}"
~/bbmap/bbduk.sh ktrimright=t k=27 hdist=1 edist=0 ref=${BBDUKADAPTER} qtrim=rl trimq=20 minlength=20 trimbyoverlap=t minoverlap=20 ordered=t qin=33 in=${input1} in2=${input2} out=results/${FILE1}.trimmed.fastq out2=results/${FILE2}.trimmed.fastq
STEP="trimmed"

#Trim QC
if [ ${QCSKIP} != 1 ]; then
	fastqc -t 8 -f fastq  results/${FILE1}.trimmed.fastq -noextract -o ${FASTQCDIR}
	fastqc -t 8 -f fastq  results/${FILE2}.trimmed.fastq -noextract -o ${FASTQCDIR}
fi

#Normalize
if [ ${NORMSKIP} != 1 ]; then
	echo -e "${GREEN}Normalizing read data${NC}"
	~/bbmap/bbnorm.sh target=40 mindepth=6 threads=12 qin=auto in=results/${FILE1}.trimmed.fastq in2=results/${FILE2}.trimmed.fastq out=results/${FILE1}.norm.fastq out2=results/${FILE2}.norm.fastq

	#Normalize QC
	if [ ${QCSKIP} != 1 ]; then
		echo -e "${GREEN}Normalized quality checks${NC}"
		fastqc -t 8 -f fastq  results/${FILE1}.norm.fastq -noextract -o ${FASTQCDIR}
		fastqc -t 8 -f fastq  results/${FILE2}.norm.fastq -noextract -o ${FASTQCDIR}
	fi
	STEP="norm"
fi


#Kaiju Taxonomic binning
echo -e "${GREEN}Kaiju taxonomic binning${NC}"
kaiju -z 20 -t ${KAIJUDBDIR}/nodes.dmp -f ${KAIJUDBDIR}/kaiju_db_viruses.fmi -i results/${FILE1}.${STEP}.fastq -j results/${FILE2}.${STEP}.fastq -o results/kaiju.out -v
#kaiju2table -t ${KAIJUDBDIR}/nodes.dmp -n ${KAIJUDBDIR}/names.dmp -r genus -o results/kaiju_summary.tsv results/kaiju.out -v -p #-l superkingdom,genus,species
kaiju-addTaxonNames -t ${KAIJUDBDIR}/nodes.dmp -n ${KAIJUDBDIR}/names.dmp -i results/kaiju.out -o results/kaiju.names.out -u -p -v

#Extract classfied viral seqs
echo -e "${GREEN}Splitting out viral sequences${NC}"
awk '$8=="Viruses;"{print $2}' results/kaiju.names.out > results/viralset.txt
seqtk subseq results/${FILE1}.${STEP}.fastq results/viralset.txt >  results/viral_R1.fastq
seqtk subseq results/${FILE2}.${STEP}.fastq results/viralset.txt >  results/viral_R2.fastq

###ASSEMBLY CAN BEGIN###
echo -e "${GREEN}Assembling via spades${NC}"
spades.py -k 21,33,55,77,99 -1 results/${FILE1}.${STEP}.fastq -2 results/${FILE2}.${STEP}.fastq -o spades_output

#Re-Kaiju annotate
echo -e "${GREEN}Kaiju annotations of SPADES output${NC}"
kaiju -z 20 -t ${KAIJUDBDIR}/nodes.dmp -f ${KAIJUDBDIR}/kaiju_db_viruses.fmi -i spades_output/contigs.fasta -o spades_output/kaiju.out -v
kaiju-addTaxonNames -t ${KAIJUDBDIR}/nodes.dmp -n ${KAIJUDBDIR}/names.dmp -i spades_output/kaiju.out -o spades_output/kaiju.names.out -u -v -r superkingdom,genus,species
${SCRIPTDIR}/fastannot.py -f spades_output/contigs.fasta -a spades_output/kaiju.names.out > spades_output/annotation.fasta
