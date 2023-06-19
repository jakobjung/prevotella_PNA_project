#!/bin/bash
main(){
    # go to project directory where the reads and reference
    # sequences are stored:
    PROJECT=../data
    echo "Start trimming"
    rename_trim_rna_libs
    echo "Trimming done. Start mapping"
    align_rna_reads_genome
    echo "Finished mapping. Start connecting all tab files"
    featureCounts -T 5 -t CDS -g locus_tag \
		  -a $PROJECT/reference_sequences/tot_seqs.gff \
		  -o $PROJECT/rna_align/counttable_bb.txt -M \
		  $PROJECT/rna_align/*.bam
   
    featureCounts -T 5 -t rRNA -g locus_tag \
		  -a $PROJECT/reference_sequences/tot_seqs.gff \
		  -o $PROJECT/rna_align/rRNA_counttable_bb.txt \
		  $PROJECT/rna_align/*.bam
    
    featureCounts -T 5 -t tRNA -g locus_tag \
		  -a $PROJECT/reference_sequences/tot_seqs.gff \
		  -o $PROJECT/rna_align/tRNA_counttable_bb.txt \
		  $PROJECT/rna_align/*.bam
}

rename_trim_rna_libs(){
    mkdir -pv $PROJECT/libs
    for NAME in $(ls $PROJECT/fastq/*.fq.gz)
    do
        echo "$NAME starts trimming nowwwwwww"
        NEWNAME=${NAME##*/}
        NEWNAME=${NEWNAME%.fq.gz}_trimmed.fastq.gz
        echo $NEWNAME
       # bbduk trims low quality bases and removes adapters:
        #bbduk.sh  in=$NAME \
	#		     ref=../data/reference_sequences/adapters.fa -Xmx4g t=30\
	#		     out=$PROJECT/libs/${NEWNAME} ktrim=r k=20 mink=4\
	#		     hdist=1 qtrim=r trimq=10 ftl=12 overwrite=t
	cutadapt --trim-n --match-read-wildcards -u 16 -n 4 -a AGATCGGAAGAGCACACGTCTG -a AAAAAAAA \
		 -a GAACTCCAGTCAC -e 0.2 --nextseq-trim 20 -m 15 \
		 -o $PROJECT/libs/${NEWNAME} $NAME 
	fastqc $PROJECT/libs/${NEWNAME}
    done
}


align_rna_reads_genome(){
    mkdir -p $PROJECT/rna_align
    DIR=$PROJECT/rna_align
    for i in $(ls $PROJECT/libs/*.fastq.gz)
    do
        NAME=${i##*/}
        NAME=${NAME%_trimmed.fastq.gz}
        echo "Starting mapping for sample: $NAME"
        #salmon quant -i $PROJECT/reference_sequences/metatranscriptome_index -l A -r $i \
	#       --fldMean 300 --fldSD 300 --validateMappings --writeUnmappedNames --incompatPrior 0.0 --recoverOrphans --softclipOverhangs -o $DIR/$NAME
	#fastqc $i
        bbmap.sh in=$i trimreaddescription=t  t=20 \
			     ref=$PROJECT/reference_sequences/tot_seqs.fasta \
			     k=12 ambig=random outm=$DIR/$NAME.sam #outu=$DIR/${NAME}_unmapped.sam
	samtools sort -O BAM -@ 40 $DIR/$NAME.sam > $DIR/$NAME.bam
	echo  "done sorting, start indexing"
	samtools index $DIR/$NAME.bam
	echo "done indexing, start bamcoverage"
	bamCoverage -b $DIR/$NAME.bam -o $DIR/$NAME.bigwig -of bigwig -bs 1
	#rm $i
    done
}

main  
