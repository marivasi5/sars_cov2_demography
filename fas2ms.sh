pop=$1
fas2ms.pl fas=${pop}_gisaid_cov2020_sequences_hc_20200402_oneline_fullseq_ns_20200402.fasta.aln outgroup="BAT|PANGOLIN" > ${pop}_gisaid_cov2020_sequences_hc_20200402_oneline_fullseq_ns_20200402.fasta.ms 2>${pop}.info

seqs=`grep [01] ${pop}_gisaid_cov2020_sequences_hc_20200402_oneline_fullseq_ns_20200402.fasta.ms  | grep -v [^01] | wc -l | awk '{print $1}'`
echo "For population $pop there are $seqs sequences"

