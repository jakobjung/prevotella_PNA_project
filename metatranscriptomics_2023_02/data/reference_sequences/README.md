```bash
sed -E "s/>.+\\[gene=([^]]+)(.*)\\[locus_tag=([^]]+).*/>b_ovatus|\\1|\\3/" b_ovatus.fasta | sed -E "s/>.+\\[locus_tag=([^]]+).*/>b_ovatus|\\1|\\1/" > database_transcriptomics_flame.fasta

sed -E "s/>.+\\[gene=([^]]+)(.*)\\[locus_tag=([^]]+).*/>b_theta|\\1|\\3/" b_theta.fasta | sed -E "s/>.+\\[locus_tag=([^]]+).*/>b_theta|\\1|\\1/" >> database_transcriptomics_flame.fasta

sed -E "s/>.+\\[gene=([^]]+)(.*)\\[locus_tag=([^]]+).*/>b_theta|\\1|\\3/" p_copri.fasta | sed -E "s/>.+\\[locus_tag=([^]]+).*/>p_copri|\\1|\\1/" >> database_transcriptomics_flame.fasta 



```



```bash
grep -P "exon\t" b_theta_wg.gff3 > b_theta_noncoding.gff
gffread b_theta_noncoding.gff -g b_theta_whole.fasta -w b_theta_nc.fasta -F

grep -P "exon\t" b_ovatus_whole.gff3 > b_ovatus_noncoding.gff
gffread b_ovatus_noncoding.gff -g b_ovatus_whole.fasta -w b_ovatus_nc.fasta -F

grep -P "exon\t" p_copri_whole.gff3 > p_copri_noncoding.gff
gffread p_copri_noncoding.gff -g p_copri_whole.fasta -w p_copri_nc.fasta -F

```



```bash
sed -E "s/^>.+gbkey=([^ ]+).+locus_tag=([^ ]+).+/>p_copri|\\1|\\2/" p_copri_nc.fasta > ncRNAs.fasta 
sed -E "s/^>.+gbkey=([^ ]+).+locus_tag=([^ ]+).+/>b_theta|\\1|\\2/" b_theta_nc.fasta >> ncRNAs.fasta 
sed -E "s/^>.+gbkey=([^ ]+).+locus_tag=([^ ]+).+/>b_ovatud|\\1|\\2/" b_ovatus_nc.fasta >> ncRNAs.fasta 

```

