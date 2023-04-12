```bash
sed -E "s/>.+\\[gene=([^]]+)(.*)\\[locus_tag=([^]]+).*/>b_ovatus|\\1|\\3/" b_ovatus.fasta | sed -E "s/>.+\\[locus_tag=([^]]+).*/>b_ovatus|\\1|\\1/" > database_transcriptomics_flame.fasta

sed -E "s/>.+\\[gene=([^]]+)(.*)\\[locus_tag=([^]]+).*/>b_theta|\\1|\\3/" b_theta.fasta | sed -E "s/>.+\\[locus_tag=([^]]+).*/>b_theta|\\1|\\1/" >> database_transcriptomics_flame.fasta

sed -E "s/>.+\\[gene=([^]]+)(.*)\\[locus_tag=([^]]+).*/>b_theta|\\1|\\3/" p_copri.fasta | sed -E "s/>.+\\[locus_tag=([^]]+).*/>p_copri|\\1|\\1/" >> database_transcriptomics_flame.fasta 



```

