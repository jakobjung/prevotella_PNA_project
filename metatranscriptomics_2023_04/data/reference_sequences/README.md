```bash
sed -E "s/.*locus_tag=([^]]+).*/>\\1/" b_thera_cds.fna > b_theta_cds.fasta


sed -E "s/.*locus_tag=([^]]+).*/>\\1/" b_theta_only_cds_original.fasta >  b_theta_cds_orig.fasta


proteinortho6.pl --project=portho_b_theta ./b_theta_cds_orig.fasta ./b_theta_cds.fasta --p=blastn



```

