# data

```bash
grep -P "(\tgene\t)|(\tpseudogene\t)|(\rRNA\t)|(\tRNA\t)" ./reference_sequences/p_copri_whole.gff3 |\
cut -f9 | \
sed -r 's/^.*gene=([^;]+).*;locus_tag=([^; ]+).*$/\1\t\2/' | \
sed -r 's/^ID=.*;locus_tag=([^;]+).*$/\1\t\1/'  > link_lt_gn.tab
```

