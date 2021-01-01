for i in *stat.rnk *ape.logFC.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx KEGG_PPAR_SIGNALING_PATHWAY.gmt \
-set_min 5 -set_max 1500 \
-rpt_label $i.KEGG_PPAR -plot_top_x 40 -nperm 1000
done
