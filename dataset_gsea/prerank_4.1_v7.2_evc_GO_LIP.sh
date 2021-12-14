for i in EvC*rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx GO_LIPID_TRANSFER_ACTIVITY.gmt \
-set_min 5 -set_max 1500 \
-scoring_scheme classic \
-rpt_label $i.classic.GO_LIP_TRANS -plot_top_x 40 -nperm 1000
done

for i in EvC*rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx GO_LIPID_TRANSFER_ACTIVITY.gmt \
-set_min 5 -set_max 1500 \
-rpt_label $i.weighted.GO_LIP_TRANS -plot_top_x 40 -nperm 1000
done
