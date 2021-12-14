for i in EvC.norm.stat.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx KEGG_PPAR_SIGNALING_PATHWAY.gmt \
-set_min 5 -set_max 1500 \
-rpt_label evc_KEGG_PPAR -plot_top_x 40 -nperm 1000
done

for i in hpm3.norm.res.A.C.stat.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx HALLMARK_TNFA_SIGNALING_VIA_NFKB.gmt \
-set_min 5 -set_max 1500 \
-rpt_label avc_TNF -plot_top_x 40 -nperm 1000    
done

for i in norm.G1.Rev.stat.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx HALLMARK_E2F_TARGETS.gmt \
-set_min 5 -set_max 1500 \
-rpt_label g1vrev_E2F -plot_top_x 40 -nperm 1000          
done

for i in muscle.stat.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx GO_SKELETAL_MUSCLE_CONTRACTION.gmt \
-set_min 5 -set_max 1500 \
-rpt_label svc_SKEMUSC -plot_top_x 40 -nperm 1000        
done

