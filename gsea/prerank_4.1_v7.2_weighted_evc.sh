for i in *stat.rnk *ape.logFC.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx /net/ostrom/data/bcc/charliew/annotationFiles/h.all.v7.2.symbols.gmt \
-set_min 5 -set_max 1500 \
-rpt_label $i.h -plot_top_x 40 -nperm 1000
done

for i in *stat.rnk *ape.logFC.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx /net/ostrom/data/bcc/charliew/annotationFiles/c8.all.v7.2.symbols.gmt \
-set_min 5 -set_max 1500 \
-rpt_label $i.c8 -plot_top_x 40 -nperm 1000
done

for i in *stat.rnk *ape.logFC.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx /net/ostrom/data/bcc/charliew/annotationFiles/c2.cgp.v7.2.symbols.gmt \
-set_min 5 -set_max 1500 \
-rpt_label $i.c2cgp -plot_top_x 40 -nperm 1000
done

for i in *stat.rnk *ape.logFC.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx /net/ostrom/data/bcc/charliew/annotationFiles/c2.cp.v7.2.symbols.gmt \
-set_min 5 -set_max 1500 \
-rpt_label $i.c2cp -plot_top_x 40 -nperm 1000
done

for i in *stat.rnk *ape.logFC.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx /net/ostrom/data/bcc/charliew/annotationFiles/c5.go.cc.v7.2.symbols.gmt \
-set_min 5 -set_max 1500 \
-rpt_label $i.c5cc -plot_top_x 40 -nperm 1000
done


for i in *stat.rnk *ape.logFC.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx /net/ostrom/data/bcc/charliew/annotationFiles/c5.go.mf.v7.2.symbols.gmt \
-set_min 5 -set_max 1500 \
-rpt_label $i.c5mf -plot_top_x 40 -nperm 1000
done

for i in *stat.rnk *ape.logFC.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx /net/ostrom/data/bcc/charliew/annotationFiles/c5.go.bp.v7.2.symbols.gmt \
-set_min 5 -set_max 1500 \
-rpt_label $i.c5bp -plot_top_x 40 -nperm 1000
done

for i in *stat.rnk *ape.logFC.rnk
do
sh /net/ostrom/data/bcc/charliew/gsea_4/GSEA_Linux_4.1.0/gsea-cli.sh GSEAPreranked -rnk $i \
-gmx /net/ostrom/data/bcc/charliew/annotationFiles/c5.hpo.v7.2.symbols.gmt \
-set_min 5 -set_max 1500 \
-rpt_label $i.c5hpo -plot_top_x 40 -nperm 1000
done
