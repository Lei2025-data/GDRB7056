#!/bin/bash
set -e
/Bio/User/yinquan/software/miniconda3/bin/Rscript /Bio/Bin/pipeline/SCENIC/v1.0/libs/RunSCENIC.r /public2/Bio/Project/GDRB7056/GDRB7056_sup_21/4.scenic/out/4.AUC/RunSCENIC.yaml

rm -f /public2/Bio/Project/GDRB7056/GDRB7056_sup_21/4.scenic/out/4.AUC/Rplots.pdf; for i in /public2/Bio/Project/GDRB7056/GDRB7056_sup_21/4.scenic/out/4.AUC/*pdf; do convert -density 300 $i ${i/pdf/png};done
cd /public2/Bio/Project/GDRB7056/GDRB7056_sup_21/4.scenic/out/4.AUC/Step3_RegulonActivity_tSNE_colByActivity; rename ' ' '' * ; for i in *pdf;do convert -density 300 $i ${i/pdf/png};done;cd -
/Bio/bin/perl /Bio/Bin/pipeline/SCENIC/v1.0/libs/cytoscape.pl /public2/Bio/Project/GDRB7056/GDRB7056_sup_21/4.scenic/out/4.AUC/tf_target.egde.tsv /public2/Bio/Project/GDRB7056/GDRB7056_sup_21/4.scenic/out/4.AUC/tf_target.node.tsv > /public2/Bio/Project/GDRB7056/GDRB7056_sup_21/4.scenic/out/4.AUC/code.js