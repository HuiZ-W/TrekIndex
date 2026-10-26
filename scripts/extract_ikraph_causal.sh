#!/bin/bash
# iKraph (Zenodo 14851275, CC BY 4.0) -> directed causal relations TSV for tools/causenet_dag.cpp.
# Keeps PubMed-derived relations with a direction (12: node1 -> node2, 21: node2 -> node1)
# of types 2 Positive_Correlation, 3 Negative_Correlation, 51 Target_Disease, 52 Biomarker_Disease.
# support = number of supporting PubMed documents, summed over relation types.
#   usage: extract_ikraph_causal.sh iKraph_full/PubMedList.json out.tsv
LC_ALL=C awk '
/"id": "/ { if(id!="") emit(); match($0,/"id": "[^"]*"/); id=substr($0,RSTART+7,RLENGTH-8); docs=0; next }
/^      \[$/ { docs++ }
function emit(   a,r,d,u,v) { split(id,a,"."); r=a[3]; d=a[5];
  if((r=="2"||r=="3"||r=="51"||r=="52") && (d=="12"||d=="21")) { if(d=="12"){u=a[1];v=a[2]}else{u=a[2];v=a[1]}; sup[u"\t"v]+=docs; typ[u"\t"v]=typ[u"\t"v] "" r ";" } }
END { if(id!="") emit(); for(k in sup) print k "\t" sup[k] "\t" typ[k] }' "$1" > "$2"
