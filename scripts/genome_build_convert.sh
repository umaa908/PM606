### convert lifeline data from hg19 to hg38
# convert to BED format
for HMO in $(cat lifeline_hmo_list.txt); do
    zcat HMO_GWAS_res_M3.mother_milk_HMO_${HMO}_ugml_invr.glm.linear.gz | awk '
    BEGIN { FS=OFS="\t" }
    NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="#CHROM" || $i=="CHROM") chr_col=i
        else if ($i=="POS" || $i=="BP") pos_col=i
        else if ($i=="ID") id_col=i
    }
    next
    }
    {
    chr = $chr_col
    if (chr !~ /^chr/) chr = "chr" chr
    print chr, $pos_col-1, $pos_col, $id_col
    }' > $HMO.hg19.bed
done

# used liftOver locally since I couldn't find it on CARC
# 6620 vars not converted, 5378355 vars successful
for f in *.hg19.bed; do
    base=${f%.hg19.bed}

    liftOver \
        "$f" \
        hg19ToHg38.over.chain.gz \
        "${base}.hg38.bed" \
        "${base}.unmapped.bed"
done

# back on CARC
# Make an ID -> hg38 coordinate map from the lifted BED
# Rebuild the original glm file with hg38 CHR/POS
for HMO in $(cat data/lifeline/lifeline_hmo_list.txt); do

  awk 'BEGIN{OFS="\t"} {chr=$1; sub(/^chr/,"",chr); print $4, chr, $3}' \
    "data/lifeline/${HMO}.hg38.bed" > "data/lifeline/${HMO}.hg38.map.tsv"

  awk '
    BEGIN { FS=OFS="\t" }
    NR==FNR {
      map[$1] = $2 OFS $3
      next
    }
    FNR==1 {
      for (i=1; i<=NF; i++) {
        if ($i=="#CHROM" || $i=="CHROM" || $i=="CHR") chr_col=i
        else if ($i=="POS" || $i=="BP") pos_col=i
        else if ($i=="ID" || $i=="SNP") id_col=i
      }
      print
      next
    }
    {
      if ($id_col in map) {
        split(map[$id_col], a, "\t")
        $chr_col = a[1]
        $pos_col = a[2]
        print
      }
    }
  ' "data/lifeline/${HMO}.hg38.map.tsv" <(zcat "data/lifeline/HMO_GWAS_res_M3.mother_milk_HMO_${HMO}_ugml_invr.glm.linear.gz") \
    > "data/lifeline/lifeline.${HMO}.hg38.glm.linear"

  gzip -f "data/lifeline/lifeline.${HMO}.hg38.glm.linear"
done

### convert child data from hg19 to hg38
# bed format
for file in data/child/child_*.tsv; do
    base=$(basename "$file" .tsv)
    sample=${base#child_}

    awk 'BEGIN { FS=OFS="\t" } NR>1 { print $1, $2-1, $2, $9 }' \
        "$file" > "${sample}.hg19.bed"
done

# used liftOver locally
# 3210 vars not converted, 5470947 vars successful
for HMO in $(cat child_hmo_list.txt); do
    awk 'BEGIN{OFS="\t"} {$1="chr"$1; print}' "${HMO}.hg19.bed" > "${HMO}.chr.hg19.bed"

    liftOver \
        "${HMO}.chr.hg19.bed" \
        ../liftover_lifeline/hg19ToHg38.over.chain.gz \
        "${HMO}.hg38.bed" \
        "${HMO}.unmapped.bed"
done

# # back on CARC
# Make an ID -> hg38 coordinate map from the lifted BED
# Rebuild original tsv file
for HMO in $(cat child_hmo_list.txt); do

  awk 'BEGIN{FS=OFS="\t"}


NR==FNR {
    chr[$4] = $1
    pos[$4] = $2 + 1   
    next
}


FNR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "chromosome") chrcol = i
        if ($i == "base_pair_location") bpcol = i
        if ($i == "rs_id") idcol = i
    }
    print
    next
}

{
    id = $idcol
    if (id in chr) {
        $chrcol = chr[id]
        $bpcol  = pos[id]
    }
    print
}
' ${HMO}.hg38.bed child_${HMO}.tsv > child_${HMO}_hg38.tsv
done
