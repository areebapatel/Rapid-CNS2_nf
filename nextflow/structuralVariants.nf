// Sniffles2 - fast baseline SV caller
process structuralVariants {
    label 'heavy'

    publishDir "${params.outDir}/sv/", mode: 'copy'

    input:
        tuple path(bam), path(bai)
        path(ref)
        val(id)

    output:
        path "${id}.sniffles2.vcf", emit: vcf

    script:
        // Tumour-only input, so somatic mode is the default. Sniffles renamed
        // --non-germline to --mosaic; 2.8.0 rejects the old flag outright.
        def mosaic = params.snifflesMosaic ? '--mosaic' : ''
        """
        sniffles --threads ${task.cpus} --allow-overwrite \
            ${mosaic} \
            --reference ${ref} \
            --input ${bam} \
            --vcf ${id}.sniffles2.vcf
        """
}

// Severus - somatic SV caller using a breakpoint graph. Better on complex
// rearrangements (chromothripsis, EGFRvIII-type events) than Sniffles alone.
// Tumour-only here: without --control-bam, a panel of normals (--PON) is the
// only somatic filter, so supply one via params.severusPON where possible.
process severus {
    label 'severus'

    publishDir "${params.outDir}/sv/severus/", mode: 'copy'

    // stageAs gives the two optional resources distinct staged names: both
    // default to the same NO_FILE placeholder, and Nextflow refuses to stage
    // two inputs with the same filename into one task directory.
    input:
        tuple path(bam), path(bai)
        val(id)
        path(vntrBed, stageAs: 'severus_vntr.bed')
        path(pon,     stageAs: 'severus_pon.tsv.gz')

    output:
        path "severus_${id}/**", emit: allOutputs
        path "${id}.severus.vcf", emit: vcf

    script:
        """
        VNTR_ARG=""
        if [ -s severus_vntr.bed ]; then VNTR_ARG="--vntr-bed severus_vntr.bed"; fi
        PON_ARG=""
        if [ -s severus_pon.tsv.gz ];  then PON_ARG="--PON severus_pon.tsv.gz"; fi

        severus \
            --target-bam ${bam} \
            --out-dir severus_${id} \
            -t ${task.cpus} \
            \${VNTR_ARG} \
            \${PON_ARG}

        # Tumour-only runs emit all_SVs/severus_all.vcf; with a PON a somatic
        # set is produced too. Prefer the somatic call set when present.
        if [ -f severus_${id}/somatic_SVs/severus_somatic.vcf ]; then
            cp severus_${id}/somatic_SVs/severus_somatic.vcf ${id}.severus.vcf
        elif [ -f severus_${id}/all_SVs/severus_all.vcf ]; then
            cp severus_${id}/all_SVs/severus_all.vcf ${id}.severus.vcf
        else
            FOUND=\$(find severus_${id} -name "*.vcf" | head -n 1)
            [ -n "\${FOUND}" ] || { echo "ERROR: Severus produced no VCF" >&2; exit 1; }
            cp "\${FOUND}" ${id}.severus.vcf
        fi
        """
}

process annotSV {
    label 'rapid_cns'

    publishDir "${params.outDir}/sv/", mode: 'copy'

    input:
        path(svVcf)
        val(annotSvAnnot)
        val(id)
        val(caller)

    output:
        path "${id}_${caller}_annotsv.tsv", emit: svAnno

    script:
        // AnnotSV defaults to GRCh37; this pipeline is hg38 only.
        """
        AnnotSV \
            -SVinputFile ${svVcf} \
            -outputFile ${id}_${caller}_annotsv.tsv \
            -outputDir . \
            -genomeBuild GRCh38 \
            -annotationsDir ${annotSvAnnot}
        """
}



// Gene-level fusion screen, run per caller.
//
// Three filters keep this specific enough to put in a report:
//   1. structural sanity - only rearrangements that can actually fuse two
//      genes (BND/TRA, or DEL/DUP/INV above minFusionLen). Without this a
//      50 bp deletion straddling a gene boundary counts as a fusion, which is
//      why an unfiltered Sniffles VCF yields more candidates than records.
//   2. evidence - read support and mapping quality where the caller reports it.
//   3. a curated allowlist (data/cns_fusions.tsv) matched on gene pair AND
//      expected structure; only those become Tier 1, i.e. reportable.
//
// Both breakends are derived generically: a BND takes its partner from the ALT
// bracket, everything else from INFO END/CHR2. Gene assignment uses whole gene
// bodies from refGene, because the panel BED is exon-level and breakpoints fall
// in introns.
process svFusions {
    label 'rapid_cns'

    publishDir "${params.outDir}/sv/fusions/", mode: 'copy'

    input:
        path(svVcf)
        path(refGene)
        path(knownFusions)
        val(id)
        val(caller)

    output:
        path "${id}_${caller}_fusions_reportable.tsv", emit: reportable
        path "${id}_${caller}_fusions_all.tsv",        emit: allCandidates
        path "${id}_${caller}_EGFRvIII.txt",           emit: egfrviii

    script:
        """
        if [ "${caller}" = "savana" ]; then
            awk -F'\\t' '/^#/ {print; next} \$8 ~ /SOURCE=SUPPLEMENTARY/ {print}' \
                ${svVcf} > input.vcf
        else
            cp ${svVcf} input.vcf
        fi

        zcat ${refGene} \
          | awk -F'\\t' 'BEGIN{OFS="\\t"} \$3 ~ /^chr[0-9XY]+\$/ {print \$3, \$5, \$6, \$13}' \
          | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct > gene_bodies.bed

        python3 - "${id}" "${caller}" "${knownFusions}" \
                 "${params.minFusionLen}" "${params.minFusionReads}" "${params.minFusionMapq}" <<'PY'
import bisect, re, sys

sample, caller, known_f = sys.argv[1], sys.argv[2], sys.argv[3]
MIN_LEN, MIN_READS, MIN_MAPQ = int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6])

genes, starts = {}, {}
for line in open("gene_bodies.bed"):
    c, s, e, g = line.rstrip("\\n").split("\\t")
    genes.setdefault(c, []).append((int(s), int(e), g))
for c, v in genes.items():
    starts[c] = [x[0] for x in v]

def genes_at(chrom, pos):
    v = genes.get(chrom)
    if not v:
        return set()
    i = bisect.bisect_right(starts[chrom], pos)
    out = set()
    for s, e, g in v[max(0, i - 40):i]:
        if s <= pos <= e:
            out.update(g.split(","))
    return out

# curated list: (geneA, geneB) -> (svtype, min_len, max_len, entity)
known = {}
for line in open(known_f):
    if line.startswith("#") or line.startswith("gene5"):
        continue
    f = line.rstrip("\\n").split("\\t")
    if len(f) < 6:
        continue
    known[(f[0], f[1])] = (f[2], int(f[3]), int(f[4]), f[5])

def match_known(ga, gb, svtype, span):
    # Best match over all gene combinations at the two breakends. A named pair
    # (both genes given) outranks a wildcard entry, so KIAA1549-BRAF reports as
    # pilocytic astrocytoma rather than the generic BRAF-fusion label.
    best = None
    for (ka, kb), (ktype, lo, hi, entity) in known.items():
        if ka == "*" and kb == "*":
            continue
        specific = "*" not in (ka, kb)
        for x in ga:
            for y in gb:
                if x == y:
                    continue
                for p, q in ((x, y), (y, x)):
                    if (ka == p or ka == "*") and (kb == q or kb == "*"):
                        if (ktype != "ANY" and ktype != svtype
                                and svtype not in ("BND", "TRA")):
                            continue
                        if lo and span and not (lo <= span <= hi):
                            continue
                        cand = (1 if specific else 2, entity, f"{p}--{q}")
                        if best is None or cand[0] < best[0]:
                            best = cand
    return best

def info_get(info, key):
    m = re.search(rf"(?:^|;){key}=([^;\\t]+)", info)
    return m.group(1) if m else None

rows = []
for line in open("input.vcf"):
    if line.startswith("#"):
        continue
    f = line.rstrip("\\n").split("\\t")
    chrom, pos, alt, info = f[0], int(f[1]), f[4], f[7]
    svtype = info_get(info, "SVTYPE") or "."

    m = re.search(r"[\\[\\]]([^\\[\\]:]+):(\\d+)[\\[\\]]", alt)
    if m:
        chrom2, pos2 = m.group(1), int(m.group(2))
    else:
        chrom2 = info_get(info, "CHR2") or chrom
        end = info_get(info, "END")
        if end is None:
            continue
        pos2 = int(end)

    span = abs(pos2 - pos) if chrom == chrom2 else 0

    # 1. structural sanity
    interchrom = chrom != chrom2
    if not (svtype in ("BND", "TRA") or interchrom or span >= MIN_LEN):
        continue

    # 2. evidence
    sup = (info_get(info, "TUMOUR_READ_SUPPORT") or info_get(info, "SUPPORT") or "")
    if not sup:
        sr = info_get(info, "SUPP_READS")          # severus: a:b:c:d:e:f
        if sr:
            try:
                sup = str(max(int(x) for x in sr.split(":")))
            except ValueError:
                sup = ""
    if sup.isdigit() and int(sup) < MIN_READS:
        continue
    mq = info_get(info, "MAPQ") or info_get(info, "ORIGIN_MAPQ_MEAN")
    if mq:
        try:
            if float(mq) < MIN_MAPQ:
                continue
        except ValueError:
            pass

    ga, gb = genes_at(chrom, pos), genes_at(chrom2, pos2)
    if not ga or not gb:
        continue

    # one row per event: a breakpoint in a gene-dense locus must not fan out
    # into a row per gene combination
    hit = match_known(ga, gb, svtype, span)
    tier = str(hit[0]) if hit else "3"
    label = hit[2] if hit else f"{sorted(ga)[0]}--{sorted(gb)[0]}"
    rows.append((label, chrom, pos, chrom2, pos2, svtype,
                 span or ".", sup or ".", mq or ".", tier,
                 hit[1] if hit else "",
                 ",".join(sorted(ga)), ",".join(sorted(gb))))

seen, uniq = set(), []
for r in rows:
    key = tuple(sorted([(r[1], r[2]), (r[3], r[4])]))
    if key not in seen:
        seen.add(key); uniq.append(r)
uniq.sort(key=lambda r: (r[9], -int(r[7]) if str(r[7]).isdigit() else 0))

hdr = "fusion\\tchrom1\\tpos1\\tchrom2\\tpos2\\tsvtype\\tspan\\treads\\tmapq\\ttier\\tentity\\tgenes1\\tgenes2\\n"
with open(f"{sample}_{caller}_fusions_all.tsv", "w") as out:
    out.write(hdr)
    for r in uniq:
        out.write("\\t".join(map(str, r)) + "\\n")

SIG = {"1": "entity-defining", "2": "potentially significant"}
rep = [r for r in uniq if r[9] in SIG]
with open(f"{sample}_{caller}_fusions_reportable.tsv", "w") as out:
    out.write(hdr.rstrip("\\n") + "\\tsignificance\\n")
    for r in rep:
        out.write("\\t".join(map(str, r)) + "\\t" + SIG[r[9]] + "\\n")
tier1 = [r for r in rep if r[9] == "1"]

print(f"{caller}: {len(uniq)} pass filters, {len(rep)} reportable, {len(tier1)} entity-defining")
for r in rep:
    print(f"  TIER{r[9]}  {r[0]}  {r[5]}  {r[6]} bp  {r[7]} reads  -> {r[10]}")
PY

        awk -F'\\t' '!/^#/ && \$1=="chr7" && \$2>=55019017 && \$2<=55211628 {
                 len = 0
                 if (match(\$8, /SVLEN=-?[0-9]+/)) len = substr(\$8, RSTART+6, RLENGTH-6) + 0
                 if (len < 0) len = -len
                 if (len >= 50000 && len <= 250000)
                     print "CANDIDATE\\t" \$1 "\\t" \$2 "\\t" \$3 "\\tSVLEN=" len
             }' input.vcf > egfrviii.hits || true
        if [ -s egfrviii.hits ]; then
            { echo "EGFRvIII (${caller}): candidate intragenic EGFR deletion(s)"; cat egfrviii.hits; } \
                > ${id}_${caller}_EGFRvIII.txt
        else
            echo "EGFRvIII (${caller}): not detected" > ${id}_${caller}_EGFRvIII.txt
        fi
        cat ${id}_${caller}_EGFRvIII.txt
        """
}
