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

// Gene-level fusion screen, run per caller. Filters on structure (BND/TRA, or
// DEL/DUP/INV above minFusionLen), evidence (reads and MAPQ where reported), and
// the curated list in data/cns_fusions.tsv; only curated matches are reportable.
// Breakends come from the ALT bracket for BND and INFO END/CHR2 otherwise. Genes
// are assigned from refGene gene bodies - the panel BED is exon-level, and
// breakpoints fall in introns. See the README for the full description.
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
        path "${id}_${caller}_EGFRvIII.tsv",           emit: egfrviii

    script:
        """
        if [ "${caller}" = "savana" ]; then
            # keep breakpoints backed by supplementary alignments, including
            # SOURCE=CIGAR/SUPPLEMENTARY; only CIGAR-only calls are dropped
            awk -F'\\t' '/^#/ {print; next} \$8 ~ /SOURCE=[^;]*SUPPLEMENTARY/ {print}' \
                ${svVcf} > input.vcf
        else
            cp ${svVcf} input.vcf
        fi

        zcat ${refGene} \
          | awk -F'\\t' 'BEGIN{OFS="\\t"} \$3 ~ /^chr[0-9XY]+\$/ {print \$3, \$5, \$6, \$13}' \
          | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct > gene_bodies.bed

        # exon bounds of the EGFR reference transcript, for the EGFRvIII check
        zcat ${refGene} \
          | awk -F'\\t' '\$2=="${params.egfrTranscript}" && \$3=="chr7" {print \$10; print \$11; exit}' \
          > egfr_exons.txt

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

def get_vaf(info, fmt, sample):
    # sniffles reports VAF, savana TUMOUR_AF (one value per breakend), severus a
    # VAF format field; fall back to DV/(DV+DR) where only counts are given.
    for key in ("VAF", "TUMOUR_AF"):
        v = info_get(info, key)
        if v:
            try:
                return float(v.split(",")[0])
            except ValueError:
                pass
    if fmt and sample:
        keys, vals = fmt.split(":"), sample.split(":")
        if "VAF" in keys and keys.index("VAF") < len(vals):
            try:
                return float(vals[keys.index("VAF")])
            except ValueError:
                pass
        if "DV" in keys and "DR" in keys:
            try:
                dv = float(vals[keys.index("DV")])
                dr = float(vals[keys.index("DR")])
                if dv + dr > 0:
                    return dv / (dv + dr)
            except (ValueError, IndexError):
                pass
    return None

# EGFRvIII is the in-frame loss of exons 2-7, so both breakends sit inside EGFR:
# the 5-prime one in intron 1, the 3-prime one in intron 7. Testing only one
# breakend flags amplicon rearrangements that delete the promoter instead, and a
# size floor misses the real thing - real deletions run to a few tens of kb, well
# under the textbook ~190 kb. Introns are read from the annotation rather than
# hard-coded: exon numbering is transcript-specific, and shorter EGFR
# transcripts lack exon 4, which shifts every downstream exon by one.
EGFR_INTRON1 = EGFR_INTRON7 = None
try:
    _f = open("egfr_exons.txt").read().split()
    _st = [int(x) for x in _f[0].rstrip(",").split(",")]
    _en = [int(x) for x in _f[1].rstrip(",").split(",")]
    if len(_st) > 8:
        EGFR_INTRON1 = (_en[0], _st[1])
        EGFR_INTRON7 = (_en[6], _st[7])
except (OSError, IndexError, ValueError):
    pass
if EGFR_INTRON1 is None:
    print("WARNING: EGFR reference transcript not found - skipping the EGFRvIII check")

def is_egfrviii(chrom, pos, chrom2, pos2, svtype):
    if EGFR_INTRON1 is None:
        return False
    if chrom != "chr7" or chrom2 != "chr7" or svtype not in ("DEL", "BND", "TRA"):
        return False
    a, b = min(pos, pos2), max(pos, pos2)
    return (EGFR_INTRON1[0] <= a <= EGFR_INTRON1[1]
            and EGFR_INTRON7[0] <= b <= EGFR_INTRON7[1])

rows, viii_rows = [], []
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

    # 1. evidence
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
    vaf = get_vaf(info, f[8] if len(f) > 8 else "", f[9] if len(f) > 9 else "")
    vaf_s = "." if vaf is None else f"{vaf * 100:.1f}%"

    # 2. EGFRvIII, before the size filter below: it is intragenic, so the gene
    # pairing further down would discard it (both breakends are in EGFR).
    # EGFRvIII is an intragenic deletion, not a gene fusion, so it is reported in
    # its own table. The caller's own SVTYPE is kept - savana encodes the event
    # as a breakend pair, and calling that a DEL would misstate what it reported.
    if is_egfrviii(chrom, pos, chrom2, pos2, svtype):
        viii_rows.append((chrom, min(pos, pos2), max(pos, pos2), svtype,
                          span or ".", sup or ".", vaf_s, mq or "."))
        continue

    # 3. structural sanity
    interchrom = chrom != chrom2
    if not (svtype in ("BND", "TRA") or interchrom or span >= MIN_LEN):
        continue

    ga, gb = genes_at(chrom, pos), genes_at(chrom2, pos2)
    if not ga or not gb:
        continue

    # one row per event: a breakpoint in a gene-dense locus must not fan out
    # into a row per gene combination
    hit = match_known(ga, gb, svtype, span)
    tier = str(hit[0]) if hit else "3"
    label = hit[2] if hit else f"{sorted(ga)[0]}--{sorted(gb)[0]}"
    rows.append((label, chrom, pos, chrom2, pos2, svtype,
                 span or ".", sup or ".", vaf_s, mq or ".", tier,
                 hit[1] if hit else "",
                 ",".join(sorted(ga)), ",".join(sorted(gb))))

seen, uniq = set(), []
for r in rows:
    key = tuple(sorted([(r[1], r[2]), (r[3], r[4])]))
    if key not in seen:
        seen.add(key); uniq.append(r)
uniq.sort(key=lambda r: (r[10], -int(r[7]) if str(r[7]).isdigit() else 0))

hdr = ("fusion\\tchrom1\\tpos1\\tchrom2\\tpos2\\tsvtype\\tspan\\treads\\tvaf\\tmapq"
       "\\ttier\\tentity\\tgenes1\\tgenes2\\n")
with open(f"{sample}_{caller}_fusions_all.tsv", "w") as out:
    out.write(hdr)
    for r in uniq:
        out.write("\\t".join(map(str, r)) + "\\n")

SIG = {"1": "entity-defining", "2": "potentially significant"}
rep = [r for r in uniq if r[10] in SIG]
with open(f"{sample}_{caller}_fusions_reportable.tsv", "w") as out:
    out.write(hdr.rstrip("\\n") + "\\tsignificance\\n")
    for r in rep:
        out.write("\\t".join(map(str, r)) + "\\t" + SIG[r[10]] + "\\n")
tier1 = [r for r in rep if r[10] == "1"]

print(f"{caller}: {len(uniq)} pass filters, {len(rep)} reportable, {len(tier1)} entity-defining")
for r in rep:
    print(f"  TIER{r[10]}  {r[0]}  {r[5]}  {r[6]} bp  {r[7]} reads  VAF {r[8]}  -> {r[11]}")

# EGFRvIII, in its own table: the same deletion reported by both mate records of
# a breakend pair collapses to one row
seen_v, viii = set(), []
for r in sorted(viii_rows, key=lambda r: r[1]):
    key = (r[1] // 1000, r[2] // 1000)
    if key not in seen_v:
        seen_v.add(key); viii.append(r)
with open(f"{sample}_{caller}_EGFRvIII.tsv", "w") as out:
    out.write("chrom\\tpos1\\tpos2\\tsvtype\\tspan\\treads\\tvaf\\tmapq\\n")
    for r in viii:
        out.write("\\t".join(map(str, r)) + "\\n")

print(f"{caller}: EGFRvIII candidates: {len(viii)}")
for r in viii:
    print(f"  EGFRvIII  {r[0]}:{r[1]}-{r[2]}  {r[4]} bp  {r[5]} reads  VAF {r[6]}")
PY

        """
}
