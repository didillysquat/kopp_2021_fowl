import sys

rtg_path = sys.argv[1]
try:
    rtg_iteration = int(sys.argv[2])
    out_file = f"rtg.stats.summary.iter_{rtg_iteration}.txt"
except IndexError:
    rtg_iteration = None
    out_file = f"rtg.stats.summary.txt"

snps = []
insertions = []
deletions = []
ti_tv = []
insertion_deletion_ratio = []
no_snps_warning_samples = []
bad_insertion_deletion_ratio_samples = []
bad_ti_tv_ratio_samples = []
sample_name = ''
with open(rtg_path, 'r') as f:
    for line in [_.rstrip() for _ in f]:
        if line.startswith("Sample Name:"):
            sample_name = line.split()[-1]
        elif line.startswith('SNPs'):
            num_snps = int(line.split()[-1])
            if num_snps > 0:
                snps.append(int(line.split()[-1]))
            else:
                no_snps_warning_samples.append(sample_name)
                snps.append(0)
            continue
        elif line.startswith('Insertions'):
            if not line.split()[-1] == '-':
                insertions.append(int(line.split()[-1]))
            else:
                insertions.append(0)
            continue
        elif line.startswith('Deletions'):
            if not line.split()[-1] == '-':
                deletions.append(int(line.split()[-1]))
            else:
                deletions.append(0)
            continue
        elif line.startswith('SNP Transitions/Transversions'):
            if not line.split()[-2] == '-':
                ti_tv.append(float(line.split()[-2]))
            else:
                bad_ti_tv_ratio_samples.append(sample_name)
                ti_tv.append(False)
            continue
        elif line.startswith('Insertion/Deletion ratio'):
            if not line.split()[-2] == '-':
                insertion_deletion_ratio.append(float(line.split()[-2]))
            else:
                bad_insertion_deletion_ratio_samples.append(sample_name)
                insertion_deletion_ratio.append(False)
            continue

with open(out_file, 'w') as f:
    f.write(f"Summary of the rtg vcfstats file:\n\n")
    f.write(f"total SNPs: {sum(snps)}\n\n")

    # Weighted by the number of SNPs in the sample
    weighted_average_insertions = sum([snp * ins for snp, ins in zip(snps, insertions)]) / sum(snps)
    f.write(f"weighted average insertions: {weighted_average_insertions:.2f}\n\n")

    # Weighted by the number of SNPs in the sample
    weighted_average_deletions = sum([snp * dele for snp, dele in zip(snps, deletions)]) / sum(snps)
    f.write(f"weighted average deletions: {weighted_average_deletions:.2f}\n\n")

    # Weighted by the number of SNPs in the sample
    # Don't count those samples that had bad ti/tv ratios
    tot_snp = 0
    weighted_count = 0
    for snp, t in zip(snps, ti_tv):
        if t:
            weighted_count += snp*t
            tot_snp += snp
    if tot_snp > 0:
        weighted_average_ti_tv = weighted_count/tot_snp
        f.write(f"Weighted average SNP Transitions / Transversions ratio: {weighted_average_ti_tv:.2f}\n")
        if bad_ti_tv_ratio_samples:
            f.write(f"WARNING: {len(bad_ti_tv_ratio_samples)} samples did not have Ti/Tv ratios counted. "
                    f"These samples were ignored in the calculation of this average.\n\n")
    else:
        f.write(f"WARNING: 0 samples had insertion/deletion ratios counted.")
    # Weighted by the number of SNPs in the sample
    tot_snp = 0
    weighted_count = 0
    for snp, indel in zip(snps, insertion_deletion_ratio):
        if indel:
            weighted_count += snp * indel
            tot_snp += snp
    if tot_snp > 0:
        weighted_average_ins_dels = weighted_count/tot_snp
        f.write(f"Weighted average SNP Insertion / Deletion ratio: {weighted_average_ins_dels:.2f}\n")
        if bad_insertion_deletion_ratio_samples:
            f.write(f"WARNING: {len(bad_insertion_deletion_ratio_samples)} samples did not have insertion/deletion "
                    f"ratios counted. "
                    f"These samples were ignored in the calculation of this average.\n\n")
    else:
        f.write(f"WARNING: 0 samples had insertion/deletion ratios counted.")

sys.exit(0)
