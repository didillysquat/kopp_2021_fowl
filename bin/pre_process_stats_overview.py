"""
A script for summarising the pre processing results
in terms of the number of reads at each step of the QC and after mapping.

As input it will take the full paths to the output directories and the path to the input meta-info.tsv.
"""

import sys
import pandas as pd
import os
import re
import ntpath

class PreProcSummary:
    def __init__(self):
        self.cwd = sys.argv[1]
        self.tsv_path = os.path.join(self.cwd, sys.argv[2])

        # These are the headers/metrics that we will pull out of the output files from the above directories and populate the output df with
        self.columns = [
            "RGSM", "RGID",	"RGLB",	"RGPL",	"RGPU", "filename_one", "filename_two",
            # Pre adapter trim stats
            "reads_pre_trim_one", "reads_pre_trim_two", "reads_pre_trim_total",
            # Post adapter trim stats
            "reads_post_trim_one", "reads_post_trim_two", "reads_post_trim_total",
            # Reads lost
            "reads_trimmed_lost_one", "reads_trimmed_lost_two", "reads_trimmed_lost_total", "reads_trimmed_lost_total_proportion",
            # Pre deduplication bam stats
            "raw_total_sequences_pre_deduplication",	"reads_mapped_pre_deduplication",	"reads_mapped_and_paired_pre_deduplication",	"reads_unmapped_pre_deduplication",	"reads_mapped_proportion_pre_deduplication",
            # Deduplication MarkDuplicatesSpark
            "unpaired_reads_examined_for_deduplication_MarkDuplicatesSpark", "read_pairs_examined_for_deduplication_MarkDuplicatesSpark", "unpaired_read_duplicated_MarkDuplicatesSpark",
            "read_pair_duplicates_MarkDuplicatesSpark", "proportion_duplication_MarkDuplicatesSpark", "estimated_library_complexity_MarkDuplicatesSpark", "proportion_library_sequenced_MarkDuplicatesSpark",
            # Deduplication EstimateLibraryComplexity
            "read_pairs_examined_EstimateLibraryComplexity", "read_pair_duplicates_EstimateLibraryComplexity", "proportion_duplication_EstimateLibraryComplexity",
            "estimated_library_complexity_EstimateLibraryComplexity", "proportion_library_sequenced_EstimateLibraryComplexity",
            # Post deduplication bam stats
            "raw_total_sequences_post_deduplication", "reads_mapped_post_deduplication", "reads_mapped_and_paired_post_deduplication", "reads_unmapped_post_deduplication", "reads_mapped_proportion_post_deduplication",
            # Coverage stats
            "average_coverage", "average_coverage_stdev",
            ]

        self.meta_info_df = pd.read_csv(self.tsv_path, sep="\t")
        self.meta_info_df.set_index("RGSM", inplace=True, drop=False)

        # Make dict that is sample name to read_one read_two
        # Meta dict
        meta_info_dict = {}
        for sample_name, ser in self.meta_info_df.iterrows():
            meta_info_dict[ser.RGSM] = [
                    ser["RGSM"], ser["RGID"], ser["RGLB"], ser["RGPL"], ser["RGPU"], ntpath.basename(ser["file_name_one"]), ntpath.basename(ser["file_name_two"]),
                    # Pre adapter trim stats
                    0,0,0,
                    # Post adapter trim stats
                    0,0,0,
                    # Reads lost
                    0,0,0,0.0,
                    # Pre deduplication bam stats
                    0,0,0,0,0.0,
                    # Deduplication MarkDuplicatesSpark
                    0,0,0,
                    0,0.0,0,0.0,
                    # Deduplication EstimateLibraryComplexity
                    0,0,0.0,
                    0,0.0,
                    # Post deduplication bam stats
                    0,0,0,0,0.0,
                    # Coverage stats
                    0.0,0.0
                ]

        self.summary_df = pd.DataFrame.from_dict(orient="index", columns=self.columns, data=meta_info_dict)

    def populate_summary_dict(self):
    
        self._pre_trim()

        self._post_trim_seq_counts()

        self._mapping_stats_pre_deduplication()

        self._markduplicatesspark_complexity()
        
        self._estimatelibrarycomplexity_complexity()
        
        self._mapping_stats_post_deduplication()     

        self._coverage_stats()

        self.summary_df.to_csv(os.path.join(self.cwd, "preprocessing_overview.tsv"), sep="\t", index=False)

    def _pre_trim(self):
        print("Collecting pre-trimming seq counts")
        pre_trim_file_list = [_ for _ in os.listdir(self.cwd) if _.endswith(".pre_trim.fastqc.html")]
        for sample in self.summary_df.index:
            print(f"\r{sample}", end="")
            # Do read 1
            extracted_seq_count = self._extract_seq_count_from_html(read_name=self.summary_df.at[sample, "filename_one"], list_of_objects=pre_trim_file_list, base_data_dir=self.cwd)
            self.summary_df.at[sample, "reads_pre_trim_one"] = extracted_seq_count

            # Do read 2
            extracted_seq_count = self._extract_seq_count_from_html(read_name=self.summary_df.at[sample, "filename_two"], list_of_objects=pre_trim_file_list, base_data_dir=self.cwd)
            self.summary_df.at[sample, "reads_pre_trim_two"] = extracted_seq_count

            # Do total
            self.summary_df.at[sample, "reads_pre_trim_total"] = self.summary_df.at[sample, "reads_pre_trim_one"] + self.summary_df.at[sample, "reads_pre_trim_two"]
        print("\n")

    def _coverage_stats(self):
        print("Collecting coverage stats")
        coverage_stats_file_list = [_ for _ in os.listdir(self.cwd) if _.endswith(".depth.txt")]
        for sample in self.summary_df.index:
            print(f"\r{sample}", end="")
            read_name=self.summary_df.at[sample, "filename_one"]
            mapping_file = self._find_item_match(seq_file_name=read_name, list_of_objects=coverage_stats_file_list)
            average_coverage = None
            coverage_stdev = None
            with open(os.path.join(self.cwd, mapping_file), "r") as f:
                for line in [_.rstrip() for _ in f]:
                    if "Average" in line:
                        average_coverage = float(re.search('([0-9]+\.[0-9]+)', line).group(0))
                    elif "Stdev" in line:
                        coverage_stdev = float(re.search('([0-9]+\.[0-9]+)', line).group(0))
            if average_coverage is None or coverage_stdev is None:
                raise RuntimeError(f"Error in extracting coverage stats for {sample}")
            self.summary_df.at[sample, "average_coverage"] = average_coverage
            self.summary_df.at[sample, "average_coverage_stdev"] = coverage_stdev
        print("\n")

    def _mapping_stats_post_deduplication(self):
        print("Collecting mapping stats post deduplication")
        post_dedup_mapping_stats_file_list = [_ for _ in os.listdir(self.cwd) if _.endswith("bam.postdedup.stats.txt")]
        for sample in self.summary_df.index:
            raw_total_sequences, reads_mapped,  reads_mapped_and_paired, reads_unmapped = self._extract_bam_stats(sample=sample, file_list=post_dedup_mapping_stats_file_list)
            self.summary_df.at[sample, "raw_total_sequences_post_deduplication"] = raw_total_sequences
            self.summary_df.at[sample, "reads_mapped_post_deduplication"] = reads_mapped
            self.summary_df.at[sample, "reads_mapped_and_paired_post_deduplication"] = reads_mapped_and_paired
            self.summary_df.at[sample, "reads_unmapped_post_deduplication"] = reads_unmapped
            self.summary_df.at[sample, "reads_mapped_proportion_post_deduplication"] = reads_mapped/raw_total_sequences
        print("\n")

    def _estimatelibrarycomplexity_complexity(self):
        print("Collecting EstimateLibraryComplexity stats")
        duplicate_stats_file_list = [_ for _ in os.listdir(self.cwd) if _.endswith(".estlibcomp.txt")]
        for sample in self.summary_df.index:
            print(f"\r{sample}", end="")
            read_name=self.summary_df.at[sample, "filename_one"]
            duplicate_file = self._find_item_match(seq_file_name=read_name, list_of_objects=duplicate_stats_file_list)
            with open(os.path.join(self.cwd, duplicate_file), "r") as f:
                lines = [_.rstrip() for _ in f]
            for i, line in enumerate(lines):
                if line.startswith("LIBRARY"):
                    components = lines[i+1].split()
                    READ_PAIRS_EXAMINED = int(components[2])
                    READ_PAIR_DUPLICATES = int(components[6])
                    PERCENT_DUPLICATION = float(components[8])
                    ESTIMATED_LIBRARY_SIZE = int(components[9])
                    break

            self.summary_df.at[sample, "read_pairs_examined_EstimateLibraryComplexity"] = READ_PAIRS_EXAMINED
            self.summary_df.at[sample, "read_pair_duplicates_EstimateLibraryComplexity"] = READ_PAIR_DUPLICATES
            self.summary_df.at[sample, "proportion_duplication_EstimateLibraryComplexity"] = PERCENT_DUPLICATION

            # NB the complexity estimates are calculated using the MarkDuplicates from GATK and so are based only on the mapped reads (paired and unpaired)
            # As such, the sequenced_library_complexity is the reads_mapped statistic from the samtools stats output that was run on the deduplicated
            # bam file.
            self.summary_df.at[sample, "estimated_library_complexity_EstimateLibraryComplexity"] = ESTIMATED_LIBRARY_SIZE
            self.summary_df.at[sample, "proportion_library_sequenced_EstimateLibraryComplexity"] = (READ_PAIRS_EXAMINED - READ_PAIR_DUPLICATES) / ESTIMATED_LIBRARY_SIZE
        print("\n")

    def _markduplicatesspark_complexity(self):
        print("Collecting MarkDuplicatesSpark complexity stats")
        mark_dup_spark_deduplicate_stats_file_list = [_ for _ in os.listdir(self.cwd) if _.endswith(".merged.deduplicated.sorted.metrics.txt")]
        for sample in self.summary_df.index:
            print(f"\r{sample}", end="")
            read_name=self.summary_df.at[sample, "filename_one"]
            duplicate_file = self._find_item_match(seq_file_name=read_name, list_of_objects=mark_dup_spark_deduplicate_stats_file_list)
            with open(os.path.join(self.cwd, duplicate_file), "r") as f:
                lines = [_.rstrip() for _ in f]
            for i, line in enumerate(lines):
                if line.startswith("LIBRARY"):
                    components = lines[i+1].split()
                    unpaired_reads_examined_for_deduplication = int(components[1])
                    read_pairs_examined_for_deduplication = int(components[2])
                    unpaired_read_duplicated = int(components[5])
                    paired_read_duplicates = int(components[6])
                    proportion_duplication = float(components[8])
                    estimated_library_complexity = int(components[9])
                    break
            
            self.summary_df.at[sample, "unpaired_reads_examined_for_deduplication_MarkDuplicatesSpark"] = unpaired_reads_examined_for_deduplication
            self.summary_df.at[sample, "read_pairs_examined_for_deduplication_MarkDuplicatesSpark"] = read_pairs_examined_for_deduplication
            self.summary_df.at[sample, "unpaired_read_duplicated_MarkDuplicatesSpark"] = unpaired_read_duplicated
            self.summary_df.at[sample, "read_pair_duplicates_MarkDuplicatesSpark"] = paired_read_duplicates
            self.summary_df.at[sample, "proportion_duplication_MarkDuplicatesSpark"] = proportion_duplication

            # NB the complexity estimates are calculated using the MarkDuplicates from GATK and so are based only on the mapped reads (paired and unpaired)
            # As such, the sequenced_library_complexity is the reads_mapped statistic from the samtools stats output that was run on the deduplicated
            # bam file.
            self.summary_df.at[sample, "estimated_library_complexity_MarkDuplicatesSpark"] = estimated_library_complexity
            # Because the sequence complexity estimate is the number of sequences (not reads!; i.e. a sequence needs two reads)
            # we will calculate the proportion of library seuquenced based on the number of unique read pairs
            sequenced_complexity = (read_pairs_examined_for_deduplication - paired_read_duplicates)
            self.summary_df.at[sample, "proportion_library_sequenced_MarkDuplicatesSpark"] = sequenced_complexity / estimated_library_complexity
        print("\n")

    def _mapping_stats_pre_deduplication(self):
        print("Collecting mapping stats pre deduplication")
        pre_dedup_mapping_stats_file_list = [_ for _ in os.listdir(self.cwd) if _.endswith("bam.prededup.stats.txt")]
        for sample in self.summary_df.index:
            raw_total_sequences, reads_mapped,  reads_mapped_and_paired, reads_unmapped = self._extract_bam_stats(sample=sample, file_list=pre_dedup_mapping_stats_file_list)

            self.summary_df.at[sample, "raw_total_sequences_pre_deduplication"] = raw_total_sequences
            self.summary_df.at[sample, "reads_mapped_pre_deduplication"] = reads_mapped
            self.summary_df.at[sample, "reads_mapped_and_paired_pre_deduplication"] = reads_mapped_and_paired
            self.summary_df.at[sample, "reads_unmapped_pre_deduplication"] = reads_unmapped
            self.summary_df.at[sample, "reads_mapped_proportion_pre_deduplication"] = reads_mapped/raw_total_sequences
        print("\n")

    def _post_trim_seq_counts(self):
        print("Collecting post-trimming seq counts")
        post_trim_file_list = [_ for _ in os.listdir(self.cwd) if _.endswith("_fastqc.html")]
        for sample in self.summary_df.index:
            print(f"\r{sample}", end="")
            # Read 1 
            read_one_paired_count = self._extract_seq_count_from_html(read_name=self.summary_df.at[sample, "filename_one"], list_of_objects=post_trim_file_list, base_data_dir=self.cwd, must_contain="1P")
            read_one_unpaired_count = self._extract_seq_count_from_html(read_name=self.summary_df.at[sample, "filename_one"], list_of_objects=post_trim_file_list, base_data_dir=self.cwd, must_contain="1U")
            self.summary_df.at[sample, "reads_post_trim_one"] = read_one_paired_count + read_one_unpaired_count

            # Read 2
            read_two_paired_count = self._extract_seq_count_from_html(read_name=self.summary_df.at[sample, "filename_one"], list_of_objects=post_trim_file_list, base_data_dir=self.cwd, must_contain="2P")
            read_two_unpaired_count = self._extract_seq_count_from_html(read_name=self.summary_df.at[sample, "filename_one"], list_of_objects=post_trim_file_list, base_data_dir=self.cwd, must_contain="2U")
            self.summary_df.at[sample, "reads_post_trim_two"] = read_two_paired_count + read_two_unpaired_count

            # Do total
            self.summary_df.at[sample, "reads_post_trim_total"] = self.summary_df.at[sample, "reads_post_trim_one"] + self.summary_df.at[sample, "reads_post_trim_two"]

            # Do trimming summary stats
            self.summary_df.at[sample, "reads_trimmed_lost_one"] = self.summary_df.at[sample, "reads_pre_trim_one"] - self.summary_df.at[sample, "reads_post_trim_one"]
            self.summary_df.at[sample, "reads_trimmed_lost_two"] = self.summary_df.at[sample, "reads_pre_trim_two"] - self.summary_df.at[sample, "reads_post_trim_two"]
            self.summary_df.at[sample, "reads_trimmed_lost_total"] = self.summary_df.at[sample, "reads_pre_trim_total"] - self.summary_df.at[sample, "reads_post_trim_total"]
            self.summary_df.at[sample, "reads_trimmed_lost_total_proportion"] = self.summary_df.at[sample, "reads_trimmed_lost_total"] / self.summary_df.at[sample, "reads_pre_trim_total"]
        print("\n")

    def _extract_bam_stats(self, sample, file_list):
        print(f"\r{sample}", end="")
        read_name=self.summary_df.at[sample, "filename_one"]
        mapping_file = self._find_item_match(seq_file_name=read_name, list_of_objects=file_list)
        raw_total_sequences = None
        reads_mapped = None
        reads_mapped_and_paired = None
        reads_unmapped = None
        with open(os.path.join(self.cwd, mapping_file), "r") as f:
            for line in [_.rstrip() for _ in f]:
                if "raw total sequences" in line:
                    raw_total_sequences = int(re.search('([0-9]+)', line).group(0))
                elif "reads mapped and paired" in line:
                    reads_mapped_and_paired = int(re.search('([0-9]+)', line).group(0))
                elif "reads mapped" in line:
                    reads_mapped = int(re.search('([0-9]+)', line).group(0))
                elif "reads unmapped" in line:
                    reads_unmapped = int(re.search('([0-9]+)', line).group(0))
        if raw_total_sequences is None or reads_mapped is None or reads_mapped_and_paired is None or reads_unmapped is None:
            raise RuntimeError(f"Error extracting the mapping stats for {sample}")
        return raw_total_sequences, reads_mapped, reads_mapped_and_paired, reads_unmapped

    def _extract_seq_count_from_html(self, read_name, list_of_objects, base_data_dir, must_contain=""):
        matched_file = self._find_item_match(seq_file_name=read_name, list_of_objects=list_of_objects, must_contain=must_contain)
        with open(os.path.join(base_data_dir, matched_file), "r") as f:
            html_line = f.read()
        return int(re.search('<td>Total Sequences</td><td>([0-9]+)</td>', html_line).group(1))
        
    def _find_item_match(self, seq_file_name, list_of_objects, must_contain=""):
        """given a seq file name, unique match to a single object in the list of objects"""
        for i in range(1, len(seq_file_name)):
            restart = False
            seq_file_name_short = ''.join(list(seq_file_name)[:i])
            matches = []
            # now look to see how many matches there are
            for file_obj_name in list_of_objects:
                if seq_file_name_short in file_obj_name and must_contain in file_obj_name:
                    matches.append(file_obj_name)
                    if len(matches) > 1:
                        restart = True
                        break
            if not restart:
                if len(matches) == 1:
                    return matches[0]
                else:
                    if len(matches) == 0:
                        # if the last character is '_', try replacing this with '.' and see if this gives a unique match
                        if seq_file_name_short[-1] == "_":
                            seq_file_name_short = seq_file_name_short[:-1] + "."
                            for file_obj_name in list_of_objects:
                                if seq_file_name_short in file_obj_name and must_contain in file_obj_name:
                                    matches.append(file_obj_name)
                            if len(matches) == 1:
                                return matches[0]
                        else:
                            raise RuntimeError("No unique match found")
        raise RuntimeError("No unique match found")

PreProcSummary().populate_summary_dict()
        
