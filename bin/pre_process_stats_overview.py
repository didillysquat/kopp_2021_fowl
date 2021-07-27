"""
A script for summarising the pre processing results
in terms of the number of reads at each step of the QC and after mapping.

As input it will take the full paths to the output directories and the path to the input meta-info.tsv.
"""

import sys
import pandas as pd
import os
import re

class PreProcSummary:
    def __init__(self):
        self.cwd = sys.argv[1]
        self.tsv_path = os.path.join(self.cwd, "file_list.tsv")

        # These are the headers/metrics that we will pull out of the output files from the above directories and populate the output df with
        self.columns = [
            "RGSM", "RGID",	"RGLB",	"RGPL",	"RGPU", "filename_one", "filename_two",
            "reads_pre_trim_one", "reads_pre_trim_two", "reads_pre_trim_total", "reads_post_trim_one",
            "reads_post_trim_two", "reads_post_trim_total", "reads_trimmed_lost_one", "reads_trimmed_lost_two", "reads_trimmed_lost_total",
            "reads_trimmed_lost_total_percent", "reads_mapped", "reads_mapped_and_paired", "reads_unmapped",
            "reads_mapped_percent", "average_coverage", "average_coverarge_stdev", "unpaired_reads_examined_for_duplication",
            "paired_reads_examined_for_duplication", "unpaired_read_duplicated", "paired_read_duplicates", "percent_duplication",
            "sequenced_library_complexity", "estimated_library_complexity", "percent_library_sequenced"
            ]

        self.meta_info_df = pd.read_csv(self.tsv_path, sep="\t")
        self.meta_info_df.set_index("file_name", inplace=True)

        # Make dict that is sample name to read_one read_two
        # Meta dict
        meta_info_dict = {}
        for file_name, ser in self.meta_info_df.iterrows():
            if "R1" in file_name:
                meta_info_dict[ser.RGSM] = [
                        ser["RGSM"], ser["RGID"], ser["RGLB"], ser["RGPL"], ser["RGPU"], file_name, file_name.replace("_R1_", "_R2_"),
                        0,0,0,0,
                        0,0,0,0,0,
                        0.0,0,0,0,
                        0.0,0.0,0.0,0,
                        0,0,0,0.0,
                        0,0,0.0
                    ]

        self.summary_df = pd.DataFrame.from_dict(orient="index", columns=self.columns, data=meta_info_dict)

    def populate_summary_dict(self):
        # Go metric by metric
        # Pre trim reads
        # <td>Total Sequences</td><td>10788937</td>
        pre_trim_file_list = [_ for _ in os.listdir(self.cwd) if _.endswith(".pre_trim.fastqc.html")]
        # Find the file that we should be working with
        print("Collecting pre-trimming seq counts")
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

        # Here the pre trim should be done.

        # Now the post
        # All of the fastqc files are named with the R1 string. The R1 and R2 are differentiated by 1U and 2U and 1P and 2P
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
            self.summary_df.at[sample, "reads_trimmed_lost_total_percent"] = self.summary_df.at[sample, "reads_trimmed_lost_total"] / self.summary_df.at[sample, "reads_pre_trim_total"]
        print("\n")

        print("Collecting mapping stats")
        mapping_stats_file_list = [_ for _ in os.listdir(self.cwd) if _.endswith("bam.stats.txt")]
        for sample in self.summary_df.index:
            print(f"\r{sample}", end="")
            read_name=self.summary_df.at[sample, "filename_one"]
            mapping_file = self._find_item_match(seq_file_name=read_name, list_of_objects=mapping_stats_file_list)
            
            reads_mapped = None
            reads_mapped_and_paired = None
            reads_unmapped = None
            with open(os.path.join(self.cwd, mapping_file), "r") as f:
                for line in [_.rstrip() for _ in f]:
                    if "reads mapped and paired" in line:
                        reads_mapped_and_paired = int(re.search('([0-9]+)', line).group(0))
                    elif "reads mapped" in line:
                        reads_mapped = int(re.search('([0-9]+)', line).group(0))
                    elif "reads unmapped" in line:
                        reads_unmapped = int(re.search('([0-9]+)', line).group(0))
            if reads_mapped is None or reads_mapped_and_paired is None or reads_unmapped is None:
                raise RuntimeError(f"Error extracting the mapping stats for {sample}")

            self.summary_df.at[sample, "reads_mapped"] = reads_mapped
            self.summary_df.at[sample, "reads_mapped_and_paired"] = reads_mapped_and_paired
            self.summary_df.at[sample, "reads_unmapped"] = reads_unmapped
            self.summary_df.at[sample, "reads_mapped_percent"] = 1 - (reads_unmapped/reads_mapped)
        print("\n")

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
            self.summary_df.at[sample, "average_coverarge_stdev"] = coverage_stdev
        print("\n")

        print("Collecting complexity stats")
        duplicate_stats_file_list = [_ for _ in os.listdir(self.cwd) if _.endswith(".metrics.txt")]
        for sample in self.summary_df.index:
            print(f"\r{sample}", end="")
            read_name=self.summary_df.at[sample, "filename_one"]
            duplicate_file = self._find_item_match(seq_file_name=read_name, list_of_objects=duplicate_stats_file_list)
            with open(os.path.join(self.cwd, duplicate_file), "r") as f:
                lines = [_.rstrip() for _ in f]
            for i, line in enumerate(lines):
                if line.startswith("LIBRARY"):
                    components = lines[i+1].split()
                    unpaired_reads_examined_for_duplication = int(components[1])
                    paired_reads_examined_for_duplication = int(components[2])
                    unpaired_read_duplicated = int(components[5])
                    paired_read_duplicates = int(components[6])
                    percent_duplication = float(components[8])
                    sequenced_library_complexity = int(components[9])
                    break

            self.summary_df.at[sample, "unpaired_reads_examined_for_duplication"] = unpaired_reads_examined_for_duplication
            self.summary_df.at[sample, "paired_reads_examined_for_duplication"] = paired_reads_examined_for_duplication
            self.summary_df.at[sample, "unpaired_read_duplicated"] = unpaired_read_duplicated
            self.summary_df.at[sample, "paired_read_duplicates"] = paired_read_duplicates
            self.summary_df.at[sample, "percent_duplication"] = percent_duplication

            self.summary_df.at[sample, "sequenced_library_complexity"] = (self.summary_df.at[sample, "reads_mapped"] + self.summary_df.at[sample, "reads_unmapped"]) - ((self.summary_df.at[sample, "reads_mapped"] + self.summary_df.at[sample, "reads_unmapped"]) * percent_duplication)
            self.summary_df.at[sample, "estimated_library_complexity"] = sequenced_library_complexity
            self.summary_df.at[sample, "percent_library_sequenced"] = self.summary_df.at[sample, "sequenced_library_complexity"] / self.summary_df.at[sample, "estimated_library_complexity"]
        print("\n")

        self.summary_df.to_csv(os.path.join(self.cwd, "preprocessing_overview.tsv"), sep="\t", index=False)

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
        
