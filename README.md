# DecodingUs Tools

A suite of future tools intended to be used with the main [DecodingUs](https://github.com/JamesKane/decodingus) project.

Disclaimer: This project is barely ALPHA release quality.  Use at your own risk!

## Installation
```shell
git clone https://github.com/JamesKane/DecodingUsTools.git
cd DecodingUsTools
cargo build --release
```

After building, the binary will be in target/release/decodingus-tools.  Add this directory to your PATH or 
copy the binary to a location already in your PATH.

## Usage

### Coverage Analysis

```shell
decodingus-tools coverage --reference <REFERENCE_FILE> <BAM_FILE> -o <cov_report.txt>
```

Results in a coverage report for each contig in the BAM. Requires the original reference.
Intended to replace `samtools coverage` and `gatk CallableLoci` with a single tool.

### Y-DNA Branch Finding

```shell
decodingus-tools find-y-branch --reference <REFERENCE_FILE> <BAM_FILE> <OUTPUT_FILE>
```

Find the closest YDNA branch

|Haplogroup|Score|Matching_SNPs|Mismatching_SNPs|Ancestral_Matches|No_Calls|Total_SNPs|Cumulative_SNPs|Depth|
|----------|------|-------------|----------------|-----------------|--------|----------|---------------|-----|
|R-FGC29071|2.75|3|0|2|1|9|1917|53|
|R-FGC29067|3.08|7|0|0|2|4|1911|52|

### mtDNA Branch Finding

```shell
decodingus-tools find-mt-branch --reference <REFERENCE_FILE> <BAM_FILE> <OUTPUT_FILE>
```

Find the closest mtDNA branch

|Haplogroup| Score |Matching_SNPs|Mismatching_SNPs|Ancestral_Matches|No_Calls|Total_SNPs|Cumulative_SNPs|Depth|
|----------|-------|-------------|----------------|-----------------|--------|----------|---------------|-----|
|U5a1b1g|1.65|1|0|0|0|1|55|15|

### Fix Surjected BAM

```shell
decodingus-tools fix-surjected-bam --reference <REFERNECE_FILE> -o <OUTPUT_BAM> <SURJECTED_BAM>
```

When a GAM file is surjected back via `vg surject` to a linear BAM, the results are mixed up compared to a traditional linear reference. This automates:
- Reheadering the @SQ details with a known reference to match the order
- Removing the PanSN-spec prefix from the individual reads
- Invoking samtools to sort the final result


### Caching
The branch finding algorithms store files in these locations:
* Linux: ~/.cache/decodingus-tools
* macOS: ~/Library/Caches/com.decodingus.decodingus-tools
* Windows: %LOCALAPPDATA%\decodingus-tools\Cache

Each tree will only be retrieved once per week if the file exists.  Remove the weekly file if you wish to force a reload.

## [![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.

