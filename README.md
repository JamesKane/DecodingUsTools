# DecodingUs Tools

A suite of future tools intended to be used with the main [DecodingUs](https://github.com/JamesKane/decodingus) project.

Disclaimer: This project is barely ALPHA release quality.  Use at your own risk!

Usage: 
```shell
decodingus-tools coverage <BAM_FILE> -o <cov_report.txt>
```

Results in a coverage report for each contig in the BAM.

```shell
decodingus-tools find-y-branch <BAM_FILE> <OUTPUT_FILE>
```

Find the closest YDNA branch

```shell
decodingus-tools find-mt-branch <BAM_FILE> <OUTPUT_FILE>
```

Find the closest mtDNA branch

### Caching
The branch finding algorithms store files in these locations:
* Linux: ~/.cache/decodingus-tools
* macOS: ~/Library/Caches/com.decodingus.decodingus-tools
* Windows: %LOCALAPPDATA%\decodingus-tools\Cache

Each tree will only be retrieved once per week if the file exists.  Remove the weekly file if you wish to force a reload.

## [![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.

