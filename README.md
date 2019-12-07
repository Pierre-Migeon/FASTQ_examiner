# FASTQ_examiner
## A tool to check validity of FASTQ files, correct if necessary, and produce basic summary charts and statistics for these

*tool in developement...*

FASTQ examiner is a tool written in python to do basic sanity checking of FASTQ files. First, files are checked for validity, wrapping, and truncation. Wrapped files are unwrapped, and any malformed or truncated entries in the FASTQ files are removed. Subsequent to these steps, summary statistics and graphs for the input files are produced. 

usage: 
	fastq_looker.py [-h] -f1 FASTQ_1 [-f2 FASTQ_2] [-i]


Diagnostic graphs produced can be useful for understanding fastq data quality or other status. For example, extreme 5' nucleotide bias in this case suggests sequencing adapters have yet to be removed:
![graph](https://user-images.githubusercontent.com/8321639/70365885-95504880-1848-11ea-9321-5fb1756d2e7f.png)


I used this tool to quickly make some random fastq files to play with during development:
M. Frampton, R. Houlston (2012) Generation of Artificial FASTQ Files to Evaluate the Performance of Next-Generation Sequencing Pipelines.
*PLoS ONE* 7 (11), http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0049110

See the following for a description of the fastq format:
Cock PJ, Fields CJ, Goto N, Heuer ML, Rice PM. The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Res. 2010;38(6):1767–1771. doi:10.1093/nar/gkp1137

