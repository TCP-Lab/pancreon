# Pancreon
Channelome profiling in PANC-1 cell line

## Dependencies
### Bash
- kerblam (>= 1.0.0-rc.1)
- x.FASTQ (>= 2.0.0)
### R
- ggplot2 (>= 3.5.0)
- dplyr (>= 1.1.4)
- r4tcpl (>= 1.5.1)
- readABF (>= 1.0.2)

## Run the analysis
```bash
kerblam data fetch
kerblam run make_geneset
kerblam run quality
kerblam run profiler
```
