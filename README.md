# Pairwise sequence alignment

A simple pairwise sequence alignment CLI.

## Usage
```
A command line interface for easily aligning two sequences.

Usage: align.exe [OPTIONS] <X> <Y>

Arguments:
  <X>  First sequence
  <Y>  Second sequence

Options:
  -g, --global                   Use global alignment, default
  -s, --semi-global              Use semi-global alignment
  -l, --local                    Use local alignment
  -n, --line-width <LINE_WIDTH>  The number of characters to show on a single line in the alignment [default: 50]
  -h, --help                     Print help
  -V, --version                  Print version
```

## Example usage
```
AKTGLSHLGYGMDV AKEGLAFLGYGMDV
```
![example result](inc/example.png)