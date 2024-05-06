# Freyja Summarizer
Creates a summary file (Freyja_aggregated_summary.csv) from Freyja aggregated lineage data (aggregated-freyja.tsv).

## Getting Started
### Getting the script
Clone the Freyja summarizer repository to your local machine:

```bash
git clone https://github.com/DrB-S/freyja_summarizer.git
```

`freyja_summarizer` can be found in the freyja_summarizer directory just created.

### Getting the dependencies

freyja_summarizer requires the following dependencies:
  - python3
    - argparse
    - logging
    - polars
    - re

### Install dependencies
```bash
pip install argparse logging polars re
```

## Running freyja_summarizer.py
```usage: freyja_summarizer.py [-h] [-i INPUT] [-o OUT] [-t TYPE] [-v]

options:
  -h, --help                show this help message and exit  
  -i INPUT, --input INPUT   input Freyja summary file name (default = 'aggregated_freyja.tsv')
 
  -o OUT, --out OUT         final file name (default = 'Freyja_aggregated_summary')
  -t TYPE, --type TYPE      file extention for final image (default = 'csv')
  -v, --version             print version and exit
  ```
