#!/usr/bin/python3

###########################################
# Freyja Summarizer v0.9.0                #
# Written by Stephen Beckstrom-Sternberg  #
# Assisted by Khalil Khoury               #
# Creates a simplified Freyja summary     #
# with lineage abundances & coverage      #
# from aggregated Freyja file             #
# using polars instead of pandas          #
###########################################

import re
import argparse
import logging
import polars as pl

logging.basicConfig(
    format="%(asctime)s - %(message)s", datefmt="%y-%b-%d %H:%M:%S", level=logging.INFO
)

"""
Parse command-line arguments or use defaults
    
    Image types: eps, pdf (default), png, svg, tiff
"""
parser = argparse.ArgumentParser()
parser.add_argument(
    "-i",
    "--input",
    type=str,
    help="input aggregated Freyja file name",
    default="aggregated-freyja.tsv",
)
parser.add_argument(
    "-o", "--out", type=str, help="final file name", default="Freyja_aggregated_summary"
)
parser.add_argument(
    "-t", "--type", type=str, help="file extension for output table", default="csv"
)
parser.add_argument(
    "-v",
    "--version",
    help="print version and exit",
    action="version",
    version="%(prog)s " + "v0.9.0",
)
args = parser.parse_args()

"""input file"""
aggregated_file = args.input
"""output file"""
outputFile = args.out + "." + args.type


"""Simplify each abundance value to 4 decimal places and remove leading zero"""
def abundance_simplify(abundance):
    if not isinstance(abundance, str):  # Check if abundance is already a string
        abundance = str(abundance)  # Convert abundance to string if it's not
    abundance = re.sub(r"1\.00(\d)+\D*", r"1.00", abundance)
    abundance = re.sub(r"0\.(\d)(\d)(\d)(\d)(\d)+\D*", r".\1\2\3\4", abundance)
    abundance = float(abundance)  # Convert abundance to float (f64)
    return abundance

"""Get simplified lineage names"""
def perform_name_replacements(lineage):
    patterns = [
        (r"^\[\('", ""),  # Pattern 1: Remove beginning "(\['"
        (r"^',", ""),  # Pattern 2: Remove beginning "'"
        (r" \[.*\)\]'", ""),  # Pattern 3: Remove "[ ... \)]'"
        (r" ..*$", ""),  # Pattern 4: Remove everything after the space
        (r"^'", ""),  # Pattern 5: Remove beginning "'"
        (r"'$", ""),  # Pattern 6: Remove ending "'"
    ]
    """Apply each pattern sequentially"""
    for pattern, replacement in patterns:
        lineage = re.sub(pattern, replacement, lineage)
    return lineage

"""Simplify coverage to 2 decimal places"""
def coverage_simplify(coverage):
    coverage = float(coverage)
    coverage = round(coverage, 2)
    return coverage

"""Process the aggregated Freyja file"""
def process_aggregated_freya_file(aggregated_file):
    print("Opening ", aggregated_file, "\n")  # input filename

    """Read aggregated Freyja file"""
    with open(aggregated_file, "r") as input_file:
        lines = input_file.readlines()

        """Initialize list and hash of lineages"""
        lineageKeys = new_lineageKeys = {}
        strainHash = new_sH = {}

    """Skip header line and parse each line"""
    for line in lines[1:]:
        sample, summarized, lineages, abundances, resid, coverage = line.split("\t")
        """Parse sample, coverage and summarized data"""
        sample = re.sub(r"_variants.tsv", "", sample)

        coverage = coverage_simplify(coverage)
        """Create sample_covg tuple of sample and coverage"""
        sample_covg = sample, coverage

        """Parse summarized data"""
        lin_abundances = summarized.split("), (")

        """Separate lineage and abundance from within lin_abundances for simplifying"""
        new_lin_abundances = []
        for i in lin_abundances:
            (lineage, abundance) = i.split(", ")  # get tuple from each i
            lineage = perform_name_replacements(lineage)
            abundance = float(abundance_simplify(abundance))
            linHash = {lineage: abundance}

            """Append linHash to new_lin_abundances, rather than replacing it"""
            new_lin_abundances.append(linHash)

            """Add to lineageKeys hash"""
            new_lineageKeys = {lineage: 1}
            lineageKeys.update(new_lineageKeys)

        """Create sorted list of lineages"""
        sortedLineages = sorted(list(lineageKeys.keys()))

        """Create new_sH row and add to strainHash"""
        new_sH = {sample_covg: new_lin_abundances}
        strainHash.update(new_sH)
    return (strainHash, sortedLineages)

def save_df(strainHash, sortedLineages):
    covg_pct = "Coverage (%)"
    print("Sorted lineage list: ", sortedLineages)
    print("Number of lineages: " + str(len(sortedLineages)))

    """Initialize an empty list to hold the rows"""
    rows = []

    """Flatten the nested data before populating the rows list"""
    for sample_covg, lin_abundances in strainHash.items():
        sample, covg = sample_covg  # Unpack the tuple

        thisRow = {"Sample": sample, covg_pct: covg}
        for abundance_dict in lin_abundances:
            for lineage, abundance in abundance_dict.items():
                thisRow[lineage] = abundance_simplify(abundance)
        rows.append(thisRow)

    """Create a DataFrame from the flattened rows"""
    df = pl.DataFrame(rows)

    """Sort the DataFrame by 'Sample'"""
    df = df.sort("Sample")

    """Get the abundance columns excluding 'Sample' and 'Coverage (%)'"""
    abundance_columns = [
        col for col in df.columns if col not in ["Sample", covg_pct]
    ]

    """Sort the abundance columns alphabetically by lineage name"""
    sorted_abundance_columns = sorted(abundance_columns, key=lambda x: x.lower())

    """Select columns in the desired order"""
    df = df.select(["Sample", covg_pct] + sorted_abundance_columns)

    return df

if __name__ == "__main__":
    """Process aggregated Freyja file"""
    (strainHash, sortedLineages) = process_aggregated_freya_file(aggregated_file)
    print("\n\nsortedLineages: " + str(sortedLineages))

    """Process strainHash - create new version of StrainHash new_strainHash, 
    with the header row containing sample, coverage, and all lineageKeys, and each data row 
    containing the abundance for each lineage, including null values"""
    df = save_df(strainHash, sortedLineages)

    """Write modified dataframe as csv"""
    df.write_csv(outputFile)
    print("Results written to ", outputFile)

    """Print dataframe to screen"""
    with pl.Config(tbl_cols=-1, tbl_rows=-1):
        print(df)
        
