#!/usr/bin/python3

###########################################
# Freyja Summarizer v0.8.3                #
# Written by Stephen Beckstrom-Sternberg  #
#  Assisted by Khalil Khoury              #
#  Creates a simplified Freyja summary    #
#   with lineage abundances & coverage    #
#   from aggregated Freyja file           #
###########################################

import re
import pandas as pd
import numpy as np

"""input file"""
aggregated_file = "aggregated-freyja.tsv"
"""output file"""
outputFile = "Freyja_aggregated_summary.tsv"

"""Simplify each abundance value to 4 decimal places and remove leading zero"""
def abundance_simplify(abundance):
    abundance = re.sub(r"1\.00(\d)+\D*", r"1.00", abundance)
    abundance = re.sub(r"0\.(\d)(\d)(\d)(\d)(\d)+\D*", r".\1\2\3\4", abundance)
#    abundance = re.sub(r"0(\d)\.(\d)(\d)(\d)+", r"\1.\2\3%", abundance)
#    abundance = re.sub(r'\)\]', r"", abundance)
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
#        print("Number of lin_abundances: ", len(lin_abundances))

        """Separate lineage and abundance from within lin_abundances for simplifying"""
        new_lin_abundances = []
        for i in lin_abundances:
            (lineage, abundance) = i.split(", ")  # get tuple from each i
            lineage = perform_name_replacements(lineage)
            abundance = abundance_simplify(abundance)

            linHash = {lineage: abundance}
            
            """Append linHash to new_lin_abundances, rather than replacing it"""
            new_lin_abundances.append(linHash)

            """Add to lineageKeys hash"""
            new_lineageKeys = {lineage: 1}
            lineageKeys.update(new_lineageKeys)
        
        """Create sorted list of lineages"""
        sortedLineages = sorted(list(lineageKeys.keys()))

        """Create new_sH row and add to strainHash"""
        new_sH = {sample_covg:new_lin_abundances}
    
#        print("new_sH: " +str(new_sH))
        strainHash.update(new_sH)

#    return (strainHash, sortedLineages)
    return (strainHash, sortedLineages)


def save_df(strainHash,sortedLineages):
    new_strainHash = {}
    ##NEW STUFF
    
    
    lineageHeader = ["Lineages" for _ in range(len(sortedLineages))]
    print("LINEAGEHEADER: " +str(lineageHeader))
#    lineageHeader = ["Lineages"]*len(sortedLineages)
#    lineageTuple = list(zip(lineageHeader, sortedLineages))

#    headerList = ["Sample","Coverage (%)"] + lineageTuple
    headerList = ["Sample","Coverage (%)"] + sortedLineages
    print("Sorted lineage list: ", sortedLineages)
    print("Number of lineages: " +str(len(sortedLineages)))
    print("\ndataframe headerList: " +str(headerList))
    
    """Create empty dataframe with column headers"""
    df  = pd.DataFrame(columns = headerList)
#    print("\ndf after creating empty dataframe with named column headers:\n", df)
    """Populate dataframe with rows from strainHash"""
    for sample_covg, lin_abundances in strainHash.items():
        sample, covg = sample_covg  # Unpack the tuple
    
        thisRow = {"Sample": sample, "Coverage (%)": covg}
        for abundance_dict in lin_abundances:
            for lineage, abundance in abundance_dict.items():
                thisRow[lineage] = abundance_simplify(abundance)    
        df = df.append(thisRow, ignore_index=True)  # add row to DataFrame
    df.replace(np.nan, 0, inplace=True)
    df.sort_values(by=["Sample"], inplace=True)
#    print("\ndf after appending data rows to the DataFrame:\n", df)
    return(new_strainHash, df)

if __name__ == "__main__":

    """Process aggregated Freyja file"""
    (strainHash,sortedLineages) = process_aggregated_freya_file(aggregated_file)
    
#    print("\n\nstrainHash after processing: " +str(strainHash))
    print("\n\nsortedLineages: " +str(sortedLineages))
    
    """Process strainHash - create new version of StrainHash new_strainHash, 
    with the header row containing sample, coverage, and all lineageKeys, and each data row 
    containing the abundance for each lineage, including null values"""

    """Create header, prepending Sample and Coverage to the sorted lineages"""
    
    (new_strainHash, df) = save_df(strainHash, sortedLineages)
    
    """Write out the modified dataframe as tsv"""
    open(outputFile,'w')
    df.to_csv(outputFile, sep="\t", index=False)
    print('Results written to ', outputFile)