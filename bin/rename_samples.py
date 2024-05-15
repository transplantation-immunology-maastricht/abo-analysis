#!/usr/bin/env python3

import pandas as pd
import sys
import re


def main(final_export_file, deobfuscation):
    """
    This function reads in a file exported from Soft and the ABO pipeline final_export file and 
    replaces SequencingAcc# with Patient# to allow export into MatchPoint
    """
    # regex pattern to capture everything before "_barcode" followed by digits
    pattern = r'(.+)_barcode\d+$'


    # Excel file from soft PCR HLA exported file
    renaming_file = pd.read_excel(deobfuscation, index_col=None, na_values=['NA'], usecols="A,C")
    renaming_file = renaming_file.rename(columns={'Acc#': 'Sample ID', 'Patient Name': 'Grid_number'})
    renaming_file['Grid_number'] = renaming_file['Grid_number'].str.split(',').str[0]
    renaming_file['Grid_number'] = renaming_file['Grid_number'].astype(str)
    renaming_file_filtered = renaming_file[~renaming_file['Grid_number'].str.match('^[a-zA-Z]')]

    # Check if "Sample ID" contains "_barcode" and split if it does

    # Python complains when same file is modified in-place:
    # .. Error message.. A value is trying to be set on a copy of a slice from a DataFrame.
    # Try using .loc[row_indexer,col_indexer] = value instead. 
    # The Solution ids to create a copy of the file before modification

    renaming_file_filtered = renaming_file_filtered.copy()
    renaming_file_filtered.loc[:, "Sample ID"] = renaming_file_filtered["Sample ID"].apply(
        lambda x: re.match(pattern, x).group(1) if re.match(pattern, x) else x
    )
        
    final_export = pd.read_csv(final_export_file, sep=",")

    # Similar solution to line 24-30 above
    final_export = final_export.copy()

    final_export.loc[:, "Sample ID"] = final_export["Sample ID"].apply(
        lambda x: re.match(pattern, x).group(1) if re.match(pattern, x) else x
    )

    # Left join final_export with samples using "Sample ID" as the key
    final_export_grid = pd.merge(final_export, renaming_file_filtered, how="left")
    col_order = ['Grid_number'] + [col for col in final_export_grid.columns if col != 'Grid_number']
    final_export_grid = final_export_grid[col_order]
    final_export_grid['Grid_number'] = final_export_grid['Grid_number'].astype(str)

    # Rename 'Sample ID' to 'SequencingAcc#'
    final_export_grid = final_export_grid.rename(columns={'Sample ID': 'SequencingAcc#'})
    
    # Rename 'Grid_number' to 'Sample ID'
    final_export_grid = final_export_grid.rename(columns={'Grid_number': 'Sample ID'})

    # Create a copy with both sequencing Acc# and Patient ID #
    final_export_grid_noAccesion = final_export_grid.copy()

    # Drop sequencing Acc# column
    final_export_grid_noAccesion.drop('SequencingAcc#', axis=1, inplace=True)

    # Write to file
    final_export_grid.to_csv("MatchPointExport_with_sequencingAcc.txt", index=False)
    final_export_grid_noAccesion.to_csv("MatchPointExport.txt", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python rename_samples.py <final_export_file> <deobfuscation>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])


sys.exit(0)
