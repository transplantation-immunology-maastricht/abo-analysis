#!/usr/bin/env python3
import os
import re
import pandas as pd
import openpyxl
from openpyxl import Workbook
from openpyxl import load_workbook

# Define input and output directories
input_dir = "/data/abo_genotypes/abo-analysis/output"

def parse_exon7(filename):    
    # Open the file for reading
    with open(filename, 'r') as f:
        lines = f.readlines()

        # Extract the required lines from the first section
        pos1 = int(lines[2].split(':')[1].strip())
        nuc1 = lines[3].split(':')[1].strip()
        blood1 = lines[4].split(':')[1].strip()
        count1 = int(lines[6].split(':')[1].strip())
        mat1, mis1, ins1, del1, a1, g1, c1, t1 = [int(x) for x in lines[8].split()]

        # Extract the required lines from the second section
        pos2 = int(lines[10].split(':')[1].strip())
        nuc2 = lines[11].split(':')[1].strip()
        blood2 = lines[12].split(':')[1].strip()
        count2 = int(lines[14].split(':')[1].strip())
        mat2, mis2, ins2, del2, a2, g2, c2, t2 = [int(x) for x in lines[16].split()]

    # Create a DataFrame from the extracted values
    df = pd.DataFrame({
        'Exon': ['7', '7'],
        'Position': [pos1, pos2],
        '#Reads': [count1, count2],
        'Mat': [mat1, mat2],
        'Mis': [mis1, mis2],
        'Ins': [ins1, ins2],
        'Del': [del1, del2],
        'A': [a1, a2],
        'G': [g1, g2],
        'C': [c1, c2],
        'T': [t1, t2]
    })

    # Remove the string "Exon 7 pos" from the Position column
    df['Position'] = df['Position'].apply(lambda x: x)

    # Add the 'Type' column based on values in the [Ins Del A G C T] columns
    def get_type(row):
        if row['Position'] == 422:
            if row['A'] > 80:
                return 'B'
            elif row['C'] > 80:
                return 'A or O'
            elif abs(row['A'] - row['C']) <= 20:
                return '(A or O) and B'
        elif row['Position'] == 429:
            if row['G'] > 80:
                return 'A or O'
            elif row['C'] > 80:
                return 'B'
            elif abs(row['G'] - row['C']) <= 20:
                return '(A or O) and B'
        return ''
    
    df['Type'] = df.apply(get_type, axis=1)

    # Reorder the columns
    df = df[['Exon', 'Position', '#Reads', 'Mat', 'Mis', 'Ins', 'Del', 'A', 'G', 'C', 'T', 'Type']]
    
    return df


def parse_exon6(filename):    
    # Open the file for reading
    with open(filename, 'r') as f:
        lines = f.readlines()

        # Extract the required lines from the first section
        pos1 = int(lines[2].split(':')[1].strip())
        count1 = int(lines[6].split(':')[1].strip())
        mat1, mis1, ins1, del1, a1, g1, c1, t1 = [int(x) for x in lines[8].split()]

    # Create a DataFrame from the extracted values
    df = pd.DataFrame({
        'Exon': ['6'],
        'Position': [pos1],
        '#Reads': [count1],
        'Mat': [mat1],
        'Mis': [mis1],
        'Ins': [ins1],
        'Del': [del1],
        'A': [a1],
        'G': [g1],
        'C': [c1],
        'T': [t1]
    })

    # Remove the string "Exon 6 pos" from the Position column
    df['Position'] = df['Position'].apply(lambda x: x)

    # Add the 'Type' column
    max_g = df['G'].max()
    max_del = df['Del'].max()
    if max_g > 80 and max_g > max_del:
        df['Type'] = 'A or B'
    elif max_del > 80 and max_del > max_g:
        df['Type'] = 'O'
    elif abs(max_g - max_del) <= 20:
        df['Type'] = 'O and (A or B)'
    else:
        df['Type'] = ''

    # Reorder the columns
    df = df[['Exon', 'Position', '#Reads', 'Mat', 'Mis', 'Ins', 'Del', 'A', 'G', 'C', 'T', 'Type']]

    return df

# Add Pheno, Geno, Expecetd columns to dataframe.

def add_phenotype_genotype(df):
    # Find the column numbers with "Type" in the column name
    type_cols = [i for i, col in enumerate(df.columns) if 'Type' in col]
    
    # Assign the column numbers to Type1, Type2, and Type3
    Type1 = df.columns[type_cols[0]]
    Type2 = df.columns[type_cols[1]]
    Type3 = df.columns[type_cols[2]]
    
    print(Type1, Type2, Type3)
    # Add the Phenotype, Genotype, and Expected columns
    df['Phenotype'] = ''
    df['Genotype'] = ''
    df['Expected'] = ''

    for i in range(len(df)):
        if ((df.loc[i, Type1] == 'A or B').all() and (df.loc[i, Type2] == 'A or O').all() and (df.loc[i, Type3] == 'A or O')).any():
            df.at[i, 'Phenotype'] = 'A'
            df.at[i, 'Genotype'] = 'AA'
            df.at[i, 'Expected'] = 'AA'
        # elif ((df.loc[i, Type1] == 'A or B') and (df.loc[i, Type2] == 'B') and (df.loc[i, Type3] == 'B')).empty().any().all().item():
        #     df.at[i, 'Phenotype'] = 'B'
        #     df.at[i, 'Genotype'] = 'BB'
        #     df.at[i, 'Expected'] = 'BB'
        # elif (df.loc[i, Type1] == 'A or B') and (df.loc[i, Type2] == '(A or O) and B') and (df.loc[i, Type3] == '(A or O) and B'):
        #     df.at[i, 'Phenotype'] = 'AB'
        #     df.at[i, 'Genotype'] = 'AB'
        #     df.at[i, 'Expected'] = 'AB'
        # elif (df.loc[i, Type1] == 'O and (A or B)') and (df.loc[i, Type2] == 'A or O') and (df.loc[i, Type3] == 'A or O'):
        #     df.at[i, 'Phenotype'] = 'A'
        #     df.at[i, 'Genotype'] = 'AO'
        #     df.at[i, 'Expected'] = 'AO'
        # elif (df.loc[i, Type1] == 'O and (A or B)') and (df.loc[i, Type2] == '(A or O) and B') and (df.loc[i, Type3] == '(A or O) and B'):
        #     df.at[i, 'Phenotype'] = 'B'
        #     df.at[i, 'Genotype'] = 'BO'
        #     df.at[i, 'Expected'] = 'BO'
        # elif (df.loc[i, Type1] == 'O') and (df.loc[i, Type2] == 'A or O') and (df.loc[i, Type3] == 'A or O'):
        #     df.at[i, 'Phenotype'] = 'O'
        #     df.at[i, 'Genotype'] = 'OO'
        #     df.at[i, 'Expected'] = 'OO'
        # elif (df.loc[i, Type1] == 'O and (A or B)') and (df.loc[i, Type2] == 'A or O') and (df.loc[i, Type3] == 'Unknown'):
        #     df.at[i, 'Phenotype'] = 'A'
        #     df.at[i, 'Genotype'] = 'AO'
        #     df.at[i, 'Expected'] = 'Unknown'
        # elif ((df.loc[i, Type1] == 'O and (A or B)') and (df.loc[i, Type2] == 'A or O') and (df.loc[i, Type3] == 'Unknown')).all():
        #     df.at[i, 'Phenotype'] = 'A'
        #     df.at[i, 'Genotype'] = 'AO'
        #     df.at[i, 'Expected'] = 'Unknown'
        else:
            df.at[i, 'Phenotype'] = 'Indeterminate'
            df.at[i, 'Genotype'] = 'Indeterminate'
            df.at[i, 'Expected'] = 'Indeterminate'
    return df

# create an empty list to store the dataframes
results = []

# Define the column headers
exon6 = ['Exon 6 pos 22'] * 10
exon7_422 = ['Exon 7 pos 422'] * 10
exon7_429 = ['Exon 7 pos 429'] * 10
max_len = max(len(exon6), len(exon7_422), len(exon7_429))
exon6 += [''] * (max_len - len(exon6))
exon7_422 += [''] * (max_len - len(exon7_422))
exon7_429 += [''] * (max_len - len(exon7_429))
header_cols = ['', ''] + exon6 + exon7_422 + exon7_429 + ['', '', '']
header_rows = ['Barcode','Sample'] + ['#Reads', 'Mat', 'Mis', 'Ins', 'Del', 'A', 'G', 'C', 'T', 'Type'] * 3 + ['Phenotype','Genotype', 'Expected']
columns = pd.MultiIndex.from_arrays([header_cols, header_rows])

# Loop through all files in the input directory
for filename in os.listdir(input_dir):
    # Check if the file is a sample folder (for HSS Perth, these files start with IMM[0-9] with _barcode[0-9])
    import re

    if os.path.isdir(os.path.join(input_dir, filename)):
        match = re.match(r"^IMM-[0-9]+-[0-9]+_barcode\d+", filename)
        if match:
            print("Processing file: " + filename)
        # process the file
        # if os.path.isdir(os.path.join(input_dir, filename)) and filename.startswith("IMM") and "barcode" in filename and "POS" not in filename:
            # Get the sample name and barcode from the folder name
            sample_name, barcode = filename.split("_")
            # print("\nProcessing Sample %s with barcode %s" %(barcode, sample_name))

            # Process exon 6 and 7 data using the files
            exon6_file = os.path.join(input_dir, filename, "exon6", "ABOPhenotype.txt")
            exon7_file = os.path.join(input_dir, filename, "exon7", "ABOPhenotype.txt")
            
            df_exon6 = parse_exon6(exon6_file)
            # print(df_exon6)
            df_exon7 = parse_exon7(exon7_file)
            # print(df_exon7)

            # create two subset dataframes, one for each row of df for exon7
            df_exon7_pos422 = df_exon7.iloc[[0]].reset_index(drop = True)
            df_exon7_pos429 = df_exon7.iloc[[1]].reset_index(drop = True)

            print(df_exon7_pos422)
            print(df_exon7_pos429)

            sample = pd.DataFrame({
                'Barcode': [barcode],
                'Sample': [sample_name]
                })

            # concatenate the dataframes side by side

            merged_df = pd.concat([sample, df_exon6, df_exon7_pos422, df_exon7_pos429], 
                                axis = 1, join = "inner").drop(['Exon', 'Position'], axis = 1)
            

            merged_df['Barcode'] = merged_df['Barcode'].str.replace('barcode', '', case = False)
            # Convert the Barcode column to numeric
            merged_df['Barcode'] = pd.to_numeric(merged_df['Barcode'], errors = 'coerce')

            merged_df = add_phenotype_genotype(merged_df)

            # print(merged_df)

            # # set the column headers to the MultiIndex
            merged_df.columns = columns

            # print(merged_df.shape)
            # print(columns.shape)
        
            
            

            # append the merged_df to the results list
            results.append(merged_df)
            print("Done adding Sample %s with barcode %s to merged data frame\n" %(barcode, sample_name))


# concatenate all the dataframes in the results list into one dataframe
final_df = pd.concat(results)

print(final_df.to_string(index = False), "\n")

# # Write the concatenated dataframe to a text file
final_df.to_csv(os.path.join(input_dir,'result.txt'), sep = '\t', index = False)

# # write the final dataframe to an excel file
final_df.to_excel(os.path.join(input_dir,'result.xlsx'), sheet_name = "ABO_Result", index = True)

# # Write the DataFrame to an Excel file with merged cells and no row_num.. 
# with pd.ExcelWriter(os.path.join(input_dir,'result.xlsx'), engine = 'openpyxl') as writer:
#     final_df.to_excel(writer, sheet_name = 'ABO_Result')
#     workbook  = writer.book
#     worksheet = writer.sheets['ABO_Result']
#     for merged in worksheet.merged_cells:
#         worksheet.unmerge_cells(str(merged))
#         worksheet.merge_cells(str(merged))
#     writer.save()

