import argparse
import pandas as pd

def extract_ko_numbers(proteomics_file: str, ko_annotation_file: str, output_file: str):
    """
    Extract KEGG Orthology (KO) numbers from a proteomics file
    and annotate them using a KO definition table.

    Parameters
    ----------
    proteomics_file : str
        Path to the proteomics Excel file.
    ko_annotation_file : str
        Path to the .txt file containing KO annotations.
    output_file : str
        Path to the output CSV file.
    """
    # 1. Read proteomics.xlsx
    df_prot = pd.read_excel(proteomics_file)
    
    # 2. Read KO annotation file
    ko_df = pd.read_csv(
        ko_annotation_file, 
        sep="\t", 
        header=None, 
        names=["proteinID", "KO", "description"], 
        usecols=[0, 1, 2]
    )

    # 3. Merge both tables based on proteinID
    merged = df_prot.merge(ko_df, on="proteinID", how="left")

    # 4. Save the result
    merged.to_csv(output_file, index=False)
    print(f"✅ File saved as {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract KO numbers from proteomics data")
    parser.add_argument("--proteomics_file", required=True, help="Input proteomics file (.xlsx)")
    parser.add_argument("--ko_annotation_file", required=True, help="KO annotation file (.txt)")
    parser.add_argument("--output_file", required=True, help="Output CSV file")

    args = parser.parse_args()

    extract_ko_numbers(args.proteomics_file, args.ko_annotation_file, args.output_file)
