import sys
import pandas as pd
from Bio import Entrez
import argparse
import time
import csv
from pandas import read_excel
from Bio import SeqIO
import ftputil

# python3 CF_isolates_get_BioSampleID.py -i /Bacterial\ Cultures\ Catalogue.xlsx -o CF_isolates_BioSampleID.csv
arguments = sys.argv

parser = argparse.ArgumentParser(description='Get RefSeq ID')
parser.add_argument("-i", "--input", required=True, help="Select input file")
parser.add_argument("-o", "--output", required=True, help="Direct the output to a file")
args = parser.parse_args()

input_file = args.input
output_file = args.output
my_sheet = 'CF isolates'

# ---------------------------------------------------------------
# Get GenBank ID - NCBI
# ---------------------------------------------------------------
def get_old_id(old_id):
    """Get Assembly ID from Accession ID"""
    path = '/genomes/archive/old_refseq/Bacteria/%s' % (old_id)
    # Access to ftp NCBI with directory name for organism from RefSeq Ftp path
    host.chdir(path)
    # Get list of the directory with files for organism: GCF_<id for organism>_<Assembly ID>(_genomic.fna.gz, _assembly_report.txt,....)
    dir_list = host.listdir(host.curdir)
    nucl_id = dir_list[0].split('.')[0]
    return nucl_id

def get_biosample_id(old_id):
    """Get nnew RefSeq ID from an old RefSeq ID"""
    from Bio import Entrez
    # provide your own mail here
    Entrez.email =  # email
    nucl_id = get_old_id(old_id)
    handle = Entrez.efetch(db='nucleotide', id=nucl_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    ref_id = record.dbxrefs
    biosample_id = str([x for x in ref_id if x.startswith('BioSample:')]).strip("['BioSample']").replace(":", "")
    return biosample_id

# ---------------------------------------------------------------
# Load excel file sheet
# ---------------------------------------------------------------
sheet_df = read_excel(input_file, sheet_name = my_sheet)
print(sheet_df.head()) # shows headers with top 5 rows
old_refseq = sheet_df['Genome Name / Sample Name'].tolist()

# ---------------------------------------------------------------
# Get GenBank ID - NCBI
# ---------------------------------------------------------------
biosample_dic = {}

# Access to ftp NCBI
host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')

for idx, old_id in enumerate(old_refseq):
    try:
        biosample_id = get_biosample_id(old_id)
        biosample_dic[old_id] = biosample_id
    except:
        print("Old Refseq ID not found: " + old_id)
        biosample_dic[old_id] = ('-')
    #print(old_id)
    print(str(idx + 1) + "/" + str(len(old_refseq)))
    time.sleep(4)  # Delays for 4 seconds

# Save GenBank IDs in a Df
entrez_df = pd.DataFrame(biosample_dic.items(), columns=['Genome Name / Sample Name','Biosample Accession'])
# Replace 'Biosample Accession' with updated IDs
sheet_df.set_index('Genome Name / Sample Name', inplace=True)
sheet_df.update(entrez_df.set_index('Genome Name / Sample Name'))
sheet_df.reset_index(inplace=True)
sheet_df.to_csv(output_file, sep="\t", quoting=csv.QUOTE_NONE, index=False, header=True)
print("Done!")


