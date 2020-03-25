import sys
import pandas as pd
from Bio import Entrez
import argparse
import time
import csv
from pandas import read_excel
from Bio import SeqIO
import ftputil
from ftplib import FTP

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
    # Get first file name
    file_name = host.nlst(path)[0].strip().split('/')[-1]
    nucl_id = file_name.split('.')[0]
    return nucl_id

def get_biosample_id(old_id):
    """Get nnew RefSeq ID from an old RefSeq ID"""
    from Bio import Entrez
    # provide your own mail here
    Entrez.email = # email
    nucl_id = get_old_id(old_id)
    handle = Entrez.efetch(db='nucleotide', id=nucl_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    ref_id = record.dbxrefs
    biosample_id = str([x for x in ref_id if x.startswith('BioSample:')]).strip("['BioSample']").replace(":", "")
    if biosample_id == '':
        biosample_id = str([x for x in ref_id if x.startswith('BioProject')]).strip("['BioProject']").replace(":", "")
    return biosample_id

# ---------------------------------------------------------------
# Load excel file sheet
# ---------------------------------------------------------------
sheet_df = read_excel(input_file, sheet_name = my_sheet)
print(sheet_df.head()) # shows headers with top 5 rows
old_refseq = sheet_df['Genome Name / Sample Name'].drop_duplicates().tolist()
# ---------------------------------------------------------------
# Get GenBank ID - NCBI
# ---------------------------------------------------------------
biosample_dic = {}

# Access to ftp NCBI
#host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
#host = FTP('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')

for idx, old_id in enumerate(old_refseq):
    try:
        host = FTP('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
        biosample_id = get_biosample_id(old_id)
        biosample_dic[old_id] = biosample_id
    except:
        print("Old Refseq ID not found: " + old_id)
        biosample_dic[old_id] = ('-')
    #print(old_id)
    print(str(idx + 1) + "/" + str(len(old_refseq)))
    time.sleep(4)  # Delays for 4 seconds

#for k in not_found:
for k, v in biosample_dic.items():
    if v == 'nan':
        print('Not found: ', k)
        try:
            host = FTP('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
            biosample_id = get_biosample_id(k)
            biosample_dic[k] = biosample_id
        except:
            print("Old Refseq ID not found: " + k)
            biosample_dic[k] = ('-')
    else:
        #print(k)
        biosample_dic[k] = v
        continue
    time.sleep(4)

# Save GenBank IDs in a Df
entrez_df = pd.DataFrame(biosample_dic.items(), columns=['Genome Name / Sample Name','Biosample Accession'])
# Replace 'Biosample Accession' with updated IDs
sheet_df.set_index('Genome Name / Sample Name', inplace=True)
sheet_df.update(entrez_df.set_index('Genome Name / Sample Name'))
sheet_df.reset_index(inplace=True)
sheet_df.to_csv(output_file, sep="\t", quoting=csv.QUOTE_NONE, index=False, header=True)
print("Done!")