import argparse
import os
import asyncio
import logging
import re
import pandas as pd
import tqdm
import gzip
import requests
import time
from pyppeteer import launch
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Downloads genomes from protein data')
parser.add_argument('-i', type=str, default=None, required=True, help='TSV table input')
parser.add_argument('--ncbi', type=str, default=None, required=True, help='NCBI assembly summary file')
#parser.add_argument('-lbs', nargs='+', type=int, default=None, required=True, help='Lower bounds. Cover, Identity.')
parser.add_argument('-o', type=str, default=None, required=True, help='Output folder name')
args = parser.parse_args()

undownloaded = []

logging.basicConfig(level=logging.INFO)
logging.getLogger('pyppeteer').setLevel(logging.WARNING)

async def fetch_genome_data(protein_acc):
    url = f"https://www.ncbi.nlm.nih.gov/ipg/?term={protein_acc}"
    genome_accessions = []
    try:
        # Launch browser
        browser = await launch(headless=True)
        page = await browser.newPage()
        await page.goto(url, {'waitUntil': 'networkidle2', 'timeout': 120000})

        # Wait for page to load
        await page.waitForSelector('#ph-ipg > div.ipg-rs > table > tbody > tr > td:nth-child(7) > a',{'timeout': 120000})

        # Extract from the td element
        data_elements = await page.querySelectorAll('#ph-ipg > div.ipg-rs > table > tbody > tr > td:nth-child(7) > a')
        for element in data_elements:
            text = await page.evaluate('(element) => element.textContent', element)
            genome_accessions.append(text)

        # Close browser
        await browser.close()

    except asyncio.TimeoutError:
        logging.error(f"Timeout while processing URL: {url}")
        undownloaded.append(protein_acc)
    except Exception as e:
        logging.error(f"An error occurred while processing URL: {url}, error: {e}")
        undownloaded.append(protein_acc)
    finally:
        if 'browser' in locals():
            await browser.close()

    return genome_accessions

def extract_accession(link):
    match = re.search(r'//www.ncbi.nlm.nih.gov/protein/([a-zA-Z0-9_.]+)\?report=genbank', link)
    return match.group(1) if match else link

print("Loading protein data...")
df = pd.read_csv(args.i)[["Query Cover", "Per. ident", "Accession  "]]
df.rename(columns={"Query Cover": "qc", "Per. ident": "pi", "Accession  ": "acc"}, inplace=True)
df['acc'] = df['acc'].apply(lambda x: extract_accession(x))

print("Loading NCBI data...")
ncbi = pd.read_csv(args.ncbi, sep='\t', header=1)[['#assembly_accession', 'gbrs_paired_asm', 'ftp_path']]

for index, row in df.iterrows():
    if float(row["qc"][:-1]) >= 70 and float(row["pi"]) >= 70:
        genome_acc = asyncio.get_event_loop().run_until_complete(fetch_genome_data(row['acc']))
        if genome_acc:
            print(f"Protein Accession: {row['acc']}, Genome Accession: {', '.join(genome_acc)}")
            # Create directory for protein accession if it doesn't exist
            protein_dir = os.path.join(args.o, row['acc'])
            os.makedirs(protein_dir, exist_ok=True)

            for gen_ac in genome_acc:
                url = ""
                if gen_ac.startswith("GCA_"):
                    if gen_ac not in ncbi['#assembly_accession'].values:
                        print(f"Genome Accession: {gen_ac}, FTP path not found")
                        continue
                    url = ncbi.loc[ncbi['#assembly_accession'] == gen_ac, 'ftp_path'].values[0]
                elif gen_ac.startswith("GCF_"):
                    if gen_ac not in ncbi['gbrs_paired_asm'].values:
                        print(f"Genome Accession: {gen_ac}, FTP path not found")
                        continue
                    url = ncbi.loc[ncbi['gbrs_paired_asm'] == gen_ac, 'ftp_path'].values[0]
                else:
                    print("Unknown accession type")
                    continue
                acc_dir = url[url.rfind('/') + 1:]
                url = f"{url}/{acc_dir}_genomic.gbff.gz"

                try:
                    response = requests.get(url, stream=True)
                    response.raise_for_status()

                    gz_file_path = os.path.join(protein_dir, f"{acc_dir}_genomic.gbff.gz")
                    gb_file_path = os.path.join(protein_dir, f"{acc_dir}_genomic.gb")

                    with open(gz_file_path, "wb") as output_handle:
                        for chunk in response.iter_content(chunk_size=8192):
                            output_handle.write(chunk)

                    # Unzip and change extension to .gb
                    with gzip.open(gz_file_path, "rt") as gzipped_file:
                        with open(gb_file_path, "w") as gb_file:
                            records = tqdm.tqdm(SeqIO.parse(gzipped_file, "genbank"), desc=acc_dir, leave=True, dynamic_ncols=True)
                            recs_written = SeqIO.write(records, gb_file, "genbank")

                    # Remove gzipped file after unzipping
                    os.remove(gz_file_path)

                    print(f"Downloaded and converted {recs_written} records for {row['acc']}")

                    time.sleep(5)

                except requests.HTTPError as e:
                    print(f"HTTP error occurred: {e}")
                except Exception as e:
                    print(f"An error occurred: {e}")
        else:
            print(f"Protein Accession: {row['acc']}, Genome Accession: Not found")


text = open(f"{args.o}/undownloaded.txt", "w")
for item in undownloaded:
    text.write(item + "\n")
text.close()
