import subprocess
import argparse
import os
import shutil
import asyncio
import time
import logging
import re
import pandas as pd
from domainator.utils import parse_seqfiles
from pyppeteer import launch

parser = argparse.ArgumentParser(description='Downloads genomes from protein data')
parser.add_argument('-i', type=str, default=None, required=True, help='TSV table input')
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
        # Launch the browser
        browser = await launch(headless=True)
        page = await browser.newPage()
        await page.goto(url, {'waitUntil': 'networkidle2', 'timeout': 120000})

        # Wait for the page to load (you can adjust the selector)
        await page.waitForSelector('#ph-ipg > div.ipg-rs > table > tbody > tr > td:nth-child(7) > a',
                                   {'timeout': 120000})

        # Extract data (adjust the selector as necessary)
        data_elements = await page.querySelectorAll('#ph-ipg > div.ipg-rs > table > tbody > tr > td:nth-child(7) > a')
        for element in data_elements:
            text = await page.evaluate('(element) => element.textContent', element)
            genome_accessions.append(text)

        # Close the browser
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

df = pd.read_csv(args.i)[["Query Cover", "Per. ident", "Accession  "]]
df.rename(columns={"Query Cover": "qc", "Per. ident": "pi", "Accession  ": "acc"}, inplace=True)
df['acc'] = df['acc'].apply(lambda x: extract_accession(x))

for index, row in df.iterrows():
    if float(row["qc"][:-1]) >= 70 and float(row["pi"]) >= 70:
        genome_acc = asyncio.get_event_loop().run_until_complete(fetch_genome_data(row['acc']))
        if genome_acc:
            print(f"Protein Accession: {row['acc']}, Genome Accession: {', '.join(genome_acc)}")
            # Create directory for protein accession if it doesn't exist
            protein_dir = os.path.join(args.o, row['acc'])
            os.makedirs(protein_dir, exist_ok=True)

            for gen_ac in genome_acc:
                command = f"datasets download genome accession {gen_ac}"
                result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

                if result.returncode == 0:
                    print("Download successful.")
                    try:
                        os.rename("ncbi_dataset.zip", f"{gen_ac}.zip")
                        shutil.move(f"{gen_ac}.zip", protein_dir)
                    except FileNotFoundError:
                        print("Downloaded file not found. Make sure the download was successful.")
                    except Exception as e:
                        print(f"An error occurred while renaming or moving the file: {e}")
                else:
                    print(f"Error downloading genome accession {gen_ac}: {result.stderr}")
        else:
            print(f"Protein Accession: {row['acc']}, Genome Accession: Not found")

time.sleep(5)

text = open(f"{args.o}/undownloaded.txt", "w")
for item in undownloaded:
    text.write(item + "\n")
text.close()

# records = parse_seqfiles(args.i)
# for i, rec in enumerate(records):  # Each rec is a protein
#     protein_acc = rec.annotations["accessions"][0]
#     genome_acc = asyncio.get_event_loop().run_until_complete(fetch_genome_data(protein_acc))
#     if genome_acc:
#         print(f"Protein Accession: {protein_acc}, Genome Accession: {', '.join(genome_acc)}")
#
#         # Create directory for protein accession if it doesn't exist
#         protein_dir = os.path.join(args.o, protein_acc)
#         os.makedirs(protein_dir, exist_ok=True)
#
#         for gen_ac in genome_acc:
#             command = f"datasets download genome accession {gen_ac}"
#             result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#
#             if result.returncode == 0:
#                 print("Download successful.")
#                 try:
#                     os.rename("ncbi_dataset.zip", f"{gen_ac}.zip")
#                     shutil.move(f"{gen_ac}.zip", protein_dir)
#                 except FileNotFoundError:
#                     print("Downloaded file not found. Make sure the download was successful.")
#                 except Exception as e:
#                     print(f"An error occurred while renaming or moving the file: {e}")
#             else:
#                 print(f"Error downloading genome accession {gen_ac}: {result.stderr}")
#     else:
#         print(f"Protein Accession: {protein_acc}, Genome Accession: Not found")
#
#     time.sleep(5)

# text = open(f"{args.o}/undownloaded.txt", "w")
# for item in undownloaded:
#     text.write(item + "\n")
# text.close()
