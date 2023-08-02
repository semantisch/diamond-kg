# This scripts takes in medical labels and turns them into context-enriched Knowledge graphs
# Using ChatGPT and a NER model

import datetime
import sys
import os
from pathlib import Path
import random
from stem import Signal
from stem.control import Controller
import pandas as pd
import requests
import json
from rich.console import Console
from rich.prompt import Prompt
from rich.style import Style
from rich import inspect
import time
from rich.progress import track
import hashlib
import requests
import concurrent
from concurrent.futures import ThreadPoolExecutor
import psutil

import asyncio

from time import sleep
from random import random
from concurrent.futures import as_completed

import unicodedata
import re

import psutil
import json
import os
from pprint import pprint
import requests
import time
import datetime
import os
from urllib.parse import urlparse
import csv
import base64
import random

import queue
import threading

import os
import openai
import subprocess

import argparse, sys # Pass arguments

import urllib3
import rdflib
from rdflib import Graph
from rdflib.namespace import RDF, FOAF, RDFS, OWL
from rdflib.term import URIRef, Literal, BNode
import warnings

import networkx as nx
from rdflib import Graph, Namespace
from rdflib.namespace import VOID
from rdflib.namespace import RDF
from rdflib.namespace import RDFS
from rdflib.namespace import OWL
from rdflib.namespace import XSD
from rdflib.namespace import DCTERMS
from rdflib.namespace import FOAF
from rdflib.namespace import QB
from rdflib import URIRef, BNode, Literal
import uuid
from bs4 import BeautifulSoup
import openai

import hashlib

import json

# import matplotlib.pyplot as plt
# import pandas as pd
#
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# sns.set()
# import spacy
# import pyLDAvis.gensim_models
# # pyLDAvis.enable_notebook()# Visualise inside a notebook
# import en_core_web_md
# from gensim.corpora.dictionary import Dictionary
# from gensim.models import LdaMulticore
# from gensim.models import CoherenceModel

import logging

## NER:

import os
from typing import Optional

# import numpy as np
# import requests
# import spacy
# import torch
# from fastapi import APIRouter, Body, HTTPException
# from fastapi.responses import JSONResponse
# from pydantic import BaseModel
# from transformers import BertForSequenceClassification, BertTokenizer

biolink_context = {
    "APO": "http://purl.obolibrary.org/obo/APO_",
    "AspGD": "http://www.aspergillusgenome.org/cgi-bin/locus.pl?dbid=",
    "BFO": "http://purl.obolibrary.org/obo/BFO_",
    "BIGG.METABOLITE": "http://identifiers.org/bigg.metabolite/",
    "BIGG.REACTION": "http://identifiers.org/bigg.reaction/",
    "BIOGRID": "http://identifiers.org/biogrid/",
    "BIOSAMPLE": "http://identifiers.org/biosample/",
    "BSPO": "http://purl.obolibrary.org/obo/BSPO_",
    "BTO": "http://purl.obolibrary.org/obo/BTO_",
    "CAID": "http://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=",
    "CAS": "http://identifiers.org/cas/",
    "CATH": "http://identifiers.org/cath/",
    "CATH.SUPERFAMILY": "http://identifiers.org/cath.superfamily/",
    "CDD": "http://identifiers.org/cdd/",
    "CHADO": "http://gmod.org/wiki/Chado/",
    "CHEBI": "http://purl.obolibrary.org/obo/CHEBI_",
    "CHEMBL.COMPOUND": "http://identifiers.org/chembl.compound/",
    "CHEMBL.MECHANISM": "https://www.ebi.ac.uk/chembl/mechanism/inspect/",
    "CHEMBL.TARGET": "http://identifiers.org/chembl.target/",
    "CID": "http://pubchem.ncbi.nlm.nih.gov/compound/",
    "CIO": "http://purl.obolibrary.org/obo/CIO_",
    "CL": "http://purl.obolibrary.org/obo/CL_",
    "CLINVAR": "http://identifiers.org/clinvar",
    "CLO": "http://purl.obolibrary.org/obo/CLO_",
    "COAR_RESOURCE": "http://purl.org/coar/resource_type/",
    "COG": "https://www.ncbi.nlm.nih.gov/research/cog-project/",
    "CPT": "https://www.ama-assn.org/practice-management/cpt/",
    "CTD": "http://ctdbase.org/",
    "CTD.CHEMICAL": "http://ctdbase.org/detail.go?type=chem&acc=",
    "CTD.DISEASE": "http://ctdbase.org/detail.go?type=disease&db=MESH&acc=",
    "CTD.GENE": "http://ctdbase.org/detail.go?type=gene&acc=",
    "ChemBank": "http://chembank.broadinstitute.org/chemistry/viewMolecule.htm?cbid=",
    "ClinVarVariant": "http://www.ncbi.nlm.nih.gov/clinvar/variation/",
    "DBSNP": "http://identifiers.org/dbsnp/",
    "DDANAT": "http://purl.obolibrary.org/obo/DDANAT_",
    "DGIdb": "https://www.dgidb.org/interaction_types",
    "DOID": "http://purl.obolibrary.org/obo/DOID_",
    "DOID-PROPERTY": "http://purl.obolibrary.org/obo/doid#",
    "DRUGBANK": "http://identifiers.org/drugbank/",
    "DrugCentral": "http://drugcentral.org/drugcard/",
    "EC": "http://www.enzyme-database.org/query.php?ec=",
    "ECO": "http://purl.obolibrary.org/obo/ECO_",
    "ECTO": "http://purl.obolibrary.org/obo/ECTO_",
    "EDAM-DATA": "http://edamontology.org/data_",
    "EDAM-FORMAT": "http://edamontology.org/format_",
    "EDAM-OPERATION": "http://edamontology.org/operation_",
    "EDAM-TOPIC": "http://edamontology.org/topic_",
    "EFO": "http://www.ebi.ac.uk/efo/EFO_",
    "EGGNOG": "http://identifiers.org/eggnog/",
    "ENSEMBL": "http://identifiers.org/ensembl/",
    "ENVO": "http://purl.obolibrary.org/obo/ENVO_",
    "ExO": "http://purl.obolibrary.org/obo/ExO_",
    "FAO": "http://purl.obolibrary.org/obo/FAO_",
    "FB": "http://identifiers.org/fb/",
    "FBcv": "http://purl.obolibrary.org/obo/FBcv_",
    "FMA": "http://purl.obolibrary.org/obo/FMA_",
    "FOODON": "http://purl.obolibrary.org/obo/FOODON_",
    "FYPO": "http://purl.obolibrary.org/obo/FYPO_",
    "GAMMA": "http://translator.renci.org/GAMMA_",
    "GENEPIO": "http://purl.obolibrary.org/obo/GENEPIO_",
    "GENO": "http://purl.obolibrary.org/obo/GENO_",
    "GO": "http://purl.obolibrary.org/obo/GO_",
    "GOLD.META": "http://identifiers.org/gold.meta/",
    "GOP": "http://purl.obolibrary.org/obo/go#",
    "GOREL": "http://purl.obolibrary.org/obo/GOREL_",
    "GSID": "https://scholar.google.com/citations?user=",
    "GTEx": "https://www.gtexportal.org/home/gene/",
    "GTOPDB": "https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId=",
    "HAMAP": "http://identifiers.org/hamap/",
    "HANCESTRO": "http://www.ebi.ac.uk/ancestro/ancestro_",
    "HCPCS": "http://purl.bioontology.org/ontology/HCPCS/",
    "HGNC": "http://identifiers.org/hgnc/",
    "HGNC.FAMILY": "http://identifiers.org/hgnc.family/",
    "HMDB": "http://identifiers.org/hmdb/",
    "HP": "http://purl.obolibrary.org/obo/HP_",
    "HsapDv": "http://purl.obolibrary.org/obo/HsapDv_",
    "IAO": "http://purl.obolibrary.org/obo/IAO_",
    "ICD10": "https://icd.codes/icd9cm/",
    "ICD9": "http://translator.ncats.nih.gov/ICD9_",
    "IDO": "http://purl.obolibrary.org/obo/IDO_",
    "INCHI": "http://identifiers.org/inchi/",
    "INCHIKEY": "http://identifiers.org/inchikey/",
    "INO": "http://purl.obolibrary.org/obo/INO_",
    "INTACT": "http://identifiers.org/intact/",
    "IUPHAR.FAMILY": "http://identifiers.org/iuphar.family/",
    "KEGG": "http://identifiers.org/kegg/",
    "KEGG.BRITE": "http://www.kegg.jp/entry/",
    "KEGG.COMPOUND": "http://identifiers.org/kegg.compound/",
    "KEGG.DGROUP": "http://www.kegg.jp/entry/",
    "KEGG.DISEASE": "http://identifiers.org/kegg.disease/",
    "KEGG.DRUG": "http://identifiers.org/kegg.drug/",
    "KEGG.ENVIRON": "http://identifiers.org/kegg.environ/",
    "KEGG.ENZYME": "http://www.kegg.jp/entry/",
    "KEGG.GENE": "http://www.kegg.jp/entry/",
    "KEGG.GLYCAN": "http://identifiers.org/kegg.glycan/",
    "KEGG.MODULE": "http://identifiers.org/kegg.module/",
    "KEGG.ORTHOLOGY": "http://identifiers.org/kegg.orthology/",
    "KEGG.PATHWAY": "https://www.kegg.jp/entry/",
    "KEGG.RCLASS": "http://www.kegg.jp/entry/",
    "KEGG.REACTION": "http://identifiers.org/kegg.reaction/",
    "LOINC": "http://loinc.org/rdf/",
    "MAXO": "http://purl.obolibrary.org/obo/MAXO_",
    "MEDDRA": "http://identifiers.org/meddra/",
    "MESH": "http://id.nlm.nih.gov/mesh/",
    "METANETX.REACTION": "https://www.metanetx.org/equa_info/",
    "MGI": "http://identifiers.org/mgi/",
    "MI": "http://purl.obolibrary.org/obo/MI_",
    "MIR": "http://identifiers.org/mir/",
    "MONDO": "http://purl.obolibrary.org/obo/MONDO_",
    "MP": "http://purl.obolibrary.org/obo/MP_",
    "MPATH": "http://purl.obolibrary.org/obo/MPATH_",
    "MSigDB": "https://www.gsea-msigdb.org/gsea/msigdb/",
    "NBO": "http://purl.obolibrary.org/obo/NBO_",
    "NBO-PROPERTY": "http://purl.obolibrary.org/obo/nbo#",
    "NCBIGene": "http://identifiers.org/ncbigene/",
    "NCBITaxon": "http://purl.obolibrary.org/obo/NCBITaxon_",
    "NCIT": "http://purl.obolibrary.org/obo/NCIT_",
    "NCIT-OBO": "http://purl.obolibrary.org/obo/ncit#",
    "NDC": "http://identifiers.org/ndc/",
    "NDDF": "http://purl.bioontology.org/ontology/NDDF/",
    "NLMID": "https://www.ncbi.nlm.nih.gov/nlmcatalog/?term=",
    "OBAN": "http://purl.org/oban/",
    "OBI": "http://purl.obolibrary.org/obo/OBI_",
    "OBOREL": "http://purl.obolibrary.org/obo/RO_",
    "OGMS": "http://purl.obolibrary.org/obo/OGMS_",
    "OIO": "http://www.geneontology.org/formats/oboInOwl#",
    "OMIM": "http://purl.obolibrary.org/obo/OMIM_",
    "OMIM.PS": "https://www.omim.org/phenotypicSeries/",
    "ORCID": "https://orcid.org/",
    "ORPHA": "http://www.orpha.net/ORDO/Orphanet_",
    "ORPHANET": "http://identifiers.org/orphanet/",
    "PANTHER.FAMILY": "http://www.pantherdb.org/panther/family.do?clsAccession=",
    "PANTHER.PATHWAY": "http://identifiers.org/panther.pathway/",
    "PATO": "http://purl.obolibrary.org/obo/PATO_",
    "PCO": "http://purl.obolibrary.org/obo/PCO_",
    "PFAM": "http://identifiers.org/pfam/",
    "PHARMGKB.PATHWAYS": "http://identifiers.org/pharmgkb.pathways/",
    "PHAROS": "http://pharos.nih.gov",
    "PIRSF": "http://identifiers.org/pirsf/",
    "PMID": "http://www.ncbi.nlm.nih.gov/pubmed/",
    "PO": "http://purl.obolibrary.org/obo/PO_",
    "PR": "http://purl.obolibrary.org/obo/PR_",
    "PRINTS": "http://identifiers.org/prints/",
    "PRODOM": "http://identifiers.org/prodom/",
    "PROSITE": "http://identifiers.org/prosite/",
    "PUBCHEM.COMPOUND": "http://identifiers.org/pubchem.compound/",
    "PUBCHEM.SUBSTANCE": "http://identifiers.org/pubchem.substance/",
    "PW": "http://purl.obolibrary.org/obo/PW_",
    "PathWhiz": "http://smpdb.ca/pathways/#",
    "PomBase": "https://www.pombase.org/gene/",
    "REACT": "http://www.reactome.org/PathwayBrowser/#/",
    "REPODB": "http://apps.chiragjpgroup.org/repoDB/",
    "RFAM": "http://identifiers.org/rfam/",
    "RGD": "http://identifiers.org/rgd/",
    "RHEA": "http://identifiers.org/rhea/",
    "RNACENTRAL": "http://identifiers.org/rnacentral/",
    "RO": "http://purl.obolibrary.org/obo/RO_",
    "RXCUI": "https://mor.nlm.nih.gov/RxNav/search?searchBy=RXCUI&searchTerm=",
    "RXNORM": "http://purl.bioontology.org/ontology/RXNORM/",
    "ResearchID": "https://publons.com/researcher/",
    "SEED.REACTION": "https://modelseed.org/biochem/reactions/",
    "SEMMEDDB": "https://skr3.nlm.nih.gov/SemMedDB",
    "SEPIO": "http://purl.obolibrary.org/obo/SEPIO_",
    "SGD": "http://identifiers.org/sgd/",
    "SIDER.DRUG": "http://identifiers.org/sider.drug/",
    "SIO": "http://semanticscience.org/resource/SIO_",
    "SMART": "http://identifiers.org/smart/",
    "SMPDB": "http://identifiers.org/smpdb/",
    "SNOMED": "http://www.snomedbrowser.com/Codes/Details/",
    "SNOMEDCT": "http://www.snomedbrowser.com/Codes/Details/",
    "SO": "http://purl.obolibrary.org/obo/SO_",
    "STATO": "http://purl.obolibrary.org/obo/STATO_",
    "STY": "http://purl.bioontology.org/ontology/STY/",
    "SUPFAM": "http://identifiers.org/supfam/",
    "ScopusID": "https://www.scopus.com/authid/detail.uri?authorId=",
    "TAXRANK": "http://purl.obolibrary.org/obo/TAXRANK_",
    "TCDB": "http://identifiers.org/tcdb/",
    "TIGRFAM": "http://identifiers.org/tigrfam/",
    "TO": "http://purl.obolibrary.org/obo/TO_",
    "UBERGRAPH": "http://translator.renci.org/ubergraph-axioms.ofn#",
    "UBERON": "http://purl.obolibrary.org/obo/UBERON_",
    "UBERON_CORE": "http://purl.obolibrary.org/obo/uberon/core#",
    "UBERON_NONAMESPACE": "http://purl.obolibrary.org/obo/core#",
    "UMLS": "http://identifiers.org/umls/",
    "UMLSSG": "https://lhncbc.nlm.nih.gov/semanticnetwork/download/sg_archive/SemGroups-v04.txt",
    "UNII": "http://identifiers.org/unii/",
    "UNIPROT.ISOFORM": "http://identifiers.org/uniprot.isoform/",
    "UO-PROPERTY": "http://purl.obolibrary.org/obo/uo#",
    "UPHENO": "http://purl.obolibrary.org/obo/UPHENO_",
    "UniProtKB": "http://identifiers.org/uniprot/",
    "VANDF": "https://www.nlm.nih.gov/research/umls/sourcereleasedocs/current/VANDF/",
    "VMC": "https://github.com/ga4gh/vr-spec/",
    "WB": "http://identifiers.org/wb/",
    "WBPhenotype": "http://purl.obolibrary.org/obo/WBPhenotype_",
    "WBVocab": "http://bio2rdf.org/wormbase_vocabulary",
    "WIKIDATA": "https://www.wikidata.org/wiki/",
    "WIKIDATA_PROPERTY": "https://www.wikidata.org/wiki/Property:",
    "WIKIPATHWAYS": "http://identifiers.org/wikipathways/",
    "WormBase": "https://www.wormbase.org/get?name=",
    "XCO": "http://purl.obolibrary.org/obo/XCO_",
    "XPO": "http://purl.obolibrary.org/obo/XPO_",
    "Xenbase": "http://www.xenbase.org/gene/showgene.do?method=display&geneId=",
    "ZFIN": "http://identifiers.org/zfin/",
    "ZP": "http://purl.obolibrary.org/obo/ZP_",
    "alliancegenome": "https://www.alliancegenome.org/",
    "apollo": "https://github.com/GMOD/Apollo",
    "biolink": "https://w3id.org/biolink/vocab/",
    "bioschemas": "https://bioschemas.org/",
    "dcat": "http://www.w3.org/ns/dcat#",
    "dcid": "https://datacommons.org/browser/",
    "dct": "http://purl.org/dc/terms/",
    "dctypes": "http://purl.org/dc/dcmitype/",
    "dictyBase": "http://dictybase.org/gene/",
    "doi": "https://doi.org/",
    "fabio": "http://purl.org/spar/fabio/",
    "faldo": "http://biohackathon.org/resource/faldo#",
    "foaf": "http://xmlns.com/foaf/0.1/",
    "foodb.compound": "http://foodb.ca/compounds/",
    "gff3": "https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#",
    "gpi": "https://github.com/geneontology/go-annotation/blob/master/specs/gpad-gpi-2-0.md#",
    "gtpo": "https://rdf.guidetopharmacology.org/ns/gtpo#",
    "interpro": "https://www.ebi.ac.uk/interpro/entry/",
    "isbn": "https://www.isbn-international.org/identifier/",
    "isni": "https://isni.org/isni/",
    "issn": "https://portal.issn.org/resource/ISSN/",
    "linkml": "https://w3id.org/linkml/",
    "medgen": "https://www.ncbi.nlm.nih.gov/medgen/",
    "metacyc.reaction": "https://identifiers.org/metacyc.reaction:",
    "mirbase": "http://identifiers.org/mirbase",
    "oboInOwl": "http://www.geneontology.org/formats/oboInOwl#",
    "oboformat": "http://www.geneontology.org/formats/oboInOwl#",
    "os": "https://github.com/cmungall/owlstar/blob/master/owlstar.ttl",
    "owl": "http://www.w3.org/2002/07/owl#",
    "pav": "http://purl.org/pav/",
    "prov": "http://www.w3.org/ns/prov#",
    "qud": "http://qudt.org/1.1/schema/qudt#",
    "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
    "schema": "http://schema.org/",
    "skos": "http://www.w3.org/2004/02/skos/core#",
    "wgs": "http://www.w3.org/2003/01/geo/wgs84_pos",
    "xsd": "http://www.w3.org/2001/XMLSchema#",
}
path = "data"

#######################################################################################

def chunked_dict(input_dict, chunk_size):
    """
    Yield successive chunks of key-value pairs from a dictionary.

    :param input_dict: Input dictionary
    :param chunk_size: The size of each chunk
    :return: Yields chunks of key-value pairs
    """
    items = list(input_dict.items())
    for i in range(0, len(items), chunk_size):
        yield dict(items[i:i + chunk_size])

def chunked_list(input_list, chunk_size):
    """
    Yield successive chunks of key-value pairs from a dictionary.

    :param input_list: Input list
    :param chunk_size: The size of each chunk
    :return: Yields chunks of values
    """
    # items = list(input_list.items())
    items = input_list
    for i in range(0, len(items), chunk_size):
        yield items[i:i + chunk_size]

def fetch_data_from_api(search_string, offset, limit):
    """
    Fetch data from the API using the given search string, offset, and limit via POST request.

    :param search_string: The string to search.
    :param offset: The offset for the results.
    :param limit: The limit for the number of results.
    :return: A list containing the parsed JSON data from the API.
    """

    print('fetch_data_from_api')

    try:
        # API URL
        base_url = "https://name-resolution-sri.renci.org/lookup"

        # Construct the full URL with query parameters
        url = f"{base_url}?string={search_string}&offset={offset}&limit={limit}"

        # Headers
        headers = {
            'accept': 'application/json'
        }

        # Making a POST request
        response = requests.post(url, headers=headers, data='')

        # Check if the request was successful
        if response.status_code == 200:
            # Return the JSON response parsed as a list
            return response.json()
        else:
            # In case of an error return an empty list
            print(f"Error: Unable to fetch data, HTTP Status code {response.status_code}")
            return []
    except Exception as e:
        return None # This None value means we didn't want to bother ;)

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r", color="black"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    # print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd) # Replace prev line with progress
    console.print(f'\r{prefix} |{bar}| {iteration}/{total} | {percent}% {suffix} ', style=color)
    if iteration == total:
        print()

def logTime():
    start_time = time.time()
    print("===========================================================================")
    print(f"Started: {time.strftime('%X')}")
    print("==================================================================")
    return start_time

def endTime(start_time):
    duration = time.time() - start_time
    print("==================================================================")
    print(f"Finished in {duration} seconds / {time.strftime('%X')}")
    print("===========================================================================")

def parse_json_from_file(file_name):
    """
    Parse JSON data from a file.

    Args:
        file_name (str): Name of the file to read the JSON data from.

    Returns:
        dict: The JSON data parsed into a Python dictionary.
    """
    try:
        # Read and parse the JSON data from the file
        with open(file_name, 'r', encoding='utf-8') as file:
            data = json.load(file)
        return data

    except json.JSONDecodeError:
        print('The file does not contain valid JSON data')

    except Exception as e:
        print(f'An error occurred: {e}')

def save_data_as_json(data, file_name):
    """
    Save a text as a UTF-8 encoded JSON file.

    Args:
        text (str): Text to be saved (must be a string representation of a dictionary).
        file_name (str): Name of the file to save the text in.
    """
    try:
        # Write the data to a JSON file with UTF-8 encoding
        with open(file_name, 'w', encoding='utf-8') as file:
            json.dump(data, file, ensure_ascii=False, indent=4)

        print(f'The file was saved successfully as {file_name}')

    except json.JSONDecodeError:
        print('The provided text is not a valid JSON string')

    except Exception as e:
        print(f'An error occurred: {e}')

def genJSON(sentences):

    print('genJSON')

    if os.path.isfile(f"{args.output}"):
        all = parse_json_from_file(f"{args.output}")
    else:
        all = {}

    newResults = {}
    newSentences = []

    for sentence in sentences:
        sentence = sentence["text"]
        if sentence in all:
            newResults[sentence] = all[sentence]
        else:
            newSentences.append(sentence)

    if len(newSentences) > 0:
        try:
            print('genJSON: if len(newSentences) > 0')

            if not args.prompt:
                print("No prompt type selected. Refer to the README: https://github.com/semantisch/diamond-kg.")
                sys.exit()
            elif args.prompt == "triples":
                 promptText = f"""
Paragraphs (separated by a full stop): {" ".join(newSentences)}
Instructions:
Inside each Paragraph, identify all triples and output a JSON dictionary with full Paragraph texts as keys and lists of triples formatted as [SUBJECT, PREDICATE, OBJECT] as values.
Only output the resulting JSON.
"""
            elif args.prompt == "freeContext":
                 promptText = f"""
Paragraphs (separated by " | "): {" | ".join(newSentences)}
Instructions:
Inside each Paragraph, identify all context words/phrases and output a JSON dictionary with full Paragraph texts as keys and lists of {{[context type]: [corresponding words/phrases as a list of values (more than one possible)]}} as values.
Other context types are self-explanatory.
Only output the resulting JSON.
"""
            elif args.prompt == "definedContext":
                 promptText = f"""
Paragraphs (separated by " | "): {" | ".join(newSentences)}
Instructions:
Inside each Paragraph, identify all context words/phrases (allowed context types: "age group", "co-morbidity", "symptom", "co-therapy", "adjunct therapy", "past therapies", "treatment duration", "conditional", "co-prescribed medication", "genetics", "temporal aspects") and output a JSON dictionary with full Paragraph texts as keys and lists of {{[context type]: [corresponding words/phrases as a list of values (more than one possible)]}} as values.
Definition of context types:
"conditional" - a statement about when the medication is appropriate to use
"target" - the condition (symptom or illness) that is intended to be treated
"co-prescribed medication"- drugs commonly prescribed together with the given drug (not therapeutic procedures! -> for that use "co-therapy")
"co-therapy" - procedures or therapies that should be applied in combination with the drug (not medications or substances! -> for that use "co-prescribed medication")
"co-morbidity" - diseases or conditions that commonly occur together (with a target condition) in the same patients
"genetics"-  particular genetic strains of a disease
"temporal aspects"-  information which explains at what life stage, disease stage, or treatment phase a drug should be administered
Other context types are self-explanatory.
Only output the resulting JSON.
"""

            # print(promptText)


            # response = openai.ChatCompletion.create(model="gpt-4", messages = [{"role": "user", "content": promptText}], max_tokens= 6144)

            response = openai.ChatCompletion.create(model="gpt-4", messages = [{"role": "user", "content": promptText}], max_tokens= maxTokens)
            # print(response["choices"][0]["message"]["content"])
            response = response["choices"][0]["message"]["content"]

            print('genJSON: response')
            print(response)

            # Parse the JSON
            data = json.loads(response)

            if args.prompt == "triples":
                for sentence in data:
                    if not sentence in all:
                        all[sentence] = []
                        for triple in data[sentence]:
                            print('genJSON: new triple')
                            newTriple = [{"label" : triple[0], "id" : fetch_data_from_api(triple[0], 0, 5)}, {"label" : triple[1], "id" : fetch_data_from_api(triple[1], 0, 5)}, {"label" : triple[2], "id" : fetch_data_from_api(triple[2], 0, 5)}]
                            print(newTriple)
                            all[sentence].append(newTriple)
                # sys.exit()
            else:
                for sentence in data:
                    for context in data[sentence]:
                        contextValues = data[sentence][context]
                        newContextValues = []
                        for contextValue in contextValues:
                            print('genJSON: new context value')
                            newContextValue = {"contextValue" : contextValue, "id" : fetch_data_from_api(contextValue, 0, 5)}
                            print(newContextValue)
                            newContextValues.append(newContextValue )
                        data[sentence][context] = newContextValues
                    newResults[sentence] = data[sentence]
                    if sentence not in all:
                        all[sentence] = newResults[sentence]

            save_data_as_json(all, f"{args.output}")

        except Exception as e:
            print(e)
            sys.exit()

    # return graph
    # print(newResults)
    return newResults

def get_hash(text):
    return hashlib.sha256(text.encode()).hexdigest()

def get_unique_hash(text, algorithm='sha256'):
    """
    Function to get a unique hash for a given text.

    :param text: The text you want to hash
    :param algorithm: The hashing algorithm you want to use (default is 'sha256')
    :return: The hexadecimal hash string
    """
    # Validate the input
    if not isinstance(text, str):
        raise ValueError("Input text must be a string")

    # Select the hashing algorithm
    hash_function = None
    if algorithm == 'sha256':
        hash_function = hashlib.sha256()
    elif algorithm == 'sha3_256':
        hash_function = hashlib.sha3_256()
    elif algorithm == 'blake2s':
        hash_function = hashlib.blake2s()
    else:
        raise ValueError("Unsupported algorithm. Supported algorithms are: 'sha256', 'sha3_256', 'blake2s'")

    # Hash the text
    hash_function.update(text.encode())

    # Return the hexadecimal hash string
    return hash_function.hexdigest()

#######################################################################################

logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)

console = Console()
danger_style = Style(color="red", blink=True, bold=True)

parser=argparse.ArgumentParser()
parser.add_argument("--input", "-i", help="Input file")
parser.add_argument("--output", "-o", help="Output file")
parser.add_argument("--apiKey", "-a", help="OpenAI API key")
parser.add_argument("--organization", "-org", help="OpenAI API organization")
parser.add_argument("--prompt", "-p", help="Prompt type")
parser.add_argument("--maxTokens", "-m", help="Maximum tokens")
parser.add_argument("--chunkSize", "-c", help="Chunk size")
args=parser.parse_args()

if not args.apiKey or not args.organization:
    print("OpenAPI key or organization undefined.")
    sys.exit()

# if not args.maxTokens:
#     print("Max tokens not defined (max value: 6144).")
#     sys.exit()

# openai.organization = "org-B7EpIIq1IYHKKXHqeoakLpOz"
openai.organization = args.organization
# openai.api_key = "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"
openai.api_key = args.apiKey

openai.Model.list()

#######################################################################################


if args.input:
    sentences = parse_json_from_file(f"{args.input}")
else:
    print('Define --input file path')
    sys.exit()

if not args.output:
    print('Define --output file path')
    sys.exit()

newResults = []
counter = 0

maxTokens = 6144
if args.maxTokens :
    maxTokens  = int(args.maxTokens)

chunkSize = 10
if args.chunkSize:
    chunkSize = int(args.chunkSize)

for chunk in chunked_list(sentences, chunkSize):
    try:
        counter = counter + 1
        print(chunk)
        # sys.exit()
        newResults.append(genJSON(chunk))
    except Exception as e:
        print(e)
    printProgressBar(counter, len(sentences), prefix = 'From all inputs:', suffix = 'Complete', length = 50)

print(newResults)

print(all)
