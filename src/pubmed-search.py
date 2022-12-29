#!/usr/bin/env python

from Bio import Entrez
import json, os, socket, sys, xmltodict

page_size = 100
page = 0

ohost, oport = os.environ["SR_OUTPUT"].split(":")
sr_output = socket.create_connection((ohost, int(oport)))
config = json.load(open(os.environ['SR_CONFIG']))
Entrez.email = config['reviewer']
query = config['current_step'].get('query')
if not query:
    print("No query specified", file=sys.stderr)
    sys.exit(1)

def get_abstract(article):
    abstract = article['PubmedArticleSet']['PubmedArticle']['MedlineCitation']['Article'].get('Abstract')
    if abstract:
        texts = abstract.get('AbstractText')
        if isinstance(texts, str):
            return texts
        abstract = ''
        if not isinstance(texts, list):
            texts = [texts]
        for x in texts:
            text = x.get('#text')
            if text:
                if abstract:
                    abstract += '\n\n'
                abstract += text
        if abstract:
            return abstract

def get_title(article):
    return article['PubmedArticleSet']['PubmedArticle']['MedlineCitation']['Article'].get('ArticleTitle')

def get_article(id):
    handle = Entrez.efetch(db="pubmed", id=id, rettype="xml", retmode="text")
    article = handle.read()
    handle.close()
    article = xmltodict.parse(article)
    article['abstract'] = get_abstract(article)
    article['title'] = get_title(article)
    return {
        "data": article,
        "type": "document",
        "url": "https://www.ncbi.nlm.nih.gov/pubmed/" + id,
    }

while True:
    handle = Entrez.esearch(db="pubmed",
                            term=query,
                            retmax=page_size,
                            retstart=page * page_size)
    record = Entrez.read(handle)
    handle.close()

    for id in record["IdList"]:
        sr_output.send(json.dumps(get_article(id)).encode())
        sr_output.send(b"\n")

    if len(record["IdList"]) < page_size:
        break
    page += 1
