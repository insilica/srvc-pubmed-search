#!/usr/bin/env python

from Bio import Entrez

Entrez.email = "john@insilica.co"

# Search for articles with the term "bronchitis" in the title
query = "bronchitis[Title]"
handle = Entrez.esearch(db="pubmed", term=query)
record = Entrez.read(handle)

# Print the number of articles found
print(f"{record['Count']} articles found")

# Print the IDs of the articles
print(record["IdList"])
