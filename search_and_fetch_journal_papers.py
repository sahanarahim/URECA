#@markdown Run this cell to load in get_papers(search_term) function
from Bio import Entrez

import os
import openai
import json
import time
import pickle 

prompt = '''
Write a very short summary about the functions of genes in this abstract. The summary must show pair-wise relationships, for example:
gene: !affects! Process
gene: !localizes to! X
gene: !interacts with! Y
gene: !enhances! Z
gene: !represses! U
gene: !synthesizes! I

Please provide only one statement per line, and ensure that each line contains exactly two actors. If a relationship involves more than two actors, please break it down into multiple separate lines.

<ENTER PROMPT HERE>

VERY SHORT, CONCISE SUMMARY CONTAINING ALL INFORMATION WITH TWO ACTORS PER LINE: 

'''

openai.api_key = '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'


def gpt_analyzer(abstract):
    print('Abstract:',  abstract)
    my_prompt = prompt.replace('<ENTER PROMPT HERE>',abstract)

    response = openai.Completion.create(
      model="text-davinci-003",
      prompt=my_prompt,
      temperature=0,
      max_tokens=256,
      top_p=1.0,
      frequency_penalty=0.0,
      presence_penalty=1
    )
    print(response['choices'][0]['text']+'\n')

    return response['choices'][0]['text']

def search_and_fetch_journal_papers(journal_name, email, year_cutoff, retmax=100000):
    Entrez.email = email
    query = f'"{journal_name}"[Journal] AND "{year_cutoff}/01/01"[Date - Publication] : "3000"[Date - Publication]'
    handle = Entrez.esearch(db='pubmed',
                            sort='relevance',
                            retmax=retmax,
                            term=query,
                            usehistory='y')
    results = Entrez.read(handle)
    handle.close()

    pmid_list = results['IdList']
    handle = Entrez.efetch(db='pubmed', id=','.join(pmid_list), retmode='xml')
    papers = Entrez.read(handle)['PubmedArticle']
    handle.close()

    paper_details = []

    for paper in papers:
        pmid = paper['MedlineCitation']['PMID']
        title = paper['MedlineCitation']['Article']['ArticleTitle']
        abstract = paper['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', None)
        abstract = ' '.join(abstract) if abstract is not None else ''
        date = paper['MedlineCitation']['Article']['ArticleDate'][0] if paper['MedlineCitation']['Article'].get('ArticleDate') else paper['MedlineCitation']['DateRevised']
        journal = journal_name
        doi = paper['PubmedData']['ArticleIdList'][-1] if paper['PubmedData'].get('ArticleIdList') else None
        author_list = paper['MedlineCitation']['Article'].get('AuthorList', [])
        authors = [f"{author.get('LastName', '')} {author.get('ForeName', '')}" for author in author_list if author.get('LastName')]

        paper_info = {
            'pmid': str(pmid),
            'date': f"{date['Year']}-{date['Month']}-{date['Day']}",
            'journal': journal,
            'title': str(title),
            'abstract': abstract,
            'doi': str(doi),
            'authors': authors
        }

        paper_details.append(paper_info)

    print(len(paper_details), journal_name)
    return paper_details

email = 'mutwil@gmail.com'
year_cutoff = 2023

papers = []

papers += search_and_fetch_journal_papers("Plant Cell", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Plant Physiol", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Plant J", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Plant Cell Physiol", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Plant Mol Biol", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Front Plant Sci", email, year_cutoff)
# papers += search_and_fetch_journal_papers("J Exp Bot", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Mol Plant", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Planta", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Plant Signal Behav", email, year_cutoff)
# papers += search_and_fetch_journal_papers("New Phytol", email, year_cutoff)
# papers += search_and_fetch_journal_papers("BMC Plant Biol", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Mol Plant Microbe Interact", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Plant Cell Environ", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Nat Plants", email, year_cutoff)
# papers += search_and_fetch_journal_papers("Physiol Plant", email, year_cutoff)
# papers += search_and_fetch_journal_papers("J Integr Plant Biol", email, year_cutoff)
# papers += search_and_fetch_journal_papers("J Plant Physiol", email, year_cutoff)

counter = 0
for i in papers:
    print(i['abstract'], str(i['pmid']))
    pubmed, abstract = str(i['pmid']), i['abstract']
    if pubmed+'.txt' not in os.listdir(os.getcwd()+'/annotations'):
        try:
            result = gpt_analyzer(abstract)
            v = open(os.getcwd()+'/annotations/%s.txt' % pubmed,'w')
            v.write(str(abstract)+'\n\n'+str(result))
            v.close()
        except:
            print('gpt failed', i)
    else:
        print(i, 'already done')
    counter +=1
    print(counter, 'out of', len(papers))