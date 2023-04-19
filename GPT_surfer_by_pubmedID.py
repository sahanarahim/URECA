#@markdown Run this cell to load in get_papers(search_term) function
from Bio import Entrez

import os
import openai
import json


# Enter your email address for Entrez authentication
Entrez.email = "mutwil@gmail.com"
openai.api_key = "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"

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

def get_papers(id_list):

  # Get the list of IDs for the search results

  # Fetch the abstracts for the IDs
  handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
  records = handle.read()

  # Split the records by abstract
  records = records.split("\n\n")

  # Create a text file to write the titles and abstracts to

  paper2abs = {}
  for record in records:
      
      pmid = None
      date = None
      journal = None
      title = None
      abstract = None
      doi = None
      author = None

      # Split the record by line
      lines = record.split("\n")
      doing_abstract, doing_title = 0,0

      # Extract the full title and abstract from the lines
      for i, line in enumerate(lines):
          #print(line)
          if line!='':
            if line.startswith('PMID- '):
              pmid = line[6:]
            if line.startswith('DP  '):
              date = line[6:]
            if line.startswith('TA'):
              journal = line[6:]
            if line.startswith('LID'):
              doi = line[6:].split()[0]
            if line.startswith('AID'):
              doi = line[6:].split()[0]

            if line.startswith('AU'):
              author = line[6:].split()[0]

            if line.startswith("TI  - "):
                    title = line[6:]
                    doing_title = 1
                    continue
            
            #print("1", line[0], line[0]==' ', doing_title)
            if line[0]==' ' and doing_title==1:
              title += line[6:]
            else:
              doing_title = 0

            if line.startswith("AB  - "):
                    abstract = line[6:]
                    doing_abstract = 1
                    continue
            
            if line[0]==' ' and doing_abstract == 1:
                    abstract += " " + line[6:]
            else:
              doing_abstract = 0            

      # Write the title and abstract to the file
      if title and abstract:
          paper2abs[title] = [pmid, abstract, date, journal, doi, author]

  return paper2abs

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

#getting genes with experimental evidence
gene2function = {}
pubmedIDs = []
for i in open('Z:/Projects/GPT_annotator/ATH_GO_GOSLIM.txt','r').readlines()[1:]:
  splitta = i.rstrip().split('\t')
  if splitta[0] not in gene2function:
    gene2function[splitta[0]] = []
  
  if len(splitta)>10:
    if splitta[9] in ['IDA','IPI', 'EXP', 'IMP','IGI']:
     if splitta[12].count(':')==2:
        gene2function[splitta[0]] += [splitta[12].split(':')[2].strip()]
        pubmedIDs+= [splitta[12].split(':')[2].strip()]


unique_ids = list(set(pubmedIDs))

found_papers = get_papers(unique_ids)

counter = 0
for i in found_papers:
    pubmed, abstract = found_papers[i][0], found_papers[i][1]
    if pubmed+'.txt' not in os.listdir('D:/GPT_annotator/annotations'):
        result = gpt_analyzer(abstract)
        try:
            v = open('D:/GPT_annotator/annotations/%s.txt' % pubmed,'w')
            v.write(str(found_papers[i])+'\n\n'+str(result))
            v.close()
        except:
            print('failed writing', i)
    else:
        print(i, 'already done')
    counter +=1
    print(counter, 'out of', len(found_papers))


