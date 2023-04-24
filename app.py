from flask import Flask, render_template, request, redirect, url_for
import pickle
import re
from Bio import Entrez


def generate_cytoscape_js(elements):
    nodes = [
        "{ data: { id: '%s' } }" % node
        for node in set(edge["source"] for edge in elements) | set(edge["target"] for edge in elements)
    ]
    edges = [
        "{ data: { id: 'edge%s', source: '%s', target: '%s', interaction: '%s' } }" % (
            i,
            edge['source'],
            edge['target'],
            edge['interaction'],
        )
        for i, edge in enumerate(elements)
    ]

    a = open('network.txt','r').read()

    nodes = ', '.join(nodes)
    edges = ', '.join(edges)
    return a.replace('NODES',nodes).replace('EDGES',edges)

def process_network(elements):    
    #remove redundancies
    edges = []
    for i in elements:
        if i not in edges:
            edges.append(i)

    #groups similar nodes that have same source+interaction
    edgeTypes = {}
    for i in elements:
        key = (i['source'], i['interaction'])
        if key in edgeTypes:
            edgeTypes[key] += [i['target']]
        else:
            edgeTypes[key] = [i['target']]
    return edges
    
def make_text(elements):    
    '''Given all edges in the KnowledgeNet, it makes the text summary'''
    pubmedLink = '<a href="https://pubmed.ncbi.nlm.nih.gov/%s" target="_blank">%s</a>'
    #Paragraph order: node by the highest degree. Sentence order in a paragraph: node,interaction type by the number of targets.  
    topicDic = {}
    nodeDegree, nodeSentenceDegree = {}, {}
    for i in elements:
        if i.id not in topicDic:
            topicDic[i.id]={}
            nodeDegree[i.id]=0
            nodeSentenceDegree[i.id]={}
        if i.inter_type not in topicDic[i.id]:
            topicDic[i.id][i.inter_type] = [[i.target, i.publication]]
            nodeSentenceDegree[i.id][i.inter_type] = 1
            nodeDegree[i.id] += 1
        else:
            topicDic[i.id][i.inter_type] += [[i.target, i.publication]]
            nodeDegree[i.id] += 1
            nodeSentenceDegree[i.id][i.inter_type] += 1

    sorted_nodes = sorted(nodeDegree, key=lambda x: nodeDegree[x], reverse=True)    
    save = []
    for i in sorted_nodes:
        sorted_sentences = sorted(nodeSentenceDegree[i], key=lambda x: nodeSentenceDegree[i][x], reverse=True)
        tempSentences = []
        for j in sorted_sentences:
            text = '<span  style="color: red;">'+i+'</span>'+' '+'<span  style="color: orange;">'+j+'</span>'+' '
            temp = {}
            tempRefs = []
            for k in topicDic[i][j]:
                if k[0] not in temp:
                    temp[k[0]] = [pubmedLink %(k[1],k[1])]
                else:
                    temp[k[0]] += [pubmedLink %(k[1],k[1])]
            
            for target in temp:
                tempRefs += [target + ' ('+', '.join(list(set(temp[target])))+')']
            tempSentences.append(text+', '.join(tempRefs))

        finishedSentence= '. '.join(tempSentences)+'.'
        save.append(finishedSentence)
    return '<br><br>'.join(save)

def find_terms(my_search, genes):   
    forSending = []
    elements = [] 

    if len(my_search)>0:
        for i in genes:
            if '!' in my_search: #exact search
                if my_search.upper().strip().replace('!','') == i.strip():
                    for j in genes[i]:
                        if j[0]!='' and j[2]!='':
                            forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type
                            elements.append({"source": j[0].replace("'","").replace('"',''), "target": j[2].replace("'","").replace('"',''), "interaction": j[1].replace("'","").replace('"','')})                
            if '?' in my_search: #exact word
                print(my_search.upper().strip().replace('?',''))
                if my_search.upper().strip().replace('?','') in i.strip().split():
                    for j in genes[i]:
                        if j[0]!='' and j[2]!='':
                            forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type
                            elements.append({"source": j[0].replace("'","").replace('"',''), "target": j[2].replace("'","").replace('"',''), "interaction": j[1].replace("'","").replace('"','')})                

            
            if my_search.upper().strip() in i.strip():
                for j in genes[i]:
                    if j[0]!='' and j[2]!='':
                        forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type
                        elements.append({"source": j[0].replace("'","").replace('"',''), "target": j[2].replace("'","").replace('"',''), "interaction": j[1].replace("'","").replace('"','')})
    return elements, forSending

app = Flask(__name__)

class Gene:
    def __init__(self, id, description, inter_type, publication):
        self.id = id
        self.target = description
        self.inter_type = inter_type
        self.publication = publication

@app.route('/')
def index():
    v = open('stats.txt','r').read().rstrip().split()

    journals, numbers = open('journal_statistics.txt','r').read().splitlines()
    piechart = open('piechart.txt','r').read()
    piechart = piechart.replace('JOURNALS', journals)
    piechart = piechart.replace('NUMBERS', numbers)

    return render_template('index.html', entities = v[1], papers = v[0], piechart_code = piechart)
    
@app.route('/author', methods=['POST'])
def author():
    try:
        my_search = request.form["author"].lower()
        replacements = {"ä": "ae", "ö": "oe", "ü": "ue", "ß": "ss", "é": "e", "ô": "o", "î": "i", "ç": "c"}
        my_search = ''.join(replacements.get(c, c) for c in my_search)
    except:
        my_search =''

    if my_search!='':
        with open('authors', 'rb') as f:
            # Load the object from the file
            papers = pickle.load(f)
           
        hits = []
        for author in papers:      
            if len(set(my_search.split())&set(author.lower().split()))==len(set(my_search.split())):
                hits+=papers[author]
                    
        # provide your email address to the Entrez API
        Entrez.email = "mutwil@gmail.com"

        # search query to find papers by an author
        search_query = my_search+"[Author]"

        # perform the search and retrieve the count of papers
        handle = Entrez.esearch(db="pubmed", term=search_query)
        record = Entrez.read(handle)
        count = record["Count"]
        
        forSending = []
        if hits!=[]:
            with open('allDic', 'rb') as file:
                genes = pickle.load(file)
            
            
            elements = []
            papers = []
            for i in genes:
                for j in genes[i]:
                    if j[3] in hits:
                        if j[0]!='' and j[2]!='':
                            papers.append(j[3])
                            forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type
                            elements.append({"source": j[0].replace("'","").replace('"',''), "target": j[2].replace("'","").replace('"',''), "interaction": j[1].replace("'","").replace('"','')})                
                        break

    
    if forSending!=[]:
        elements = process_network(elements)
        cytoscape_js_code = generate_cytoscape_js(elements)
        warning = ''
        summaryText = make_text(forSending)
        
        if len(elements)>400:
            warning = 'The network might be too large to be displayed, so click on "Layout Options", select the edge types that you are interested in and click "Recalculate layout".'

        return render_template('author.html', genes=forSending, cytoscape_js_code=cytoscape_js_code, ncbi_count=count, author= my_search, connectome_count=len(set(papers)), warning=warning, summaryText=summaryText)
    else:
        return render_template('not_found.html')
        
@app.route('/title', methods=['POST'])
def title():

    my_search = request.form['title'].lower()
    pmids = []
    for i in my_search.split(';'):
        pmids+=i.split()

            
    forSending = []
    print(my_search)
    if pmids!=[]:
        
        with open('titles', 'rb') as f:
            # Load the object from the file
            papers = pickle.load(f)

        hits = list(set(pmids)&set(papers))
        
        
        if hits!=[]:
            with open('allDic', 'rb') as file:
                genes = pickle.load(file)
            
            
            elements = [] 
            for i in genes:
                for j in genes[i]:
                    if j[3] in hits:
                        if j[0]!='' and j[2]!='':
                            forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type 
                            elements.append({"source": j[0].replace("'","").replace('"',''), "target": j[2].replace("'","").replace('"',''), "interaction": j[1].replace("'","").replace('"','')})                
                        break

    if forSending!=[]:
        elements = process_network(elements)
        cytoscape_js_code = generate_cytoscape_js(elements)
        summaryText = make_text(forSending)
        return render_template('gene.html', genes=forSending, cytoscape_js_code=cytoscape_js_code, number_papers = len(hits), search_term = my_search, summaryText=summaryText)
    else:
        return render_template('not_found.html')

@app.route('/search', methods=['POST'])
def search():
    with open('allDic', 'rb') as file:
        genes = pickle.load(file)

    try:
        my_search = request.form['gene_id']
    except:
        my_search='cesa'
         
    if len(my_search)>0:
        split_search = my_search.split(';')
        
        forSending = []
        elements = []

        for term in split_search:
            results = find_terms(term, genes)
            elements += results[0]
            forSending += results[1]
            
        papers = []
        for i in forSending:
            papers+=[i.publication]
            
        summaryText = make_text(forSending)
        warning = ''
        if len(elements)>400:
            warning = 'The network might be too large to be displayed, so click on "Layout Options",  select the edge types that you are interested in and click "Recalculate layout".'
        elements = process_network(elements)
        cytoscape_js_code = generate_cytoscape_js(elements)

    if forSending!=[]:
        return render_template('gene.html', genes=forSending, cytoscape_js_code=cytoscape_js_code, search_term = my_search, number_papers = len(set(papers)), warning = warning, summaryText=summaryText)
    else:
        return render_template('not_found.html')

@app.route('/help')
def help():
    return render_template('help.html')

if __name__ == '__main__':
    import os
    items, edges = [],0
    for i in os.listdir(os.getcwd()+'/annotations/'):
        a = open(os.getcwd()+'/annotations/'+i,'r').read()
        
        if len(a)>0:
            findings = a.split('\n\n')[1].split('\n')
            
            for j in findings:
                if j.count('!')==2:
                    splitta = j.split('!')
        
                    agentA, _, agentB = splitta
                    agentA = agentA.split(':')[0].upper()
                    agentB = agentB.strip().upper()
                    edges+=1
                    items+=[agentA]
                    items+=[agentB]
                
    
    v = open('stats.txt','w')
    v.write(str(len(os.listdir(os.getcwd()+'/annotations/')))+'\t'+str(len(set(items))))
    v.close()
    app.run(debug=True)