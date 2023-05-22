from flask import Flask, render_template, request, redirect, url_for
import pickle
import re
from Bio import Entrez
import networkx as nx
import json
import os

# Importing custom utilities 
from utils import edges as e

'''
# == DELETE AFTER ==
connections = {}
def find_connections(f):
    try: 
        for j in open(f'./annotations/{f}').readlines()[2:]:
            n1, edge, n2 = list(map(lambda x : x.strip().replace(':', '').title(), j.split("!")))
            key = f"{n1}%?%{edge}"
            
            if key not in connections:
                connections[key] = []
            connections[key].append(n2)
    except:
        pass

for f in os.listdir("./annotations"):
    find_connections(f)
connection_keys = sorted(connections, key = lambda x : len(connections[x]), reverse = True)
'''
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
    elements = list(set(elements))    
    if len(elements)>500:  # only perform node filtration if the number of relationship is larger than 500
        #using Networkx multiDiGraph to model input data as Graph
        G =nx.MultiDiGraph()
        for i in elements:
            G.add_edge(i[0], i[1], relation = i[2])
        
        G,ref =nodeDegreeFilter(G)
        
        return graphConverter(G,ref)
    else:
        return edgeConverter(elements)

def nodeDegreeFilter(graph):
    nodesToKeep = []
    totalDegree = 0
    for i in sorted(dict(graph.degree), key=dict(graph.degree).get, reverse=True):
        currDegree = graph.degree(i)
        if totalDegree <= 500 or nodesToKeep==[]: #must keep at least one nodes
            nodesToKeep.append(i)
            totalDegree += currDegree
        if totalDegree > 500:
            break
    lst = set()
    for j in nodesToKeep:
        lst.update(list(graph.in_edges(j))+list(graph.out_edges(j)))
    ref ={}
    for i in lst:  # create a hash map for the entities to be kept
        ref[i] = ref.get(i,0)+1
    return graph, ref
            
def graphConverter(graph, ref):
    updatedElements = []
    for k,v in graph.adjacency():
        source = str(k).replace("'","").replace('"','')
        if v:
            for i,j in v.items():
                target = str(i).replace("'","").replace('"','')
                for p,q in j.items():
                    type = str(q['relation']).replace("'","").replace('"','')
                    if (source, target) in ref:
                        updatedElements.append({"source": source, "target": target, "interaction": type})
    return updatedElements 
            
def edgeConverter(elements): #Convert Edges to default dictionary format 
    updatedElements = []
    for i in elements:
        updatedElements.append({"source": str(i[0]).replace("'","").replace('"',''), "target": str(i[1]).replace("'","").replace('"',''), "interaction": str(i[2]).replace("'","").replace('"','')})
    return updatedElements 

def make_text(elements):    
    '''Given all edges in the KnowledgeNet, it makes the text summary'''
    pubmedLink = '<span class="pubmed-link" style="color:blue" data-pubmed-id="%s" data-source="%s" data-typa="%s" data-target="%s">%s</span>'
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
    # print(sorted_nodes)
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
                    temp[k[0]] = [pubmedLink %(k[1],i,j,k[0],k[1])]
                else:
                    temp[k[0]] += [pubmedLink %(k[1],i,j,k[0],k[1])]
            
            for target in temp:
                tempRefs += [target + ' ('+', '.join(list(set(temp[target])))+')']
            tempSentences.append(text+', '.join(tempRefs))

        finishedSentence= '. '.join(tempSentences)+'.'
        save.append(finishedSentence)
    return '<br><br>'.join(save)

def find_terms_helper(gene, genes):
    forSending = []
    elements =[]
    for j in genes[gene]:
        if j[0]!='' and j[2]!='':
            forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type
            elements.append((j[0].replace("'","").replace('"',''),  j[2].replace("'","").replace('"',''), j[1].replace("'","").replace('"','')))
    return elements, forSending

def find_terms(my_search, genes):   
    if my_search:
        if '!' in my_search: #exact search -Specific word
            for i in genes:
                if my_search.upper().strip().replace('!','') == i.strip():
                    elements, forSending = find_terms_helper(i,genes)
                    return list(set(elements)), forSending
        
        if '?' in my_search: # include any term that contains the substring of the query 
            forSending = []
            elements = [] 
            for i in genes:
                if my_search.upper().strip().replace('?','') in i.strip():
                    outputOne, outputTwo = find_terms_helper(i,genes)
                    elements.extend(outputOne)
                    forSending.extend(outputTwo)
            return list(set(elements)), forSending
        
        if '&' in my_search: # include any related terms 
            forSending = []
            elements = []
            for i in genes:
                substring = my_search.upper().strip().replace('&','')
                string = i.strip()
                if substring in string:
                    if not is_alphanumeric_helper(string, substring):
                        outputOne, outputTwo = find_terms_helper(i,genes)
                        elements.extend(outputOne)
                        forSending.extend(outputTwo)
            return list(set(elements)), forSending
        
        if '@' in my_search: # include any genes with similar gene alias 
            with open('geneAlias', 'rb') as file:
                adjMat = pickle.load(file)
            try:
                terms = adjMat[my_search.upper().strip().replace('@','')]
            except:
                terms = []
            forSending = []
            elements = []
            for i in terms:
                for j in genes:
                    if i.upper().strip() in j.strip().split():
                        outputOne, outputTwo = find_terms_helper(j,genes)
                        elements.extend(outputOne)
                        forSending.extend(outputTwo)
            return list(set(elements)), forSending
        
        forSending = []
        elements = []
        for i in genes:
            if my_search.upper().strip() in i.strip().split():  # default search - terms that contain the specific query word
                outputOne, outputTwo = find_terms_helper(i,genes)
                elements.extend(outputOne)
                forSending.extend(outputTwo)
        return list(set(elements)), forSending

def is_alphanumeric_helper (string, substring):
    patternLeft = r'[a-zA-Z0-9]'+ re.escape(substring)
    patternRight = re.escape(substring) + r'[a-zA-Z0-9]'
    matchesLeft = re.findall(patternLeft, string)
    matchesRight = re.findall(patternRight, string)
    
    return (len(matchesLeft) > 0) or (len(matchesRight) > 0)

app = Flask(__name__)

class Gene:
    def __init__(self, id, description, inter_type, publication):
        self.id = id
        self.target = description
        self.inter_type = inter_type
        self.publication = publication
    
    def __repr__(self):
        return str(self.__dict__)
    
    def __str__(self):
        return str(self.__dict__)
    
    def getElements(self):
        return (self.id, self.target, self.inter_type)

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
            with open('allDic2', 'rb') as file:
                genes = pickle.load(file)
            
            
            elements = []
            papers = []
            for i in genes:
                for j in genes[i]:
                    if j[3] in hits:
                        if j[0]!='' and j[2]!='':
                            papers.append(j[3])
                            forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type
                            elements.append((j[0].replace("'","").replace('"',''),  j[2].replace("'","").replace('"',''), j[1].replace("'","").replace('"','')))                
                        break

    
    if forSending!=[]:
        updatedElements = process_network(elements)
        cytoscape_js_code = generate_cytoscape_js(updatedElements)
        warning = ''
        summaryText = make_text(forSending)
        
        if len(elements)>500:
            warning = 'The network might be too large to be displayed, so click on "Layout Options", select the edge types that you are interested in and click "Recalculate layout".'

        return render_template('author.html', genes=forSending, cytoscape_js_code=cytoscape_js_code, ncbi_count=count, author= my_search, connectome_count=len(set(papers)), warning=warning, summary=summaryText)
    else:
        return render_template('not_found.html')
        
@app.route('/title', methods=['POST'])
def title():
    remove_redundant_edges = request.form.get('redundancy')
    my_search = request.form['title'].lower()
    pmids = []
    for i in my_search.split(';'):
        pmids+=i.split()

            
    forSending = []
    if pmids!=[]:
        
        with open('titles', 'rb') as f:
            # Load the object from the file
            papers = pickle.load(f)

        hits = list(set(pmids)&set(papers))
        
        
        if hits!=[]:
            with open('allDic2', 'rb') as file:
                genes = pickle.load(file)
            
            
            elements = [] 
            for i in genes:
                for j in genes[i]:
                    if j[3] in hits:
                        if j[0]!='' and j[2]!='':
                            forSending.append(Gene(j[0], j[2], j[1], j[3])) #source, target, type 
                            elements.append((j[0].replace("'","").replace('"',''),  j[2].replace("'","").replace('"',''), j[1].replace("'","").replace('"','')))                
                        break

    if forSending!=[]:
        updatedElements = process_network(elements)
        cytoscape_js_code = generate_cytoscape_js(updatedElements)
        summaryText = make_text(forSending)
        return render_template('gene.html', genes=forSending, cytoscape_js_code=cytoscape_js_code, number_papers = len(hits), search_term = my_search, summary=summaryText)
    else:
        return render_template('not_found.html')

@app.route('/search', methods=['POST'])
def search():
    with open('allDic2', 'rb') as file:
        genes = pickle.load(file)

    try:
        my_search = request.form['gene_id']
        remove_redundant_edges = request.form.get('redundancy') is not None
        if remove_redundant_edges:
            print("removing passive edges")
    except:
        my_search = 'cesa'
         
    if len(my_search) > 0:
        split_search = my_search.split(';')
        forSending = []
        elements = []

        for term in split_search:
            results = find_terms(term, genes)#, remove_redundant_edges)
            elements += results[0]
            forSending += results[1]
            
        # remove redundancies
        elements = list(set(elements))
        
        warning = ''
        if len(elements)>500:
            warning = 'The network might be too large to be displayed, so click on "Layout Options",  select the edge types that you are interested in and click "Recalculate layout".'
        
        updatedElements = process_network(elements)
        cytoscape_js_code = generate_cytoscape_js(updatedElements)

        papers = []
        for i in forSending:
            papers+=[i.publication]
            
        summaryText = make_text(forSending)
        
    if forSending!=[]:
        return render_template('gene.html', genes=forSending, cytoscape_js_code=cytoscape_js_code, search_term = my_search, number_papers = len(set(papers)), warning = warning, summary=summaryText)
        
    else:
        return render_template('not_found.html')

@app.route('/help')
def help():
    '''
    Renders the help template
    '''
    return render_template('help.html')

if __name__ == '__main__':
    # import os
    # items, edges = [],0
    # for i in os.listdir(os.getcwd()+'/annotations/'):
    #     a = open(os.getcwd()+'/annotations/'+i,'r',encoding='ISO-8859-1').read()
        
    #     if len(a)>0:
    #         findings = a.split('\n\n')[1].split('\n')
            
    #         for j in findings:
    #             if j.count('!')==2:
    #                 splitta = j.split('!')
        
    #                 agentA, _, agentB = splitta
    #                 agentA = agentA.split(':')[0].upper()
    #                 agentB = agentB.strip().upper()
    #                 edges+=1
    #                 items+=[agentA]
    #                 items+=[agentB]
                
    
    # v = open('stats.txt','w')
    # v.write(str(len(os.listdir(os.getcwd()+'/annotations/')))+'\t'+str(len(set(items))))
    # v.close()
    app.run()