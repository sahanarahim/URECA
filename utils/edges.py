'''
This is a dictionary that maps active edges to passive ones (i.e., ones in passive voice)
This will be used to remove all instances of passive edges.

UPDATE 04/05/2023 - I've gotten started on this dictionary, but it's still very much a work in 
progress.
'''
PASSIVE_EDGES = {
    ('ENHANCED BY') : 'ENHANCES',
    ('INCREASED BY') : 'INCREASES',
    ('UPREGULATED BY', 'UP-REGULATED BY') : 'UP-REGULATES',
    ('REPRESSED BY') : 'REPRESSES',
    ('DOWNREGULATED BY', 'DOWN-REGULATED BY') : 'DOWN-REGULATES',
    ('INHIBITED BY') : 'INHIBTS',
    ('ACTIVATED BY') : 'ACTIVATES',
    ('INFLUENCED BY') : 'INFLUENCES',
    ('CONTROLLED BY') : 'CONTROLLED',
    ('MODULATED BY') : 'MODULATES',
    ('OVEREXPRESSED IN') : 'OVEREXPRESSES',
    ('TARGETED BY') : 'TARGETS',
    ('MODULATED BY') : 'MODULATES',
    ('CAUSED BY') : 'CAUSES'
}

GROUPINGS = {
    'CLUSTER 19' : ['temperatures', 'susceptible', 'basic', 'zea', 'tolerance', 'domestication', 'locus', 
                    'controlled', 'areas', 'transcriptome', 'coldinduced', 'genes', 'bzip', 'differential', 
                    'represses', 'could', 'major', 'zipper', 'maize', 'transcription'],
    'CLUSTER 4' : ['infections', 'sources', 'plants', 'signatures', 'ca2', 'far', 'receptors', 'free', 
                   'new', 'cytosolic', 'diverse', 'launch', 'signaling', 'branches', 'identification', 
                   'elevation', 'immunity', 'among', 'potential', 'homeostasis'],
    'CLUSTER 12' : ['sequences', 'plants', 'modules', 'technologies', 'time', 'coordinate', 
                    'biology', 'extraordinary', 'fully', 'include', 'gene', 'assemblies', 'compare', 
                    'highquality', 'improvement', 'advances', 'enhancers', 'characteristics', 
                    'identification', 'lag'],
    'CLUSTER 16' : ['functional', 'establish', 'diversification', 'find', 
                    'phenotypes', 'overall', 'posttranscriptional', 'unrecognized', 
                    'connection', 'gene', 'signaling', 'show', 'precursor', 'mrnas', 
                    'findings', 'mutants', 'mutant', 'constitute', 'rice', 'describe'],
    'CLUSTER ???' : []
}


def search_passive_edges(edge: str) -> str:
    '''
    Given an 'edge', this function tries to search for an edge inside to_search.

    If a key is found, then the upper-cased version of that key is returned.  
    Otherwise, a 'None' is returned.
    '''
    for i in list(PASSIVE_EDGES.items()):
        if edge.upper() in i[0]:
            return i[1].upper()
    return None

def simplify_edge(edge: str, to_search: dict) -> str:
    '''
    Determines the simplified edge that an edge should be transformed into given its value.  This function is 
    heavily dependent on the SIMPLIFICATIONS dictionary inside edges.py, so as time goes on, do add more groupings and their associated 
    edges!

    If no grouping is found, then the original edge is returned
    '''
    for pairs in list(to_search.items()):
        if edge.higher() in pairs[1]:
            return pairs[0]
    return edge
