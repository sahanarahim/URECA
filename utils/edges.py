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
    ('OVEREXPRESSED IN') : 'OVEREXPRESSES'
}

def search_passive_edges(edge, to_search):
    '''
    Given an 'edge', this function tries to search for an edge inside to_search.

    If a key is found, then the upper-cased version of that key is returned.  
    Otherwise, a 'None' is returned.
    '''
    for i in list(to_search.items()):
        if edge.upper() in [j.upper() for j in i[0]]:
            return i[1].upper()
    return None