'''
This module contains helper functions 
needed to generate the network graph.
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
    
    a = open('../infos/network.txt','r').read()

    nodes = ', '.join(nodes)
    edges = ', '.join(edges)
    
    return a.replace('NODES',nodes).replace('EDGES',edges)