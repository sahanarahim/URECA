from flask import Flask, render_template, request, redirect, url_for
import pickle

def generate_cytoscape_js(elements):
    nodes = [
        "{ data: { id: '%s' } }" % node.replace("'", r"\'")
        for node in set(edge["source"] for edge in elements) | set(edge["target"] for edge in elements)
    ]
    edges = [
        "{ data: { id: 'edge%s', source: '%s', target: '%s', interaction: '%s' } }" % (
            i,
            edge['source'].replace("'", r"\'"),
            edge['target'].replace("'", r"\'"),
            edge['interaction'].replace("'", r"\'"),
        )
        for i, edge in enumerate(elements)
    ]

    script = f"""
    document.addEventListener('DOMContentLoaded', function () {{
      const cy = cytoscape({{
        container: document.getElementById('cy'),
    
        elements: [
          // Nodes
          {', '.join(nodes)},
    
          // Edges
          {', '.join(edges)}
        ],
    
        style: [
          {{
            selector: 'node',
            style: {{
              'background-color': '#666',
              'label': 'data(id)',
              'font-size': '6px',
              'text-halign': 'center',
              'text-valign': 'center',
              'text-wrap': 'wrap',
              'text-max-width': 100
            }},
          }},
    
          {{
            selector: 'edge',
            style: {{
              'width': 3,
              'line-color': '#ccc',
              'target-arrow-color': '#ccc',
              'target-arrow-shape': 'triangle',
              'curve-style': 'bezier',
              'label': 'data(interaction)',
              'font-size': '6px',
              'text-rotation': 'autorotate',
              'text-margin-x': '5px',
              'text-margin-y': '-5px'
            }},
          }},
        ],
    
        layout: {{
          name: 'cose'
        }}
      }});
    }}
    );
    """
    return script





app = Flask(__name__)

class Gene:
    def __init__(self, id, description, inter_type, publication, frequent):
        self.id = id
        self.target = description
        self.inter_type = inter_type
        self.publication = publication
        self.frequent = frequent

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/search', methods=['POST'])
def search():
    with open('D:/GPT_annotator/allDic', 'rb') as file:
        genes = pickle.load(file)

    try:
        my_search = request.form['gene_id']
    except:
        my_search='cesa'
        
    forSending = []
    elements = []
    if len(my_search)>0:
        for i in genes:
            if my_search.upper() in i:
                for j in genes[i]:
                    for k in genes[i][j]:
                        if k[0]==1:
                            forSending.append(Gene(i, j, k[1], k[2], 'Yes'))
                            elements.append({"source": j, "target": i, "interaction": k[1]})
                        else:
                            forSending.append(Gene(i, j, k[1], k[2], 'No'))
                            elements.append({"source": j, "target": i, "interaction": k[1]})

            
    cytoscape_js_code = generate_cytoscape_js(elements)
    if forSending!=[]:
        return render_template('gene.html', genes=forSending, cytoscape_js_code=cytoscape_js_code)
    else:
        return render_template('not_found.html')



if __name__ == '__main__':
    app.run(debug=True)
