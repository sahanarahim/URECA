from flask import Flask, render_template

# == IMPORTING THE ROUTES  ==
from routes.author_search import author_search
from routes.term_searches import term_searches
from routes.title_searches import title_searches

app = Flask(__name__)

# == REGISTERING THE ROUTES ==
app.register_blueprint(author_search)
app.register_blueprint(term_searches)
app.register_blueprint(title_searches)

# == Fetching the 'index', the 'help', and the 'features' page ==
'''
I left these routes here as they're trivial enough for now.
'''

@app.route('/', methods = ['GET'])
def index():
    v = open('stats.txt','r').read().rstrip().split()
    return render_template('index.html', entities = v[1], papers = v[0])

@app.route('/help', methods = ['GET'])
def help():
    '''
    Returns the help page.
    '''
    return render_template('help.html')

@app.route('/features', methods = ['GET'])
def features():
    '''
    Renders the features page.
    '''
    journals, numbers = open('journal_statistics.txt','r').read().splitlines()
    piechart = open('piechart.txt','r').read().replace('JOURNALS', journals).replace('NUMBERS', numbers)
    return render_template('features.html', piechart_code = piechart)

if __name__ == '__main__':
    app.run()