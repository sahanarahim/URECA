import unittest
import os 
import sys 

# Just to circumvent that ImportError: "attempted relative 
# import with no known parent package" error 
# (a bit hacky though)
sys.path.append(os.getcwd())

# Importing user modules:
import utils.edges as e
from app import app

class BackEndTestCases(unittest.TestCase):
    def setUp(self):
        app.testing = True
        self.app = app
        self.example_values = [('Tomato levels', 'INCREASED BY', 'Ketchup'), 
                               ('Mustard', 'CAUSES', 'indigestion'),
                               ('Chocolate cravings', 'RepResseD BY', 'wilLpower'),
                               ('GeNOType', 'MODulated by', 'stuff')]
        self.example_groupings = ['temperatuREs', 'neW', 'seQUences', 'FUnctional', 'blEh', 'fully']
        
    def test_search_passive_edges_should_return_right_edge(self):
        answers = ['INCREASES', None, 'REPRESSES', 'MODULATES']
        for i in range(len(self.example_values)):
            self.assertEqual(e.search_passive_edges(self.example_values[i][1]), answers[i])
    
    def test_find_grouping_should_return_right_group(self):
        answers = ['CLUSTER 19', 'CLUSTER 4', 'CLUSTER 12', 'CLUSTER 16', '???', 'CLUSTER 12']
        for i in range(len(self.example_groupings)):
            self.assertEqual(e.find_grouping(self.example_groupings[i]), answers[i])

if __name__ == '__main__':
    unittest.main()
        