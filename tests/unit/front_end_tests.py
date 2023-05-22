import unittest
import sys 
import os

# == HACKY, BUT IT WORKS ==
# Just to circumvent that ImportError: "attempted relative 
# import with no known parent package" error 
sys.path.append(os.getcwd())

from app import app

class FrontEndTests(unittest.TestCase): 
    def setUp(self):
        app.testing = True
        self.app = app.test_client()

    def test_home_page_should_be_accessible_and_have_right_contents(self):
        res = self.app.get('/')
        self.assertEqual(res.status_code, 200)
        self.assertIn(b"PlantConnectome", res.data)
        self.assertIn(b"Enter a gene, molecule, compartment, stress, cell type, organ, or any other related term", res.data)
        self.assertIn(b"Enter author name (without special characters, e.g., accents)", res.data)
        self.assertIn(b"Enter paper pubmed ID, separated by semicolon", res.data)
        self.assertIn(b"Search Tips:", res.data)
        self.assertIn(b"Future:", res.data)
        self.assertIn(b"Processed journals and abstracts", res.data)
    
    def test_help_page_should_be_accessible_and_have_right_contents(self):
        res = self.app.get('/help')
        self.assertEqual(res.status_code, 200)
        self.assertIn(b'How was the database constructed?', res.data)
        self.assertIn(b'How can the database be searched?', res.data)
        self.assertIn(b'Which papers were analyzed?', res.data)
        self.assertIn(b'How to use the KnowledgeNetwork viewer?', res.data)
        self.assertIn(b'Can I see the paper abstracts that revealed a given association?', res.data),
        self.assertIn(b'How accurate is GPT in extracting the information from abstracts?', res.data)

if __name__ == '__main__':
    unittest.main()