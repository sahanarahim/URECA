import unittest 
from app import app

class FrontEndTestCases(unittest.TestCase):
    '''
    These are some unit test cases that Kevin has came up with to ensure that everything on 
    the website that should appear appears.  Like the class name suggests, this class 
    basically contains test cases for the front-end aspect of the database - so, it attempts 
    to verify that the page appears as it should be regardless if somebody changes the
    application's templates!

    Feel free to add in more test cases as time goes on, but DO NOT change that "setUp" method!
    '''
    def setUp(self):
        app.testing = True 
        self.app = app.test_client()

    def test_home_page_should_be_accessible(self):
        res = self.app.get('/')
        self.assertEqual(res.status_code, 200)
    
    def test_home_page_should_contain_right_contents(self):
        res = self.app.get('/')
        self.assertIn(b"PlantConnectome", res.data)
        self.assertIn(b"Enter a gene, molecule, compartment, stress, cell type, organ, or any other related term", res.data)
        self.assertIn(b"Enter author name (without special characters, e.g., accents)", res.data)
        self.assertIn(b"Enter paper pubmed ID, separated by semicolon", res.data)
        self.assertIn(b"Search Tips:", res.data)
        self.assertIn(b"Future:", res.data)
        self.assertIn(b"Processed journals and abstracts", res.data)
        pass
    
    def test_help_page_should_be_accessible(self):
        res = self.app.get('/help')
        self.assertEqual(res.status_code, 200)
    
    def test_help_page_should_contain_right_contents(self):
        res = self.app.get('/help')
        self.assertIn(b'How was the database constructed?', res.data)
        self.assertIn(b'How can the database be searched?', res.data)
        self.assertIn(b'Which papers were analyzed?', res.data)
        self.assertIn(b'How to use the KnowledgeNetwork viewer?', res.data)
        self.assertIn(b'Can I see the paper abstracts that revealed a given association?', res.data),
        self.assertIn(b'How accurate is GPT in extracting the information from abstracts?', res.data)
        pass



if __name__ == '__main__':
    unittest.main()