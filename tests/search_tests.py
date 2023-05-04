import unittest 
from ..app import app 

class SearchingTests(unittest.TestCase):
    def setUp(self):
        self.app = app 
    

if __name__ == '__main__':
    unittest.main()