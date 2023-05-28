import unittest
import os 
import sys 

# Just to circumvent that ImportError: "attempted relative 
# import with no known parent package" error 
# (a bit hacky though)
sys.path.append(os.getcwd())

# Importing user modules:
import utils.search as e
from app import app

class BackEndTestCases(unittest.TestCase):
    pass

if __name__ == '__main__':
    unittest.main()
        