from flask_frozen import Freezer
from app import app

f = Freezer(app)

f.freeze()