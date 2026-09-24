# conftest.py  (project root, next to app.py)
import sys
import os

# make sure the root and create_graph are importable
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "create_graph"))