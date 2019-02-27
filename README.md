
# Pygssw
A simple Python wrapper around [gssw](https://github.com/vgteam/gssw), enabling sequence to graph alignment from Python. Only tested with Python 3.

## Install
Pygssw can installed with pip:
```
pip3 install pygssw
```

# Usage
```python
from pygssw import align
nodes = [1, 2, 3] 
edges = [(1, 2), (2, 3)]
node_sequences = ["AAA", "CCC", "TTT"]  # Same order as nodes list
read = "AAACTCTTT"

alignment, score = align(nodes, node_sequences, edges, read)
print(alignment, score)
```

Note: It seems that GSSW only supports dense sorted nodes, and no edges going from a higher node id to a lower. 
Pygssw attempts to converting all node IDs before calling GSSW. After aligning, it converts the node IDs back.



## Installing from source
This guide assumes linux.
* Clone this repository
* Install requirements:
```bash
apt install swig
apt install python3-dev
```
* Install:
```bash
cd pygssw   # Folder pygssw inside the repo
swig -python gssw.i
gcc -fpic  -c -O3 -msse4 gssw.c gssw_wrap.c -I/usr/include/python3.6m
ld -shared gssw.o gssw_wrap.o -o _gssw.so
cd ..
pip3 install -e .
```
NOTE: Change /usr/include/python3.6m to python-dev executable (requires intallation of python-dev, i.e. `apt-get install python3-dev`)


