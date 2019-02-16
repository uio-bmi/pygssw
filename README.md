
# pygssw
Simple python wrapper around GSSW.

## Install
(Assumes linux)
Install requirements first:
```bash
apt install swig
apt install python3-dev

```
Then install:
```bash
cd pygssw   # Folder pygssw inside the repo
swig -python gssw.i
gcc -fpic  -c -O3 -msse4 gssw.c gssw_wrap.c -I/usr/include/python3.6m
ld -shared gssw.o gssw_wrap.o -o _gssw.so
cd ..
pip3 install -e .
```
NOTE: Change /usr/include/python3.6m to python-dev executable (requires intallation of python-dev, i.e. `apt-get install python3-dev`)

# Usage
```python
from pygssw.align import align
nodes = [1, 2, 3] 
edges = [(1, 2), (2, 3)]
node_sequences = ["AAA", "CCC", "TTT"]  # Same order as nodes list
read = "AAACTCTTT"

alignment, score = align(nodes, node_sequences, edges, read)
print(alignment, score)
```



