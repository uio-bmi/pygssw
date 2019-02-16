
# pygssw
Simple python wrapper around GSSW.

## Install
```bash
cd pygsswswig -python gssw.i
swig -python gssw.i
gcc -fpic  -c -O3 -msse4 gssw.c gssw_wrap.c -I/usr/include/python3.6m
ld -shared gssw.o gssw_wrap.o -o _gssw.so
```

Change /usr/include/python3.6m to python-dev executable (requires intallation of python-dev, i.e. `apt-get install python3.dev`)

