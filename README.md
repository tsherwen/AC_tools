# AC_tools - Atmospheric Chemistry (AC) tools

This repository contains a portable collection of functions used for 
working with global/regional chemical transport model (CTM) ouput and
observations.

The module is setup to be used as a submodule, with collections of
functions held in module files that can be imported in entirely or 
seperately.
 
e.g. 

```python
import AC_tools as AC

AC.get_GC_output()
```


Example analysis code for GEOS-Chem using AC_tools is available in the 
scripts folder.

For more infomation, please visit the wiki: https://github.com/tsherwen/AC_tools/wiki


INSTALL
=======

```bash
mkdir -p $HOME/python
cd $HOME/python
git clone --recursive https://github.com/tsherwen/AC_tools/
```
