# Testing `TRUSS`

I used AI to write a little test script based on what is in the design documentation markdown file, specifically using the sample input and output. This script can test any executable (including other scripts)

## Using `test_truss.py`

Using it is as simple as running the test program against an executable

```
cd TRUSS/src
gfortran -o trusstran truss.f90
../test/test_truss.py trusstran
```
