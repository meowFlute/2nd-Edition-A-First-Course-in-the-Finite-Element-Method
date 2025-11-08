# A First Course in the Finite Element Method, 2nd Edition

## What is this repo?
Every winter I buy a few textbooks to explore. Winter of 2025 into 2026 one of the ones I've chosen is a 1992 copy of _A First Course in the Finite Element Method_, which happens to be the 2nd Edition **AND IT COMES WITH A 5.25 INCH FLOPPY DISK** so figuring out a way to read that will be part of the fun.

### What is on the disk?
Back in 1992 the programming language of science was Fortran, and as a result the book is all about Fortran. The floppy disk is Fortran programs and the book references them extensively.

## Loose plan
I've never really used Fortran, but I do have experience in C and python. I plan on creating modern versions of the programs from the book in all 3 languages. Maybe after this I jump into the NASTRAN-related projects or something, this should be an approachable jumping-off point

## Repo Layout

### `/` the root level of the repo
This level will contain any necessary build files. I tend to use either vanilla makefiles or CMake to build C and I'm guessing I'll do the same with Fortran. They'll probably autogenerate a build folders that I ignore explicitly in the .gitignore file

### The `PROGRAM` folders
This level will also contain 6 folders representing the programs that shipped on the 5.25 inch floppy:
- `TRUSS` - can be used to solve two- and three-dimensional truss problems (Chapter 4)
- `PFRAME` - can be used to solve large two-dimensional rigid plane frame and grid problems (Chapter 6)
- `CSFEP` - can be used to solve large plan stress/strain problems (Chapter 8)
- `HEAT` - can be used to solve two-dimensional heat-transfer (without mass transport) problems (Chapter 13)
- `FLUID` - can be used to solve two-dimensional fluid-flow (without rotation, compressibility, and only in steady state) problems (Chapter 14) 
- `DFRAME` - can be used to solve time-dependent plane frame problems (Chapter 16)

In each of these folders, you will find the following folder structure

#### `src/` the source folder
This is where the Python, C, and Fortran files will reside which can be compiled and/or run to execute the program and solve problems

#### `test/` the test folder
There are tons of fully-worked-out example problems in this book, and I plan to use these to create test cases. Will try to use these to make it so that each language has identical functionality for each program to what is described in the book.

#### `doc/` the documentation folder
I don't plan on making super thorough documentation (buy the book!) but some markdown files on usage, design notes, etc. tend to be useful to me so I can figure out what the heck I was doing when I come back to this in several years to admire stuff I forgot about.

#### `build/` the build folder (not checked in, but generated during build)
When you run the build process the compiled and linked executables will be put here along with any intermediate/generated files
