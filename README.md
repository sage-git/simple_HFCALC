# Simple Quantum Chemstry Calculation

Refactoring of the sample program of "MODERN QUANTUM CHEMISTRY" (Szabo and Ostlund)

## About

Simply changed the sample program to fortran 90/95 style.
* Replaced the COMMON blocks with modules
* divied into some separated source files
* introduced variable declearations due to `implicit none`
* Eliminated line numbers

## MEMO

Original codes are found in the book or here:

http://www.ccl.net/cca/software/SOURCES/FORTRAN/szabo/index.html

## DEMO

```
scons
./szabo_test > result.log
```

When debug mode is needed:
```
scons debug=1
```

If you want to use the intrinsic error function of gfortran:

```
scons my_erf=0
```

With the gfortran error function, results will differ a little bit.