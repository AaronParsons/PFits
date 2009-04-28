#ifndef FITSIO_WRAP_H
#define FITSIO_WRAP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Python.h>
#include "numpy/arrayobject.h"
#include "fitsio.h"

// Python object that holds handle to a FITS file
typedef struct {
    PyObject_HEAD
    fitsfile *fptr;
} FITSObject;

extern PyTypeObject FITSType;

// Python object that holds handle to an HDU 
typedef struct{
    PyObject_HEAD
    fitsfile *fptr;
    int hdunum;
    int hdutype;
    PyObject *parent;
    PyObject *cards;
} HDUObject;

extern PyTypeObject HDUType;
static PyObject *HDUObject_NEWC(PyTypeObject *, fitsfile *, int, FITSObject *);


#endif
