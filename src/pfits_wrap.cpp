#include "pfits_wrap.h"

#define QUOTE(a) # a
#define CHKSTATUS(status,err_text,rv) \
    if (status != 0) { \
        fits_get_errstatus(status, err_text); \
        PyErr_Format(PyExc_ValueError, "CFITSIO: %s", err_text); \
        return rv;}
#define CHKNULL(var) \
    if (var == NULL) { \
        PyErr_Format(PyExc_RuntimeError, "Failed to create %s", QUOTE(var)); \
        return NULL;}

/*____*___*_____*____*******************************************************
|  ___|_ _|_   _/ ___| 
| |_   | |  | | \___ \ 
|  _|  | |  | |  ___) |
|_|   |___| |_| |____/ 
**************************************************************************/

static PyObject *FITSObject_new(PyTypeObject *type,
        PyObject *args, PyObject *kwds) {
    FITSObject *self;
    self = (FITSObject *) type->tp_alloc(type, 0);
    if (self != NULL) {
        self->fptr = NULL;
    }
    return (PyObject *) self;
}

// Deallocate memory when Python object is deleted
static void FITSObject_dealloc(FITSObject* self) {
    int status = 0;
    if (self->fptr != NULL) fits_close_file(self->fptr, &status);
    self->ob_type->tp_free((PyObject*)self);
}

// Initialize object (__init__)
static int FITSObject_init(FITSObject *self, PyObject *args) {
    int status = 0;
    char *filename, readwrite='r', err_text[FLEN_ERRMSG];
    if (!PyArg_ParseTuple(args, "s|c", &filename, &readwrite)) return -1;
    switch (readwrite) {
        case 'r':
            fits_open_file(&(self->fptr), filename, READONLY, &status);
            break;
        case 'w':
            fits_open_file(&(self->fptr), filename, READWRITE, &status);
            break;
        case 'c':
           fits_create_file(&(self->fptr), filename, &status);
            break;
        default:
            PyErr_Format(PyExc_ValueError, "mode must be 'r','c', or 'w'");
            return -1;
    }
    CHKSTATUS(status,err_text,-1);
    return 0;
}

/*        _     _       _         
 __ _ ___| |_  | |_  __| |_  _ ___
/ _` / -_)  _| | ' \/ _` | || (_-<
\__, \___|\__|_|_||_\__,_|\_,_/__/
|___/       |___|                  */
// Return a list of hdus in this FITS file
PyObject *FITSObject_get_hdus(FITSObject *self) {
    int status = 0, hdunum, i;
    char err_text[FLEN_ERRMSG];
    PyObject *rv;
    HDUObject *hdu;
    fits_get_num_hdus(self->fptr, &hdunum, &status);
    CHKSTATUS(status,err_text,NULL);
    rv = PyTuple_New(hdunum);
    CHKNULL(rv);
    for (i=0; i < hdunum; i++) {
        hdu = (HDUObject *)HDUObject_NEWC(&HDUType, self->fptr, i+1, self);
        if (hdu == NULL) {
            Py_DECREF(rv);
            CHKNULL(hdu);
        }
        // HDUs can't exist without the FITS file, so incref for every HDU
        Py_INCREF(self);
        fits_movabs_hdu(hdu->fptr, hdu->hdunum, &(hdu->hdutype), &status);
        PyTuple_SET_ITEM(rv, i, (PyObject *)hdu);
    }
    if (status != 0) {
        Py_DECREF(rv);   // This will also DECREF all consituent hdus
        CHKSTATUS(status,err_text,NULL);
    }
    return rv;
}

/*__ _ _                                  _           
 / _(_) |_ ___ __ __ ___ _ __ _ _ __ _ __(_)_ _  __ _ 
|  _| |  _(_-< \ V  V / '_/ _` | '_ \ '_ \ | ' \/ _` |
|_| |_|\__/__/  \_/\_/|_| \__,_| .__/ .__/_|_||_\__, |
                               |_|  |_|         |___/  */
// Bind methods to object
static PyMethodDef FITSObject_methods[] = {
    {"get_hdus", (PyCFunction)FITSObject_get_hdus, METH_NOARGS,
        "get_hdus()\nReturn the number of FITSs in the file."},
    {NULL}  // Sentinel
};

PyTypeObject FITSType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "pfits.FITS", /*tp_name*/
    sizeof(FITSObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)FITSObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "This class provides the basic interfaces to FITS files.  FITS(filename)",       /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    FITSObject_methods,     /* tp_methods */
    0,                      /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)FITSObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    FITSObject_new,       /* tp_new */
};
/****_*___***_***_*******************************************************
| | | |  _ \| | | |
| |_| | | | | | | |
|  _  | |_| | |_| |
|_| |_|____/ \___/ 
*************************************************************************/

static PyObject *HDUObject_NEWC(PyTypeObject *type, 
        fitsfile *fptr, int hdunum, FITSObject *parent) {
    HDUObject *self;
    self = (HDUObject *) type->tp_alloc(type, 0);
    if (self != NULL) {
        self->fptr = fptr;
        self->hdunum = hdunum;
        self->parent = (PyObject *)parent;
        self->cards = NULL;
    }
    return (PyObject *) self;
}

// Allocate memory for Python object 
static PyObject *HDUObject_new(PyTypeObject *type, 
        PyObject *args, PyObject *kwds) {
    return HDUObject_NEWC(type, NULL, 0, NULL);
}

// Deallocate memory when Python object is deleted
static void HDUObject_dealloc(HDUObject* self) {
    if (self->parent != NULL) Py_DECREF(self->parent);
    if (self->cards != NULL) Py_DECREF(self->cards);
    self->ob_type->tp_free((PyObject*)self);
}

// Initialize object (__init__)
static int HDUObject_init(HDUObject *self) {
    PyErr_Format(PyExc_RuntimeError, "HDUs must be instantiated through a FITS file");
    return -1;
}

/*          _    _                  
 __ _ ___| |_ | |_ _  _ _ __  ___ 
/ _` / -_)  _||  _| || | '_ \/ -_)
\__, \___|\__|_\__|\_, | .__/\___|
|___/       |___|  |__/|_|          */
// Return the type of ths HDU (IMAGE_HDU, ASCII_TBL, BINARY_TBL)
PyObject *HDUObject_get_type(HDUObject *self, void *closure) {
    switch (self->hdutype) {
        case IMAGE_HDU: return PyString_FromString("IMAGE_HDU");
        case ASCII_TBL: return PyString_FromString("ASCII_TBL");
        case BINARY_TBL: return PyString_FromString("BINARY_TBL");
        default:
            PyErr_Format(PyExc_ValueError, "Invalid HDU type");
            return NULL;
    }
}

/*         _ _    _                  _    
| |__ _  _(_) |__| |  __ __ _ _ _ __| |___
| '_ \ || | | / _` | / _/ _` | '_/ _` (_-<
|_.__/\_,_|_|_\__,_|_\__\__,_|_| \__,_/__/
                  |___|                     */
// Generate a list of all header cards + attach to HDUObject->cards
int HDUObject_build_cards(HDUObject *self) {
    int status=0, i, nkeys, namelen;
    char card[FLEN_CARD], name[FLEN_KEYWORD], value[FLEN_VALUE];
    char comment[FLEN_COMMENT], err_text[FLEN_ERRMSG];

    if (self->cards != NULL) return 0;
    fits_movabs_hdu(self->fptr, self->hdunum, NULL, &status);
    fits_get_hdrspace(self->fptr, &nkeys, NULL, &status);
    CHKSTATUS(status,err_text,-1);
    self->cards = PyTuple_New(nkeys);
    for (i=0; i < nkeys; i++) {
        fits_read_record(self->fptr, i+1, card, &status);
        if (status != 0) card[0] = '\0';
        fits_get_keyname(card, name, &namelen, &status);
        fits_parse_value(card, value, comment, &status);
        if (status != 0) {
            PyTuple_SET_ITEM(self->cards, i, Py_BuildValue("s",card));
            status = 0;
        } else {
            PyTuple_SET_ITEM(self->cards, i, 
                Py_BuildValue("s#ss",name,namelen,value,comment));
        }
    }
    return 0;
}

/*        _                    _    
 __ _ ___| |_   __ __ _ _ _ __| |___
/ _` / -_)  _| / _/ _` | '_/ _` (_-<
\__, \___|\__|_\__\__,_|_| \__,_/__/
|___/       |___|                     */
// Return list of cards to the outside world
PyObject *HDUObject_get_cards(HDUObject *self, void *closure) {
    // Build cards if we don't have it already (build_cards won't build it 2x)
    if (HDUObject_build_cards(self) == -1) return NULL;
    // Return a copy of the card list we have built
    Py_INCREF(self->cards);
    return self->cards;
}

/*               
| |_____ _  _ ___
| / / -_) || (_-<
|_\_\___|\_, /__/
         |__/     */
// Return list of card titles
PyObject *HDUObject_keys(HDUObject *self) {
    PyObject *keys, *key;
    int nkeys, i;
    
    // Build cards if we don't have it already (build_cards won't build it 2x)
    if (HDUObject_build_cards(self) == -1) return NULL;
    nkeys = (int)PyTuple_Size(self->cards);
    keys = PyList_New(0);
    for (i=0; i < nkeys; i++) {
        key = PyTuple_GetItem(self->cards, (Py_ssize_t) i);
        if (PyTuple_Check(key) && PyTuple_Size(key) > 0) {
            key = PyTuple_GetItem(key, 0);
            if (PyString_Check(key) && PyString_Size(key) > 0) {
                Py_INCREF(key);
                PyList_Append(keys, key);
            }
        }
    }
    return keys;        
}

/*              _            
| |_  __ _ ___ | |_____ _  _ 
| ' \/ _` (_-< | / / -_) || |
|_||_\__,_/__/_|_\_\___|\_, |
            |___|       |__/  */
// Return whether or not key is a keyword in header
PyObject * HDUObject_has_key(HDUObject *self, PyObject *key) {
    PyObject *rv;
    int status=0;
    char value[FLEN_VALUE];

    // Make sure our key is a string of proper length
    if (!PyString_Check(key) || PyString_Size(key) >= FLEN_KEYWORD) {
        rv = Py_False;
    } else {
        if (HDUObject_build_cards(self) == -1) return NULL;
        fits_movabs_hdu(self->fptr, self->hdunum, NULL, &status);
        fits_read_key(self->fptr, TSTRING, PyString_AsString(key),
            value, NULL, &status);
        rv = Py_True;
        if (status != 0) rv = Py_False;
    }
    Py_INCREF(rv);
    return rv;
}

/*                    _           
 _ __  __ _ _ __ _ __(_)_ _  __ _ 
| '  \/ _` | '_ \ '_ \ | ' \/ _` |
|_|_|_\__,_| .__/ .__/_|_||_\__, |
           |_|  |_|         |___/  */
// Get the length of the header
Py_ssize_t HDUObject_len(HDUObject *self) {
    if (HDUObject_build_cards(self) == -1) return NULL;
    return PyTuple_Size(self->cards);
}

// Get an item (like hdu['NAXIS'])
PyObject * HDUObject_getitem(HDUObject *self, PyObject *key) {
    int status=0;
    char value[FLEN_VALUE], err_text[FLEN_ERRMSG];

    // Make sure our key is a string of proper length
    if (!PyString_Check(key) || PyString_Size(key) >= FLEN_KEYWORD) {
        PyErr_Format(PyExc_KeyError, "keyword not found in header");
        return NULL;
    }
    if (HDUObject_build_cards(self) == -1) return NULL;
    fits_movabs_hdu(self->fptr, self->hdunum, NULL, &status);
    fits_read_key(self->fptr, TSTRING, PyString_AsString(key),
        value, NULL, &status);
    // Special status check to return a KeyError
    if (status != 0) {
        fits_get_errstatus(status, err_text);
        PyErr_Format(PyExc_KeyError, err_text);
        return NULL;
    }
    return PyString_FromString(value);
}

// Build the mapping that gets tied to the HDU object initialization
static PyMappingMethods HDUObject_as_mapping = {
    (lenfunc)HDUObject_len,         /* mp_length */
    (binaryfunc)HDUObject_getitem,  /* mp_subscript */
    0,                              /* mp_ass_subscript */
};

/*  _       
 __| |_ _ _ 
(_-<  _| '_|
/__/\__|_|   */
// Return entire header as a string for str(hdu)
PyObject * HDUObject_str(HDUObject *self) {
    PyObject *rv;
    int status=0, ncards;
    char *hdrstr, err_text[FLEN_ERRMSG];

    fits_movabs_hdu(self->fptr, self->hdunum, NULL, &status);
    // Memory gets allocated for hdrstr in fits_hdr2str
    fits_hdr2str(self->fptr, 0, NULL, 0, &hdrstr, &ncards, &status);
    CHKSTATUS(status,err_text,NULL);
    rv = PyString_FromString(hdrstr);
    free(hdrstr);
    return rv;
}

/*            _       _                  _        _      _        
(_)_ __  __ _| |_  __| |_  _    __ _ ___| |_   __| |__ _| |_ __ _ 
| | '  \/ _` | ' \/ _` | || |  / _` / -_)  _| / _` / _` |  _/ _` |
|_|_|_|_\__, |_||_\__,_|\_,_|__\__, \___|\__|_\__,_\__,_|\__\__,_|
        |___/              |___|___/       |___|                   */

inline int get_numpy_dtype(int f_dtype, int *n_dtype) {
    switch (f_dtype) {
        case TBIT: case TLOGICAL: n_dtype[0] = NPY_BOOL;      return 0;   
        case TSBYTE:              n_dtype[0] = NPY_BYTE;      return 0;
        case TBYTE:               n_dtype[0] = NPY_UBYTE;     return 0;
        case TSHORT:              n_dtype[0] = NPY_SHORT;     return 0;
        case TUSHORT:             n_dtype[0] = NPY_USHORT;    return 0;
        case TINT:                n_dtype[0] = NPY_INT;       return 0;
        case TUINT:               n_dtype[0] = NPY_UINT;      return 0;
        case TLONG:               n_dtype[0] = NPY_LONG;      return 0;
        case TULONG:              n_dtype[0] = NPY_ULONG;     return 0;
        case TFLOAT:              n_dtype[0] = NPY_FLOAT;     return 0;
        case TLONGLONG:           n_dtype[0] = NPY_LONGLONG;  return 0;;
        case TDOUBLE:             n_dtype[0] = NPY_DOUBLE;    return 0;
        case TCOMPLEX:            n_dtype[0] = NPY_CFLOAT;    return 0;
        case TDBLCOMPLEX:         n_dtype[0] = NPY_CDOUBLE;   return 0;
        case TSTRING:             n_dtype[0] = NPY_STRING;    return 0;
        default:
            PyErr_Format(PyExc_ValueError,"Unknown FITS data type %d",f_dtype);
            return -1;
    }
}

// Read image from an HDU of type IMAGE_HDU
// Return a numpy array of the entire contents
PyObject *imghdu_get_data(HDUObject *self) {
    PyArrayObject *rv;
    int i, bitpix, status=0, ndim, f_dtype, n_dtype, anynul;
    long dims[NPY_MAXDIMS], fpixel[NPY_MAXDIMS], nelements=1;
    npy_intp npy_dims[NPY_MAXDIMS];
    char err_text[FLEN_ERRMSG], nulval[8];

    fits_movabs_hdu(self->fptr, self->hdunum, NULL, &status);
    fits_get_img_param(self->fptr, NPY_MAXDIMS, &bitpix, &ndim, dims, &status);
    // Get the data type of image *after* scaling etc.
    fits_get_img_equivtype(self->fptr, &bitpix, &status);
    CHKSTATUS(status,err_text,NULL);
    if (ndim == 0) { // Get out of here if there be no data
        Py_INCREF(Py_None);
        return Py_None;
    }
    for (i=0; i < ndim; i++) {
        if (dims[i] <= 0) {
            Py_INCREF(Py_None);
            return Py_None;
        }
        npy_dims[i] = dims[i];
        nelements *= dims[i];
        fpixel[i] = 1;
    }
    switch (bitpix) {
        case BYTE_IMG:
            f_dtype = TBYTE;
            ((char *)nulval)[0] = 0;
            break;
        case SHORT_IMG:
            f_dtype = TSHORT;
            ((short *)nulval)[0] = -1; 
            break;
        case LONG_IMG:
            f_dtype = TLONG;
            ((long *)nulval)[0] = -1;
            break;
        case LONGLONG_IMG:
            f_dtype = TLONGLONG;
            ((long long *)nulval)[0] = -1;
            break;
        case FLOAT_IMG:
            f_dtype = TFLOAT;
            ((float *)nulval)[0] = NAN;
            break;
        case DOUBLE_IMG:
            f_dtype = TDOUBLE;
            ((double *)nulval)[0] = NAN; 
            break;
        default:
            PyErr_Format(PyExc_ValueError, 
                "Unrecognized FITS image data format %d", bitpix);
            return NULL;
    }
    get_numpy_dtype(f_dtype,&n_dtype);
    rv = (PyArrayObject *) PyArray_SimpleNew(ndim, npy_dims, n_dtype);
    CHKNULL(rv);
    fits_read_pix(self->fptr, f_dtype, fpixel, nelements, nulval, 
        rv->data, &anynul, &status);
    if (status != 0) {
        Py_DECREF(rv);
        CHKSTATUS(status,err_text,NULL);
    }
    return PyArray_Return(rv);
}

/*      _             _        _      _        
 __ ___| |   __ _ ___| |_   __| |__ _| |_ __ _ 
/ _/ _ \ |  / _` / -_)  _| / _` / _` |  _/ _` |
\__\___/_|__\__, \___|\__|_\__,_\__,_|\__\__,_|
        |___|___/       |___|                   */

#define MKNULVAL(nulval,type,nvals,ii,val) \
    nulval = (char *) malloc(nvals*sizeof(type));\
    for (ii=0; ii < nvals; ii++) { ((type *)nulval)[ii] = val; }

char *malloc_nulval(int n_dtype, long nvals) {
    char *nulval=NULL;
    long ii;
    switch (n_dtype) {
        case NPY_BOOL: case NPY_BYTE: case NPY_UBYTE:
            MKNULVAL(nulval,char,nvals,ii,'\x00'); break;
        case NPY_SHORT: case NPY_USHORT:
            MKNULVAL(nulval,short,nvals,ii,0); break;
        case NPY_INT: case NPY_UINT:
            MKNULVAL(nulval,int,nvals,ii,0); break;
        case NPY_LONG: case NPY_ULONG:
            MKNULVAL(nulval,long,nvals,ii,0); break;
        case NPY_FLOAT:
            MKNULVAL(nulval,float,nvals,ii,NAN); break;
        case NPY_LONGLONG:
            MKNULVAL(nulval,long long,nvals,ii,0); break;
        case NPY_DOUBLE:
            MKNULVAL(nulval,double,nvals,ii,NAN); break;
        case NPY_CFLOAT:
            MKNULVAL(nulval,float,2*nvals,ii,NAN); break;
        case NPY_CDOUBLE:
            MKNULVAL(nulval,double,2*nvals,ii,NAN); break;
        case NPY_STRING:
            MKNULVAL(nulval,char,nvals,ii,' '); nulval[nvals-1] = '\0'; break;
    }
    return nulval;
}

// Return a numpy array containing the data of a table column
PyObject *col_get_data(HDUObject *self, int colnum, long nrows) {
    PyArrayObject *rv;
    PyObject *rvtup;
    int status=0, f_dtype, n_dtype, anynul, ndim;
    long repeat, offset, width, ii, jj, naxes[NPY_MAXDIMS];
    npy_intp npy_dims[NPY_MAXDIMS];
    char err_text[FLEN_ERRMSG], *nulval, **buf=NULL;

    // Get the data type of this column
    fits_get_eqcoltype(self->fptr,colnum+1,&f_dtype,&repeat,&width,&status);
    fits_read_tdim(self->fptr, colnum+1, NPY_MAXDIMS, &ndim, naxes, &status);
    CHKSTATUS(status,err_text,NULL);

    //printf("colnum:%d,dtype:%d\n",colnum+1,f_dtype);
    // Check if this is a variable length column
    if (f_dtype < 0) {
        f_dtype = -f_dtype;
        // Match column data type with numpy data type
        if (get_numpy_dtype(f_dtype, &n_dtype) != 0) return NULL;
        ndim = 1;
        // Numpy doesn't do var length arrays, so we'll make a tuple instead
        rvtup = PyTuple_New(nrows);
        for (ii=0; ii < nrows; ii++) {
            // Get length of this particular row (offset is unnecessary)
            fits_read_descript(self->fptr,colnum+1,ii+1,
                &repeat,&offset,&status);
            //printf("row:%ld,repeat:%ld,offset:%ld\n", ii+1,repeat, offset);
            if (status != 0) {
                Py_DECREF(rvtup);
                CHKSTATUS(status,err_text,NULL)
            };
            if (f_dtype == TSTRING) {
                nulval = malloc_nulval(n_dtype, width+1);
                if (nulval == NULL) { Py_DECREF(rvtup); CHKNULL(nulval); }
                npy_dims[0] = 1;
                rv = (PyArrayObject *) PyArray_New(&PyArray_Type,ndim,npy_dims,
                    n_dtype, NULL, NULL, width+1, 0, NULL);
                if (rv == NULL) { free(nulval); Py_DECREF(rvtup); CHKNULL(rv); }
                // Fill buf with string terminators b/c fits_read_col 
                // stops writing at first term, but numpy backtracks from
                // end to determine str size.
                for (jj=0; jj < width+1; jj++) { (rv->data)[jj] = '\0'; }
                // Read data into the string buffer
                fits_read_col(self->fptr, TSTRING, colnum+1, ii+1, 1,
                    1, nulval, &(rv->data), &anynul, &status);
                free(nulval);
            } else {
                nulval = malloc_nulval(n_dtype, repeat);
                if (nulval == NULL) { Py_DECREF(rvtup); CHKNULL(nulval); }
                npy_dims[0] = repeat;
                rv = (PyArrayObject *) PyArray_SimpleNew(ndim,npy_dims,n_dtype);
                if (rv == NULL) { free(nulval); Py_DECREF(rvtup); CHKNULL(rv); }
                fits_read_col(self->fptr, f_dtype, colnum+1, ii+1, 1, 
                    repeat, nulval, rv->data, &anynul, &status);
                free(nulval);
            }
            if (status != 0) {
                Py_DECREF(rvtup); Py_DECREF(rv);
                CHKSTATUS(status,err_text,NULL)
            };
            PyTuple_SET_ITEM(rvtup, ii, PyArray_Return(rv));
        }
        return rvtup;
    } else {
        // Match column data type with numpy data type
        if (get_numpy_dtype(f_dtype, &n_dtype) != 0) return NULL;

        if (f_dtype == TSTRING) {
            ndim == 1;
            npy_dims[0] = nrows;
            //printf("width:%d\n", width);
            rv = (PyArrayObject *) PyArray_New(&PyArray_Type, ndim, npy_dims, 
                n_dtype, NULL, NULL, width+1, 0, NULL);
            CHKNULL(rv);
            buf = (char **) malloc(nrows*sizeof(char *));
            if (buf == NULL) { Py_DECREF(rv); CHKNULL(buf); }
            nulval = malloc_nulval(n_dtype, width+1);
            if (nulval == NULL) { free(buf); Py_DECREF(rv); CHKNULL(nulval); }
            for (ii=0; ii < nrows; ii++) {
                buf[ii] = rv->data + ii*rv->strides[0];
                // Fill buf with string terminators b/c fits_read_col 
                // stops writing at first term, but numpy backtracks from
                // end to determine str size.
                for (jj=0; jj < width+1; jj++) { buf[ii][jj] = '\0'; }
            }
            // Read data into the string buffer
            fits_read_col(self->fptr, TSTRING, colnum+1, 1, 1,
                nrows, nulval, buf, &anynul, &status);
            free(buf); free(nulval);
        } else {
            // Final axis is implicitly nrows
            // If first dim is 1, don't bother making it an axis
            if (ndim == 1 && naxes[0] == 1) {
                npy_dims[0] = nrows;
            } else {
                for (ii=0; ii < ndim; ii++) { npy_dims[ii] = naxes[ii]; }
                npy_dims[ndim++] = nrows;
            }
            //printf("ndim:%d,n_dtype:%d\n", ndim, n_dtype);
            //for (ii=0; ii<ndim; ii++) {printf("dim%d:%d\n",ii,npy_dims[ii]);}
            rv = (PyArrayObject *) PyArray_SimpleNew(ndim, npy_dims, n_dtype);
            CHKNULL(rv);
            nulval = malloc_nulval(n_dtype, repeat);
            if (nulval == NULL) { Py_DECREF(rv); CHKNULL(nulval); }
            // Read the data from the FITS column
            fits_read_col(self->fptr, f_dtype, colnum+1, 1, 1,
                (long) PyArray_SIZE(rv), nulval, rv->data, &anynul, &status);
            free(nulval);
        }
        if (status != 0) {
            Py_DECREF(rv);
            CHKSTATUS(status,err_text,NULL);
        }
        return PyArray_Return(rv);
    }
}
 
/*   _    _ _       _                  _        _      _        
| |_| |__| | |_  __| |_  _    __ _ ___| |_   __| |__ _| |_ __ _ 
|  _| '_ \ | ' \/ _` | || |  / _` / -_)  _| / _` / _` |  _/ _` |
 \__|_.__/_|_||_\__,_|\_,_|__\__, \___|\__|_\__,_\__,_|\__\__,_|
                         |___|___/       |___|                   */

// Read all columns from an HDU of type BINARY_TBL
// Return a tuple of the contents of each column
PyObject *tblhdu_get_data(HDUObject *self) {
    PyObject *rv, *item;
    int status=0, ncols, colnum;
    long nrows;
    char err_text[FLEN_ERRMSG], colname[FLEN_VALUE], keyname[FLEN_KEYWORD];

    fits_get_num_rows(self->fptr, &nrows, &status);
    fits_get_num_cols(self->fptr, &ncols, &status);
    CHKSTATUS(status,err_text,NULL);
    rv = PyDict_New(); CHKNULL(rv);

    // Loop over all columns collecting data
    for (colnum=0; colnum < ncols; colnum++) {
        // Lookup the column name for the record array
        sprintf(keyname, "TTYPE%d", colnum+1);
        fits_read_key(self->fptr, TSTRING, keyname, colname, NULL, &status);
        if (status != 0) { Py_DECREF(rv); CHKSTATUS(status,err_text,NULL); }
        item = col_get_data(self,colnum,nrows);
        if (item == NULL) { Py_DECREF(rv); return NULL; }
        PyDict_SetItemString(rv, colname, item);
    }
    return rv;
}

/*        _        _      _        
 __ _ ___| |_   __| |__ _| |_ __ _ 
/ _` / -_)  _| / _` / _` |  _/ _` |
\__, \___|\__|_\__,_\__,_|\__\__,_|
|___/       |___|                   */
// Get the data in this HDU (IMAGE_HDU, BINARY_TBL, or ASCII_TBL)
PyObject *HDUObject_get_data(HDUObject *self) {
    int status=0;
    char err_text[FLEN_ERRMSG];

    fits_movabs_hdu(self->fptr, self->hdunum, NULL, &status);
    CHKSTATUS(status,err_text,NULL);
    // Figure out which HDU Type this is and act accordingly
    switch (self->hdutype) {
        // In an IMAGE_HDU, read image into a numpy array
        case IMAGE_HDU: return imghdu_get_data(self);
        // In a BINARY_TBL or ASCII_TBL, read data into a numpy array, unless 
        // it's string data, in which case return a tuple of strings.
        case BINARY_TBL: case ASCII_TBL: return tblhdu_get_data(self);
        default:
            PyErr_Format(PyExc_ValueError, "Unsupported HDU type");
            return NULL;
    }
}
            
/*       _                                  _           
| |_  __| |_  _  __ __ ___ _ __ _ _ __ _ __(_)_ _  __ _ 
| ' \/ _` | || | \ V  V / '_/ _` | '_ \ '_ \ | ' \/ _` |
|_||_\__,_|\_,_|  \_/\_/|_| \__,_| .__/ .__/_|_||_\__, |
                                 |_|  |_|         |___/  */
// Bind methods to object
static PyMethodDef HDUObject_methods[] = {
    {"keys", (PyCFunction)HDUObject_keys, METH_NOARGS,
        "keys()\nReturn a list of keywords in this HDU's header."},
    {"has_key", (PyCFunction)HDUObject_has_key, METH_O,
        "has_key(key)\nReturn if key is a keyword in this HDU's header."},
    {"get_data", (PyCFunction)HDUObject_get_data, METH_NOARGS,
        "get_data()\nReturn this HDU's data as a tuple with an entry for the data in each column.  Each column may contain a numpy array or a tuple of strings."},
    {NULL}  // Sentinel
};

static PyGetSetDef HDUObject_getset[] = {
    {"type", (getter)HDUObject_get_type},
    {"cards", (getter)HDUObject_get_cards},
    {NULL}  // Sentinel
};

PyTypeObject HDUType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "pfits.HDU", /*tp_name*/
    sizeof(HDUObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)HDUObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    &HDUObject_as_mapping,     /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    (reprfunc)HDUObject_str,   /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,        /*tp_flags*/
    "This class provides the basic interfaces to HDUs in FITS files.  HDU(filename)",       /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    HDUObject_methods,     /* tp_methods */
    0,                      /* tp_members */
    HDUObject_getset,       /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)HDUObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    HDUObject_new,       /* tp_new */
};

/*_**__***********_*******_**************************************************
|  \/  | ___   __| |_   _| | ___ 
| |\/| |/ _ \ / _` | | | | |/ _ \
| |  | | (_) | (_| | |_| | |  __/
|_|  |_|\___/ \__,_|\__,_|_|\___|
***************************************************************************/
// Module methods
static PyMethodDef pfits_methods[] = {
    {NULL}  /* Sentinel */
};

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

// Module init
PyMODINIT_FUNC initpfits(void) {
    PyObject* m;
    FITSType.tp_new = PyType_GenericNew;
    HDUType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&FITSType) < 0 || PyType_Ready(&HDUType) < 0) return;
    m = Py_InitModule3("pfits", pfits_methods,
    "This is a hand-written Python wrapper.");
    import_array();
    Py_INCREF(&FITSType);
    PyModule_AddObject(m, "FITS", (PyObject *)&FITSType);
    Py_INCREF(&HDUType);
    PyModule_AddObject(m, "HDU", (PyObject *)&HDUType);
}


