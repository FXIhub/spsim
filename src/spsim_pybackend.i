%module spsim_pybackend
%{
  /* Includes the header in the wrapper code */
#include <spimage.h>
#include <spimage/image.h>
#include <spimage/image_util.h>
#include "../include/spsim.h"
#include "../include/config.h"
#include "../include/diffraction.h"
#include "../include/molecule.h"
#include <numpy/arrayobject.h>
  //#include <numpy.i>
  //#include <stdio.h>
  %}

%init %{
        import_array();
%}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------
// COPIED FROM LIBSPIMAGE SPIMAGE_PYBACKEND.I

/* Parse the header file to generate wrappers */
%typemap(out) Complex {
    $result = PyComplex_FromDoubles($1.re,$1.im);
}
/* Parse the header file to generate wrappers */
%typemap(in) Complex {
    $1.re = PyComplex_RealAsDouble($input);
    $1.im = PyComplex_ImagAsDouble($input);
}


%typemap(out) sp_i3matrix * {
  if($1->z == 1){
    /* Swap the order of the dimensions so we can plot things easily in imshow.
       This is bound to create confusion at some point in the future. */
    npy_intp dims[2] = {$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(2, dims, PyArray_INT, $1->data);
  }else{
    npy_intp dims[3] = {$1->z,$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(3, dims, PyArray_INT, $1->data);
  }
}

%typemap(in) sp_i3matrix * {
  /* Make sure input is a NumPy array and an Int */
  PyArrayObject *arr;
  if(PyArray_Check($input) == 0 || PyArray_ISINTEGER($input) == 0){
    PyErr_SetString( PyExc_TypeError, "not an int array" );
    return NULL;
  } 
  arr = (PyArrayObject *)($input);
  //  if ((arr = (PyArrayObject *) PyArray_ContiguousFromObject($input,PyArray_INT, 0, 0)) == NULL) return NULL;
  npy_intp * dim = PyArray_DIMS(arr);
  $1->x = dim[0];
  $1->y = dim[1];
  $1->z = dim[2];
  $1->data = (int *)arr->data;
  $input = (PyObject *)arr;
}

%typemap(out) sp_c3matrix * {
  if($1->z == 1){
    /* Swap the order of the dimensions so we can plot things easily in imshow.
       This is bound to create confusion at some point in the future. */
    npy_intp dims[2] = {$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(2, dims, PyArray_CFLOAT, $1->data);
  }else{
    npy_intp dims[3] = {$1->z,$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(3, dims, PyArray_CFLOAT, $1->data);
  }
}

%typemap(out) sp_3matrix * {
  if($1->z == 1){
    /* Swap the order of the dimensions so we can plot things easily in imshow.
       This is bound to create confusion at some point in the future. */
    npy_intp dims[2] = {$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(2, dims, PyArray_FLOAT, $1->data);
  }else{
    npy_intp dims[3] = {$1->z,$1->y,$1->x};
    $result = PyArray_SimpleNewFromData(3, dims, PyArray_FLOAT, $1->data);
  }
}

%typemap(out) float[3] {
  /*$result = PyTuple_New(3);
  PyTuple_SetItem($result, 0, PyFloat_FromDouble($1[0]));
  PyTuple_SetItem($result, 1, PyFloat_FromDouble($1[1]));
  PyTuple_SetItem($result, 2, PyFloat_FromDouble($1[2]));*/
  npy_intp dims[1] = {3};
  $result = PyArray_SimpleNewFromData(1, dims, PyArray_FLOAT, $1);
}

   // I need to comment this out without knowing why this has been added originally...
   /*
%typemap(in) int * {
  PyArrayObject *arr;
  if(PyArray_Check($input) == 0 || PyArray_ISINTEGER($input) == 0 
     || PyArray_TYPE($input) != NPY_INT){
    PyErr_SetString( PyExc_TypeError, "int * argument not an int32 numpy array" );
    return NULL;
  } 
  arr = (PyArrayObject *)($input);
  $1 = (int *)arr->data;
}
   */

// ----------------------------------------------------------------------------------------------------------------------------------------------------------



// For garbage collection, not quite sure if that already works
/*
%typemap(newfree) Image * {
   sp_image_free($1);
}
%newobject make_image;
%newobject make_cimage;
*/

/*
%typemap(out) Diffraction_Pattern_Py * %{
npy_intp dims_pattern[1] = {$1->HKL_list_size};
npy_intp dims_vec[1] = {3};
npy_intp dims_hkl[2] = {$1->HKL_list_size,3};
npy_intp dims_rot[2] = {0,0}; 
if (!$1->rot) { 
$result = Py_BuildValue("{sOsOsOsOsO}",
			  "ints", PyArray_SimpleNewFromData(1, dims_pattern, PyArray_FLOAT, $1->ints),
			  "F", PyArray_SimpleNewFromData(1, dims_pattern, PyArray_CFLOAT, $1->F),
			  "max_s", PyArray_SimpleNewFromData(1, dims_vec, PyArray_FLOAT, $1->max_s),
			  "d_s", PyArray_SimpleNewFromData(1, dims_vec, PyArray_FLOAT, $1->d_s),
			  "HKL_list", PyArray_SimpleNewFromData(1, dims_hkl, PyArray_FLOAT, $1->HKL_list));
} else {
dims_rot[0] = $1->rot->cols;
dims_rot[1] = $1->rot->rows;
$result = Py_BuildValue("{sOsOsOsOsOsO}",
			  "ints", PyArray_SimpleNewFromData(1, dims_pattern, PyArray_FLOAT, $1->ints),
			  "F", PyArray_SimpleNewFromData(1, dims_pattern, PyArray_CFLOAT, $1->F),
			  "max_s", PyArray_SimpleNewFromData(1, dims_vec, PyArray_FLOAT, $1->max_s),
			  "d_s", PyArray_SimpleNewFromData(1, dims_vec, PyArray_FLOAT, $1->d_s),
			  "HKL_list", PyArray_SimpleNewFromData(1, dims_hkl, PyArray_FLOAT, $1->HKL_list),
			  "rot", PyArray_SimpleNewFromData(2, dims_rot, PyArray_FLOAT, $1->rot->data));
}
%}
*/

%include "spimage.h"
%include "spimage/image.h"
%include "spimage/image_util.h"
%include "../include/spsim.h"
%include "../include/config.h"
%include "../include/diffraction.h"
%include "../include/molecule.h"

