// -*- c++ -*-

%define utilLib_DOCSTRING
"Swig-exposed classes for eotest"
%enddef

%feature("autodoc", "1");
%module(package="utilLib", docstring=utilLib_DOCSTRING) utilLib

%include "lsst/p_lsstSwig.i"
%lsst_exceptions()

%{
#define PY_ARRAY_UNIQUE_SYMBOL LSST_UTIL_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"

#include "lsst/pex/logging.h"
#include "lsst/afw.h"
#include "lsst/afw/table.h"
#include "lsst/meas/algorithms.h"

#include "lsst/eotest/ImageTools.h"
%}

%import "lsst/meas/algorithms/algorithmsLib.i"

%include "ndarray.i"
%init %{
    import_array();
%}

%include "lsst/eotest/ImageTools.h"
