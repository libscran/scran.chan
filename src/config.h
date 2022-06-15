#ifndef CONFIG_H
#define CONFIG_H

#include "Rcpp.h"

// Must be included before any knncolle Annoy include,
// so as to get rid of stderr usage in R libraries.
#define __ERROR_PRINTER_OVERRIDE__  REprintf 

#endif
