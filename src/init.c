#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _ABAEnrichment_binom_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ABAEnrichment_binom_randset(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ABAEnrichment_conti_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ABAEnrichment_conti_randset(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ABAEnrichment_hyper_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ABAEnrichment_hyper_randset(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _ABAEnrichment_unlock_environment(SEXP);
extern SEXP _ABAEnrichment_wilcox_category_test(SEXP, SEXP, SEXP, SEXP);
extern SEXP _ABAEnrichment_wilcox_randset(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_ABAEnrichment_binom_category_test",  (DL_FUNC) &_ABAEnrichment_binom_category_test,  4},
    {"_ABAEnrichment_binom_randset",        (DL_FUNC) &_ABAEnrichment_binom_randset,        5},
    {"_ABAEnrichment_conti_category_test",  (DL_FUNC) &_ABAEnrichment_conti_category_test,  4},
    {"_ABAEnrichment_conti_randset",        (DL_FUNC) &_ABAEnrichment_conti_randset,        5},
    {"_ABAEnrichment_hyper_category_test",  (DL_FUNC) &_ABAEnrichment_hyper_category_test,  4},
    {"_ABAEnrichment_hyper_randset",        (DL_FUNC) &_ABAEnrichment_hyper_randset,        6},
    {"_ABAEnrichment_unlock_environment",   (DL_FUNC) &_ABAEnrichment_unlock_environment,   1},
    {"_ABAEnrichment_wilcox_category_test", (DL_FUNC) &_ABAEnrichment_wilcox_category_test, 4},
    {"_ABAEnrichment_wilcox_randset",       (DL_FUNC) &_ABAEnrichment_wilcox_randset,       5},
    {NULL, NULL, 0}
};

void R_init_ABAEnrichment(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
