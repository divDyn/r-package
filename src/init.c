#include <R.h>                                                         
#include <Rinternals.h>                                                
#include <stdlib.h> // for NULL                                        
#include <R_ext/Rdynload.h>                                            
                                                                       
/* FIXME:                                                              
   Check these declarations against the C/Fortran source code.         
*/                                                                     
                                                                       
/* .Call calls */                                                      
extern SEXP _divDyn_Counts(SEXP, SEXP);                                
extern SEXP _divDyn_CRbinwise(SEXP, SEXP);         
extern SEXP _divDyn_seqduplicated(SEXP);                     
                                                                       
static const R_CallMethodDef CallEntries[] = {                         
    {"_divDyn_Counts",    (DL_FUNC) &_divDyn_Counts,    2},            
    {"_divDyn_CRbinwise", (DL_FUNC) &_divDyn_CRbinwise, 2},            
    {"_divDyn_seqduplicated", (DL_FUNC) &_divDyn_seqduplicated, 1},        
    {NULL, NULL, 0}                                                    
};                                                                     
                                                                       
void R_init_divDyn(DllInfo *dll)                                       
{                                                                      
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);            
    R_useDynamicSymbols(dll, FALSE);                                   
}                                                                      

