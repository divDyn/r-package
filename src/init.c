#include <R.h>                                                         
#include <Rinternals.h>                                                
#include <stdlib.h> // for NULL                                        
#include <R_ext/Rdynload.h>                                            
                                                                       
/* FIXME:                                                              
   Check these declarations against the C/Fortran source code.         
*/                                                                     
                                                                       
/* .Call calls */                                                      
extern SEXP _divDyn_Counts(SEXP, SEXP);                                
extern SEXP _divDyn_CRbinwise(SEXP, SEXP,SEXP);
extern SEXP _divDyn_seqduplicated(SEXP); 
extern SEXP _divDyn_fillLogical(SEXP, SEXP);
extern SEXP _divDyn_fillNumeric(SEXP, SEXP, SEXP);     
extern SEXP _divDyn_fillCharacter(SEXP, SEXP);                      
                                                                       
static const R_CallMethodDef CallEntries[] = {                         
    {"_divDyn_Counts",    (DL_FUNC) &_divDyn_Counts,    2},            
    {"_divDyn_CRbinwise", (DL_FUNC) &_divDyn_CRbinwise, 3},
    {"_divDyn_seqduplicated", (DL_FUNC) &_divDyn_seqduplicated, 1}, 
    {"_divDyn_fillLogical", (DL_FUNC) &_divDyn_fillLogical, 2},
    {"_divDyn_fillNumeric", (DL_FUNC) &_divDyn_fillNumeric, 3},          
    {"_divDyn_fillCharacter", (DL_FUNC) &_divDyn_fillCharacter, 2},          
             
    {NULL, NULL, 0}                                                    
};                                                                     
                                                                       
void R_init_divDyn(DllInfo *dll)                                       
{                                                                      
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);            
    R_useDynamicSymbols(dll, FALSE);                                   
}                                                                      

