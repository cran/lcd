#include <queue>
#include <vector>
using namespace std;


int
is_linked_c(double *amat, int p, int i, int j)
{
    queue<int> q;
    vector<int> un;
    int k;
    for(int l = 0; l < p; l++)
	un.push_back(1);
    q.push(i);
    while(!q.empty()){
	k = q.front();
	q.pop();
	if(k == j){
	    return 1;
	} else {
	    for(int l = 0; l < p; l++)
		if(amat[k*p + l] == 1.0 && amat[k + p*l] == 1.0 && un[l] == 1){
		    q.push(l);
		    un[l] = 0;
		}
	}
    }
    return 0;
}


extern "C" {

#include <R.h>
#include <Rinternals.h>    

SEXP
is_linked(SEXP amat, SEXP iR, SEXP jR)
{
    int p = INTEGER(getAttrib(amat, R_DimSymbol))[0];
    int i = INTEGER(iR)[0] - 1;
    int j = INTEGER(jR)[0] - 1;
    int res;
    SEXP result;
    PROTECT_INDEX pi;
    PROTECT_WITH_INDEX(amat, &pi);
    REPROTECT(amat = coerceVector(amat, REALSXP), pi);
    res = is_linked_c(REAL(amat), p, i, j);
    PROTECT(result = allocVector(INTSXP,1));
    INTEGER(result)[0] = res;
    UNPROTECT(2);
    return result;
}




}// end of extern "C"



    
	    
    
    
