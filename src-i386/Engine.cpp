#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Counts(NumericVector tax, NumericVector bin)
{
	int nOccPoints = bin.size();
	int nBins = max(bin)+1;
	int nTax = max(tax)+1;
	
	LogicalMatrix TotalMatrix (nTax, nBins);
	
	//fill the matrix
	for(int i=0;i<nOccPoints;i++){
		TotalMatrix(tax(i),bin(i)) = 1;
	}
	
	// start the counting
	
	NumericVector t1(nBins);
	NumericVector t2u(nBins);
	NumericVector t2d(nBins);
	NumericVector t3(nBins);
	NumericVector tP(nBins);
	
	NumericVector tGFu(nBins);
	NumericVector tGFd(nBins);
	
	NumericVector singleton(nBins);
	NumericVector tThrough(nBins);
	NumericVector tExtNoSing(nBins);
	NumericVector tOriNoSing(nBins);
	NumericVector divSIB(nBins);
	
	// go through
	int FAD;
	int LAD;
	
	// through all taxa
	for(int j=0;j<nTax;j++){
		FAD = 0;
		LAD = 0;
		
		// through all time slices
		for(int i=0; i<nBins;i++){
			
			// the FAD
			if((FAD==0) & (TotalMatrix(j,i)==1)){
				FAD = i;
				
			}
			
			// the LAD
			if((TotalMatrix(j,i))==1){
				LAD = i;
				divSIB(i)++;
				
			}
			
			// Alroy counts
			// only count these if it can be meaningful
			if((i>0) & (i<nBins)){
				// t1 taxa
				if((TotalMatrix(j,i-1)==0) & (TotalMatrix(j,i+1)==0) & (TotalMatrix(j,i)==1)){
					t1(i)++;

				}
				
				//t3
				if((TotalMatrix(j,i-1)==1) & (TotalMatrix(j,i+1)==1) & (TotalMatrix(j,i)==1)){
					t3(i)++;
				}
				
				//tP
				if((TotalMatrix(j,i-1)==1) & (TotalMatrix(j,i+1)==1) & (TotalMatrix(j,i)==0)){
					tP(i)++;
				}
				
				
				
			}
			//t2d
			if(i>0){
				
				if((TotalMatrix(j,i-1)==1) & (TotalMatrix(j,i)==1)){
					t2d(i)++;
				}
			}
			
			//t2u
			if(i<(nBins-1)){
				
				if((TotalMatrix(j,i+1)==1) & (TotalMatrix(j,i)==1)){
					t2u(i)++;
				}
			}
			
			//tGFu
			if((i>0)& (i<(nBins-2))){	
				if((TotalMatrix(j,i+2)==1) & (TotalMatrix(j,i+1)==0)& (TotalMatrix(j,i-1)==1)){
					tGFu(i)++;
				}
			
			}
			//tGFd
			if((i>1)& (i<(nBins-1))){	
				if((TotalMatrix(j,i-2)==1) & (TotalMatrix(j,i-1)==0)& (TotalMatrix(j,i+1)==1)){
					tGFd(i)++;
				}
			
			}
			
		}
		
		// now that we know the FAD and the LAD
		// singleton
		if(FAD == LAD){
			singleton(FAD)++;
			
		}else{
			tOriNoSing(FAD)++;
			tExtNoSing(LAD)++;
	
			for(int k=(FAD+1);k<LAD;k++){
				tThrough(k)++;
			}
		
			
		}
		
		
		
	}
	
	NumericMatrix endMatrix(nBins,12);
	endMatrix(_,0) = t1;
	endMatrix(_,1) = t2d;
	endMatrix(_,2) = t2u;
	endMatrix(_,3) = t3;
	endMatrix(_,4) = tP;
	endMatrix(_,5) = tGFd;
	endMatrix(_,6) = tGFu;
	endMatrix(_,7) = singleton;
	endMatrix(_,8) = tOriNoSing;
	endMatrix(_,9) = tExtNoSing;
	endMatrix(_,10) = tThrough;
	endMatrix(_,11) = divSIB;

	return endMatrix;
	
}

