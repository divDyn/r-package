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
	
	NumericVector s1d(nBins);
	NumericVector s2d(nBins);
	NumericVector s3d(nBins);
	
	NumericVector s1u(nBins);
	NumericVector s2u(nBins);
	NumericVector s3u(nBins);
	
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
			if((i>0) & (i<(nBins-1))){
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
			
			//s1d
			if((i>0)& (i<(nBins-2))){	
				if((TotalMatrix(j,i-1)==1) & (TotalMatrix(j,i)==1) & (TotalMatrix(j,i+1)==0) & (TotalMatrix(j,i+2)==0)){
					s1d(i)++;
				}
			
			}
			
			//s2d
			if((i>0)& (i<(nBins-2))){	
				if((TotalMatrix(j,i-1)==1) & (TotalMatrix(j,i)==0) & (TotalMatrix(j,i+1)==1) & (TotalMatrix(j,i+2)==0)){
					s2d(i)++;
				}
			
			}
			
			//s3d
			if((i>0)& (i<(nBins-2))){	
				if((TotalMatrix(j,i-1)==1) & (TotalMatrix(j,i)==0) & (TotalMatrix(j,i+1)==0) & (TotalMatrix(j,i+2)==1)){
					s3d(i)++;
				}
			
			}
			
			//s1u
			if((i>1)& (i<(nBins-1))){	
				if((TotalMatrix(j,i-2)==0) & (TotalMatrix(j,i-1)==0) & (TotalMatrix(j,i)==1) & (TotalMatrix(j,i+1)==1)){
					s1u(i)++;
				}
			
			}
			
			//s2u
			if((i>1)& (i<(nBins-1))){	
				if((TotalMatrix(j,i-2)==0) & (TotalMatrix(j,i-1)==1) & (TotalMatrix(j,i)==0) & (TotalMatrix(j,i+1)==1)){
					s2u(i)++;
				}
			
			}
			
			//s3u
			if((i>1)& (i<(nBins-1))){	
				if((TotalMatrix(j,i-2)==1) & (TotalMatrix(j,i-1)==0) & (TotalMatrix(j,i)==0) & (TotalMatrix(j,i+1)==1)){
					s3u(i)++;
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
	
	NumericMatrix endMatrix(nBins,18);
	endMatrix(_,0) = t1;
	endMatrix(_,1) = t2d;
	endMatrix(_,2) = t2u;
	endMatrix(_,3) = t3;
	endMatrix(_,4) = tP;
	endMatrix(_,5) = tGFd;
	endMatrix(_,6) = tGFu;
	endMatrix(_,7) = s1d;
	endMatrix(_,8) = s2d;
	endMatrix(_,9) = s3d;
	endMatrix(_,10) = s1u;
	endMatrix(_,11) = s2u;
	endMatrix(_,12) = s3u;
	endMatrix(_,13) = singleton;
	endMatrix(_,14) = tOriNoSing;
	endMatrix(_,15) = tExtNoSing;
	endMatrix(_,16) = tThrough;
	endMatrix(_,17) = divSIB;

	return endMatrix;
	
}


int randWrapper(const int n) {
	return floor(unif_rand()*n);
}

Rcpp::NumericVector randomShuffle(Rcpp::NumericVector a) {

    // clone a into b to leave a alone
    Rcpp::NumericVector b = Rcpp::clone(a);

    std::random_shuffle(b.begin(), b.end(), randWrapper);

    return b;
}

// [[Rcpp::export]]
NumericMatrix CRbinwise(NumericVector binVar, int quota){
	// do a single loop for finding the  different bins
	int n = binVar.size();
	NumericVector indexVector(n);
	
	// define storage matrix
	NumericMatrix storage(n, 3);
	
	// initialize
	int changeCount=0;
	storage(0,0) = 0;
	
	// initialize the first bin number
	storage(0, 2) = binVar(0);
	
	for(int i=0; i<n; i++){
		// store the values for later use
		indexVector(i) =i;
		
		// check relationship to previous
		if(i>0){
			if(binVar(i)!=binVar(i-1)){
				storage(changeCount,1) = i-1; 
				// increment
				changeCount++;
				
				// the next one starts
				storage(changeCount, 0) = i;
				
				// which bin are you talking about?
				storage(changeCount, 2) = binVar(i);
			}
		}
		
		// if you reached the last row
		if(i==(n-1)){
			storage(changeCount, 1) = i;
		
		}
	
	}


	// define two integer variables
	int cStart =0;
	int cEnd = 0;
	
	int counter=0;
	
	NumericVector tempVect;
	NumericMatrix endMatrix (n, 2);
	
	//number of entries saved
	int endValue = 0;
	int vectLength = 0;
	
	// the last row
	changeCount++;
	
	// do the random sampling and store stuff in a random variable
	// how many bins are there?
	for(int j=0; j<changeCount;j++){
		// starting and ending points
		cStart= storage(j, 0);
		cEnd = storage(j, 1);
		
		
		// what will be the length of this vector? is it larger than the quota?
		vectLength = cEnd-cStart+1; 
		
		if(vectLength>=quota){
			
			// vector of indices
			NumericVector tempVect(vectLength);
			counter=0;
			for(int i=cStart;i<(cEnd+1);i++){
				tempVect(counter)= i;
				counter++;				
				
			}
			
			// shuffle vector in some way
			tempVect=randomShuffle(tempVect);
			
			for(int i=0; i<quota; i++){
				endMatrix(endValue,0)= tempVect(i);
				endMatrix(endValue,1) = storage(j, 2);
				endValue++;			
				
			}
		}else{
			// return a placeholder value ? 
			endMatrix(endValue,0) = -9;
			endMatrix(endValue,1) = storage(j, 2);
			
			
			// for which bin
			endValue++;
			
			
		}
	}
	
	// clean up the endMAtrix
	NumericMatrix realEnd(endValue, 2);
	
	for(int i=0; i<endValue;i++){
		realEnd(i,_) = endMatrix(i,_);
		
	}

	return realEnd;
}
