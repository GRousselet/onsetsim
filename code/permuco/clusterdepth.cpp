#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
IntegerMatrix get_clusterdepth_head(IntegerMatrix cluster, String border){

  double acc = 0;


  // initialize the result matrix
  IntegerMatrix res(cluster.nrow(),cluster.ncol());

  for(int rowi = 0; rowi < res.nrow(); rowi++){
    acc = 0;
    for(int coli = 0; coli < res.ncol(); coli++){
      if(cluster(rowi,coli)==0){
        acc=0;
      }else if(cluster(rowi,coli)>0){
        acc = acc+1;
      }
      res(rowi,coli) = acc;
    }
  }

  // handle border

  if(border == "rm"){
    for(int rowi = 0; rowi < res.nrow(); rowi++){
      int coli = 0;
      while((res(rowi,coli)>=(coli+1))&(coli<res.ncol())){
        res(rowi,coli)=0;
        coli = coli+1;
      }
    }
  }


  if(border == "reverse"){
    for(int rowi = 0; rowi < res.nrow(); rowi++){
      if(res(rowi,0)>0){

        //length of border cluster
        int lcli=0;
        while((res(rowi,lcli)>0) & (lcli<res.ncol())){
          lcli = lcli+1;
        }


        for(int coli = 0; coli < lcli; coli++){
          res(rowi,coli) = lcli-coli;
        }
      }
    }
  }




  return res;

}

//[[Rcpp::export]]
NumericMatrix depth_distribution_head(NumericMatrix distribution, IntegerMatrix head){
  int max = 0;
  for(int rowi = 0; rowi < head.nrow(); rowi++){
    for(int coli = 0; coli < head.ncol(); coli++){
      if(max < head(rowi,coli)){
        max = head(rowi,coli);
      }
    }
  }

  NumericMatrix res(head.nrow(),max);


  for(int rowi = 0; rowi < head.nrow(); rowi++){
    for(int coli = 0; coli < head.ncol(); coli++){
      if(head(rowi,coli)>0){
        if(res(rowi,head(rowi,coli)-1) < distribution(rowi,coli)){
          res(rowi,head(rowi,coli)-1) = distribution(rowi,coli);
        }
      }
    }
  }




  return res;

}








