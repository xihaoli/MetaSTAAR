// [[Rcpp::depends(RcppArmadillo)]]

#include <STAAR.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace STAAR;
using namespace Rcpp;

// [[Rcpp::export]]
List Indiv_Score_Test_meta(arma::vec U, arma::vec V)
{
	int i;

	// number of markers
	int p = U.n_elem;

	arma::vec pvalue;
	pvalue.zeros(p);

	arma::vec Uscore_se;
	Uscore_se.zeros(p);

	double test_stat = 0;

	for(i = 0; i < p; i++)
	{
		Uscore_se(i) = sqrt(V(i));

		if (V(i) == 0)
		{
			pvalue(i) = 1;
		}
		else
		{
			test_stat = pow(U(i),2)/V(i);
			pvalue(i) = R::pchisq(test_stat,1,false,false);
		}

	}

	return List::create(Named("Uscore") = U, Named("Uscore_se") = Uscore_se, Named("pvalue") = pvalue);
}

