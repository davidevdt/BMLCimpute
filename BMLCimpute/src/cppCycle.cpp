	/* 
	
	Package BMLCimpute
	Function multilevelLCMI 

	This RCPP function provides the Gibbs Sampler used to perform Multiple Imputation through the Bayesian Multilevel Latent Class Models. 

	In the first part, functions for random numbers and random variables generators are defined. 

	In the second part, functions used for the multiple imputation (likelihood and posterior memberships computation, allocation of units across the latent classes, etc...) are introduced.

	In the third part, the main function ('cppCycle'), which is the function called in R by 'multilevelLCMI', is executed.        
	
	*/
		
	#include <Rcpp.h>
	#include <stdio.h>
	#include <iostream>
	#include <stdlib.h>
	#include <math.h>
	#include <time.h> 
		
	using namespace Rcpp;


	#define ZEROFLOAT 1.0E-300
	#define MAXNEG -1.0E300

	/*

		RANDOM VARIABLES GENERATORS
				
	*/
	
	// Uniform 1
	double ranf(void){
		double r;
		r = ((double) (rand() + 1.0)) / (RAND_MAX + 2);
		return r;
	}
	
	// Uniform 2
	double rand_(void){
		double r = 0;
		while( r == 0 || r == 1 ){
			r = ((double) rand()) / RAND_MAX;
		}
		return r;
	}

	// Multinomial
	int multdev(int K_, NumericVector pr){
		int i = -1;
		double cp = 0, r;
		r = rand_();
		do{
			i++;
			cp += pr(i);
		} while( r > cp && i < K_ - 1 );
		return i;	
	}

	// Exponential
	float expodev(void){
		static float q[8] = {
			0.6931472, 0.9333737, 0.9888778, 0.9984959, 0.9998293, 0.9999833, 0.9999986, 1.0
		};
		static long i;
		static float sexpo, a, u, ustar, umin;
		static float *q1 = q;
		a = 0.0;
		u = ranf();
		goto S30;
		S20:
			a += *q1;
		S30:
			u += u;
		if( u <= 1.0 ) goto S20;
		u -= 1.0;
		if( u > *q1 ) goto S60;
		sexpo = a + u;
		return sexpo;
		S60:
			i = 1;
		ustar = ranf();
		umin = ustar;	
		S70:
			ustar = ranf();
		if( ustar < umin ) umin = ustar;
		i += 1;
		if( u > *(q+i-1) ) goto S70;
		sexpo = a + umin**q1;
		return sexpo;
	}


	// Gamma
	double gamdev(double ia){
		double am, e, s, v1, v2, x, y,p,b;
		double zero = 1.0E-30;

		if(ia == 0) return zero;
		else if( ia < 1){
			b = 1.0+0.3678794 * ia;
			S130:
				p = b * (ranf());
			if (p >= 1.0) goto S140;
			x = exp(log(p) / ia);
			if (expodev() < x) goto S130;
			return x;
			S140:
				x = -log((b-p) / ia);
			if (expodev() < (1.0 - ia)*log(x)) goto S130;
			return x;
		} 
		else if( ia == 1 ){
			x = -log(ranf());
		} else {
			do{
				do{
					do{
						v1 = 2.0*ranf()-1.0;
						v2 = 2.0*ranf()-1.0;	
					} while( v1*v1+v2*v2 > 1.0 );
					y= v2 / (v1 + zero);
					am = (ia - 1);
					s = sqrt(2.0*am+1.0);
					x = s*y+am;
				} while( x <= 0.0 );
				if( fabs(am*log(x / am) - s*y) < 1000)
					e = (1.0+y*y) * exp(am*log(x/am)-s*y);
				else 
					e = 0;
			} while ( ranf() > e );
		}
		return x;	
	}

	// Dirichlet
	void dirdev(int n, Rcpp::NumericVector a, Rcpp::NumericVector p){
		int i;
		double ptot = 0.0;

		for( i=0; i<n; i++ ) {
			p(i) = gamdev(a(i));
			ptot += p(i);
		}
		for( i=0; i<n; i++ )
			p(i) /= ptot;
	}
	
	/*

		Multiple Imputation Functions 
				
	*/
	
	
	// Log-likelihood and posterior memberships calculations 
	double countLoglik(int L, int K, int J, int T, List pclasslev1, List piresp, List datJ, List R, NumericMatrix piK, NumericMatrix jointlev2, NumericVector piW, NumericVector grouploglik, double loglik, IntegerVector nj, List npattern, bool doVar2, List piresp2, IntegerMatrix Y2, IntegerMatrix R2, int T2, IntegerVector npattern2, IntegerVector pat2, NumericMatrix condp1, double scale, IntegerVector npatj, List prc10_, List pres_, NumericMatrix prc10, NumericMatrix prp2, NumericMatrix pres, NumericMatrix tt_, NumericMatrix tmp, IntegerMatrix dat, IntegerMatrix Rj, double largest1, double tot1, double tot2, int i, int j, int l, int t, int k, int t2){			
		
		loglik = 0.0;

		
		
		for( j=0; j<J; j++ ){
			npatj = as<IntegerVector>(npattern[j]);
			grouploglik(j) = 0.0;
			prc10_ = pclasslev1(j);
			for( l=0; l<L; l++ ){
				if( npattern2(j) > 0 ){
					jointlev2(l,j) = log(piW(l));
				}
				prc10 = as<NumericMatrix>(prc10_(l));
				for( i=0; i<nj(j); i++){
					if( npatj(i) > 0 ){
						for( k = 0; k < K; k++ ){
							prc10(k, i) = log(piK(l, k));
						}
					}
				}
				prc10_(l) = prc10;
			}	
			pclasslev1(j) = prc10_;
		}			
		
		
		for( j=0; j<J; j++ ){
			npatj = as<IntegerVector>(npattern[j]);
			prc10_ = pclasslev1(j);
			dat = as<IntegerMatrix>(datJ(j));
			Rj = as<IntegerMatrix>(R(j));	
				
			if( doVar2 ){
				if( npattern2(j) > 0 ){
					for( t2 = 0; t2 < T2; t2++ ){				
						prp2 = as<NumericMatrix>(piresp2(t2));
						if( R2(j,t2) == 0 ){
							for( l=0; l<L; l++ ){
								jointlev2(l,j) += log(prp2(l,Y2(j,t2) - 1));
							}
						}else{
							for( l=0; l<L; l++ ){
								jointlev2(l,j) += 0.0;
							}
						}									
					}
				}else{
					for( l=0; l<L; l++){
						jointlev2(l,j) = jointlev2(l,pat2(j));
					}
				}						
			}				
				
			for( l=0; l<L; l++ ){					
				prc10 = as<NumericMatrix>(prc10_(l));	
				for( t=0; t<T; t++ ){
					pres_ = piresp(t);
					pres = as<NumericMatrix>(pres_(l));
					for( i=0; i<nj(j); i++ ){
						if( npatj(i) > 0 ){
							if( Rj(i,t) == 0 ){
								for( k=0; k<K; k++ ){
									if( pres(k,(dat(i,t)-1)) <= ZEROFLOAT ){
										prc10(k,i) += log(ZEROFLOAT);
									}else{
										prc10(k,i) += log(pres(k,(dat(i,t) - 1)));
									}									
								}
							}
						}
					}
				}
				prc10_(l) = prc10;
			}
			pclasslev1(j) = prc10_;
		}
		
		
		for( j=0; j<J; j++ ){
			npatj = as<IntegerVector>(npattern[j]);
			prc10_ = pclasslev1(j);
			tmp = NumericMatrix(L,nj(j));
			tot2 = 0.0;
			for( l=0; l<L; l++ ){				
				prc10 = as<NumericMatrix>(prc10_(l));
				for( i=0; i<nj(j); i++ ){
					if( npatj(i) > 0 ){
						tot1 = 0.0;
						largest1 = max(prc10(_,i));
						prc10(_,i) = exp(prc10(_,i) - largest1);
						tot1 = sum(prc10(_,i));
						prc10(_,i) = prc10(_,i) / tot1; 
						tmp(l,i) += npatj(i) * (log(tot1) + largest1);
					}
				} 
				prc10_(l) = prc10;
				condp1(l,j) = sum(tmp(l,_)); 
				tot2 += exp((jointlev2(l,j) + (condp1(l,j))) * scale);	
			}
			grouploglik(j) += log(tot2);
			pclasslev1(j) = prc10_;	
		}	
		loglik = sum(grouploglik);	
		return loglik;		
	}
		
		
	
	// Allocate level-2 and level-1 units to latent classes by using the posterior memberships just calculated 
	void allocation(int L, int K, int J, IntegerVector nj, int T, List pclasslev1, NumericMatrix jointlev2, NumericVector nlatw, NumericMatrix nlatpk, List datJ, List R, List nresp, IntegerVector ncat, List z, IntegerVector w, int b, double prilev2, double prilev1, double priresp, List pat, List npattern, List nlat2, double prioresp2, bool doVar2, int T2, IntegerVector ncat2, IntegerMatrix R2, IntegerMatrix Y2, IntegerVector sel, NumericMatrix condp1, bool saveAllocations, IntegerVector npatj, IntegerVector zj, IntegerVector patj, NumericVector probsb, List prc10_, List nresp_, NumericMatrix prc10, NumericMatrix nlt2, NumericMatrix nrsp, IntegerMatrix dat, IntegerMatrix Rj, int zz, int j, int l, int i, int ii, int k, int t, int s, int t2, int s2){
		
		// Initialize counts to prior hyperparameters
		
		saveAllocations = (is_true(any(sel==(b))) && b>0);		

		
		if( doVar2 ){
			for( t2=0; t2<T2; t2++ ){
				nlt2 = as<NumericMatrix>(nlat2(t2));
				for( l=0; l<L; l++ ){
					for( s2=0; s2<ncat2(t2); s2++ ){
						nlt2(l,s2) = prioresp2;
					}
				}
				nlat2(t2) = nlt2;						
			}	
		}		
		
		for( l=0; l<L; l++ ){
			nlatw(l) = prilev2;
			for( k=0; k<K; k++){
				nlatpk(l,k) = prilev1;
				for( t=0; t<T; t++ ){
					nresp_ = nresp(t);
					nrsp = as<NumericMatrix>(nresp_(l));
					for( s=0; s<ncat(t); s++ ){
						nrsp(k,s) = priresp;
					}				
					nresp_(l) = nrsp;
					nresp(t) = nresp_;
				}
			}			
		}		
		
		
		// Sampling step
		
		for( j=0; j<J; j++ ){
			dat = as<IntegerMatrix>(datJ(j));
			Rj = as<IntegerMatrix>(R(j));
			prc10_ = pclasslev1(j);
			
			probsb = jointlev2(_,j) + condp1(_,j);
			probsb = exp(probsb - max(probsb));
			probsb = probsb/sum(probsb);
			w(j) = multdev(L,probsb);
				
			nlatw(w(j)) += 1.0;	
				
			if(doVar2){
				for( t2=0; t2<T2; t2++ ){
					if( R2(j,t2)==0 ){
						nlt2 = as<NumericMatrix>(nlat2(t2));
						nlt2(w(j),(Y2(j,t2)-1)) += 1.0;
						nlat2(t2) = nlt2;
					}
				}
			}				

			npatj = as<IntegerVector>(npattern[j]);		
			prc10 = as<NumericMatrix>(prc10_(w(j)));
			
			for( i=0; i<nj(j); i++ ){
				if( npatj(i) > 0 ){
					if( !saveAllocations ){
						for( ii=0; ii<npatj(i); ii++ ){
							zz = multdev(K,prc10(_,i));
							nlatpk(w(j),zz) += 1.0;
							for( t=0; t<T; t++ ){
								if( Rj(i,t) == 0 ){
									nresp_ = nresp(t);
									nrsp = as<NumericMatrix>(nresp_(w(j)));
									nrsp(zz,(dat(i,t) - 1)) += 1.0;
									nresp_(w(j)) = nrsp;
									nresp(t) = nresp_;
								}
							}						
						}				
					}else{
						zj = as<IntegerVector>(z(j));	
						patj = as<IntegerVector>(pat[j]);
						for( ii=0; ii<nj(j); ii++ ){
							if( patj(ii) == i ){
								zj(ii) = multdev(K, prc10(_,i));
								nlatpk(w(j),zj(ii)) += 1.0;
								for( t=0; t<T; t++ ){
									if( Rj(i,t) == 0 ){
										nresp_ = nresp(t);
										nrsp = as<NumericMatrix>(nresp_(w(j)));
										nrsp(zj(ii),(dat(i,t)-1)) += 1.0;
										nresp_(w(j)) = nrsp;
										nresp(t) = nresp_;
									}
								}					
							}else{
								continue;
							}
						}
						z(j) = zj;
					}				
				}			
			}				
		}
		
	}
	
	// Draw parameter values from conditional posteriors 
	void posterior(int L, int K, NumericVector nlatw, NumericMatrix nlatpk, List nresp, NumericVector piW, NumericMatrix piK, List piresp, IntegerVector ncat, int T, bool doVar2, int T2, List nlat2, List piresp2, IntegerVector ncat2, NumericVector pix, NumericVector pp2, NumericVector pp, List pres_, List nresp_, NumericMatrix prp2, NumericMatrix pres, NumericMatrix nlt2, NumericMatrix nrsp, int l, int k, int t, int t2){
		
		
		dirdev(L, nlatw, piW);

		
		for(l=0; l<L; l++){
			pix = NumericVector(K);
			dirdev(K,nlatpk(l,_),pix);
			piK(l,_) = pix;
		}
		
		if( doVar2 ){
			for( t2=0; t2<T2; t2++ ){
				prp2 = as<NumericMatrix>(piresp2(t2));
				nlt2 = as<NumericMatrix>(nlat2(t2));
				for( l=0; l<L; l++ ){
					pp2 = NumericVector(ncat2(t2));
					dirdev(ncat2(t2),nlt2(l,_),pp2);
					prp2(l,_) = pp2;
				}
				piresp2(t2) = prp2;
			}
		}
		

		for( t=0; t<T; t++ ){
			nresp_ = nresp(t);
			pres_ = piresp(t);
			for( l=0; l<L; l++ ){
				nrsp = as<NumericMatrix>(nresp_(l));
				pres = as<NumericMatrix>(pres_(l));
				for( k=0; k<K; k++ ){
					pp = NumericVector(ncat(t));
					dirdev(ncat(t),nrsp(k,_),pp);
					pres(k,_) = pp;
				}
				pres_(l) = pres;
			}
			piresp(t) = pres_;
		}
	}
	
	
	// Perform the imputations 
	void imputation(IntegerVector w, List z, List piresp, List imp, int J, IntegerVector nj, List R, int T, int nimp, IntegerVector ncat, bool doVar2, List imp2, int T2, IntegerMatrix R2, List piresp2, IntegerVector ncat2, IntegerVector zj, List pres_, List imp_, NumericMatrix prp2, NumericMatrix pres, IntegerMatrix Rj, IntegerMatrix imp2_, IntegerMatrix impj, int ind, int j, int t, int i, int t2){		
		
		if( doVar2 ){
			for( t2=0; t2<T2; t2++ ){
				if( is_true(any(R2(_,t2) == 1)) ){
					imp2_ = as<IntegerMatrix>(imp2(t2));
					prp2 = as<NumericMatrix>(piresp2(t2));
					ind = 0;
					for( j=0; j<J; j++ ){
						if( R2(j,t2) == 1){							
							imp2_(ind,nimp) = (multdev(ncat2(t2),prp2(w(j),_))) + 1; 
							ind += 1;
						}
					}
					imp2(t2) = imp2_;
				}
			}		
		}		
		
		
		for( j=0; j<J; j++ ){
			Rj = as<IntegerMatrix>(R(j));
			zj = as<IntegerVector>(z(j));
			imp_ = imp(j);				
			for( t=0; t<T; t++ ){
				if( is_true(any(Rj(_,t) == 1)) ){
					impj = as<IntegerMatrix>(imp_(t));
					pres_ = piresp(t);
					pres = as<NumericMatrix>(pres_(w(j)));
					ind = 0;
					for( i=0; i<nj(j); i++ ){
						if( Rj(i,t) == 1 ){
							impj(ind,nimp) = (multdev(ncat(t),pres(zj(i),_))) + 1;
							ind += 1;
						}
					}
					imp_(t) = impj;
					imp(j) = imp_;	
				}				
			}				
		}		
	}	

	// Observed patterns (level-1)
	Rcpp::List countPatterns(Rcpp::List R, Rcpp::List datJ, Rcpp::List npattern, Rcpp::List pat, IntegerVector nj, int J, int T){

		for( int j=0; j<J; j++ ){
			int nrec = nj(j);
			SEXP a_ = R[j];
			IntegerMatrix Rj_(a_);
			SEXP b_ = datJ[j];
			IntegerMatrix Y(b_);
			SEXP c_ = npattern[j];
			IntegerVector npatj_(c_);
			SEXP d_ = pat[j];
			IntegerVector patj_(d_);
			
			for(int i=0;i<nrec;i++){
				int r = 0;
				int found = 0;
				while( !found && r < i ){
					found = 1;
					int k = 0;
					while( found && k < T ){
						if( Rj_(r,k) == 1 ){
							if( Rj_(i,k) == 1 ){
								k++;	
							}else{
								found = 0;
							}
						}else{
							if( Rj_(i,k) == 1 ){
								found = 0;
							}else if( Rj_(i,k) == 0 && Y(r,k) != Y(i,k) ){
								found = 0;
							}else{
								k++;
							}
						}
					}
					if( !found ){
						r++;
					}else{
						npatj_(r)++;
						patj_(i) = r;	
					}
				}
				if( !found ){
					npatj_(i)++;
					patj_(i) = i;
				}
			}
			
			npattern[j] = clone(npatj_);
			pat[j] = clone(patj_);		
		}			
		return Rcpp::List::create(npattern,pat);					
	}
	
	
	
	// Observed patterns (level-2)
	Rcpp::List countPattern2(IntegerMatrix R2, IntegerMatrix Y2, IntegerVector npattern2, IntegerVector pat2, int J, int T2){

		for( int i=0; i<J; i++ ){
			int r = 0;
			int found = 0;
			while( !found && r < i ){
				found = 1;
				int k = 0;
				while( found && k < T2 ){
					if( R2(r,k) == 1 ){
						if( R2(i,k) == 1 ){
							k++;	
						}else{
							found = 0;
						}
					}else{
						if( R2(i,k) == 1 ){
							found=0;
						}else if( R2(i,k) == 0 && Y2(r,k) != Y2(i,k) ){
							found = 0;
						}else{
							k++;
						}
					}
				}
				if( !found ){
					r++;
				}else{
					npattern2(r)++;
					pat2(i) = r;	
				}
			}
			if( !found ){
				npattern2(i)++;
				pat2(i) = i;
			}		
		}			
		return Rcpp::List::create(npattern2,pat2);					
	}
	
	// Function that stores the current parameter value in order to calculate posterior means and standard deviations of the model parameters
	double store(double loglik, double loglikm, NumericVector piW, NumericVector piWm, NumericVector piWs, NumericMatrix piK,  NumericMatrix piKm, NumericMatrix piKs, List piresp, List pirespm, List piresps, List piresp2, List piresp2m, List piresp2s, IntegerVector ncat, IntegerVector ncat2, int L, int K, int T, int T2, bool doVar2, List pres_, List pr1m__, List pr1s__, NumericMatrix prp2, NumericMatrix pres, NumericMatrix pr1m, NumericMatrix pr1s, int l, int k, int t, int c1, int t2, int c2){
		
		for(l=0;l<L;l++){
			piWm(l) += piW(l);
			piWs(l) += (piW(l) * piW(l));
			for( k=0; k<K; k++ ){
				piKm(l,k) += piK(l,k);
				piKs(l,k) += (piK(l,k) * piK(l,k));
			}			
		}

		for( t=0; t<T; t++ ){
			pres_ = piresp(t);
			pr1m__ = pirespm(t);
			pr1s__ = piresps(t);
			for( l=0; l<L; l++ ){
				pres = as<NumericMatrix>(pres_(l));
				pr1m = as<NumericMatrix>(pr1m__(l));
				pr1s = as<NumericMatrix>(pr1s__(l));
				for( k=0; k<K; k++ ){
					for( c1=0; c1<ncat(t); c1++ ){
						pr1m(k,c1) += pres(k,c1);
						pr1s(k,c1) += (pres(k,c1) * pres(k,c1));
					}
				}
				pr1m__(l) = pr1m;
				pr1s__(l) = pr1s;
			}
			pirespm(t) = pr1m__;	
			piresps(t) = pr1s__;
		}
		
		if(doVar2){			

			NumericMatrix prp2; 
			NumericMatrix pr2m;
			NumericMatrix pr2s;
			
			for( t2=0; t2<T2; t2++ ){
				prp2 = as<NumericMatrix>(piresp2(t2));
				pr2m = as<NumericMatrix>(piresp2m(t2));
				pr2s = as<NumericMatrix>(piresp2s(t2));
				for( l=0; l<L; l++ ){
					for( c2=0; c2<ncat2(t2); c2++ ){
						pr2m(l,c2) += prp2(l,c2);
						pr2s(l,c2) += (prp2(l,c2) * prp2(l,c2));
					}
				}
				piresp2m(t2) = pr2m;
				piresp2s(t2) = pr2s;
			}
		}
		
		loglikm += loglik;
		return loglikm;		
	}
	
	
	
	

	

	
	// [[Rcpp::export]]
	
	
	List cppCycle(int L, int K, int J, int it1, int it2, int it3, int itprint, double prilev2, double prilev1, double priresp, double priresp2, bool estimates, bool doVar2, IntegerVector sel, int T2, IntegerVector ncat2, IntegerMatrix Y2, IntegerMatrix R2, List nlat2, List piresp2, IntegerVector nj, int T, IntegerVector ncat, List datJ, List R, List npattern, List pat, IntegerVector w, List z, List nresp, List piresp, List pclasslev1, NumericMatrix jointlev2, NumericVector grouploglik, double loglik, NumericVector nlatw, NumericMatrix nlatpk, NumericMatrix piK, NumericVector piW, List imp, List imp2, int nimp, bool count, NumericVector piWm, NumericMatrix piKm, List pirespm, List piresp2m, NumericVector piWs, NumericMatrix piKs, List piresps, List piresp2s, double f, bool plotloglik, NumericVector pllk, NumericMatrix condp1, IntegerVector npattern2, IntegerVector pat2, double scale){
		
		
	// cppCycle : main function called by 'multilevelLCMI'; it runs the Gibbs sampler and calls the functions to store the parameter estimates (if 'estimates = TRUE'). Furthermore, log-likelihoods and distributions of the latent classes are calculated. 
	
	
		
		int b = 0;
		int l_ = 0; 
		int k = 0;
		int ha1 = 0;
		int ha2 = 0;
		int i = 0;
		int j = 0;
		int l = 0;
		int t = 0;
		int t2 = 0; 
		int ii = 0; 
		int s = 0; 
		int s2 = 0;
		int c1 = 0; 
		int c2 = 0;
		
		int zz = 0;
		int ind = 0;
		
		double maxlik = MAXNEG;
		double loglikm = 0.0;
		double DIC = 0.0;
		
		double largest1 = 0.0;
		double tot1 = 0.0;
		double tot2 = 0.0;		
		
		IntegerVector countIT(L);	
		IntegerVector countL(L);
		
		IntegerMatrix countK(L,K);
		
		bool saveAllocations = FALSE;			
		
		IntegerVector npatj;			
		IntegerVector zj;				
		IntegerVector patj; 			
		
		NumericVector probsb(L);		
		NumericVector pix;				
		NumericVector pp2;				
		NumericVector pp;					
		
		List prc10_;					
		List pres_;						
		List nresp_;					
		List imp_;						
		List pr1m__;    				
		List pr1s__;	
		
		NumericMatrix prc10;			
		NumericMatrix prp2;				
		NumericMatrix pres;				
		NumericMatrix tt_;				
		NumericMatrix tmp;				
		NumericMatrix nlt2;				
		NumericMatrix nrsp;				
		NumericMatrix pr1m;				
		NumericMatrix pr1s;				
		
		IntegerMatrix dat;				
		IntegerMatrix Rj;				
		IntegerMatrix imp2_;			
		IntegerMatrix impj; 
		
		RNGScope scope;		
		
		countPatterns(R, datJ, npattern, pat, nj, J, T);
		if( doVar2 ){
			countPattern2(R2, Y2, npattern2, pat2, J, T2);
		}
		
		// Burn-in
		for( b=1; b<=it1; b++ ){

			loglik = countLoglik(L, K, J, T, pclasslev1, piresp, datJ, R, piK, jointlev2, piW, grouploglik, loglik, nj, npattern, doVar2, piresp2, Y2, R2, T2, npattern2, pat2, condp1, scale, npatj, prc10_, pres_, prc10, prp2, pres, tt_, tmp, dat, Rj, largest1, tot1, tot2, i, j, l, t, k, t2);


			if(plotloglik){
				pllk(b-1) = loglik;
				if( loglik > maxlik ){
					maxlik = loglik;
				}
			}
			
			allocation(L, K, J, nj, T, pclasslev1, jointlev2, nlatw, nlatpk, datJ, R, nresp, ncat, z, w, -1, prilev2, prilev1, priresp, pat, npattern, nlat2, priresp2, doVar2, T2, ncat2, R2, Y2, sel, condp1, saveAllocations, npatj, zj, patj, probsb, prc10_, nresp_, prc10, nlt2, nrsp, dat, Rj, zz, j, l, i, ii, k, t, s, t2, s2);


			posterior(L, K, nlatw, nlatpk, nresp, piW, piK, piresp, ncat, T, doVar2, T2, nlat2, piresp2, ncat2, pix, pp2, pp, pres_, nresp_, prp2, pres, nlt2, nrsp, l, k, t, t2);

			if( itprint > 0 && b % (itprint) == 0 ){
				if( scale == 1.0 ){
					Rcout << "Burn-in, iteration " << b << ", loglik = " << loglik << std::endl; 
				}else{
					Rcout << "Burn-in, iteration " << b << ", loglik = " << loglik << "^(1/" << scale << ")" << std::endl;
				}
			}	
		}
		
		// Posterior computation 
		for( b=1; b<=it2; b++ ){

			loglik = countLoglik(L, K, J, T, pclasslev1, piresp, datJ, R, piK, jointlev2, piW, grouploglik, loglik, nj, npattern, doVar2, piresp2, Y2, R2, T2, npattern2, pat2, condp1, scale, npatj, prc10_, pres_, prc10, prp2, pres, tt_, tmp, dat, Rj, largest1, tot1, tot2, i, j, l, t, k, t2);

			if( plotloglik ){
				pllk((b-1) + it1) = loglik;
				if(loglik > maxlik){
					maxlik = loglik;
				}
			}

			allocation(L, K, J, nj, T, pclasslev1, jointlev2, nlatw, nlatpk, datJ, R, nresp, ncat, z, w, b, prilev2, prilev1, priresp, pat, npattern, nlat2, priresp2, doVar2, T2, ncat2, R2, Y2, sel, condp1, saveAllocations, npatj, zj, patj, probsb, prc10_, nresp_, prc10, nlt2, nrsp, dat, Rj, zz, j, l, i, ii, k, t, s, t2, s2);

			if( count ){
				ha1 = sum(nlatw > prilev2) - 1;
				countL(ha1) += 1;
				for( l_=0; l_<L; l_++ ){
					if( nlatw(l_) > prilev2 ){
						countIT(l_) += 1;
						ha2 = sum(nlatpk(l_,_) > prilev1) - 1;
						countK(l_,ha2) += 1; 
					}
				}
			}
		
			posterior(L, K, nlatw, nlatpk, nresp, piW, piK, piresp, ncat, T, doVar2, T2, nlat2, piresp2, ncat2, pix, pp2, pp, pres_, nresp_, prp2, pres, nlt2, nrsp, l, k, t, t2);

				
			if( is_true(any(sel == (b))) ){
				imputation(w, z, piresp, imp, J, nj, R, T, nimp, ncat, doVar2, imp2, T2, R2, piresp2, ncat2, zj, pres_, imp_, prp2, pres, Rj, imp2_, impj, ind, j, t, i, t2);
				nimp += 1;
			}

			if(estimates && (b%it3) == 0){
				loglikm = store(loglik, loglikm, piW, piWm, piWs, piK, piKm, piKs, piresp, pirespm, piresps, piresp2, piresp2m, piresp2s, ncat, ncat2, L, K, T, T2, doVar2, pres_, pr1m__, pr1s__, prp2, pres, pr1m, pr1s, l, k, t, c1, t2, c2);
			}

			if( itprint > 0 && b % (itprint) == 0){
				if( scale == 1.0 ){
					Rcout << "Posterior computation, iteration " << b << ", loglik = " << loglik << std::endl;
				}else{
					Rcout << "Posterior computation, iteration " << b << ", loglik = " << loglik << "^(1/" << scale << ")" <<  std::endl;
				}			
			}		
		}
		
		// Compute Posterior Means 
		if( estimates ){
			int t1;
			double mloglik;	
			
			for( l=0; l<L; l++ ){
				piWm(l) *= f;
				
				if( piWm(l) < 1.0E-300 ){
					piWm(l) = 1.0E-300;
				}
				piWs(l) *= f;
				piWs(l) -= (piWm(l) * piWm(l));
				if( piWs(l) > 0 ){
					piWs(l) = sqrt(piWs(l));
				}else{
					piWs(l) = 1.0E-300;
				}			
				for( k=0; k<K; k++ ){
					piKm(l,k) *= f;
					if( piKm(l,k) < 1.0E-300 ){
						piKm(l,k) = 1.0E-300;
					}
					piKs(l,k) *= f;
					piKs(l,k) -= (piKm(l,k) * piKm(l,k));
					if( piKs(l,k) > 0 ){
						piKs(l,k) = sqrt(piKs(l,k));
					}else{
						piKs(l,k) = 1.0E-300;
					}	
				}				
			}
			
			for( t1=0; t1<T; t1++ ){
				pr1m__ = pirespm(t1);
				pr1s__ = piresps(t1);
				for( l=0; l<L; l++ ){
					pr1m = as<NumericMatrix>(pr1m__(l));
					pr1s = as<NumericMatrix>(pr1s__(l));
					for( k=0; k<K; k++ ){
						for( c1=0; c1<ncat(t1); c1++ ){
							pr1m(k,c1) *= f;
							if( pr1m(k,c1) < 1.0E-300 ){
								pr1m(k,c1) = 1.0E-300;
							}
							pr1s(k,c1) *= f;
							pr1s(k,c1) -= (pr1m(k,c1) * pr1m(k,c1));
							if( pr1s(k,c1) > 0 ){
								pr1s(k,c1) = sqrt(pr1s(k,c1));
							}else{
								pr1s(k,c1) = 1.0E-300;
							}
						}
					}
					pr1m__(l) = pr1m;
					pr1s__(l) = pr1s;
				}
				pirespm(t1) = pr1m__;
				piresps(t1) = pr1s__;			
			}			
			
			if( doVar2 ){
				NumericMatrix pr2m;
				NumericMatrix pr2s;
				
				for( t2=0; t2<T2; t2++ ){
					pr2m = as<NumericMatrix>(piresp2m(t2));
					pr2s = as<NumericMatrix>(piresp2s(t2));
					for( l=0; l<L; l++ ){
						for( c2=0; c2<ncat2(t2); c2++ ){
							pr2m(l,c2) *= f;
							if( pr2m(l,c2) < 1.0E-300 ){
								pr2m(l,c2) = 1.0E-300;
							}
							pr2s(l,c2) *= f;
							pr2s(l,c2) -= (pr2m(l,c2) * pr2m(l,c2));
							if( pr2s(l,c2) > 0 ){
								pr2s(l,c2) = sqrt(pr2s(l,c2));
							}else{
								pr2s(l,c2) = 1.0E-300;
							}
						}						
					}
					piresp2m(t2) = pr2m;
					piresp2s(t2) = pr2s;
				}			
			}
			loglikm *= f;
			mloglik = countLoglik(L, K, J, T, pclasslev1, pirespm, datJ, R, piKm, jointlev2, piWm, grouploglik, loglik, nj, npattern,doVar2, piresp2m, Y2, R2, T2, npattern2, pat2, condp1, scale, npatj, prc10_, pres_, prc10, prp2, pres, tt_, tmp, dat, Rj, largest1, tot1, tot2, i, j, l, t, k, t2);
			
			DIC = -2*((2*loglikm) - mloglik);		
		}		

		
		return List::create(
			Rcpp::Named("piWm") = piWm,
			Rcpp::Named("piKm") = piKm,
			Rcpp::Named("pirespm") = pirespm,
			Rcpp::Named("DIC") = DIC,
			Rcpp::Named("imp") = imp,
			Rcpp::Named("piWs") = piWs,
			Rcpp::Named("piKs") = piKs,
			Rcpp::Named("piresps") = piresps,
			Rcpp::Named("imp2") = imp2,
			Rcpp::Named("piresp2m") = piresp2m,
			Rcpp::Named("piresp2s") = piresp2s,
			Rcpp::Named("pllk") = pllk,
			Rcpp::Named("countL") = countL,
			Rcpp::Named("countK") = countK, 
			Rcpp::Named("maxlik") = maxlik,
			Rcpp::Named("countIT") = countIT
		);	
	}
