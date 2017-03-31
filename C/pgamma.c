double pchisq(double x, double df, int lower_tail, int log_p) {
    return pgamma(x, df/2., 2., lower_tail, log_p);
}

double pgamma(double x, double alph, double scale, int lower_tail, int log_p) {
    x /= scale;
    return pgamma_raw(x, alph, lower_tail, log_p);
}

double pgamma_raw(double x, double alph, int lower_tail, int log_p) {
    double res;

    if (x < 1) {
	    res = pgamma_smallx (x, alph, lower_tail, log_p);
    } else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
    	double sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
	    double d = dpois_wrap (alph, x, log_p);
	    if (!lower_tail)
	        res = log_p ? R_Log1_Exp (d + sum) : 1 - d * sum;
	    else
	        res = log_p ? sum + d : sum * d;
    } else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
    	double sum;
	    double d = dpois_wrap (alph, x, log_p);
    	if (alph < 1) {
    	    if (x * DBL_EPSILON > 1 - alph)
    		    sum = R_D__1;
	        else {
		        double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
        		sum = log_p ? log (f) : f;
	        }
	    } else {
	        sum = pd_lower_series (x, alph - 1);
	        sum = log_p ? log1p (sum) : 1 + sum;
	    }
   	    if (!lower_tail)
	        res = log_p ? sum + d : sum * d;
	    else
	        res = log_p ? R_Log1_Exp (d + sum) : 1 - d * sum;
    } else {
	    res = ppois_asymp (alph - 1, x, !lower_tail, log_p);
    }

    if (!log_p && res < DBL_MIN / DBL_EPSILON) {
	    return exp (pgamma_raw (x, alph, lower_tail, 1));
    } else
	    return res;
}

static double dpois_wrap (double x_plus_1, double lambda, int give_log) {
    if (!R_FINITE(lambda))
	    return R_D__0;
    if (x_plus_1 > 1)
	    return dpois_raw (x_plus_1 - 1, lambda, give_log);
    if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
	    return R_D_exp(-lambda - lgammafn(x_plus_1));
    else {
	    double d = dpois_raw (x_plus_1, lambda, give_log);
	    return give_log ? d + log (x_plus_1 / lambda) : d * (x_plus_1 / lambda);
    }
}

static double pd_lower_series (double lambda, double y) {
    double term = 1, sum = 0;
    while (y >= 1 && term > sum * DBL_EPSILON) {
    	term *= y / lambda;
    	sum += term;
    	y--;
    }

    if (y != floor (y)) {
    	double f;
    	f = pd_lower_cf (y, lambda + 1 - y);
    	sum += term * f;
    }
    return sum;
}

static double pd_lower_cf (double y, double d) {
    double f= 0.0 /* -Wall */, of, f0;
    double i, c2, c3, c4,  a1, b1,  a2, b2;

#define	NEEDED_SCALE				\
	  (b2 > scalefactor) {			\
	    a1 /= scalefactor;			\
	    b1 /= scalefactor;			\
	    a2 /= scalefactor;			\
	    b2 /= scalefactor;			\
	}

#define max_it 200000

#ifdef DEBUG_p
    REprintf("pd_lower_cf(y=%.14g, d=%.14g)", y, d);
#endif
    if (y == 0) return 0;

    f0 = y/d;
    /* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */
    if(fabs(y - 1) < fabs(d) * DBL_EPSILON) { /* includes y < d = Inf */
#ifdef DEBUG_p
	REprintf(" very small 'y' -> returning (y/d)\n");
#endif
	return (f0);
    }

    if(f0 > 1.) f0 = 1.;
    c2 = y;
    c4 = d; /* original (y,d), *not* potentially scaled ones!*/

    a1 = 0; b1 = 1;
    a2 = y; b2 = d;

    while NEEDED_SCALE

    i = 0; of = -1.; /* far away */
    while (i < max_it) {

	i++;	c2--;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
	a1 = c4 * a2 + c3 * a1;
	b1 = c4 * b2 + c3 * b1;

	i++;	c2--;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
	a2 = c4 * a1 + c3 * a2;
	b2 = c4 * b1 + c3 * b2;

	if NEEDED_SCALE

	if (b2 != 0) {
	    f = a2 / b2;
 	    /* convergence check: relative; "absolute" for very small f : */
	    if (fabs (f - of) <= DBL_EPSILON * fmax2(f0, fabs(f))) {
#ifdef DEBUG_p
		REprintf(" %g iter.\n", i);
#endif
		return f;
	    }
	    of = f;
	}
    }

    MATHLIB_WARNING(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n",
		    f);
    return f;/* should not happen ... */
} /* pd_lower_cf() */

