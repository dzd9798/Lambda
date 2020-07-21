package com.dzd.lambda.lambda;

// thanks to reference: https://github.com/errcw/gaussian/blob/master/lib/gaussian.js
// but finally we use common-math API : http://commons.apache.org/proper/commons-math/userguide/distribution.html
public class MyNormalDistribution {

    static double erfc(double x) {
        double z = Math.abs(x);
        double t = 1 / (1 + z / 2);
        double r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 +
                   t * (0.37409196 + t * (0.09678418 + t * (-0.18628806 +
                   t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 +
                   t * (-0.82215223 + t * 0.17087277)))))))));

        return x >= 0 ? r : 2 - r;
    }

    // Probability density function
    static double pdf(double x,double mean,double standardDeviation) {
        double m = standardDeviation * Math.sqrt(2 * Math.PI);
        double e = Math.exp(-Math.pow(x - mean, 2) / (2 * standardDeviation * standardDeviation));
        return e / m;
    }

    // Cumulative density function
    static double cdf(double x,double mean,double standardDeviation) {
        return 0.5 * erfc(-(x - mean) / (standardDeviation * Math.sqrt(2)));
    }

}
