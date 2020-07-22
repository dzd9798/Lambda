package com.dzd.lambda.lambda;

import Jama.Matrix;

public class Bootstrap {
    private Matrix ahat = null;
    private Matrix L = null;
    private Matrix afixed = null;

    public Bootstrap(Matrix ahat,Matrix L){
        this.ahat = ahat.copy();
        this.L = L.copy();

        int n = ahat.getRowDimension();
        afixed = new Matrix(n,1);
        double[] afcond = new double[n];
        afcond[n-1] = ahat.get(n-1,0);
        afixed.set(n-1,0,Math.round(afcond[n-1]));

        Matrix S = new Matrix(1,n);

        for(int i = n - 1;i >= 1;i--){
            S.setMatrix(new int[]{0},0,i-1,
                    S.getMatrix(new int[]{0},0,i-1).plus(L.getMatrix(new int[]{i},0,i-1).times(afixed.get(i,0) - afcond[i])));

            afcond[i-1] = ahat.get(i-1,0) + S.get(0,i-1);

            afixed.set(i-1,0,Math.round(afcond[i-1]));
        }
    }

    public Matrix getafixed(){
        if(afixed != null){
            return afixed.copy();
        }else{
            return null;
        }
    }
}
