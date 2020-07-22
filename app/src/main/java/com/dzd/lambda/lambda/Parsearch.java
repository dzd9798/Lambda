package com.dzd.lambda.lambda;

import org.apache.commons.math3.distribution.NormalDistribution;

import Jama.Matrix;

public class Parsearch {
    private Matrix zhat = null;
    private Matrix Qzhat = null;
    private Matrix Z = null;
    private Matrix L = null;
    private Matrix D = null;
    private double P0;
    private int ncands;

    private Matrix zpar = null;
    private double[] sqnorm = null;
    private Matrix Qzpar = null;
    private Matrix Zpar = null;
    private double Ps = Double.NaN;
    private int nfixed = 0;
    private Matrix zfixed = null;

    public Parsearch(Matrix zhat,Matrix Qzhat,Matrix Z,Matrix L,Matrix D,double P0,int ncands){
        this.zhat = zhat.copy();
        this.Qzhat = Qzhat.copy();
        this.Z = Z.copy();
        this.L = L.copy();
        this.D = D.copy();
        this.P0 = P0;
        this.ncands = ncands;

        parsearch();
    }

    public Matrix getzpar(){
        if(zpar != null){
            return zpar.copy();
        }else{
            return null;
        }
    }

    public double[] getSqnorm(){
        if(sqnorm != null){
            return sqnorm.clone();
        }else{
            return null;
        }
    }

    public Matrix getQzpar(){
        if(Qzpar != null){
            return Qzpar.copy();
        }else{
            return null;
        }
    }

    public Matrix getZpar(){
        if(Zpar != null){
            return Zpar.copy();
        }else{
            return null;
        }
    }

    public double getPs(){
        return Ps;
    }

    public int getNfixed(){
        return nfixed;
    }

    public Matrix getzfixed(){
        if(zfixed != null){
            return zfixed.copy();
        }else{
            return null;
        }
    }

    private void parsearch(){
        int n = Qzhat.getRowDimension();
        Ps = 1;
        NormalDistribution normalDistribution = new NormalDistribution();
        for (int i = 0;i < D.getColumnDimension();i++){
            double cdf = normalDistribution.probability(-Double.MAX_VALUE, 0.5/Math.sqrt(D.get(i,i)));
            Ps *= (2 * cdf - 1);
        }
        int k = 1;
        while(Ps < P0 && k < n){
            k += 1;
            Ps = 1;
            for (int i = k - 1;i < D.getColumnDimension();i++){
                double cdf = normalDistribution.probability(-Double.MAX_VALUE, 0.5/Math.sqrt(D.get(i,i)));
                Ps *= (2 * cdf - 1);
            }
        }

        if(Ps > P0){
            Ssearch ssearch = new Ssearch(zhat.getMatrix(k-1,n-1,new int[]{0}),
                    L.getMatrix(k-1,n-1,k-1,n-1),D.getMatrix(k-1,n-1,k-1,n-1),ncands);
            zpar = ssearch.getafixed();
            sqnorm = ssearch.getSqnorm();

            Qzpar = Qzhat.getMatrix(k-1,n-1,k-1,n-1).copy();
            Zpar = Z.getMatrix(0,n-1,k-1,n-1);

            Matrix QP = Qzhat.getMatrix(0,k-2,k-1,n-1).times(Qzhat.getMatrix(k-1,n-1,k-1,n-1).inverse());
            if (k == 1){
                zfixed = zpar.copy();
            }else{
                zfixed = new Matrix(n,ncands);
                for(int i = 1;i <= ncands;i++){
                    zfixed.setMatrix(0,k-2,new int[]{i-1},
                            zhat.getMatrix(0,k-2,new int[]{0}).minus
                                    (QP.times(zhat.getMatrix(k-1,n-1,new int[]{0}).minus(zpar.getMatrix(0,n-k,new int[]{i-1}))))
                    );
                }
                zfixed.setMatrix(k-1,n-1,0,ncands-1, zpar);
            }
            nfixed = n - k + 1;
        }else{
            zpar = null;
            Qzpar = null;
            Zpar = null;
            sqnorm = null;
            Ps = Double.NaN;
            zfixed = zhat;
            nfixed = 0;
        }
    }

}
