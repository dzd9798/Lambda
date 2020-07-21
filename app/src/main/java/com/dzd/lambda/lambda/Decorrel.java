package com.dzd.lambda.lambda;

import Jama.Matrix;
import Jama.EigenvalueDecomposition;

//Decorrelate a covariance matrix of ambiguities
public class Decorrel {
    private Matrix Qahat = null;
    private Matrix ahat = null;
    private Matrix Qzhat = null;
    private Matrix Z = null;
    private Matrix L = null; //LtDL decomposition of Qzhat
    private Matrix D = null; //LtDL decomposition of Qzhat
    private Matrix zhat = null;
    private Matrix iZt = null;
    private Ldldecom ldldecom = null;

    public Decorrel(Matrix Qahat,Matrix ahat){
//        if(Qahat == null || ahat == null) return;
        this.Qahat = Qahat.copy();
        this.ahat = ahat.copy();

        if(!isColRowDimensionEqual() || !isSymmetric() || !isPositiveDefinite())
            return;

        this.ldldecom = new Ldldecom(Qahat);
        this.L = ldldecom.getL();
        this.D = ldldecom.getD();

        if(this.L == null || this.D == null) return;

        decorrel();
    }

    public Matrix getQzhat(){
        if(Qzhat != null){
            return Qzhat.copy();
        }else{
            return null;
        }
    }

    public Matrix getZ(){
        if(Z != null){
            return Z.copy();
        }else{
            return null;
        }
    }

    public Matrix getL(){
        if(L != null){
            return L.copy();
        }else{
            return null;
        }
    }

    public Matrix getD(){
        if(D != null){
            return D.copy();
        }else{
            return null;
        }
    }

    public Matrix getzhat(){
        if(zhat != null){
            return zhat.copy();
        }else{
            return null;
        }
    }

    public Matrix getiZt(){
        if(iZt != null){
            return iZt.copy();
        }else{
            return null;
        }
    }

    private boolean isColRowDimensionEqual(){
        return Qahat.getColumnDimension() == Qahat.getRowDimension();
    }

    private boolean isSymmetric(){
        Matrix Qahat_tran = Qahat.transpose();
        Matrix diffMatrix = Qahat_tran.minus(Qahat);
        double[][] diffMatrixValue = diffMatrix.getArray();
        for(double[] row:diffMatrixValue){
            for(double ele:row){
                if(Math.abs(ele) > 1E-6) return false;
            }
        }
        return true;
    }

    private boolean isPositiveDefinite(){
        EigenvalueDecomposition eigObj = Qahat.eig();
        for(double eigValue:eigObj.getRealEigenvalues()){
            if(eigValue < 0) return false;
        }
        return true;
    }

    private void decorrel(){
        int n = Qahat.getColumnDimension();
        iZt = Matrix.identity(n,n);
        int i1 = n - 1;
        boolean sw = true;
        while(sw){
            int i = n;
            sw = false;
            while(!sw && i > 1){
                i = i - 1;
                if(i <= i1){
                    for (int j = i + 1;j <= n;j++){
                        int mu;
//                        mu =  (int) Math.round(L.get(j-1,i-1));  // -1.5 ---> -1  0.5 ---> 1
                        if(L.get(j-1,i-1) >= 0){  // -1.5 --->-2  0.5 --->1  more results will be the same as Lambda tools on matlab but sometimes it may different(when matlab 0.499999999--->0)
                            mu = (int) Math.floor(Math.round(L.get(j-1,i-1)*1E9d)/1E9d + 0.5); // prevent sometimes 0.499999 ---> 0
                        }else{
                            mu = (int) Math.ceil(Math.round(L.get(j-1,i-1)*1E9d)/1E9d - 0.5);
                        }

                        if(mu != 0){
                            L.setMatrix(j-1,n-1,new int[]{i-1},
                                    L.getMatrix(j-1,n-1,new int[]{i-1}).minus(L.getMatrix(j-1,n-1,new int[]{j-1}).times((double)mu)));
                            iZt.setMatrix(0,n-1,new int[]{j-1},
                                    iZt.getMatrix(0,n-1,new int[]{j-1}).plus(iZt.getMatrix(0,n-1,new int[]{i-1}).times((double)mu)));
                        }
                    }
                }

                double delta = D.get(i-1,i-1) + L.get(i,i-1)*L.get(i,i-1) * D.get(i,i);
                if(delta < D.get(i,i)) {
                    double lambda = D.get(i, i) * L.get(i, i - 1) / delta;
                    double eta = D.get(i - 1, i - 1) / delta;
                    D.set(i - 1, i - 1, eta * D.get(i, i));
                    D.set(i, i, delta);

                    L.setMatrix(i - 1, i, 0, i - 2,
                            new Matrix(new double[][]{{-L.get(i, i - 1), 1}, {eta, lambda}}).times(L.getMatrix(i - 1, i, 0, i - 2)));
                    L.set(i, i - 1, lambda);

                    Matrix tmpMatrix;
                    if(n - i - 2 >= 0) {
                        tmpMatrix = new Matrix(n - i - 1, 2);
                        tmpMatrix.setMatrix(0, n - i - 2, new int[]{0},
                                L.getMatrix(i + 1, n - 1, new int[]{i}));
                        tmpMatrix.setMatrix(0, n - i - 2, new int[]{1},
                                L.getMatrix(i + 1, n - 1, new int[]{i - 1}));
                        L.setMatrix(i + 1, n - 1, i - 1, i,
                                tmpMatrix);
                    }
                    tmpMatrix = new Matrix(n,2);
                    tmpMatrix.setMatrix(0, n - 1, new int[]{0},
                            iZt.getMatrix(0, n - 1, new int[]{i}).copy());
                    tmpMatrix.setMatrix(0, n - 1, new int[]{1},
                            iZt.getMatrix(0, n - 1, new int[]{i - 1}).copy());
                    iZt.setMatrix(0, n - 1, i-1 , i,
                            tmpMatrix);

                    i1 = i;
                    sw = true;
                }
            }
        }

        Z = iZt.transpose().inverse();
        for (int i = 0;i < Z.getRowDimension();i++){
            for (int j = 0;j < Z.getColumnDimension();j++){
                Z.set(i,j,Math.round(Z.get(i,j)));
            }
        }
        Qzhat = Z.transpose().times(Qahat).times(Z);
        if(ahat != null){
            zhat = Z.transpose().times(ahat);
        }
    }

}
