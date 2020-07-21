package com.dzd.lambda.lambda;
import Jama.Matrix;

//Integer ambiguity vector search by employing the search-and-shrink technique.
public class Ssearch {
    private Matrix ahat = null;
    private Matrix L = null;
    private Matrix D = null;
    private int ncands = 0;
    private Matrix afixed = null;
    private double[] sqnorm = null;
    public Ssearch(Matrix ahat,Matrix L,Matrix D,int ncands){
        if(ahat.getColumnDimension() != 1 || ahat.getRowDimension() != D.getRowDimension())
            return;
        this.ahat = ahat.copy();
        this.L = L.copy();
        this.D = D.copy();
        this.ncands = ncands;

        ssearch();
    }

    public Matrix getafixed(){
        if(afixed != null){
            return afixed.copy();
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

    private void ssearch(){
        int n = ahat.getRowDimension();
        afixed = new Matrix(n,ncands);
        sqnorm = new double[ncands];

        double Chi2 = Double.MAX_VALUE;
        double[] dist = new double[n];
        boolean endsearch = false;
        int count = 0;

        double[] acond = new double[n];
        acond[n-1] = ahat.get(n-1,0);
        double[] zcond = new double[n];
        zcond[n-1] = Math.round(acond[n-1]);
        double left = acond[n-1] - zcond[n-1];
        double[] step = new double[n];
        step[n-1] = Math.signum(left);
        if(step[n-1] == 0.0) step[n-1] = 1.0;

        int imax = ncands;

        Matrix S = new Matrix(n,n);
        int k = n;

        while (!endsearch){
            double newdist = dist[k-1] + left*left/D.get(k-1,k-1);
            if(newdist < Chi2){
                if (k != 1){
                    k = k - 1;
                    dist[k-1] = newdist;
                    S.setMatrix(new int[]{k-1},0,k-1,
                            S.getMatrix(new int[]{k},0,k-1).plus(L.getMatrix(new int[]{k},0,k-1).times(zcond[k]-acond[k])));
                    acond[k-1] = ahat.get(k-1,0) + S.get(k-1,k-1);
                    zcond[k-1] = Math.round(acond[k-1]);
                    left = acond[k-1] - zcond[k-1];
                    step[k-1] = Math.signum(left);
                    if(step[k-1] == 0.0) step[k-1] = 1.0;
                }else{
                    if(count < ncands - 1){
                        count = count + 1;
                        afixed.setMatrix(0,n-1,new int[]{count-1},
                                new Matrix(zcond,n));
                        sqnorm[count - 1] = newdist;
                    }else{
                        afixed.setMatrix(0,n-1,new int[]{imax-1},
                                new Matrix(zcond,n));
                        sqnorm[imax-1] = newdist;
                        Chi2 = sqnorm[0];
                        imax = 1;
                        for(int i = 1;i < sqnorm.length;i++){
                            if(sqnorm[i] > Chi2){
                                Chi2 = sqnorm[i];
                                imax = i + 1;
                            }
                        }
                    }
                    zcond[0] = zcond[0] + step[0];
                    left = acond[0] - zcond[0];
                    step[0] = -step[0] - Math.signum(step[0]);
                }
            }else{
                if(k == n){
                    endsearch = true;
                }else{
                    k += 1;
                    zcond[k-1] = zcond[k-1] + step[k-1];
                    left = acond[k-1] - zcond[k-1];
                    step[k-1] = -step[k-1] - Math.signum(step[k-1]);
                }
            }
        }
        int[] order = arraySort(sqnorm);
        afixed.setMatrix(0,n-1,0,ncands-1,
                afixed.getMatrix(0,n-1,order));
    }

    private int[] arraySort(double[] arr) {
        double temp;
        int index;
        int k=arr.length;
        int[]Index= new int[k];
        for(int i=0;i<k;i++) {
            Index[i]=i;
        }

        for(int i=0;i<arr.length;i++) {
            for(int j=0;j<arr.length-i-1;j++) {
                if(arr[j]>arr[j+1]) {
                    temp = arr[j];
                    arr[j] = arr[j+1];
                    arr[j+1] = temp;

                    index=Index[j];
                    Index[j] = Index[j+1];
                    Index[j+1] = index;
                }
            }
        }
        return Index;
    }

}
