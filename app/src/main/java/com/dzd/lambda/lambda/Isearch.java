package com.dzd.lambda.lambda;

import org.apache.commons.math3.special.Gamma;
import java.io.PrintWriter;
import java.io.StringWriter;
import Jama.Matrix;

public class Isearch {
    private Matrix ahat = null;
    private Matrix L = null;
    private Matrix D = null;
    private int ncands = 0;

    private Matrix afixed = null;
    private double[] sqnorm = null;

    public Isearch(Matrix ahat,Matrix L,Matrix D,int ncands){
        this.ahat = ahat;
        this.L = L;
        this.D = D;
        this.ncands = ncands;

        isearch();
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

    private void isearch(){
        double Chi2 = chistart(1);

        Matrix Linv = L.inverse();
        Matrix Dinv = D.inverse();

        int n = ahat.getRowDimension();

        double[] right = new double[n+1];
        right[n] = Chi2;
        double[] left = new double[n+1];
        double[] dq = new double[n];
        for(int i = 0;i < n-1;i++){
            dq[i] = Dinv.get(i+1,i+1) / Dinv.get(i,i);
        }
        dq[n-1] = 1 / Dinv.get(n-1,n-1);

        boolean cand_n = false;
        boolean c_stop = false;
        boolean endsearch = false;

        int ncan = 0;

        int i = n + 1;
        int iold = i;

        afixed = new Matrix(n,ncands);
        sqnorm = new double[ncands];

        double[] lef = new double[n];
        Matrix distl = new Matrix(n,1);
        double[] endd = new double[n];
        double reach = 0,delta = 0;

        while(!endsearch){
            i -= 1;

            if (iold <= i){
                lef[i-1] = lef[i-1] + Linv.get(i,i-1);
            }else{
                lef[i-1] = 0;
                for(int j = i+1;j <= n;j++){
                    lef[i-1] = lef[i-1] + Linv.get(j-1,i-1) * distl.get(j-1,0);
                }
            }

            iold = i;
            right[i-1] = (right[i]-left[i]) * dq[i-1];
            reach = Math.sqrt(right[i-1]);
            delta = ahat.get(i-1,0) - reach - lef[i-1];
            distl.set(i-1,0,Math.ceil(delta) - ahat.get(i-1,0));

            if(distl.get(i-1,0) > reach - lef[i-1]){
                cand_n = false;
                c_stop = false;

                while(!c_stop && i < n){
                    i += 1;
                    if(distl.get(i-1,0) < endd[i-1]){
                        distl.set(i-1,0,distl.get(i-1,0)+1);
                        left[i-1] = Math.pow(distl.get(i-1,0) + lef[i-1],2);
                        c_stop = true;
                        if(i == n) cand_n = true;
                    }
                }

                if(i == n && !cand_n) endsearch = true;
            }else{
                endd[i-1] = reach - lef[i-1] - 1;
                left[i-1] = Math.pow(distl.get(i-1,0) + lef[i-1],2);
            }

            if (i == 1){
                double t = Chi2 - (right[0] - left[0]) * Dinv.get(0,0);
                endd[0] = endd[0] + 1;

                StringWriter stringWriter;
                PrintWriter printWriter;
                stringWriter = new StringWriter();
                printWriter = new PrintWriter(stringWriter);
                distl.print(printWriter, 3, 4);

                while(distl.get(0,0) <= endd[0]){
                    if(ncan < ncands){
                        ncan += 1;
                        afixed.setMatrix(0,n-1,new int[]{ncan-1},distl.plus(ahat));
                        sqnorm[ncan-1] = t;
                    }else{
                        double[] sortnorm = sqnorm.clone();
                        int[] isort = arraySort(sortnorm);
                        double maxnorm = sortnorm[ncands-1];
                        int ipos = isort[ncands-1];
                        if (t < maxnorm){
                            afixed.setMatrix(0,n-1,new int[]{ipos},distl.plus(ahat));
                            sqnorm[ipos] = t;
                        }
                    }

                    t += (2 * (distl.get(0,0)+lef[0]) + 1) * Dinv.get(0,0);
                    distl.set(0,0,distl.get(0,0)+1);

                }

                cand_n = false;
                c_stop = false;

                while(!c_stop && i < n){
                    i += 1;
                    if (distl.get(i-1,0) < endd[i-1]){
                        distl.set(i-1,0,distl.get(i-1,0) + 1);
                        left[i-1] = Math.pow(distl.get(i-1,0)+lef[i-1],2);
                        c_stop = true;
                        if(i == n) cand_n = true;
                    }
                }

                if(i == n && !cand_n) endsearch = true;
            }

        }

        int[] order = arraySort(sqnorm);
        afixed.setMatrix(0,n-1,0,ncands-1,
                afixed.getMatrix(0,n-1,order));
        for (int x = 0;x < afixed.getRowDimension();x++){
            for(int y = 0;y < afixed.getColumnDimension();y++){
                afixed.set(x,y,Math.round(afixed.get(x,y)));
            }
        }
    }

    private double chistart(double factor){
        double Chi2 = 0;
        int n = ahat.getRowDimension();
        if (ncands <= n+1){
            double[] Chi = new double[n+1];
            Matrix iQ = L.inverse().times(D.inverse()).times(L.transpose().inverse());
            for(int k = n;k >= 0;k--){
                Matrix afloat = ahat.copy();
                Matrix afixed = ahat.copy();

                for(int i = n;i >= 1;i--){
                    double dw = 0;
                    for(int j = n;j >= i;j--){
                        dw += L.get(j-1,i-1) * (afloat.get(j-1,0) - afixed.get(j-1,0));
                    }

                    afloat.set(i-1,0,afloat.get(i-1,0)-dw);
                    if(i != k){
                        afixed.set(i-1,0,Math.round(afloat.get(i-1,0)));
                    }else{
                        double tmp = Math.round(afloat.get(i-1,0));
                        afixed.set(i-1,0,tmp + Math.signum(afloat.get(i-1,0)-tmp));
                    }
                }


                Matrix tmpMatrix = ahat.minus(afixed).transpose().times(iQ).times(ahat.minus(afixed));
                Chi[n-k] = tmpMatrix.get(0,0);
            }

            arraySort(Chi);
            Chi2 = Chi[ncands-1] + 1e-6;
        }else{
            double Vn = (2.0/n) * (Math.pow(Math.PI,n/2.0) / Gamma.gamma(n/2.0));
            double prodD = 1;
            for (int i = 0;i < D.getColumnDimension();i++){
                prodD *= D.get(i,i);
            }
            Chi2 = factor * Math.pow(((double)ncands) / Math.sqrt(prodD * Vn),2.0/n);
        }

        return Chi2;
    }

    //升序排列
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
