package com.dzd.lambda.lambda;

import android.content.Context;
import android.util.Log;

import com.dzd.lambda.R;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.json.JSONArray;
import org.json.JSONObject;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Locale;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

public class Lambda {
    private Context context = null;
    private Matrix ahat = null;
    private Matrix Qahat = null;
    private int method;
    private double P0;
    boolean FFRT = false;
    private int ncands;

    private Matrix afixed = null;
    private double[] sqnorm = null;
    private double Ps;
    private Matrix Qzhat = null;
    private Matrix Z = null;
    private int nfixed;
    private double mu;

    public static final int ILS_SSEARCH = 1,ILS_ISEARCH = 2,IR = 3,IB = 4,PAR = 5,ILS_SSEARCH_RT = 6;

    public Lambda(Context context,Matrix ahat, Matrix Qahat,int method,Object...objects){
        this.context = context;
        this.ahat = ahat.copy();
        this.Qahat = Qahat.copy();
        this.method = method;
        //only one unique solution is available when apply IR and IB method
        if(method == IR || method == IB){
            ncands = 1;
        }else{
            ncands = 2;
        }

        if(method == PAR){
            P0 = 0.995;
        }

        if(method == ILS_SSEARCH_RT){
            P0 = 0.001;
            FFRT = true;
        }else{
            mu = 1;
        }

        int i = 0;
        if(objects != null) {
            while (i < objects.length){
                if(objects[i] instanceof String){
                    switch ( ((String)objects[i]).toUpperCase() ){
                        case "P0":
                            if( i + 1 >= objects.length || !(objects[i+1] instanceof Double)){
                                Log.w("Lambda","not valid input arg");
                                return;
                            }
                            P0 = (Double) objects[i + 1];
                            if(method == ILS_SSEARCH_RT && (P0 != 0.01 && P0 != 0.001) ){
                                Log.w("Lambda","Fixed failure rate must be either 0.01 or 0.001");
                                return;
                            }else if(method == PAR && P0 > 1.0){
                                Log.w("Lambda","User-defined success rate P0 cannot be larger than 1");
                                return;
                            }
                            i += 2;
                            break;
                        case "MU":
                            if( i + 1 >= objects.length || !(objects[i+1] instanceof Double)){
                                Log.w("Lambda","not valid input arg");
                                return;
                            }
                            mu = (Double) objects[i + 1];
                            if(method == ILS_SSEARCH_RT && mu > 1){
                                Log.w("Lambda","MU must be between 0 and 1");
                                return;
                            }
                            FFRT = false;
                            i += 2;
                            break;
                        case "NCANDS":
                            if( i + 1 >= objects.length || !(objects[i+1] instanceof Integer)){
                                Log.w("Lambda","not valid input arg");
                                return;
                            }
                            if(method == ILS_SSEARCH || method == ILS_ISEARCH || method == PAR){
                                ncands = (Integer) objects[i + 1];
                            }
                            i += 2;
                            break;
                        default:
                            Log.w("Lambda","not valid input arg");
                            return;
                    }
                }else{
                    Log.w("Lambda","not valid input arg");
                    return;
                }
            }
        }

        applyMethod();
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

    public double getPs(){
        return Ps;
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

    public int getNfixed(){
        return nfixed;
    }

    public double getMu(){
        return mu;
    }

    private void applyMethod(){
        int n = Qahat.getRowDimension();
        Matrix zfixed = new Matrix(n,ncands);
        nfixed = n;
        sqnorm = null;

        if(!isColRowDimensionEqual() || !isSymmetric() || !isPositiveDefinite()){
            Log.w("Lambda","Qahat not valid");
            return;
        }
        if(ahat.getColumnDimension() != 1){
            Log.w("Lambda","ahat not valid");
            return;
        }

        Matrix incr = new Matrix(n,1);
        for(int i = 0;i < n;i++){
            incr.set(i,0,ahat.get(i,0) - (ahat.get(i,0) % 1) );
            ahat.set(i,0,ahat.get(i,0) % 1 );
        }

        //Compute Z matrix based on the decomposition  Q=L^T*D*L
        Decorrel decorrel = new Decorrel(Qahat,ahat);
        Qzhat = decorrel.getQzhat();
        Z = decorrel.getZ();
        Matrix L = decorrel.getL();
        Matrix D = decorrel.getD();
        Matrix zhat = decorrel.getzhat();
        Matrix iZt = decorrel.getiZt();

        //Compute the bootstrapped success rate
        Ps = 1;
        NormalDistribution normalDistribution = new NormalDistribution();
        for (int i = 0;i < D.getColumnDimension();i++){
            double cdf = normalDistribution.probability(-Double.MAX_VALUE, 0.5/Math.sqrt(D.get(i,i)));
            Ps *= (2 * cdf - 1);
        }
        Ssearch ssearch;
        switch (method){
            case ILS_SSEARCH:
                ssearch = new Ssearch(zhat,L,D,ncands);
                zfixed = ssearch.getafixed();
                sqnorm = ssearch.getSqnorm();
                break;
            case ILS_ISEARCH:

                break;
            case IR:
                zfixed = new Matrix(n,1);
                for (int i = 0;i < n;i++){
                    zfixed.set(i,0,Math.round(zhat.get(i,0)));
                }
                break;
            case IB:
                Bootstrap bootstrap = new Bootstrap(zhat,L);
                zfixed = bootstrap.getafixed();
                break;
            case PAR:
                Parsearch parsearch = new Parsearch(zhat,Qzhat,Z,L,D,P0,ncands);
                Matrix zpar = parsearch.getzpar();
                sqnorm = parsearch.getSqnorm();
                Qzhat = parsearch.getQzpar(); //此处变量是否为Qzpar更好?
                Z = parsearch.getZpar();
                Ps = parsearch.getPs();
                nfixed = parsearch.getNfixed();
                zfixed = parsearch.getzfixed();
                break;
            case ILS_SSEARCH_RT:
                ssearch = new Ssearch(zhat,L,D,ncands);
                zfixed = ssearch.getafixed();
                sqnorm = ssearch.getSqnorm();
                if(FFRT){
                    if(1 - Ps > P0){
                        mu = ratioinv(this.context,P0,1-Ps,n);
                    }else{
                        mu = 1;
                    }
                }
                if(sqnorm[0]/sqnorm[1] > mu){
                    zfixed = zhat;
                    nfixed = 0;
                    ncands = 1;
                }
                break;
            default:
                Log.w("Lambda","impossible!");
                break;
        }

        afixed = iZt.times(zfixed);
        for(int i = 0;i < ncands;i++){
            afixed.setMatrix(0,n - 1,new int[]{i},
                    afixed.getMatrix(0,n - 1,new int[]{i}).plus(incr));
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

    // Determine appropriate threshold value MU for Ratio Test with Fixed
    // Failure rate Pf.
    // Use tabulated values of MU depending on the ILS failure rate and
    // the number of float ambiguities
    static public double ratioinv(Context context, double Pf_FIX, double PfILS, int n){
        int kPf = (int) Math.round(Pf_FIX*1000);
        if(kPf != 1 && kPf != 10) return Double.NaN;
        if(n < 1) return Double.NaN;
        double[] PfList = null;
        double[] muList = null;
        double[][] table_double = null;
        try {
            InputStream ratioTableInputStream = context.getResources().openRawResource(R.raw.ratiotab);
            String ratioTableStr = readTextInputStream(ratioTableInputStream);

            JSONObject jsonObject = new JSONObject(ratioTableStr);
            JSONArray jsonArray = jsonObject.getJSONArray(String.format(Locale.CHINESE,"table%d",kPf));

            if(n > jsonArray.getJSONArray(0).length() - 1){
                n = jsonArray.getJSONArray(0).length() - 1;
            }

            PfList = new double[jsonArray.length()];
            muList = new double[jsonArray.length()];
            for(int i = 0;i < jsonArray.length();i++){
                PfList[i] = jsonArray.getJSONArray(i).getDouble(0);
                muList[i] = jsonArray.getJSONArray(i).getDouble(n);
            }
        }catch (Exception e){
            e.printStackTrace();
            return Double.NaN;
        }

        UnivariateInterpolator interpolator = new LinearInterpolator();
        UnivariateFunction function = interpolator.interpolate(PfList, muList);

        return function.value(PfILS);
    }

    static private String readTextInputStream(InputStream is) throws Exception {
        InputStreamReader reader = new InputStreamReader(is);
        BufferedReader bufferedReader = new BufferedReader(reader);
        StringBuffer buffer = new StringBuffer("");
        String str;
        while ((str = bufferedReader.readLine()) != null) {
            buffer.append(str);
            buffer.append("\n");
        }
        return buffer.toString();
    }

}
