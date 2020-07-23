package com.dzd.lambda;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.util.Log;

import com.dzd.lambda.lambda.Lambda;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Arrays;

import Jama.*;

public class MainActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

//        double[][] Qahat_double = {{7.32,2.11,3.63},{2.11,7.91,6.42},{3.63,6.42,9.02}};
//        double[][] Qahat_double = {{7,2,3,3},{2,5,6,1},{3,6,9,5},{3,1,5,10}};
//        double[][] Qahat_double = {{8,2.4,3.2,3.5},{2.4,6.7,6.6,0.89},{3.2,6.6,12,5.2},{3.5,0.89,5.2,10}};
//        double[][] Qahat_double = {{1,2,3},{4,5,6},{7,8,9}};
        double[][] Qahat_double = {
                {2.35783652853986,-1.75616901096215,1.60912498519793,-1.28505369062676,2.09152141568610},
                {-1.75616901094439,3.41319590839312,-1.59072493092294,1.87266817417825,-2.30379315998065},
                {1.60912498519570,-1.59072493093218,1.17156603850701,-1.04626215578946,1.56267941773988},
                {-1.28505369061680,1.87266817418159,-1.04626215578464,1.11315269912044,-1.50251045379102},
                {2.09152141567291,-2.30379315999071,1.56267941773415,-1.50251045379411,2.21967191211453}
        };
        Matrix Qahat = new Matrix(Qahat_double);

        double[][] ahat_double = {
                {1.00001},
                {2.00001},
                {3.00001},
                {4.00001},
                {5.0001}
        };
//        double[][] ahat_double = {{1.1},{2.2},{3.3},{4.4}};
//        double[][] ahat_double = {{1.1},{2.2},{3.3}};
        Matrix ahat = new Matrix(ahat_double);

        StringWriter stringWriter;
        PrintWriter printWriter;

//        Ldldecom ldldecom = new Ldldecom(Qahat);
//        Decorrel decorrel = new Decorrel(Qahat,ahat);
//        stringWriter = new StringWriter();
//        printWriter = new PrintWriter(stringWriter);
//        if(ldldecom.getL() != null && ldldecom.getD() != null) {
//            ldldecom.getL().print(printWriter, 3, 4);
//            ldldecom.getD().print(printWriter, 3, 4);
//            Log.d("Mainactivity", "L and D Matrix" + "\n" + stringWriter.toString());
//        }else{
//            Log.d("Mainactivity", "null Matrix");
//        }
//
//        stringWriter = new StringWriter();
//        printWriter = new PrintWriter(stringWriter);
//        if(decorrel.getQzhat() != null){
//            decorrel.getQzhat().print(printWriter,3,4);
//            decorrel.getZ().print(printWriter,3,4);
//            decorrel.getL().print(printWriter,3,4);
//            decorrel.getD().print(printWriter,3,4);
//            decorrel.getzhat().print(printWriter,3,4);
//            decorrel.getiZt().print(printWriter,3,4);
//            Log.d("Mainactivity", "decorrel" + "\n" + stringWriter.toString());
//        }else{
//            Log.d("Mainactivity", "null Matrix");
//        }


        Lambda lambda = new Lambda(MainActivity.this,ahat,Qahat,Lambda.ILS_ISEARCH);
        stringWriter = new StringWriter();
        printWriter = new PrintWriter(stringWriter);
        lambda.getafixed().print(printWriter, 3, 4);
        lambda.getQzhat().print(printWriter,3,4);
        lambda.getZ().print(printWriter,3,4);
        Log.d("MainActivity","Lambda Method Matrix\n" + stringWriter.toString());
        Log.d("MainActivity","Lambda Method sqnorm\n" + Arrays.toString(lambda.getSqnorm()));
        Log.d("MainActivity","Lambda Method other\n" + lambda.getPs() + " " + lambda.getNfixed() + " " + lambda.getMu());


    }

}
