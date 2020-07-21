package com.dzd.lambda;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.util.Log;

import com.dzd.lambda.lambda.Decorrel;
import com.dzd.lambda.lambda.Lambda;
import com.dzd.lambda.lambda.Ldldecom;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Arrays;

import Jama.*;
import org.apache.commons.math3.distribution.NormalDistribution;

public class MainActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        double[][] Qahat_double = {{7.32,2.11,3.63},{2.11,7.91,6.42},{3.63,6.42,9.02}};
//        double[][] Qahat_double = {{7,2,3,3},{2,5,6,1},{3,6,9,5},{3,1,5,10}};
//        double[][] Qahat_double = {{8,2.4,3.2,3.5},{2.4,6.7,6.6,0.89},{3.2,6.6,12,5.2},{3.5,0.89,5.2,10}};
//        double[][] Qahat_double = {{1,2,3},{4,5,6},{7,8,9}};
        Matrix Qahat = new Matrix(Qahat_double);

//        double[][] ahat_double = {{1.1},{2.2},{3.3},{4.4}};
        double[][] ahat_double = {{1.1},{2.2},{3.3}};
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


        Lambda lambda = new Lambda(ahat,Qahat,Lambda.ILS_SSEARCH);
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
