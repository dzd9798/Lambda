package com.dzd.lambda;

import androidx.appcompat.app.AppCompatActivity;

import android.os.Bundle;
import android.util.Log;

import com.dzd.lambda.lambda.Ldldecom;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;

import Jama.*;

public class MainActivity extends AppCompatActivity {

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

//        double[][] Qahat_double = {{7,2,3},{2,5,6},{3,6,9}};
        double[][] Qahat_double = {{7,2,3,3},{2,5,6,1},{3,6,9,5},{3,1,5,10}};
//        double[][] Qahat_double = {{1,2,3},{4,5,6},{7,8,9}};
        Matrix Qahat = new Matrix(Qahat_double);

        Ldldecom ldldecom = new Ldldecom(Qahat);
        StringWriter stringWriter = new StringWriter();
        PrintWriter printWriter = new PrintWriter(stringWriter);
        if(ldldecom.getL() != null && ldldecom.getD() != null) {
            ldldecom.getL().print(printWriter, 3, 4);
            ldldecom.getD().print(printWriter, 3, 4);
            Log.d("Mainactivity", "L and D Matrix" + "\n" + stringWriter.toString());
        }else{
            Log.d("Mainactivity", "null Matrix");
        }
    }
}
