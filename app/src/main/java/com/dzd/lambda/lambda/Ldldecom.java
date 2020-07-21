package com.dzd.lambda.lambda;
import Jama.Matrix;

//LtDL decomposition
public class Ldldecom {
	private Matrix matrix = null;
	private Matrix L = null;
	private Matrix D = null;

	public Ldldecom(Matrix matrix) {
//		if (matrix == null) return;
		this.matrix = matrix.copy();
		LDLT();
		if (!isPositiveDefinite()){
			this.L = null;
			this.D = null;
		}
	}

	public Matrix getL() {
		if(L != null) {
			return L.copy();
		}else{
			return null;
		}
	}

	public Matrix getD() {
		if(D != null) {
			return D.copy();
		}else{
			return null;
		}
	}

	private boolean isPositiveDefinite(){
		if (D == null) return false;
		for(int i = 0;i < D.getColumnDimension();i++){
			if(D.get(i,i) < 1E-10 || Double.isNaN(D.get(i,i))) return false;
		}
		return true;
	}

	private void LDLT() {
		int m = matrix.getRowDimension();
		int n = matrix.getColumnDimension();
		if (m != n) {return;}
		D = new Matrix(n,n);
		L = new Matrix(n,n);
		Matrix Qahat = matrix.copy();
		for(int i = n-1;i >= 0;i--) {
			D.set(i, i, Qahat.get(i, i));
			L.setMatrix(new int[]{i},0,i,
				Qahat.getMatrix(new int[]{i},0,i).times(1.0/Math.sqrt(Qahat.get(i,i))));
			for(int j = 0;j <= i - 1;++j){
				Qahat.setMatrix(new int[]{j},0,j,
					Qahat.getMatrix(new int[]{j},0,j).minus(L.getMatrix(new int[]{i},0,j).times(L.get(i,j))));
			}
			L.setMatrix(new int[]{i},0,i,
				L.getMatrix(new int[]{i},0,i).times(1/L.get(i,i)));
		}
	}

}
