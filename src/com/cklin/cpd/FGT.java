package com.cklin.cpd;

public class FGT {
	//output set
	//xc	The K center points of the training set (d x K)
	//A_k	Polynomial coefficient (pd x K), where pd = nchoosek(p + d - 1 , d) = prod(1:p + d - 1)/(prod(1:p - 1)*prod(1:d))

	public double[][] out_xc, out_Ak;

	public FGT(double[][] x, double[] w, double sigma, double e, int K, int p){
		int d = x.length;
		int Nx = x[0].length;
		K = Math.min(Nx, K);
		if(K > Nx){
			throw new RuntimeException("K must be <= Nx");
		}
		// K = (int)sqrt(Nx); (default)
		int pd = nchoosek(p + d - 1 , d);
		double[][] xc = new double[d][K];
		double[][] A_k = new double[pd][K];
		double[] C_k = new double[pd];
		double[] dist_C = new double[Nx];
		double[] dx = new double[d];
		double[] prods = new double[pd];
		int[] indxc = new int[K];
		int[] indx = new int[Nx];
		//int[] xhead = new int[K];
		int[] xboxsz = new int[K];
		int[] heads = new int[d+1];
		int[] cinds = new int[pd];

		Kcenter(x , d , Nx , K ,
				xc , indxc , indx , xboxsz ,
				dist_C);

		Compute_C_k(d , p ,
				C_k ,
				heads , cinds);

		Compute_A_k(x , w , xc , C_k , sigma , d , Nx , p , K , pd ,
				A_k ,
				indx  , dx , prods , heads );

		out_xc = xc;
		out_Ak = A_k;
	}

	private void Kcenter(double[][] x, int d, int Nx, int K,
			double[][] xc, int[] indxc, int[] indx, int[] xboxsz,
			double[] dist_C) {
		double temp ;
		int i, j, ibase, ind, nd;
		int indxc_ind = 0;
		int x_ind , x_j;
		/* randomly pick one node as the first center. */
		/*	srand( (unsigned)time( NULL ) ); */
		/*	ind      = rand() % Nx; */

		ind = 1;
		indxc[indxc_ind++] = ind;
		x_j = 0;
		x_ind = ind*d;

		for (j = 0 ; j < Nx ; x_j += d , j++){
			dist_C[j] = (j==ind) ? 0.0 : ddist(x, x_j , x_ind , d);
			indx[j]   = 0;
		}

		for(i = 1 ; i < K ; i++){
			ind      = idmax(dist_C , Nx);
			indxc[indxc_ind++] = ind;
			x_j      = 0;
			x_ind    = ind*d;

			for (j = 0 ; j < Nx ; x_j += d, j++){
				temp = (j==ind) ? 0.0 : ddist(x, x_j , x_ind , d);
				if (temp < dist_C[j]){
					dist_C[j] = temp;
					indx[j]   = i;
				}
			}
		}

		for (i = 0 ; i < K ; i++){
			xboxsz[i] = 0;
		}

		for (i = 0; i < d*K; i++){
			setVal(xc, i, 0.0);
		}

		for (i = 0 , nd = 0 ; i < Nx ; i++ , nd += d){
			xboxsz[indx[i]]++;
			ibase = indx[i]*d;
			for (j = 0 ; j < d; j++){
				setVal(xc, j + ibase,
						getVal(xc, j + ibase) + getVal(x, j + nd));
				//xc[j + ibase] += x[j + nd];
			}
		}

		for (i = 0 , ibase = 0 ; i < K ; i++ , ibase += d){
			temp = 1.0/xboxsz[i];
			for (j = 0; j < d; j++){
				//xc[j + ibase] *= temp;
				setVal(xc, j + ibase,
						getVal(xc, j + ibase) * temp);

			}
		}
	}

	private double ddist(double[][] x, int x_j, int x_ind, int d) {
		int i;
		double t, s = 0.0;
		for (i = 0 ; i < d ; i++){
			//t  = (x[i] - y[i]);
			t = getVal(x, x_j + i) - getVal(x, x_ind + i);
			s += (t * t);
		}
		return s;
	}

	private int idmax(double[] x, int N) {
		int i , k = 0;
		double t = -1.0;
		for (i = 0 ; i < N ; i++ ){
			if( t < x[i] ){
				t = x[i];
				k = i;
			}
		}
		return k;
	}

	private void Compute_C_k(int d, double p,
			double[] C_k, int[] heads, int[] cinds) {
		int i, j, k, t, tail, head;
		for(i=0; i<d; i++)
			heads[i] = 0;
		heads[d] = Integer.MAX_VALUE;
		cinds[0] = 0;
		C_k[0] = 1.0;

		for (k=1 , t=1, tail = 1 ; k < p ; k++ , tail=t){
			for (i = 0; i < d; i++){
				head     = heads[i];
				heads[i] = t;
				for ( j = head ; j < tail ; j++ , t++){
					cinds[t] = (j < heads[i+1]) ? cinds[j] + 1 : 1;
					C_k[t]   = 2.0 * C_k[j];
					C_k[t]  /= (double) cinds[t];
				}
			}
		}
	}

	private void Compute_A_k(double[][] x, double[] w, double[][] xc, double[] C_k,
			double sigma, int d, int Nx, double p, int K, int pd, double[][] A_k,
			int[] indx, double[] dx, double[] prods, int[] heads) {
		int n , i , k , t , tail , j , head , ind;
		int nbase , ix2c , ix2cbase;
		double sum , ctesigma = 1.0/(sigma) , temp , temp1;

		for (i = 0; i < pd*K; i++){
			setVal(A_k, i, 0.0);
			//A_k[i] = 0.0;
		}

		for (n = 0 ; n < Nx ; n++){
			nbase    = n*d;
			ix2c     = indx[n];
			ix2cbase = ix2c*d;
			ind      = ix2c*pd;
			temp     = w[n];
			sum      = 0.0;

			for (i = 0 ; i < d ; i++){
				dx[i]    = (getVal(x, i + nbase) - getVal(xc, i + ix2cbase))*ctesigma;
				sum     += dx[i] * dx[i];
				heads[i] = 0;
			}

			prods[0] = Math.exp(-sum);


			for (k = 1 , t = 1 , tail = 1 ; k < p ; k++ , tail = t){
				for (i = 0 ; i < d; i++){
					head     = heads[i];
					heads[i] = t;
					temp1    = dx[i];
					for ( j = head; j < tail ; j++, t++){
						prods[t] = temp1 * prods[j];
					}
				}
			}

			for (i = 0 ; i < pd ; i++){
				//A_k[i + ind] += temp*prods[i];
				setVal(A_k, i+ind,
						getVal(A_k, i+ind) + temp*prods[i]);
			}
		}

		for (k = 0 ; k < K ; k++){
			ind  = k*pd;
			for (i = 0 ; i < pd ; i++){
				//A_k[i + ind] *= C_k[i];
				setVal(A_k, i+ind,
						getVal(A_k, i+ind) * C_k[i]);
			}
		}

	}

	private static int nchoosek(int n , int k) {
		int i , n_k = n - k , nchsk = 1;
		if (k < n_k){
			k   = n_k;
			n_k = n - k;
		}
		for ( i = 1 ; i <= n_k ; i++){
			nchsk *= (++k);
			nchsk /= i;
		}
		return nchsk;
	}

	// y:	Test point (d x Ny)
	public double[] predict(double[][] y, double sigma, double e){
		return predict(y, this.out_xc, this.out_Ak, sigma, e);
	}
	public double[] predict(
			double[][] y, double[][] xc, double[][] A_k,
			double sigma, double e){
		int Ny = y[0].length;
		int d = y.length;
		int K = xc[0].length;
		if(xc.length != d)
			throw new IllegalArgumentException(
					"xc must be (d x K)");
		int pd = A_k.length;
		if(A_k[0].length != K)
			throw new IllegalArgumentException(
					"A_k must be (pd , K) where pd = nchoosek(p + d - 1 , d)");
		//default value
		//double sigma = 1.0 , e = 10.0;
		return predict(y, xc, A_k, Ny, sigma, K, e, d, pd);
	}

	private double[] predict(
			double[][] y, double[][] xc,
			double[][] A_k, int Ny, double sigma,
			int K, double e, int d, int pd){

		int p, i, j, m, k, t, tail, kn, head;
		int mbase , xbase, ind;
		double sum2 , ctesigma = 1.0/(sigma) , temp , temp1;

		p = invnchoosek(d, pd);

		double[] v = new double[Ny];
		double[] dy = new double[d];
		double[] prods = new double[pd];
		int[] heads = new int[d + 1];
		for (m=0 ; m < Ny ; m++){
			temp    = 0.0;
			mbase   = m*d;
			for (kn = 0 ; kn < K ; kn++){
				xbase = kn*d;
				ind   = kn*pd;
				sum2  = 0.0;
				for (i = 0 ; i < d ; i++){
					//dy[i]    = (y[i + mbase] - xc[i + xbase])*ctesigma;
					dy[i] = (getVal(y, i + mbase) - getVal(xc, i + xbase))*ctesigma;
					//dy[i] = (y[i][m] - xc[i][kn])*ctesigma;
					sum2    += dy[i] * dy[i];
					heads[i] = 0;
				}

				if (sum2 > e) continue; /* skip to next kn */

				prods[0] = Math.exp(-sum2);

				for (k=1, t=1, tail=1 ; k < p ; k++ , tail=t){
					for (i = 0 ; i < d; i++){
						head     = heads[i];
						heads[i] = t;
						temp1    = dy[i];
						for (j = head ; j < tail ; j++ , t++){
							prods[t] = temp1 * prods[j];
						}
					}
				}

				for (i = 0 ; i < pd ; i++){
					temp += getVal(A_k, i + ind)*prods[i];
					//temp += A_k[i][kn]*prods[i];
				}
			}
			v[m] = temp;
		}
		return v;
	}
	private static int invnchoosek(int d, int cnk) {
		int i, cted=1, ctep, cte, p;
		for(i = 2 ; i <= d ; i++){
			cted *=i;
		}
		cte  = cnk*cted;
		p    = 2;
		ctep = p;
		for (i = p + 1 ; i < p + d ; i++){
			ctep *=i ;
		}
		while(ctep != cte){
			ctep = ((p+d)*ctep)/p;
			p++;
		}
		return p;
	}


	private static double getVal(double[][] A, int index){
		int n = A.length;
		int row = index / n;
		int col = index % n;
		return A[col][row];
	}

	private static void setVal(double[][] A, int index, double val){
		int n = A.length;
		int row = index / n;
		int col = index % n;
		A[col][row] = val;
	}
	/*
	private static double getVal2(double[][] A, int index){
		int n = A[0].length;
		int row = index / n;
		int col = index % n;
		return A[row][col];
	}
	private static void setVal2(double[][] A, int index, double val){
		int n = A[0].length;
		int row = index / n;
		int col = index % n;
		A[row][col] = val;
	}
	*/
}
