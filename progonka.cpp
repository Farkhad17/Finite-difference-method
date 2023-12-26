#include <iostream>
#include <vector>
#include <math.h>

std::vector <double> progonka_1(int n, std::vector <double> u, double h);
std::vector <double> progonka_2(int n, std::vector <double> u, double h, double alpha);
std::vector <std::vector <double>> inv_matr(int n, std::vector <std::vector<double>> matr);

double f(double x){
	return sin(x);
}

int main(){
	int n = 11, i;
	double alpha = 1, default_value = 0, h = 1.0/(n - 1);	
	std::vector <double> u(n, default_value);
	std::vector <double> u_n(n, default_value);
	std::vector <double> u_nk(n, default_value);
	
	u_n = progonka_1(n, u, h);
	
	for(i = 0; i < n; i++)std::cout << "u[" << i << "] = " << u_n[i] << std::endl; 
	std::cout << std::endl; 
	
	u_nk = progonka_2(n, u_n, h, alpha);
	
	for(i = 0; i < n; i++)std::cout << "u[" << i << "] = " << u_nk[i] << std::endl;
	std::cout << std::endl; 
	
	
	return 0;
}

std::vector <double> progonka_1(int n, std::vector <double> u, double h){
	int i;
	u[0] = 0;
	u[1] = -f(0)*pow(h, 2) + pow(u[0], 3)*pow(h, 2) + 2*u[0];
	for(i = 0; i < n - 1; i++){
		u[i+1] = -f(h*i)*pow(h, 2) + pow(u[i], 3)*pow(h, 2) + 2*u[i] - u[i-1];
	}
	
	return u;
}

std::vector <double> progonka_2(int n, std::vector <double> u, double h, double alpha){
	int i, j, k1 = 0;
	double default_value = 0, sum = 0, eps = 1e-6, sum1;	
	std::vector <std::vector<double>> A(n, std::vector<double>(n, default_value));	
	std::vector <std::vector<double>> inv_A(n, std::vector<double>(n, default_value));	
	std::vector <double> u_n(n, default_value);
	std::vector <double> u_nk(n, default_value);
	
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(i == j)A[i][j] = 2/pow(h, 2);
			else if((i == j + 1) || (i == j - 1))A[i][j] = -1/pow(h, 2);
			else A[i][j] = 0;
		}
	}
	
	for(j = 0; j < n; j++){
		A[0][j] = 0;
		A[n-1][j] = 0;
	}
	A[0][0] = A[n-1][n-1] = 1;
		
	inv_A = inv_matr(n, A);
	
	/*
	for(int k = 0; k < n; k++){
		for(int l = 0; l < n; l++)std::cout << inv_A[k][l] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/
	
	for(i = 0; i < n; i++){
		u_n[i] = f(i*h) - pow((i*h), 3);
	}	
	
	/*
	for(i = 0; i < n; i++){
		std::cout << u_n[i] << " ";
	}
	std::cout << std::endl;
	*/
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			u_nk[i] += inv_A[i][j]*u_n[j];
		}	
	}
	
	/*
	for(i = 0; i < n; i++)std::cout << u_nk[i] << " ";
	std::cout << std::endl;
	*/
	
	while(k1 < 100){
		for(i = 0; i < n; i++){
			u_n[i] = f(i*h) - pow(u_nk[i], 3);
		}	
	
		for(i = 0; i < n; i++){
			u_nk[i] = 0;
			for(j = 0; j < n; j++){
				u_nk[i] += inv_A[i][j]*u_n[j];
			}
		}
		
		for(i = 0; i < n - 1; i++){
			sum += h*pow((u_nk[i+1] - u_nk[i]),2);
			sum1 = sum;
		}
		//std::cout << "sum = " << sum << std::endl;
		if(sum < eps)break;
		else {
			sum = 0;
			for(i = 0; i < n; i++){
				for(j = 0; j < n; j++){
					u[i] += inv_A[i][j]*u_nk[j];
				}
			}
		}
		k1++;
	}
	k1 = 0; 
	
	std::cout << "sum = " << sum1 << std::endl << std::endl;
	
	return u_nk;
}

std::vector <std::vector <double>> inv_matr(int n, std::vector <std::vector<double>> matr){
	
	int i, j, m, max_i;
	double max_a, b, del,default_value;
	default_value = 0.0;
	std::vector <std::vector<double>> ematr(n, std::vector<double>(n, default_value));
	
	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(i == j)ematr[i][j] = 1;
			else ematr[i][j] = 0;
		}	
	}
	

	for(m = 0; m < n - 1; m++){		
		max_a = matr[m][m];
		
		for(i = m + 1; i < n; i++){
			if(fabs(max_a) < fabs(matr[i][m])){max_i = i; max_a = matr[i][m];}
		}
		
		if(max_i != m){
			for(j = 0; j < n; j++){
				b = matr[m][j]; matr[m][j] = matr[max_i][j]; matr[max_i][j] = b;
				b = ematr[m][j]; ematr[m][j] = ematr[max_i][j]; ematr[max_i][j] = b;
			}
		}
		
		for(i = m + 1; i < n; i++){
			del = matr[i][m]/matr[m][m];
			for(j = 0; j < n; j++){
				matr[i][j] -= del*matr[m][j];
				ematr[i][j] -= del*ematr[m][j];
			}
		}
	}
	
	for(i = 0; i < n; i++){
		del = matr[i][i];
		for(j = 0; j < n; j++){
			matr[i][j] /= del;
			ematr[i][j] /= del;
		}
	}
		
	for(m = n - 1; m > 0; m--){
		for(i = 0; i < m; i++){
			del = matr[i][m]/matr[m][m];
			for(j = 0; j < n; j++){
				matr[i][j] -= del*matr[m][j];
				ematr[i][j] -= del*ematr[m][j];
			}
		}
	}
	
	
	for(i = 0; i < n; i++){
		del = matr[i][i];
		for(j = 0; j < n; j++){
			matr[i][j] /= del;
			ematr[i][j] /= del;
		}
	}
	
	/*
	for(int k = 0; k < n; k++){
		for(int l = 0; l < n; l++)std::cout << matr[k][l] << " ";
		std::cout << "\t|";
		for(int l = 0; l < n; l++)std::cout << ematr[k][l] << " ";	
		std::cout << std::endl;
	}
	std::cout << std::endl;
	*/
	
	
	return ematr;
}


