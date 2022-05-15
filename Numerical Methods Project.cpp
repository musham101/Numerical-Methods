#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

// EQUATION (Integration include)
#define f(x) x * x  - 9 * x + 1

// derivative of EQUATION (edit if you changed the EQUATION using Newton Raphson Method)
#define df(x) 2 * x - 9

// EQUATIONS for linear equation (Jacobi Method & Gauss Seidel Method)
#define f1(x,y,z)  (17-y+2*z)/20
#define f2(x,y,z)  (-18-3*x+z)/20
#define f3(x,y,z)  (25-2*x+3*y)/20

// EQUATION (dy/dx) for ODE (Eular Method, Runga Kutta Method)
#define dy(x,y)  x + 2 * y

class Weedles_Rule{
	double start, end, h;
	public:
		void setA(double a){
			this->start = a;
		}
		void setB(double b){
			this->end = b;
		}
		void setInterval(double h){
			this->h = h;
		}
		
		void calulate()
		{
			int n = ((end - start)/h) + 1;
			double sum_of_2s = 0;
			double sum_of_3s = 0;
			double remaining_sum = 0;
			double *X_point = new double[n];
			double *X_value = new double[n];
			int times = n / 6;
			int current = 0;
			int s = 0;
			
			for(int j = 0; j < n; j++){
				X_point[j] = start + (j * h);
				X_value[j] = f(X_point[j]);
			}
			
			for(int i = 1; i <= times; i++){
				current = s;
				for(s = current; s < i * 6; s++){
					if(s % 2 == 0){
						sum_of_2s = sum_of_2s + *(X_value + s);
					}
					else if(s % 3 == 0){
						sum_of_3s = sum_of_3s + *(X_value + s);
					}
					else{
						remaining_sum = remaining_sum + *(X_value + s);
					}
				}
				
			}
			for(int s = 0; s < n; s++){
				if(s == 0){
					cout << "\nX = [" << X_point[s] << ", ";
				}
				else if(s == n - 1){
					cout << X_point[s] << "]";
				}
				else{
					cout << X_point[s] << ", ";
				}
			}
			cout << endl;
			for(int s = 0; s < n; s++){
				if(s == 0){
					cout << "\nF(X) = [" << X_value[s] << ", ";
				}
				else if(s == n - 1){
					cout << X_value[s] << "]";
				}
				else{
					cout << X_value[s] << ", ";
				}
			}
			cout << endl;
			double FX = (0.3 * h) * ((sum_of_2s) + (6 * sum_of_3s) + (5 * remaining_sum));
			cout << "\nIntegration of F(X) = " << FX << endl;
		}
};
	

class Booles_Rule{
	double start, end, h;
	public:
		void setA(double a){
			this->start = a;
		}
		void setB(double b){
			this->end = b;
		}
		void setInterval(double h){
			this->h = h;
		}
		
		void calulate()
		{
			int n = ((end - start)/h) + 1;
			double sum_of_4s = 0;
			double sum_of_2s = 0;
			double remaining_sum = 0;
			double *X_point = new double[n];
			double *X_value = new double[n];
			for(int s = 0; s < n; s++){
				X_point[s] = start + (s * h);
				X_value[s] = f(X_point[s]);
				if(s > 0 and s < n - 2 and s % 4 == 0){
					sum_of_4s = sum_of_4s + *(X_value + s);
				}
				else if(s > 0 and s < n - 2 and s % 4 != 0 and s % 2 == 0){
					sum_of_2s = sum_of_2s + *(X_value + s);
				}
				else{
					remaining_sum = remaining_sum + *(X_value + s);
				}
			}
			for(int s = 0; s < n; s++){
				if(s == 0){
					cout << "\nX = [" << X_point[s] << ", ";
				}
				else if(s == n - 1){
					cout << X_point[s] << "]";
				}
				else{
					cout << X_point[s] << ", ";
				}
			}
			cout << endl;
			for(int s = 0; s < n; s++){
				if(s == 0){
					cout << "\nF(X) = [" << X_value[s] << ", ";
				}
				else if(s == n - 1){
					cout << X_value[s] << "]";
				}
				else{
					cout << X_value[s] << ", ";
				}
			}
			cout << endl;
			double FX = (0.04444 * h) * ((7 * (X_value[0] + X_value[n - 1])) + (32 * remaining_sum) + (12 * sum_of_2s) + (14 * sum_of_4s));
			cout << "\nIntegration of F(X) = " << FX << endl;
		}
};

class Runga_Kutta{
	public:
		void calculate(){
			float x0, y0, xn, h, yn, k1, k2, k3, k4, k;
			int i;
			cout<<"Enter initial value of 'X': ";
			cin>> x0;
			cout<<"Enter initial value of 'Y': ";
			cin >> y0;
			cout<<"Enter value of 'X' at which you want to calculate 'Y': ";
			cin >> xn;
			cout << "Enter the interval size 'h': ";
			cin >> h;
			int n = (xn - x0)/h;
			cout<<"\nx0\ty0\tyn\n";
			cout<<"------------------\n";
			for(i=0; i < n; i++){
				k1 = h * (dy(x0, y0));
				k2 = h * (dy((x0+h/2), (y0+k1/2)));
				k3 = h * (dy((x0+h/2), (y0+k2/2)));
				k4 = h * (dy((x0+h), (y0+k3)));
				k = (k1+2*k2+2*k3+k4)/6;
				yn = y0 + k;
				cout<< x0<<"\t"<< y0<<"\t"<< yn<< endl;
				x0 = x0+h;
				y0 = yn;
			}
			cout<<"\nY("<< xn << ") = " << yn;
	    }
};

class Lagrange{
	int n;
	
	public:
		void setN(int n){
			this->n = n;
		}
		
		void calculate(){
			float *x = new float[n], *y = new float[n], xp, yp = 0, p;
			cout << "\nEnter the values of X\n";
			for(int s = 0; s < n; s++){
				cout << "X" << s << ": ";
				cin >> x[s];
			}
			cout << "\nEnter the values of Y\n";
			for(int s = 0; s < n; s++){
				cout << "Y" << s << ": ";
				cin >> y[s];
			}
			cout << "Enter the value at with you want to find Y: ";
			cin >> xp;
			for(int i = 1; i <= n; i++){
				p=1;
				for(int j = 1; j <= n; j++){
					if(i != j){
						p = p* (xp - x[j]) / (x[i] - x[j]);
					}
				}
				yp = yp + p * y[i];
			}
			cout<< endl<<"F("<< xp<< ") = "<< yp;
	    }
		
};

class Newton_Divided{
	int n;
	
	float proterm(int i, float value, float x[]){
		float pro = 1;
		for (int j = 0; j < i; j++) {
			pro = pro * (value - x[j]);
		}
		return pro;
    }
    
    void dividedDiffTable(float x[], float **y, int n){
    	for (int i = 1; i < n; i++) {
    		for (int j = 0; j < n - i; j++) {
    			y[j][i] = (y[j][i - 1] - y[j + 1][i - 1]) / (x[j] - x[i + j]);
			}
		}
	}
	
	float applyFormula(float value, float x[], float **y, int n){
		float sum = y[0][0];
		for (int i = 1; i < n; i++) {
			sum = sum + (proterm(i, value, x) * y[0][i]);
		}
		return sum;
	}
	
	void printDiffTable(float **y,int n){
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - i; j++) {
				cout << setw(15) << y[i][j];
			}
		cout << endl;
		}
	}
	
	public:
		
		void setN(int n){
			this->n = n;
		}
	
		void calculate(){
			float value;
			float *x = new float[n];
			float **y = new float* [n];
			for(int s = 0; s < n; s++){
				y[s] = new float[n];
			}
			cout << "Enter the values of X\n";
			for(int s = 0; s < n; s++){
				cout << "X" << s << ": ";
				cin >> x[s];
			}
			cout << "\nEnter the values of Y\n";
			for(int s = 0; s < n; s++){
				cout << "Y" << s << ": ";
				cin >> y[s][0];
			}
			cout << "Enter the value at with you want to find Y: ";
			cin >> value;
			dividedDiffTable(x, y, n);
			printDiffTable(y,n);
			cout << "\nF(" << value << ") = " << applyFormula(value, x, y, n) << endl;
		}
		
};

class Simpson_3by8{
	double start, end, h;
	public:
		void setA(double a){
			this->start = a;
		}
		void setB(double b){
			this->end = b;
		}
		void setInterval(double h){
			this->h = h;
		}
		
		void calulate()
		{
			int n = ((end - start)/h) + 1;
			double sum_of_3s = 0;
			double remaining_sum = 0;
			double *X_point = new double[n];
			double *X_value = new double[n];
			for(int s = 0; s < n; s++){
				X_point[s] = start + (s * h);
				X_value[s] = f(X_point[s]);
				if(s > 0 and s < n - 2 and s % 3 == 0){
					sum_of_3s = sum_of_3s + *(X_value + s);
				}
				else if(s > 0 and s < n - 2 and s % 3 != 0){
					remaining_sum = remaining_sum + *(X_value + s);
				}
			}
			for(int s = 0; s < n; s++){
				if(s == 0){
					cout << "\nX = [" << X_point[s] << ", ";
				}
				else if(s == n - 1){
					cout << X_point[s] << "]";
				}
				else{
					cout << X_point[s] << ", ";
				}
			}
			cout << endl;
			for(int s = 0; s < n; s++){
				if(s == 0){
					cout << "\nF(X) = [" << X_value[s] << ", ";
				}
				else if(s == n - 1){
					cout << X_value[s] << "]";
				}
				else{
					cout << X_value[s] << ", ";
				}
			}
			cout << endl;
			double FX = (0.375 * h) * ((X_value[0] + X_value[n - 1]) + (3 * remaining_sum) + (2 * sum_of_3s));
			cout << "\nIntegration of F(X) = " << FX << endl;
		}
};

class Simpson_1by3{
	double start, end, h;
	public:
		void setA(double a){
			this->start = a;
		}
		void setB(double b){
			this->end = b;
		}
		void setInterval(double h){
			this->h = h;
		}
		
		void calulate()
		{
			int n = ((end - start)/h) + 1;
			double sum_of_2s = 0;
			double remaining_sum = 0;
			double *X_point = new double[n];
			double *X_value = new double[n];
			for(int s = 0; s < n; s++){
				X_point[s] = start + (s * h);
				X_value[s] = f(X_point[s]);
				if(s > 0 and s < n - 2 and s % 2 == 0){
					sum_of_2s = sum_of_2s + *(X_value + s);
				}
				else if(s > 0 and s < n - 2 and s % 2 != 0){
					remaining_sum = remaining_sum + *(X_value + s);
				}
			}
			for(int s = 0; s < n; s++){
				if(s == 0){
					cout << "\nX = [" << X_point[s] << ", ";
				}
				else if(s == n - 1){
					cout << X_point[s] << "]";
				}
				else{
					cout << X_point[s] << ", ";
				}
			}
			cout << endl;
			for(int s = 0; s < n; s++){
				if(s == 0){
					cout << "\nF(X) = [" << X_value[s] << ", ";
				}
				else if(s == n - 1){
					cout << X_value[s] << "]";
				}
				else{
					cout << X_value[s] << ", ";
				}
			}
			cout << endl;
			double FX = (0.33333 * h) * ((X_value[0] + X_value[n - 1]) + (4 * remaining_sum) + (2 * sum_of_2s));
			cout << "\nIntegration of F(X) = " << FX << endl;
		}
};

class Euler_Method{
	double x, y, h, final;
	public:
		void setX(double x){
			this->x = x;
		}
		void setY(double y){
			this->y = y;
		}
		void setInterval(double h){
			this->h = h;
		}
		void setFinal(double final){
			this->final = final;
		}
		void calculate(){
			int n = (final - x)/ h;
			int i = 0;
			while(i <= n){
				if(i > 0){
					x = x + h;
				}
				cout << i + 1 << setw(15) << x << setw(15) << y <<"\n";
				if(i == n){
					cout << "\nY(" << final << ") = " << y;
					break;
				}
				double DF = dy(x, y);                      
				y = y + (h * DF);    
				i++;
				
			}
		}
};

class Newton_Backward{
	int n;
	float u_cal(float u, int n){
		float temp = u;
		for (int i = 1; i < n; i++){
			temp = temp * (u - i);
		}
		return temp;
	}
	int fact(int n){
		int f = 1;
		for (int i = 2; i <= n; i++){
			f *= i;
		}
		return f;
	}
	public:
	void setN(int n){
		this->n = n;
	}
	
	void calculate(){
		float value;
		float *x = new float[n];
		float **y = new float* [n];
		for(int s = 0; s < n; s++){
			y[s] = new float[n];
		}
		cout << "Enter the values of X\n";
		for(int s = 0; s < n; s++){
			cout << "X" << s << ": ";
			cin >> x[s];
		}
		cout << "\nEnter the values of Y\n";
		for(int s = 0; s < n; s++){
			cout << "Y" << s << ": ";
			cin >> y[s][0];
		}
		cout << "Enter the value at with you want to find Y: ";
		cin >> value;
		
		for (int i = 1; i < n; i++) {
		    for (int j = n - 1; j >= i; j--){
		    	y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
		    }
		}
		
		cout << endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j <= i; j++){
				cout << setw(4) << y[i][j] << "\t";
		    }
        cout << endl;
        }
        
        float sum = y[n - 1][0];
		float u = (value - x[n - 1]) / (x[1] - x[0]);
		for (int i = 1; i < n; i++) {
			sum = sum + (u_cal(u, i) * y[n - 1][i]) / fact(i);
        }
        
        cout << "\n Value at " << value << " is " << sum << endl;
		
	}
};

class Newton_Forward{
	int n;
	float u_cal(float u, int n){
		float temp = u;
		for (int i = 1; i < n; i++){
			temp = temp * (u - i);
		}
		return temp;
	}
	int fact(int n){
		int f = 1;
		for (int i = 2; i <= n; i++){
			f *= i;
		}
		return f;
	}
	public:
	void setN(int n){
		this->n = n;
	}
	
	void calculate(){
		float value;
		float *x = new float[n];
		float **y = new float* [n];
		for(int s = 0; s < n; s++){
			y[s] = new float[n];
		}
		cout << "Enter the values of X\n";
		for(int s = 0; s < n; s++){
			cout << "X" << s << ": ";
			cin >> x[s];
		}
		cout << "\nEnter the values of Y\n";
		for(int s = 0; s < n; s++){
			cout << "Y" << s << ": ";
			cin >> y[s][0];
		}
		cout << "Enter the value at with you want to find Y: ";
		cin >> value;
		for (int i = 1; i < n; i++) {
			for (int j = 0; j < n - i; j++){
				y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
			}
		}
		
		cout << endl;
		for (int i = 0; i < n; i++) {
			cout << setw(4) << x[i]<< "\t";
			for (int j = 0; j < n - i; j++){
				cout << setw(4) << y[i][j] << "\t";
			} 
			cout << endl;
		}
		
		float sum = y[0][0];
		float u = (value - x[0]) / (x[1] - x[0]);
		for (int i = 1; i < n; i++) {
			sum = sum + (u_cal(u, i) * y[0][i]) / fact(i);
		}
		
		cout << "\n F(" << value << ") = " << sum << endl;
	}
};

class Gauss_Seidel{
	float x0, y0, z0, e;
    int i = 1;
    
    public:
    void setX(float x){
    	this->x0 = x;
	}
	void setY(float y){
    	this->y0 = y;
	}
	void setZ(float z){
    	this->z0 = z;
	}
	void setTolerance(float e){
		this->e = e;
	}
    
    void calculate(){
    	float x1, y1, z1, e1, e2, e3;
    	cout << "#" << setw(15) << "X" << setw(15) << "Y" << setw(15) << "Z\n";
    	do{
    		cout << i << setw(15) << x1 << setw(15) << y1 << setw(15) << z1  << endl;
    		
    		x1 = f1(x0,y0,z0);
            y1 = f2(x1,y0,z0);
            z1 = f3(x1,y1,z0);
            
            e1 = fabs(x0-x1);
			e2 = fabs(y0-y1);
			e3 = fabs(z0-z1);
			
			i++;
			
			x0 = x1;
			y0 = y1;
			z0 = z1;
		}while(e1>e and e2>e and e3>e);
		
		cout << "\nSolution: \n";
		cout << "\tX = " << x1 << endl;
		cout << "\tY = " << y1 << endl;
		cout << "\tZ = " << z1 << endl;
	}
};

class Jacobi{
	float x0, y0, z0, e;
    int i = 1;
    
    public:
    void setX(float x){
    	this->x0 = x;
	}
	void setY(float y){
    	this->y0 = y;
	}
	void setZ(float z){
    	this->z0 = z;
	}
	void setTolerance(float e){
		this->e = e;
	}
    
    void calculate(){
    	float x1, y1, z1, e1, e2, e3;
    	cout << "#" << setw(15) << "X" << setw(15) << "Y" << setw(15) << "Z\n";
    	do{
    		cout << i << setw(15) << x1 << setw(15) << y1 << setw(15) << z1  << endl;
    		
    		x1 = f1(x0,y0,z0);
            y1 = f2(x0,y0,z0);
            z1 = f3(x0,y0,z0);
            
            e1 = fabs(x0-x1);
			e2 = fabs(y0-y1);
			e3 = fabs(z0-z1);
			
			i++;
			
			x0 = x1;
			y0 = y1;
			z0 = z1;
		}while(e1>e and e2>e and e3>e);
		
		cout << "\nSolution: \n";
		cout << "\tX = " << x1 << endl;
		cout << "\tY = " << y1 << endl;
		cout << "\tZ = " << z1 << endl;
	}
};

class Trapezoidal{
	double start, end, h;
	public:
		void setA(double a){
			this->start = a;
		}
		void setB(double b){
			this->end = b;
		}
		void setInterval(double h){
			this->h = h;
		}
		
		void calulate()
		{
			int n = ((end - start)/h) + 1;
			double sum = 0;
			double *X_point = new double[n];
			double *X_value = new double[n];
			for(int s = 0; s < n; s++){
				X_point[s] = start + (s * h);
				X_value[s] = f(X_point[s]);
				if(s > 0 and s < n - 2){
					sum = sum + *(X_value + s);
				}
			}
			for(int s = 0; s < n; s++){
				if(s == 0){
					cout << "\nX = [" << X_point[s] << ", ";
				}
				else if(s == n - 1){
					cout << X_point[s] << "]";
				}
				else{
					cout << X_point[s] << ", ";
				}
			}
			cout << endl;
			for(int s = 0; s < n; s++){
				if(s == 0){
					cout << "\nF(X) = [" << X_value[s] << ", ";
				}
				else if(s == n - 1){
					cout << X_value[s] << "]";
				}
				else{
					cout << X_value[s] << ", ";
				}
			}
			cout << endl;
			double FX = (0.5 * h) * ((X_value[0] + X_value[n - 1]) + 2 * sum);
			cout << "\nIntegration of F(X) = " << FX << endl;
		}
};


class Secant{
	double a, b, tolerance;
	int max_itrations;
	public:
		void setA(double a){
			this->a = a;
		}
		void setB(double b){
			this->b = b;
		}
		void setTolerance(double tolerance){
			this->tolerance = tolerance;
		}
		void setMax_itrations(int max_itrations){
			this->max_itrations = max_itrations;
		}
		void calulate()
		{
			int i = 1;
			cout << "#" << setw(15) << "a" << setw(15) << "b" <<setw(15) << "c" << setw(15) << "f(a)" << setw(15) << "f(b)" << setw(15) << "f(c)\n";
			while(i <= max_itrations){
				double FA = f(a); double FB = f(b);
				double c = (a * FB - b * FA) / (FB - FA);
				double FC = f(c);
				cout << i << setw(15) << a << setw(15) << b << setw(15) << c << setw(15) << FA << setw(15) << FB << setw(15) << FC << endl;
				if(FC == 0 or fabs(FC) < tolerance){
					cout << "X = " << c << endl;
					return;
				}
				else{
					a = b;
					b = c;
				}
				i++;
			}
			cout << "\t-----Failed to get root for the given number of itrations\n";
		}
};

class Regula_Falsi{
	double a, b, tolerance;
	int max_itrations;
	public:
		void setA(double a){
			this->a = a;
		}
		void setB(double b){
			this->b = b;
		}
		void setTolerance(double tolerance){
			this->tolerance = tolerance;
		}
		void setMax_itrations(int max_itrations){
			this->max_itrations = max_itrations;
		}
		void calulate()
		{
			int i = 1;
			cout << "#" << setw(15) << "a" << setw(15) << "b" <<setw(15) << "c" << setw(15) << "f(a)" << setw(15) << "f(b)" << setw(15) << "f(c)\n";
			while(i <= max_itrations){
				double FA = f(a); double FB = f(b);
				double c = (a * FB - b * FA) / (FB - FA);
				double FC = f(c);
				cout << i << setw(15) << a << setw(15) << b << setw(15) << c << setw(15) << FA << setw(15) << FB << setw(15) << FC << endl;
				if(FC == 0 or fabs(FC) < tolerance){
					cout << "X = " << c << endl;
					return;
				}
				if(FA * FC < 0){
					b = c;
				}
				else{
					a = c;
				}
				i++;
			}
			cout << "\t-----Failed to get root for the given number of itrations\n";
		}
};

class Newton_Raphson{
	double a, tolerance;
	int max_itrations;
	public:
		void setA(double a){
			this->a = a;
		}
		void setTolerance(double tolerance){
			this->tolerance = tolerance;
		}
		void setMax_itrations(int max_itrations){
			this->max_itrations = max_itrations;
		}
		void calculate(){
			int i = 1;
			cout << "#" << setw(15) << "a" << setw(15) << "f(a)" << setw(15) << "f'(a)" << setw(15) << "c" << setw(15) << "f(c)\n";
			while(i <= max_itrations){
				double FA = f(a); double D_FA = df(a);
				double c = a - (FA / D_FA);
				double FC = f(c);
				cout << i << setw(15) << a << setw(15) << FA << setw(15) << D_FA << setw(15) << c << setw(15) << FC << endl;
				if(FC == 0 or fabs(FC) < tolerance){
					cout << "X = " << c << endl;
					return;
				}
				else{
					a = c;
				}
				i++;
			}
			cout << "\t-----Failed to get root for the given number of itrations\n";
		}
};

class BiSection{
	double a, b, tolerance;
	int max_itrations;
	public:
		void setA(double a){
			this->a = a;
		}
		void setB(double b){
			this->b = b;
		}
		void setTolerance(double tolerance){
			this->tolerance = tolerance;
		}
		void setMax_itrations(int max_itrations){
			this->max_itrations = max_itrations;
		}
		void calulate()
		{
			int i = 1;
			cout << "#" << setw(15) << "a" << setw(15) << "b" <<setw(15) << "c" << setw(15) << "f(a)" << setw(15) << "f(b)" << setw(15) << "f(c)\n";
			while(i <= max_itrations){
				double FA = f(a); double FB = f(b);
				double c = (a + b) / 2;
				double FC = f(c);
				cout << i << setw(15) << a << setw(15) << b << setw(15) << c << setw(15) << FA << setw(15) << FB << setw(15) << FC << endl;
				if(FC == 0 or fabs(FC) < tolerance){
					cout << "X = " << c << endl;
					return;
				}
				if(FA * FC < 0){
					b = c;
				}
				else{
					a = c;
				}
				i++;
			}
			cout << "\t-----Failed to get root for the given number of itrations\n";
		}
		
};

int main(){
	while(true){
		string op;
		system("cls");
		cout << "---------[NM Project]---------\n";
	    cout << "Select the method you want to use\n";
    	cout << "- Press 1 for biseection Method\n";
    	cout << "- Press 2 for Newton Raphson Method\n";
    	cout << "- Press 3 for Regula Falsi Method\n";
    	cout << "- Press 4 for Secant Method\n";
    	cout << "- Press 5 for Trapezoidal Method\n";
    	cout << "- Press 6 for Jacobi Method\n";
    	cout << "- Press 7 for Gauss Seidel Method\n";
    	cout << "- Press 8 for Newton Forward Method\n";
    	cout << "- Press 9 for Newton Backward Method\n";
    	cout << "- Press 10 for Euler Method\n";
    	cout << "- Press 11 for Simpson 1/3 Method\n";
    	cout << "- Press 12 for Simpson 3/8 Method\n";
    	cout << "- Press 13 for Newton Divided Difference\n";
    	cout << "- Press 14 for Lagrange Method\n";
    	cout << "- Press 15 for Runga Kutta\n";
    	cout << "- Press 16 for Booles's Rule (n = 4)\n";
    	cout << "- Press 17 for Weedle's Rule (n = 6)\n";
    	cout << "- Press 0 to exit\n";
    	cout << "INPUT: ";
    	cin >> op;
    	if(op == "1"){
    		double a, b, tolerance;
			int max_itrations;
    		BiSection obj;
    		while(true){
    			system("cls");
    			cout << "---------[Bi-Section Method]---------\n";
    			cout << "Enter value of 'a': ";
        		cin >> a;
        		obj.setA(a);
        		cout << "Enter value of 'b': ";
        		cin >> b;
        		obj.setB(b);
        		cout << "Enter value of 'Tolerance': ";
        		cin >> tolerance;
        		obj.setTolerance(tolerance);
        		cout << "Enter the max number of itrations: ";
        		cin >> max_itrations;
        		obj.setMax_itrations(max_itrations);
        		obj.calulate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
    		
		}
		else if(op == "2"){
    		double a, b, tolerance;
			int max_itrations;
    		Newton_Raphson obj;
    		while(true){
    			system("cls");
    			cout << "---------[Newton Raphson Method]---------\n";
    			cout << "Enter value of 'a': ";
        		cin >> a;
        		obj.setA(a);
        		cout << "Enter value of 'Tolerance': ";
        		cin >> tolerance;
        		obj.setTolerance(tolerance);
        		cout << "Enter the max number of itrations: ";
        		cin >> max_itrations;
        		obj.setMax_itrations(max_itrations);
        		obj.calculate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
    		
		}
		else if(op == "3"){
    		double a, b, tolerance;
			int max_itrations;
    		Regula_Falsi obj;
    		while(true){
    			system("cls");
    			cout << "---------[Regula Falsi Method]---------\n";
    			cout << "Enter value of 'a': ";
        		cin >> a;
        		obj.setA(a);
        		cout << "Enter value of 'b': ";
        		cin >> b;
        		obj.setB(b);
        		cout << "Enter value of 'Tolerance': ";
        		cin >> tolerance;
        		obj.setTolerance(tolerance);
        		cout << "Enter the max number of itrations: ";
        		cin >> max_itrations;
        		obj.setMax_itrations(max_itrations);
        		obj.calulate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
    		
		}
		else if(op == "4"){
    		double a, b, tolerance;
			int max_itrations;
    		Secant obj;
    		while(true){
    			system("cls");
    			cout << "---------[Secant Method]---------\n";
    			cout << "Enter value of 'a': ";
        		cin >> a;
        		obj.setA(a);
        		cout << "Enter value of 'b': ";
        		cin >> b;
        		obj.setB(b);
        		cout << "Enter value of 'Tolerance': ";
        		cin >> tolerance;
        		obj.setTolerance(tolerance);
        		cout << "Enter the max number of itrations: ";
        		cin >> max_itrations;
        		obj.setMax_itrations(max_itrations);
        		obj.calulate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
    		
		}
		else if(op == "5"){
    		double a, b, tolerance;
    		Trapezoidal obj;
    		while(true){
    			system("cls");
    			cout << "---------[Trapezoidal Method]---------\n";
    			cout << "Enter value of 'Starting point': ";
        		cin >> a;
        		obj.setA(a);
        		cout << "Enter value of 'Ending point': ";
        		cin >> b;
        		obj.setB(b);
        		cout << "Enter value of 'Size of one interval': ";
        		cin >> tolerance;
        		obj.setInterval(tolerance);
        		obj.calulate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
    		
		}
		else if(op == "6"){
    		double a, b, c, e;
    		Jacobi obj;
    		while(true){
    			system("cls");
    			cout << "---------[Jacobi Method]---------\n";
    			cout << "Enter value of initial 'X' (Zero is good): ";
        		cin >> a;
        		obj.setX(a);
        		cout << "Enter value of initial 'X' (Zero is good): ";
        		cin >> b;
        		obj.setY(b);
        		cout << "Enter value of initial 'Z' (Zero is good)': ";
        		cin >> c;
        		obj.setZ(c);
        		cout << "Enter value of 'Tolerance': ";
        		cin >> e;
        		obj.setTolerance(e);
        		obj.calculate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
		}
		else if(op == "7"){
    		double a, b, c, e;
    		Gauss_Seidel obj;
    		while(true){
    			system("cls");
    			cout << "---------[Gauss Seidel Method]---------\n";
    			cout << "Enter value of initial 'X' (Zero is good): ";
        		cin >> a;
        		obj.setX(a);
        		cout << "Enter value of initial 'X' (Zero is good): ";
        		cin >> b;
        		obj.setY(b);
        		cout << "Enter value of initial 'Z' (Zero is good)': ";
        		cin >> c;
        		obj.setZ(c);
        		cout << "Enter value of 'Tolerance': ";
        		cin >> e;
        		obj.setTolerance(e);
        		obj.calculate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
		}
		else if(op == "8"){
    		double n;
    		Newton_Forward obj;
    		while(true){
    			system("cls");
    			cout << "---------[Newton Forward Method]---------\n";
    			cout << "Enter the total number of points: ";
        		cin >> n;
        		obj.setN(n);
        		obj.calculate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
		}
		else if(op == "9"){
    		double n;
    		Newton_Backward obj;
    		while(true){
    			system("cls");
    			cout << "---------[Newton Backward Method]---------\n";
    			cout << "Enter the total number of points: ";
        		cin >> n;
        		obj.setN(n);
        		obj.calculate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
		}
		else if(op == "10"){
    		double x, y, h, final;
    		Euler_Method obj;
    		while(true){
    			system("cls");
    			cout << "---------[Euler Method]---------\n";
    			cout << "Enter the initial 'X': ";
        		cin >> x;
        		obj.setX(x);
        		cout << "Enter the initial 'Y': ";
        		cin >> y;
        		obj.setY(y);
        		cout << "Enter the interval 'h': ";
        		cin >> h;
        		obj.setInterval(h);
        		cout << "Enter the 'X' at which you want 'Y': ";
        		cin >> final;
        		obj.setFinal(final);
        		obj.calculate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
		}
		else if(op == "11"){
    		double a, b, tolerance;
    		Simpson_1by3 obj;
    		while(true){
    			system("cls");
    			cout << "---------[Simpson 1/3 Method]---------\n";
    			cout << "Enter value of 'Starting point': ";
        		cin >> a;
        		obj.setA(a);
        		cout << "Enter value of 'Ending point': ";
        		cin >> b;
        		obj.setB(b);
        		cout << "Enter value of 'Size of one interval': ";
        		cin >> tolerance;
        		obj.setInterval(tolerance);
        		if(int((b - a)/tolerance) % 2 != 0){
        			cout << "n = " << (b - a)/tolerance << " which is not a multiple of 3\n";
        			system("pause");
        			break;
				}
        		obj.calulate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
    		
		}
		else if(op == "12"){
    		double a, b, tolerance;
    		Simpson_3by8 obj;
    		while(true){
    			system("cls");
    			cout << "---------[Simpson 3/8 Method]---------\n";
    			cout << "Enter value of 'Starting point': ";
        		cin >> a;
        		obj.setA(a);
        		cout << "Enter value of 'Ending point': ";
        		cin >> b;
        		obj.setB(b);
        		cout << "Enter value of 'Size of one interval': ";
        		cin >> tolerance;
        		obj.setInterval(tolerance);
        		if(int((b - a)/tolerance) % 3 != 0){
        			cout << "n = " << (b - a)/tolerance << " which is not a multiple of 3\n";
        			system("pause");
        			break;
				}
        		obj.calulate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
		}
		else if(op == "13"){
    		int n;
    		Newton_Divided obj;
    		while(true){
    			system("cls");
    			cout << "---------[Newton Divided Difference Method]---------\n";
    			cout << "Enter the total number of points: ";
        		cin >> n;
        		obj.setN(n);
        		obj.calculate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
		}
		else if(op == "14"){
    		int n;
    		Lagrange obj;
    		while(true){
    			system("cls");
    			cout << "---------[Lagraneg Method]---------\n";
    			cout << "Enter the total number of points: ";
        		cin >> n;
        		obj.setN(n);
        		obj.calculate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
		}
		else if(op == "15"){
			Runga_Kutta obj;
			while(true){
    			system("cls");
    			cout << "---------[Runga Kutta Method]---------\n";
        		obj.calculate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
		}
		else if(op == "16"){
    		double a, b, tolerance;
    		Booles_Rule obj;
    		while(true){
    			system("cls");
    			cout << "---------[Boole's Rule (n = 4)' Method]---------\n";
    			cout << "Enter value of 'Starting point': ";
        		cin >> a;
        		obj.setA(a);
        		cout << "Enter value of 'Ending point': ";
        		cin >> b;
        		obj.setB(b);
        		cout << "Enter value of 'Size of one interval': ";
        		cin >> tolerance;
        		obj.setInterval(tolerance);
        		if(int((b - a)/tolerance) % 4 != 0){
        			cout << "n = " << (b - a)/tolerance << " which is not a multiple of 4\n";
        			system("pause");
        			break;
				}
        		obj.calulate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
    		
		}
		else if(op == "17"){
    		double a, b, tolerance;
    		Weedles_Rule obj;
    		while(true){
    			system("cls");
    			cout << "---------[Weedle's Rule (n = 6)' Method]---------\n";
    			cout << "Enter value of 'Starting point': ";
        		cin >> a;
        		obj.setA(a);
        		cout << "Enter value of 'Ending point': ";
        		cin >> b;
        		obj.setB(b);
        		cout << "Enter value of 'Size of one interval': ";
        		cin >> tolerance;
        		obj.setInterval(tolerance);
        		if(int((b - a)/tolerance) % 6 != 0){
        			cout << "n = " << (b - a)/tolerance << " which is not a multiple of 4\n";
        			system("pause");
        			break;
				}
        		obj.calulate();
        		char op1;
        		cout << "\n- Press Y to Try again\n";
        		cout << "- Press any key to go back\n";
        		cout << "INPUT: ";
        		cin >> op1;
        		if(op1 == 'Y' or op1 == 'y'){
        			continue;
				}
				else{
					break;
		    	}
			}
    		
		}
		else if(op == "0"){
			cout << "\t [Program Ended]\n";
			break;
		}
		else{
			cout << "invalid Input\n";
			system("pause");
		}
	}
	
}