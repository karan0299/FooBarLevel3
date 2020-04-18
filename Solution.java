/* package whatever; // don't place package name! */

import java.util.ArrayList;
import java.util.Collections;
/* Name of the class has to be "Main" only if the class is public. */
class Ideone
{
	public static void main (String[] args) throws java.lang.Exception
	{
		int[][] m = {{0, 2, 1, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0,0}, {0, 0,3,0, 4}, {0, 0, 0, 0, 0}};
		
		int sums[] = new int[m.length];
        ArrayList<Integer> dest = new ArrayList<Integer>();
        ArrayList<Integer> inter = new ArrayList<Integer>();
        
            
        for(int i=0;i<m.length;i++){
            int sum = 0;
            for(int j =0;j<m[i].length;j++){
                sum = sum + m[i][j];
            }
            if(sum == 0){
                dest.add(i);
            }
            else {
                inter.add(i);
            }
            sums[i] = sum;
        }
        
        if(dest.get(0) < inter.get(inter.size()-1)) {
        	rearrange(m,inter,dest);
        }
        
        for(int i=0;i<m.length;i++){    
			for(int j=0;j<m.length;j++){ System.out.print(m[i][j] + ", "); }
			System.out.println();
        }
        
        Fraction[][] Ninv = new Fraction[inter.size()][inter.size()];
        Fraction[][] R = new Fraction[inter.size()][dest.size()];
 
        
        for(int i=0;i<m.length;i++){
            int sum = 0;
            for(int j =0;j<m[i].length;j++){
                sum = sum + m[i][j];
            }
            sums[i] = sum;
        }
        
        for(int i=0;i<inter.size();i++){
            for(int j =0;j<inter.size();j++){
            	Fraction q = new Fraction(m[i][j],sums[i]);
            	Fraction q1 = new Fraction(0,1);
            	if(i == j) {
            		q1.setNumerator(1);
            	}
            	Ninv[i][j] = q1.subtract(q);
            	Ninv[i][j].reduce();
            }
        }
        
        for(int i=0;i<inter.size();i++){
            for(int j =0;j<dest.size();j++){
            	int k = j+inter.size();
            	R[i][j] = new Fraction(m[i][k],sums[i]);
            	R[i][j].reduce();
            }
        }
        
        Fraction[][] N = new Fraction[inter.size()][inter.size()];
        inverse(Ninv,N);
        

        Fraction[][] ans = new Fraction[inter.size()][dest.size()];
        
        int[] fans = new int[dest.size()+1];
        int[] deno = new int[dest.size()];
        
        for(int i=0;i<inter.size();i++){    
			for(int j=0;j<dest.size();j++){    
				ans[i][j]= new Fraction(0,1);      
				for(int k=0;k<inter.size();k++)      
				{      
					
					ans[i][j]= ans[i][j].add(N[i][k].multiply(R[k][j]));
					ans[i][j].reduce(); 
					if(i == 0) { deno[j] = ans[i][j].getDenominator(); }
				}
			}
        }
        
        
        int lcm = getlcm(deno);
        
        for(int i=0;i<dest.size();i++){ 
        	fans[i] = ans[0][i].getNumerator() * (lcm/ans[0][i].getDenominator());
        	System.out.println(fans[i]);
        }
        
        fans[dest.size()] = lcm;
        
        System.out.println(fans[dest.size()]);
	}
	
	static void rearrange(int[][] m,ArrayList<Integer> inter, ArrayList<Integer> dest){
		int nj = 0;
		while(dest.get(nj) < inter.get(inter.size()-1)){
			nj++;
		}
		nj--;
		while(nj>=0&&dest.get(nj) < inter.get(inter.size()-1)){
			for(int i=0;i<m.length;i++){
				int temp = m[dest.get(nj)][i];
				m[dest.get(nj)][i] = m[inter.get(inter.size()-1)][i];
				m[inter.get(inter.size()-1)][i] = temp;
			}
			for(int i=0;i<m.length;i++){
				int temp1 = m[i][dest.get(nj)];
				m[i][dest.get(nj)] = m[i][inter.get(inter.size()-1)];
				m[i][inter.get(inter.size()-1)] = temp1;
			}
			int temp3 = dest.get(nj);
			dest.set(nj,inter.get(inter.size()-1));
			inter.set(inter.size()-1,temp3);
			Collections.sort(dest);
			Collections.sort(inter);
			nj--;
		}
	}
	
	static void getCofactor(Fraction[][] A, Fraction[][] temp, int p, int q, int n) 
	{ 
	    int i = 0, j = 0; 
	  
	    // Looping for each element of the matrix 
	    for (int row = 0; row < n; row++) 
	    { 
	        for (int col = 0; col < n; col++) 
	        { 
	            //  Copying into temporary matrix only those element 
	            //  which are not in given row and column 
	            if (row != p && col != q) 
	            { 
	                temp[i][j++] = A[row][col]; 
	  
	                // Row is filled, so increase row index and 
	                // reset col index 
	                if (j == n - 1) 
	                { 
	                    j = 0; 
	                    i++; 
	                } 
	            } 
	        } 
	    } 
	} 
	
	
  
	/* Recursive function for finding determinant of matrix. 
	   n is current dimension of A[][]. */
	static Fraction determinant(Fraction[][] A, int n) 
	{ 
	    Fraction D = new Fraction(0,1); // Initialize result 
	  
	    //  Base case : if matrix contains single element 
	    if (n == 1) 
	        return A[0][0]; 
	  
	    Fraction[][] temp = new Fraction[A.length][A.length]; // To store cofactors 
	  
	    int sign = 1;  // To store sign multiplier 
	  
	     // Iterate for each element of first row 
	    for (int f = 0; f < n; f++) 
	    { 
	        // Getting Cofactor of A[0][f] 
	        getCofactor(A, temp, 0, f, n); 
	        D = D.add((A[0][f].multiply(determinant(temp, n - 1))).multiply(new Fraction(sign,1))); 
			D.reduce();
	        // terms are to be added with alternate sign 
	        sign = -sign; 
	    } 
	  
	    return D; 
	} 
	  
	// Function to get adjoint of A[N][N] in adj[N][N]. 
	static void adjoint(Fraction[][] A, Fraction[][] adj) 
	{ 
	    if (A.length == 1) 
	    { 
	        adj[0][0] = new Fraction(1,1); 
	        return; 
	    } 
	  
	    // temp is used to store cofactors of A[][] 
	    int sign = 1;
	    Fraction[][] temp = new Fraction[A.length][A.length]; 
	  
	    for (int i=0; i<A.length; i++) 
	    { 
	        for (int j=0; j<A.length; j++) 
	        { 
	            // Get cofactor of A[i][j] 
	            getCofactor(A, temp, i, j, A.length); 
	  
	            // sign of adj[j][i] positive if sum of row 
	            // and column indexes is even. 
	            sign = ((i+j)%2==0)? 1: -1; 
	  
	            // Interchanging rows and columns to get the 
	            // transpose of the cofactor matrix 
	            adj[j][i] = determinant(temp, A.length-1);
	            adj[j][i] = adj[j][i].multiply(new Fraction(sign,1));
	            adj[j][i].reduce();
	        } 
	    } 
	} 
	  
	// Function to calculate and store inverse, returns false if 
	// matrix is singular 
	static void inverse(Fraction[][] A, Fraction[][] inverse) 
	{ 
	    // Find determinant of A[][] 
	    // Find adjoint 
	    Fraction[][] adj = new Fraction[A.length][A.length]; 
	    adjoint(A, adj); 
		Fraction det = determinant(A, A.length);
	    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
	    for (int i=0; i<A.length; i++) 
	        for (int j=0; j<A.length; j++){ 
	            inverse[i][j] = adj[i][j].divide(det); 
	            inverse[i][j].reduce();
	        }
	}
	
	static int getlcm(int[] element_array) 
    { 
        int lcm_of_array_elements = 1; 
        int divisor = 2; 
          
        while (true) { 
            int counter = 0; 
            boolean divisible = false; 
              
            for (int i = 0; i < element_array.length; i++) { 
  
                // lcm_of_array_elements (n1, n2, ... 0) = 0. 
                // For negative number we convert into 
                // positive and calculate lcm_of_array_elements. 
  
                if (element_array[i] == 0) { 
                    return 0; 
                } 
                else if (element_array[i] < 0) { 
                    element_array[i] = element_array[i] * (-1); 
                } 
                if (element_array[i] == 1) { 
                    counter++; 
                } 
  
                // Divide element_array by devisor if complete 
                // division i.e. without remainder then replace 
                // number with quotient; used for find next factor 
                if (element_array[i] % divisor == 0) { 
                    divisible = true; 
                    element_array[i] = element_array[i] / divisor; 
                } 
            } 
  
            // If divisor able to completely divide any number 
            // from array multiply with lcm_of_array_elements 
            // and store into lcm_of_array_elements and continue 
            // to same divisor for next factor finding. 
            // else increment divisor 
            if (divisible) { 
                lcm_of_array_elements = lcm_of_array_elements * divisor; 
            } 
            else { 
                divisor++; 
            } 
  
            // Check if all element_array is 1 indicate  
            // we found all factors and terminate while loop. 
            if (counter == element_array.length) { 
                return lcm_of_array_elements; 
            } 
        } 
    } 
}


class Fraction {
 
    int numerator;
    int denominator;
 
    /**
    * Constructor
    * 
    * @param numr
    * @param denr
    */
    public Fraction(int numr, int denr) {
	numerator = numr;
	denominator = denr;
	reduce();
    }
 
    public int getNumerator() {
	return numerator;
    }
 
    public void setNumerator(int numerator) {
	this.numerator = numerator;
    }
 
    public int getDenominator() {
	return denominator;
    }
 
    public void setDenominator(int denominator) {
	this.denominator = denominator;
    }
 
    /**
    * Calculates gcd of two numbers
    * 
    * @param numerator
    * @param denominator
    * @return
    */
    public int calculateGCD(int numerator, int denominator) {
    	if (numerator % denominator == 0) {
                 return denominator;
            }
    	return calculateGCD(denominator, numerator % denominator);
	}
 
    /**
    * Reduce the fraction to lowest form
    */
    void reduce() {
    	int gcd = calculateGCD(numerator, denominator);
    	numerator /= gcd;
    	denominator /= gcd;
    }
 
    /**
    * Adds two fractions
    * 
    * @param fractionTwo
    * @return
    */
    public Fraction add(Fraction fractionTwo) {
    	int numer = (numerator * fractionTwo.getDenominator()) + 
                                (fractionTwo.getNumerator() * denominator);
    	int denr = denominator * fractionTwo.getDenominator();
    	return new Fraction(numer, denr);
    }
 
    /**
    * Subtracts two fractions
    * 
    * @param fractionTwo
    * @return
    */
    public Fraction subtract(Fraction fractionTwo) {
        int newNumerator = (numerator * fractionTwo.denominator) - 
                                 (fractionTwo.numerator * denominator);
    	int newDenominator = denominator * fractionTwo.denominator;
    	Fraction result = new Fraction(newNumerator, newDenominator);
    	return result;
    }
 
    /**
    * Multiples two functions
    * 
    * @param fractionTwo
    * @return
    */
    public Fraction multiply(Fraction fractionTwo) {
    	int newNumerator = numerator * fractionTwo.numerator;
    	int newDenominator = denominator * fractionTwo.denominator;
    	Fraction result = new Fraction(newNumerator, newDenominator);
    	return result;
    }
 
    /**
    * Divides two fractions
    * 
    * @param fractionTwo
    * @return
    */
    public Fraction divide(Fraction fractionTwo) {
    	int newNumerator = numerator * fractionTwo.getDenominator();
    	int newDenominator = denominator * fractionTwo.numerator;
    	Fraction result = new Fraction(newNumerator, newDenominator);
    	return result;
    }
}
