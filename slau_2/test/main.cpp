#include <iostream>
#include "vector.hpp"
#include "matrix.hpp"
#include "csr_matrix.hpp"


int main() {
    try {

        std::cout << "Vector Test" << std::endl;
        Vector v1(3);
        v1[0] = 1.0; v1[1] = 2.0; v1[2] = 3.0;

        Vector v2 = v1 * 2.0;
        std::cout << "Vector * 2.0: " << v2[0] << " " << v2[1] << " " << v2[2] << std::endl;

        double dot = v1 * v1;
        std::cout << "(v1*v1): expect 14:  " << dot << std::endl;


        std::cout << "\n Matrix Test " << std::endl;
        Matrix m(2, 3);
        m(0, 0) = 1; m(0, 1) = 2; m(0, 2) = 3;
        m(1, 0) = 4; m(1, 1) = 5; m(1, 2) = 6;

        Vector v_res = m * v1;
        std::cout << "Matrix * Vector: (expect [14, 32]) " << v_res[0] << " " << v_res[1] << std::endl;


        std::cout << "\nCSR Test " << std::endl;
        CSRMatrix cm(3, 3);
    
    
    	cm.add_element(0, 0, 5.0);
    	cm.add_element(0, 2, 3.0);
    	cm.add_element(1, 1, 8.0);
    	cm.add_element(2, 0, 2.0);
    	cm.add_element(2, 2, 7.0);
    
  
    	cm.finalize();
    

        Vector v_res2 = cm * v1;
        std::cout << "CSR * Vector: (expect [14, 16,23]) " << v_res2[0] << " " << v_res2[1] << " "<<v_res2[2]<< std::endl;


    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

