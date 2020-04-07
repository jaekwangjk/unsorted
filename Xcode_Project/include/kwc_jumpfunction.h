//
// kwc_jumpfunctoin.h
// myproject
//
// Class KWC Optimizer design J(\theta) of KWC to match given GB data set.
// The latter is distributed to
// a)  f(\eta)
// b)  \laplace \eta
// c)   g(\eta)|\grad \theta|
// using Newton's Method
//
// Created by Jaekwang Kim on 3/11/20.
// Copyright © 2020 Jaekwang Kim. All rights reserved.
//

#ifndef kwc_jumpfunction_h
#define kwc_jumpfunction_h


#include <iomanip>   //std::setw

//Simple Functions. These are declared following the class declaration
double compute_eta_cusp_with_J(double J,int max_iterator);
double compute_KWC_energy(double eta_cusp, double J);

class KWC_Optimizer
{
public:
   
    int n_data;
    string file_name;
    
    KWC_Optimizer(int n, string s)
    {
        std::cout <<"KWC Optimizer Class constructed" << std::endl;
        n_data=n;
        file_name = s;
        cout <<"File name "<< file_name <<std::endl;
    }
    
    void run();
    
private:
   
    
};

void KWC_Optimizer::run()
{
    
    double max_iterator=10000;
    
    //input
    double W_data[n_data]; // Given data
    
    //output
    double J[n_data]; //Objective Function
    double tilt_angle[n_data];; // unit [radian]
    double W_kwc[n_data];
    double eta_cusp[n_data];
    
    //read data
    ifstream readFile;
    readFile.open(file_name);
    
    if (readFile.is_open())
    {
        for(int i=0; i< n_data; i ++)
        {
            readFile >> tilt_angle[i] >> W_data[i];
        }
    }else
    {
        cout << "file cannot be open" <<endl; exit(1);
    }
    
    readFile.close();
    
    
    
    double F0; // error between W_exp and W_kwc
    double F1;
  
    double W_kwc_step0;
    double W_kwc_step1;
    
    double J_step0;
    double J_step1;
    
    double eta_cusp_t=0.0;
    
    
    //Energy is zero when misorientation is zero
    eta_cusp[0]=1.0;
    W_kwc[0]=0;
    J[0]=0;
    
    //Loop over data point
    for(int i=1; i<n_data; i++)
    {
        //[1] First two points guess, choices are arbitrary
        J_step0=W_data[i];
        J_step1=0.5 * W_data[i];
        
        eta_cusp_t = compute_eta_cusp_with_J(J_step0 , max_iterator);
        W_kwc_step0=compute_KWC_energy(eta_cusp_t,J_step0);
        F0= W_data[i]-W_kwc_step0;
        
        eta_cusp_t = compute_eta_cusp_with_J(J_step1 , max_iterator);
        W_kwc_step1=compute_KWC_energy(eta_cusp_t,J_step1);
        F1= W_data[i]-W_kwc_step1;
        
        //[2] Execute Newton iteration
        double small_tol =1e-9; //Newtonain solver tolerance
        
        do
        {
            //update J
            double J_temp = J_step1;
            //Compute next step
            J_step1 = J_step1-(J_step0-J_step1)/(F0-F1) * F1;
            J_step0 = J_temp;
            
            //update W_kwc
            W_kwc_step0 = W_kwc_step1;
            eta_cusp_t = compute_eta_cusp_with_J(J_step1 , max_iterator);
            W_kwc_step1=compute_KWC_energy(eta_cusp_t,J_step1);
            
            F0= W_data[i]-W_kwc_step0;
            F1= W_data[i]-W_kwc_step1;
            
        }while(fabs(F1)> small_tol );
        
        //save solutions
        eta_cusp[i]=eta_cusp_t;
        W_kwc[i]=W_kwc_step1;
        J[i]=J_step1;
        
    }//End:Loop over data point
    
    
    //Confirm
    /*
    for(int i=0; i<n_data; i++)
    {
        std::cout << "( tilt_angle, W_data , W_kwc, J[theta] )=" <<
        "("<< std::setw(5) << tilt_angle[i]*180/M_PI<<", "
        << std::setw(10)<< W_data[i] <<", "
        << std::setw(10)<<W_kwc[i]<<", "
        << std::setw(10)<< J[i]<<")"
        <<std::endl;
    }
    */
    
    
    std::ofstream pipe;
    pipe.open("GB_KWC_tilt.txt");
    
    pipe    << std::setw(15)  << "tiltAngle"
    << std::setw(15)  << "Jtheta"
    << std::setw(15)  << "W_data" << std::endl;
    
    for(int i=0; i<n_data; i++)
    {
        pipe << std::setw(15) << tilt_angle[i]
        << std::setw(15) << J[i]
        << std::setw(15) << W_data[i]
        << std::endl;
    }
}

//Class declaration ends here

//Simple functions

double compute_eta_cusp_with_J(double J,int max_iterator)
{
    //find eta from J(\theta) using jump condition
    //2*(1-eta_cusp) = s * g(1-eta_cusp) * J
    //It uses a Newton's method to iteratively solve
    double s=1.0;
    double f=FLT_MAX; //Error: It is a minimizing function
    double grad_f=FLT_MAX;
    double eta_t=0.9;
    double tol=1e-10;  // Newtonian solver criteria
    
    double eta_cusp=0.0;
    
    for(int l=0 ; l<max_iterator; l++)
    {
        f= 2.0*(1.0-eta_t) + s * log(1-eta_t) * J;
        grad_f = -2.0 -  s*J/(1-eta_t);
        
        //Newtonian Iteration
        eta_t = eta_t - f/grad_f;
        
        if(fabs(f)<tol)
        {eta_cusp=eta_t;
            break;}
    }
    
    if(fabs(f) > tol)
        cout << "Warning Eta_cusp value may not be found correctly " <<endl;
    
    return eta_cusp;
}

//KWC Toal Energy
double compute_KWC_energy(double eta_cusp, double J)
{
    double s=1.0;
    double W_kwc= (1.0-eta_cusp)*(1.0-eta_cusp)-s*log(1-eta_cusp)*J;
    return W_kwc;
}

#endif /* kwc_jumpfunctoin_h */
