//
//  material.h
//  myproject
//
//  Created by Jaekwang Kim on 3/4/20.
//  Copyright © 2020 Jaekwang Kim. All rights reserved.
//

#ifndef material_h
#define material_h



//Module for defining Misorientation
class Material {
    public :
    
    virtual ~Material(){
        cout << "Base Material Class Destructor" << endl;
        // This destructor will be vailid for all the drived class
    }
    
    //virtual function! will be re-defined as
    virtual double Compute_CoupledEnergy (double jump_x,
                                          double jump_y,
                                          double jump_z)
    {
        return 0.0;
    };
};



class Simple_Material : public Material{
    //The real implementation
    //Linear to Jump value and no material symmetry is considered.
public:
    double Compute_CoupledEnergy(double jump_x,
                                 double jump_y,
                                 double jump_z)
    {
        double s=1.0;
        
        double Misorientation =  sqrt(jump_x * jump_x
                                      + jump_y * jump_y
                                      + jump_z * jump_z);
        
        return s * Misorientation;
    }
};



#endif /* material_h */
