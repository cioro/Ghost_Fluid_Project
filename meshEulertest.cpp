#include"Mesh2D.hpp"
#include"Euler2D.hpp"
#include<vector>
#include<fstream>
#include<iostream>
#include<assert.h>
#include<string>
#include<blitz/array.h>

/*
//Reflective boundary function
Euler::W_state Refl_Left_Bound(Euler::W_state w){
 
     
  Euler::W_state w_result;
  
  w_result.rho =   w.rho;
  w_result.u   = - w.u;
  w_result.P   =   w.P;

  return w_result;
}

Euler::W_state Refl_Right_Bound(Euler::W_state w ){
 
 Euler::W_state w_result;  

  w_result.rho =   w.rho;
  w_result.u   = - w.u;
  w_result.P   =   w.P;

  return w_result;
  }*/

//Transmissive boundary function
Euler::W_state Trans_Left_Bound(Euler::W_state w){
 
 Euler::W_state w_result;
   
  w_result.rho =  w.rho;
  w_result.u   =  w.u;
  w_result.v   =  w.v;
  w_result.P   =  w.P;

  return w_result;
}


Euler::W_state Trans_Right_Bound(Euler::W_state w){
 
  
  Euler::W_state w_result;
  
  w_result.rho =   w.rho;
  w_result.u   =   w.u;
  w_result.v   =   w.v;
  w_result.P   =   w.P;

  return w_result;
  }
//Initial condition function

Euler::U_state Sod_xdir(double x, double y){

  const double x_0 = 0.3;
  //Should do this assert(x_0 < x_max && x_0 > x_min); to avoid errors
  //but I would have to add to more args to the fcn.

  Euler e;//To allow use of CfromP and the creation of wl and ul. Not good design.
  
  if (x<=x_0){
    
    const Euler::W_state wL(1.0,0.0,0.0,1.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(0.125,0.0,0.0,0.1);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }
}

Euler::U_state Sod_ydir(double x,double y){
  
  const double y_0 = 0.5;
  
  Euler e;

 if (y<=y_0){
    
   const Euler::W_state wL(1.0,0.0,0.0,1.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (y>y_0){
    const Euler::W_state wR(0.125,0.0,0.0,0.1);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }
}
/*
Euler::U_state initial_test3(double x){
  
  const double x_0 = 0.5;
  
  Euler e;

 if (x<=x_0){
    
    const Euler::W_state wL(1.0,0.0,1000.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(1.0,0.0,0.01);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


Euler::U_state initial_test4(double x){
  
  const double x_0 = 0.4;
  
  Euler e;

 if (x<=x_0){
    
    const Euler::W_state wL(5.99924,19.5975,460.894);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(5.99242,-6.19633,46.0950);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


Euler::U_state initial_test5(double x){
  
  const double x_0 = 0.8;
  
  Euler e;

 if (x<=x_0){
    
    const Euler::W_state wL(1.0,-19.59745,1000.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(1.0,-19.59745,0.01);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


//Exact riemann solver function

Euler::U_state Exact_solver(double x){

  Euler::U_state u_empty;

  return u_empty;
}
*/

int main(int argc, char* argv[]){
  //Parameters of the problem
 
  
  double x_min = 0, x_max = 1.0; //domain length
  double y_min = 0, y_max = 1.0;
  double cfl = 0.9;
  int ncells = 100;
  /*
  if(argc == 4){
    ncells = atof(argv[3]);
  }else{
  ncells = 100;
  } */
 int nGhost = 1;

 Mesh m(ncells, x_min, x_max,y_min,y_max,cfl, Sod_ydir, Trans_Left_Bound, Trans_Right_Bound, nGhost);
  //Mesh n(ncells, x_min, x_max, cfl, initial_test1, Trans_Left_Bound, Trans_Right_Bound, nGhost);
  // double dt_m;
  // double dt_n;
  //double T_max = 0.2;
 
 //Print mesh to file
  std::string file_init_u = "initial_test_mesh";
  std::string file_init_w = "initial_test_w_";
  m.save_u_state(file_init_u);
  m.save_w_state(file_init_w);
  // n.save_w_state(file_init_);

  //Need to test modified applyBC()
  m.applyBC();
  std::cout << "Testing boundary conditions. Output points on each of the boudaries"<<"\n";
  std::cout << "Left boundary with transmissive conditions" << "\n";
  std::cout << m.Bdata(0,5).rho << "\n";
  std::cout << "Right boundary with transmissive conditions" << "\n";
  std::cout << m.Bdata(nGhost,5).rho << "\n";
  std::cout << "Top boundary with transmissive conditions" << "\n";
  std::cout << m.Bdata(5,nGhost).rho << "\n";
  std::cout << "Bottom boundary with transmissive conditions" << "\n";
  std::cout << m.Bdata(5,0).rho << "\n";
  //Need to test calculate_dt();

  std::cout << "Testing calculate_dt()" << "\n";
  double dt;
  dt = m.Calculate_dt();
  std::cout << "The time step is : " << dt << "\n";

 
  /*

  std::string filename_HLLC;
  std::string filename_WAF;
   std::string limiter;
  limiter = argv[2];
  if(argv[1]==std::string("test1")){
    T_max = 0.2;
    m.reset(initial_test1);
    n.reset(initial_test1);
    filename_HLLC = "Test1_HLLC_";
    filename_WAF = "Test1_WAF_";
    std::cout <<"Performing test1" << "\n";
  }   
  if(argv[1]== std::string("test2")){ 
    T_max = 0.15;
    m.reset(initial_test2);
    n.reset(initial_test2);
    filename_HLLC = "Test2_HLLC_";
    filename_WAF = "Test2_WAF_";
    std::cout <<"Performing test2" << "\n";
  }
  if(argv[1]== std::string("test3")){ 
    T_max = 0.012;
    m.reset(initial_test3);
    n.reset(initial_test3);  
    filename_HLLC = "Test3_HLLC_";
    filename_WAF = "Test3_WAF_";

    std::cout <<"Performing test3" << "\n";
  }
  if(argv[1]== std::string("test4")){ 
    T_max = 0.035;
    m.reset(initial_test4);
    n.reset(initial_test4);
    filename_HLLC = "Test4_HLLC_";
    filename_WAF = "Test4_WAF_";
    std::cout <<"Performing test4" << "\n";
  }
  if(argv[1]== std::string("test5")){ 
    T_max = 0.035;
    m.reset(initial_test5);
    n.reset(initial_test5);
    filename_HLLC = "Test5_HLLC_";
    filename_WAF = "Test5_WAF_";
    std::cout <<"Performing test5" << "\n";
  }

  //Initialise mesh with reflective BC

    
  //Initialise flux vector
  std::vector<Euler::U_state> flux_HLLC(m.ncells+1);
  std::vector<Euler::U_state> flux_WAF(n.ncells+1);
  //Main loop over domain for HLLC 

  std::cout << "HLLC loop starts" << "\n";
  for(double i = 0; i<T_max; i+=dt_m){

    m.applyBC();
    dt_m = m.Calculate_dt();
    flux_HLLC = HLLC(m);
    Mesh_update(m,flux_HLLC,dt_m);

    //Debugging messages
    std::cout <<"dt " << dt_m <<" The time is "<< i << "\n";   

    for(int j = 29; j < 33; j++){
      std::cout << flux_HLLC[j].rho << "\t" << flux_HLLC[j].momentum << "\t" << flux_HLLC[j].energy << "\n";
       }
    std::cout << "\n "; 

    //std::string snapshot = "Snap";
    // m.save_w_state(snapshot, Exact_solver); 
    
    m.time++;
  }

  m.save_w_state(filename_HLLC,Exact_solver);

  //Loop over whole domain for WAF
  std::cout << "The WAF loops start" << "\n";
  for(double j = 0; j<T_max; j+=dt_n){  


    n.applyBC();
    dt_n = n.Calculate_dt();
    flux_WAF = WAF(n,dt_m,limiter);
    Mesh_update(n,flux_WAF,dt_n);  

 //Debugging messages
    std::cout <<" dt " << dt_n <<" The time is "<< j << "\n";   

    for(int j = 29; j < 33; j++){
      std::cout << flux_WAF[j].rho << "\t" << flux_WAF[j].momentum << "\t" << flux_WAF[j].energy << "\n";
       }
    std::cout <<n.time << "\n";
    std::cout << "\n "; 
  
    n.time++;
  }

  n.save_w_state(filename_WAF,Exact_solver);
*/
}
