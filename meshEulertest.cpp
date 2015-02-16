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
 int nGhost = 2;

 Mesh m(ncells, x_min, x_max,y_min,y_max,cfl, Sod_xdir, Trans_Left_Bound, Trans_Right_Bound, nGhost);
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

  blitz::Array<Euler::U_state,1> flux(m.ncells+1);
  std::string limiter = "minmod";
  std::string sweep = "x-sweep";
  //This selects all the columns of the 5th row
  blitz::Array<Euler::U_state,1> input = m.Bdata(blitz::Range::all(),5); ; 
  input.resize(m.ncells+2*m.nGhost);
  flux = WAF_1D(input,dt,m.dx,m.ncells,m.nGhost,limiter,sweep); 

 for(int i = 0; i < m.ncells + 1; i++){
    std::cout <<"This is flux : " << i << "\n";
    //input(i).print();
    flux(i).print();
  }

/*
//--------------------Potential 2D WAF FUNCTION -------------------------------------
 double dt_dx = dt/m.dx;
 std::string x_sweep_output = "x_sweep_output";
 std::string y_sweep_output = "y_sweep_output";
 int f_idx = 0;
 for(int i = m.nGhost; i < m.ncells+m.nGhost; i++){
   std::cout << "This is the first row calculation " << i <<"\n";
   //Selects all the colums of the i th row.
    input = m.Bdata(blitz::Range::all(),i);
    //Calculate x-fluxes
    flux = WAF_1D(input,dt,m.dx,m.ncells,m.nGhost,limiter,sweep); 
    //Update solution U_dt+1/2
    for(int col= m.nGhost; col < m.ncells + m.nGhost; col++){
    m.Bdata(col,i).rho = m.Bdata(col,i).rho - dt_dx*(flux(f_idx+1).rho-flux(f_idx).rho);
    m.Bdata(col,i).moment_u = m.Bdata(col,i).moment_u - dt_dx*(flux(f_idx+1).moment_u - flux(f_idx).moment_u);
    m.Bdata(col,i).moment_v = m.Bdata(col,i).moment_v - dt_dx*(flux(f_idx+1).moment_v - flux(f_idx).moment_v);
    m.Bdata(col,i).energy = m.Bdata(col,i).energy - dt_dx*(flux(f_idx+1).energy-flux(f_idx).energy);
    f_idx++;
    }
    f_idx = 0;
    m.save_u_state(x_sweep_output);
    
    sweep = std::string("y-sweep");
    input = m.Bdata(i, blitz::Range::all());
    flux = WAF_1D(input,dt,m.dx,m.ncells,m.nGhost,limiter,sweep); 
    for(int row= m.nGhost; row < m.ncells + m.nGhost; row++){
    m.Bdata(i,row).rho = m.Bdata(i,row).rho - dt_dx*(flux(f_idx+1).rho-flux(f_idx).rho);
    m.Bdata(i,row).moment_u = m.Bdata(i,row).moment_u - dt_dx*(flux(f_idx+1).moment_u - flux(f_idx).moment_u);
    m.Bdata(i,row).moment_v = m.Bdata(i,row).moment_v - dt_dx*(flux(f_idx+1).moment_v - flux(f_idx).moment_v);
    m.Bdata(i,row).energy = m.Bdata(i,row).energy - dt_dx*(flux(f_idx+1).energy-flux(f_idx).energy);
    f_idx++;
    }
    f_idx=0;
     m.save_u_state(y_sweep_output);
 }  
 //---------------------------------------------------------------------------
 */
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

  //Time Loop.
  std::cout << "The time loop starts" << "\n";
  for(double j = 0; j<T_max; j+=dt){  


    //m.applyBC();
    //dt = m.Calculate_dt();
    //m.time++;
    //Calculate flux and update. X-sweep,Y-sweep
    //m.applyBC();
    //dt=m.Calculate_dt();
    //m.time++;
    //Calculate flux and update. Y-sweep, X-sweep

    //Register output to file to take snaps of time evolution of the solution. 

 }

  m.save_w_state(filename_WAF,Exact_solver);
*/
}
