#include"Mesh2D.hpp"
#include"Euler2D.hpp"
#include<vector>
#include<fstream>
#include<iostream>
#include<assert.h>
#include<string>
#include<blitz/array.h>


//Reflective boundary function
Euler::W_state Refl_Left_Bound(Euler::W_state w){
 
     
  Euler::W_state w_result;
  
  w_result.rho =   w.rho;
  w_result.u   = - w.u;
  w_result.v   = - w.v;
  w_result.P   =   w.P;

  return w_result;
}

Euler::W_state Refl_Right_Bound(Euler::W_state w ){
 
 Euler::W_state w_result;  

  w_result.rho =   w.rho;
  w_result.u   = - w.u;
  w_result.v   = - w.v;
  w_result.P   =   w.P;

  return w_result;
  }

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

  const double x_0 = 0.5;
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

Euler::U_state cylinder(double x,double y){
  
  const double x_0 = 1.0;
  const double y_0 = 1.0;
  const double r   = 0.4;
 
  Euler e;

  if ( ((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0)) <= r){
    
    const Euler::W_state wL(1.0,0.0,0.0,1.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }else  if ( ((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0)) > r  ){
    const Euler::W_state wR(0.125,0.0,0.0,0.01);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside initial cylindrical explosion fcn"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


Euler::U_state diagonal(double x,double y){
  
  const double x_0 = 0.4;
  
  Euler e;

  if (y < (1-x)){
    
    const Euler::W_state wL(1.0,0.0,0.0,1.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }else if (y >= (1-x)){
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


Euler::U_state collision_perpendicular_to_x_axis(double x,double y){
  
  const double x_0 = 0.3;
  const double x_1 = 0.7;
  
  Euler e;

 if (x<=x_0){
    
   const Euler::W_state wL(1.3764,0.394,0.0,1.5698);
   Euler::U_state uL = e.CfromP(wL);
   return uL;
 }else if(x_0 < x && x < x_1){
   const Euler::W_state wM(1.0,0.0,0.0,1.0);
   Euler::U_state uM = e.CfromP(wR);
   return uM;
 }else if (x >= x_1){
   const Euler::W_state wR(1.3764,-0.394,0.0,1.5698);
   Euler::U_state uR = e.CfromP(wR);
   return uR;
 }else{
   std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
   Euler::U_state u;
   return u;
 }

}

Euler::U_state collision_perpendicular_to_y_axis(double x,double y){
  
  const double y_0 = 0.3;
  const double y_1 = 0.7;
  
  Euler e;

 if (y<=y_0){
    
   const Euler::W_state wL(1.3764,0.394,0.0,1.5698);
   Euler::U_state uL = e.CfromP(wL);
   return uL;
 }else if(y_0 < y && y < y_1){
   const Euler::W_state wM(1.0,0.0,0.0,1.0);
   Euler::U_state uM = e.CfromP(wR);
   return uM;
 }else if (y >= y_1){
   const Euler::W_state wR(1.3764,-0.394,0.0,1.5698);
   Euler::U_state uR = e.CfromP(wR);
   return uR;
 }else{
   std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
   Euler::U_state u;
   return u;
 }

}

Euler::U_state collision_diagonal(double x,double y){
  
  const double x_0 = 0.3;
  const double x_1 = 0.7;
  
  Euler e;

  if (y <= (0.3-x)  ){
    
   const Euler::W_state wL(1.3764,0.394,0.0,1.5698);
   Euler::U_state uL = e.CfromP(wL);
   return uL;
  }else if((0.3-x) > y && y < (0.7-x)){
   const Euler::W_state wM(1.0,0.0,0.0,1.0);
   Euler::U_state uM = e.CfromP(wR);
   return uM;
  }else if (y >= (0.7-x)){
   const Euler::W_state wR(1.3764,-0.394,0.0,1.5698);
   Euler::U_state uR = e.CfromP(wR);
   return uR;
 }else{
   std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
   Euler::U_state u;
   return u;
 }

}

Euler::U_state collision_angle(double x,double y){
  
  const double x_0 = 0.3;
  const double x_1 = 0.7;
  
  Euler e;

  if (y <= (0.3-1.5*x)  ){
    
   const Euler::W_state wL(1.3764,0.394,0.0,1.5698);
   Euler::U_state uL = e.CfromP(wL);
   return uL;
  }else if((0.3-1.5*x) > y && y < (0.7-1.5*x)){
   const Euler::W_state wM(1.0,0.0,0.0,1.0);
   Euler::U_state uM = e.CfromP(wR);
   return uM;
  }else if (y >= (0.7-1.5*x)){
   const Euler::W_state wR(1.3764,-0.394,0.0,1.5698);
   Euler::U_state uR = e.CfromP(wR);
   return uR;
 }else{
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
 
  
  double x_min = 0.0, x_max = 1.0; //domain length
  double y_min = 0.0, y_max = 1.0;
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
  /* 
 for(int col = 0; col  < m.ncells+ 2*m.nGhost; col ++){

    std::cout << m.Bdata(col,102).rho << "\n";
  
    }*/
  
  double dt;
  
  double T_max = 0.2;
  std::string Snap = "Sod_x_dir_snap_";
  std::string slice_x_axis = "x_axis_slice_";
  std::string slice_y_axis = "y_axis_slice_";
  
  for(double t = 0; t<T_max; t+=dt){
    dt = m.Calculate_dt();
    m.applyBC();
    flux_and_update(m,dt,std::string("XY"));
    m.time++;
    dt = m.Calculate_dt();
    m.applyBC();
    flux_and_update(m,dt,std::string("YX"));
    m.time++;
   
    m.save_u_state(Snap);
    m.slice_x_axis(slice_x_axis);
    m.slice_y_axis(slice_y_axis);
  }
  std::string file_output = "Sod_x_dir_";
  m.save_u_state(file_output);

 /*
//--------------------Potential 2D WAF FUNCTION -------------------------------------
 double dt_dx = dt/m.dx;
 std::string x_sweep_output = "x_sweep_output";
 std::string y_sweep_output = "y_sweep_output";
 int f_idx = 0;
 for(int i = m.nGhost; i < m.ncells+m.nGhost; i++){
   std::cout << "This is the row  " << i <<"\n";
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

 } 
 
 m.applyBC();
 
 f_idx=0;
 sweep = std::string("y-sweep");
  for(int i = m.nGhost; i < m.ncells+m.nGhost; i++){
    std::cout << "This is the column  " << i <<"\n";
    //Selects the ith column all of the rows.
    input = m.Bdata(i, blitz::Range::all());
    std::cout <<"The data of the solution at the end of the array is : " << "\n";
    input(101).print();
    input(102).print();
    flux = WAF_1D(input,dt,m.dy,m.ncells,m.nGhost,limiter,sweep); 
   
    for(int row = m.nGhost; row < m.ncells + m.nGhost; row++){
      m.Bdata(i,row).rho = m.Bdata(i,row).rho - dt_dx*(flux(f_idx+1).rho-flux(f_idx).rho);
      m.Bdata(i,row).moment_u = m.Bdata(i,row).moment_u - dt_dx*(flux(f_idx+1).moment_u - flux(f_idx).moment_u);
      m.Bdata(i,row).moment_v = m.Bdata(i,row).moment_v - dt_dx*(flux(f_idx+1).moment_v - flux(f_idx).moment_v);
      m.Bdata(i,row).energy = m.Bdata(i,row).energy - dt_dx*(flux(f_idx+1).energy-flux(f_idx).energy);
      f_idx++;
    }
    f_idx=0;
    m.save_u_state(y_sweep_output);
    
 }  
  */

 //---------------------------------------------------------------------------

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
