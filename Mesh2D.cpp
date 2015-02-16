#include<vector>
#include"Mesh2D.hpp"
#include<blitz/array.h>
#include<fstream>
#include<cmath>
#include<iostream> 
#include<string>
#include<sstream>
#include<algorithm>
#include"Euler2D.hpp"

//Initialises the parameters of the grid and fills the vector of primitives with the initial conditions
Mesh::Mesh(){};
Mesh::Mesh(int Ancells, double Ax_min, double Ax_max,double Ay_min,double Ay_max,double Acfl, Euler::U_state (*f)(double x,double y), Euler::W_state (*b1)(Euler::W_state  w), Euler::W_state (*b2)(Euler::W_state w), int AnGhost) : ncells(Ancells), x_min(Ax_min),x_max(Ax_max),y_min(Ay_min),y_max(Ay_max),cfl(Acfl), time(0), boundary1(b1), boundary2(b2), nGhost(AnGhost)
 { 
   ptr_euler = new Euler();
   dx = (x_max-x_min)/(double)ncells;
   dy = (y_max-y_min)/(double)ncells;

   Bdata.resize(ncells + 2*nGhost,ncells + 2*nGhost);
   xaxis.resize(ncells + 2*nGhost);
   yaxis.resize(ncells + 2*nGhost);
   
   int x_counter = 0;//Used to correctly calculate value of xaxis (value starts at 0 not nGhost,
   int y_counter = 0;//but first element is nGhost; Difference between array index and represented value
   

   //fill in xaxis
   for(int i = nGhost; i <(ncells+nGhost); i++){
     xaxis(i) = x_min + x_counter*dx;
     x_counter++;
     
     for(int j = nGhost; j<(ncells+nGhost); j++){
       yaxis(j) = y_min + y_counter*dy;
       y_counter++;
       Bdata(i,j) = f(xaxis(i),yaxis(j));
		     }
     y_counter=0;
     
     }
 }

//Destructor
Mesh::~Mesh()
{
  delete ptr_euler;
  /* delete data;
  delete axis;
  delete data2;*/
}
/*
void Mesh::reset(Euler::U_state (*f)(double x)){

  std::vector<Euler::U_state>::iterator itdata = data.begin()+nGhost;
   std::vector<Euler::U_state>::iterator itdata2= data2.begin()+nGhost;
   std::vector<double>::iterator itaxis= axis.begin()+nGhost;

  
   for(int i = 0; i<ncells; i++){
     (*itaxis) = x_min + i*dx;
     (*itdata) = f(*itaxis);
  
     itaxis++;
     itdata++;
    
     }

     };
*/
//Prints vector of conserved variables to screen
void Mesh::print()const
{
  
  //std::vector<Euler::U_state>::const_iterator itdata = data.begin()+nGhost;
  //std::vector<Euler::U_state>::const_iterator itdata2= data2.begin()+nGhost;
  // std::vector<double>::const_iterator itaxis= axis.begin()+nGhost;
  

  for (int i = nGhost; i<(ncells+nGhost); i++){
    for (int j = nGhost; j<(ncells+nGhost); j++){
      std::cout << xaxis(i) << "\t" << yaxis(j) << "\t" << Bdata(i,j).rho <<"\t"<< Bdata(i,j).moment_u <<"\t" << Bdata(i,j).moment_v <<"\t" << Bdata(i,j).energy << "\n";
    }
    std::cout <<"\n";
  }
     //std::cout << *itaxis << "\t";
     //(*itdata).print();
     //itaxis++;
     //itdata++;
   
  
}

//Print to a file the 1D vector of conserved variables and the exact solution
void Mesh::save_u_state(std::string filename)const
{
  std::string dir = "data/";
  std::stringstream ss;
  ss << dir << filename << time;
  std::string tmppath = ss.str();
  
  std::cout << "CREATING FILE \n";
  FILE * outfile = fopen(tmppath.c_str(),"w");
  //std::vector<Euler::U_state>::const_iterator itdata= data.begin()+nGhost;
  //std::vector<double>::const_iterator itaxis= axis.begin()+nGhost;
  for(int i=nGhost; i<ncells+nGhost; i++)
    {
      for(int j=nGhost; j<ncells+nGhost; j++){
	fprintf(outfile, "%.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n", xaxis(i), yaxis(j),Bdata(i,j).rho, Bdata(i,j).moment_u,Bdata(i,j).moment_v,Bdata(i,j).energy);
      }
      fprintf(outfile,"\n");
  }
  fclose(outfile);
  
}

//Prints to a file the 1D vector of primitive variables and the exact solution 
void Mesh::save_w_state(std::string filename)const
{
  std::string dir = "data/";
  std::stringstream ss;
  ss << dir << filename << time;
  std::string tmppath = ss.str();
  
  std::cout << "CREATING FILE \n";
  FILE * outfile = fopen(tmppath.c_str(),"w");
  //std::vector<Euler::U_state>::const_iterator itdata= data.begin()+nGhost;
  //std::vector<double>::const_iterator itaxis= axis.begin()+nGhost;
 
  Euler::W_state w_print;

for(int i=nGhost; i<ncells+nGhost; i++)
    {
      for(int j=nGhost; j<ncells+nGhost; j++){
	w_print = ptr_euler->PfromC(Bdata(i,j));

	fprintf(outfile, "%.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \n", xaxis(i), yaxis(j),w_print.rho,w_print.u,w_print.v,w_print.P);
      }
      fprintf(outfile,"\n");
      
    }
      fclose(outfile);
  
}

//Implements the boundary conditions. The actual boundary condition function should be in the main file
void Mesh::applyBC(){

  //Left Boundary of Square Mesh
  //Loop over Ghost left columns
  for(int x_dir = 0; x_dir < nGhost; x_dir++){

    //Loop over each row
    for(int y_dir = nGhost; y_dir<ncells+nGhost; y_dir++){

      Euler::W_state w_Left_End = ptr_euler -> PfromC(Bdata(x_dir+nGhost,y_dir));
      Euler::W_state w_BC_Left = boundary1(w_Left_End);
      Bdata(nGhost-1-x_dir,y_dir) = ptr_euler-> CfromP(w_BC_Left);

    }
  }
  //Right Boundary of Square Mesh
  //Loop over Ghost right columns
  for(int x_dir = 0; x_dir < nGhost; x_dir++ ){
    //Loop over each row
    for (int y_dir = nGhost; y_dir<ncells+nGhost; y_dir++){
      Euler::W_state w_Right_End = ptr_euler -> PfromC(Bdata ((nGhost+ncells-1)-x_dir,y_dir));
	Euler::W_state w_BC_Right = boundary2(w_Right_End);
	Bdata(x_dir+(nGhost+ncells),y_dir)= ptr_euler -> CfromP(w_BC_Right);
      }
  }


  //Top Boundary of Square Mesh
  //Loop over Ghost top rows
  for(int y_dir = 0; y_dir < nGhost; y_dir++){

    //Loop over each column
    for(int x_dir = nGhost; x_dir<ncells+nGhost; x_dir++){

      Euler::W_state w_Left_End = ptr_euler -> PfromC(Bdata(x_dir,y_dir+nGhost));
      Euler::W_state w_BC_Left = boundary1(w_Left_End);
      Bdata(x_dir,nGhost-1-y_dir) = ptr_euler-> CfromP(w_BC_Left);
      
      //data[nGhost-1-j]= ptr_euler -> CfromP(w_BC_Left);

    }
  }

  //Bottom Boundary of Square Mesh
  //Loop over Ghost bottom rows
  for(int y_dir = 0; y_dir < nGhost; y_dir++ ){
    //Loop over each column
    for (int x_dir = nGhost; x_dir < ncells+nGhost; x_dir++){
	Euler::W_state w_Right_End = ptr_euler -> PfromC(Bdata(x_dir,(ncells+nGhost)-y_dir));
	Euler::W_state w_BC_Right = boundary2(w_Right_End);
	Bdata(x_dir,(nGhost + ncells)+y_dir)= ptr_euler -> CfromP(w_BC_Right);
      }
  }




}

//Calculates the adquate size of the time step dt. See page 183 from Toro(ed.2009)
double Mesh::Calculate_dt(){

  double speed_x=0.0;
  double speedtemp_x=0.0;
  double speed_y=0.0;
  double speedtemp_y=0.0;
  double min_coef = 0.0;//This is min of (S_x/dx,S_y/dy)-The minimum of the wave speed over the space step.
  for(int i=nGhost; i<nGhost+ncells;i++){
    for(int j = nGhost; j<nGhost+ncells;j++){

      Euler::W_state w = ptr_euler->PfromC(Bdata(i,j));
      speedtemp_x = ptr_euler->a(w) + fabs(w.u);
      speedtemp_y = ptr_euler->a(w) + fabs(w.v);
    
      if(fabs(speedtemp_x) > fabs(speed_x)){
	speed_x = speedtemp_x;
      }
      if(fabs(speedtemp_y) > fabs(speed_y)){
	speed_x = speedtemp_y;
      }          

    }


  }

  min_coef = std::min((dx/speed_x),(dy/speed_y));
  double dt;
  //If time < 5 then dt = 0.2
  
 if(time < 10){
    double  cfl_init = 0.2;
    dt= cfl_init*min_coef;
    
    std::cout << "Inside calculate dt function, cfl = " << cfl_init << "\n"; 

    return dt;

    }
  std::cout << "Inside calculate dt function, cfl = " << cfl << "\n";

  dt=(cfl*dx)*min_coef;
  return dt;
}
/*
//Flux  and update function for 2D Euler Solver. Using the WAF TVD Scheme. See Toro chp.16
void flux_and_update(Mesh &m, double dt, std::string limiter, std::string sweep_order){

//Logic to choose sweep order
//x-sweep
   //WAF matrix sweep loop. Loop over all the rows
   //1D WAF x_dir fluxes
//n+1/2 update
//y-sweep
     //WAF matrix sweep loop. Loop over all the columns
     //1D WAF y_dir fluxes
//n+1 update

}




//HLLC flux calculator (this is a free function not a member function of class Mesh)
//Check Toro(ed.2009) p.331 for summary of HLLC method
std::vector<Euler::U_state> HLLC(Mesh &m){
 
  //Total vector of fluxes
  std::vector<Euler::U_state> flux(m.ncells+1);
  
  double gamma = m.ptr_euler->gamma;
  //Variables used

  double P_star,P_L,P_R; //Presures
  double rho_star,rho_L,rho_R;//densities
  double u_L,u_R;//particle/gas speed in cell

  double a_L, a_R;// sound speed in cell
  double S_L, S_R,S_star;// wave speed in cell left wave, right wave, contact wave

  Euler::U_state U_star_L, U_star_R; // Star states of conserved var
  
  double P_pvrs;
  double rho_bar, a_bar;//average density and average sound speed of left and right cells

  Euler::W_state w_temp_left;
  Euler::W_state w_temp_right;

  double q_R,q_L;

  Euler::U_state U_state_L;
  Euler::U_state U_state_R;
  Euler::W_state W_L;
  Euler::W_state W_R;
  Euler::U_state U_state_L_star;
  Euler::U_state U_state_R_star;

  double star_coef_left; // The coeficient in eq. 10.73 from Toro(ed.2009); (S_k-u_k)/(S_k-u_star_k);
  double star_coef_right; // The coeficient in eq. 10.73 from Toro(ed.2009);
  
  std::vector<Euler::U_state>::iterator itflux = flux.begin();
  
  // m.data[m.nGhost].print();
  // m.data[0].print();

  //Loop over whole domain
  for(int i = m.nGhost-1; i < m.ncells+m.nGhost; i++){
    
    //Select U_state and initialise W_state

    U_state_L = m.data[i];
    U_state_R = m.data[i+1];

    W_L = m.ptr_euler->PfromC(U_state_L);
    W_R = m.ptr_euler->PfromC(U_state_R);


    //---------Pressure estimate-------------------------

    P_L = W_L.P;
    P_R = W_R.P;

    u_L = W_L.u;
    u_R = W_R.u;

    rho_L = W_L.rho;
    rho_R = W_R.rho;
   
    a_L =m.ptr_euler->a(W_L);
    a_R = m.ptr_euler->a(W_R);

    // std::cout <<"Inside the HLLC function" << "\n";
        
    rho_bar = 0.5*(rho_L + rho_R);
    a_bar = 0.5*(a_L + a_R);

    P_pvrs = 0.5*(P_L + P_R)-0.5*(u_R-u_L)*rho_bar*a_bar;

    P_star = std::max(0.0,P_pvrs);
      
    //---------------------------------------------------
   
    //----------Wave speed estimate----------------------
    
    //Calculate q_R

    if(P_star <= P_R){
      q_R = 1.0;
    }
    else{
      q_R = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_R)-1));
    }
        
    //Calculate q_L
    if(P_star <= P_L){
      q_L = 1.0;
    }
    else{
      q_L = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_L)-1));
    }
 
    //Calculate S_R and S_L
    S_L = u_L -a_L*q_L;
    S_R = u_R + a_R*q_R;
    
    double numerator = P_R - P_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R); 
    double denominator = rho_L*(S_L-u_L)-rho_R*(S_R-u_R);

    S_star = numerator/denominator; 

    //---------------------------------------------------

    //HLLC flux

    if ( 0.0 <= S_L){
      
      Euler::U_state F_L;
      F_L = m.ptr_euler->flux(U_state_L);
      *itflux = F_L;
    }

    if( (S_L <= 0.0) && (S_star >= 0.0)){

      star_coef_left = rho_L*((S_L-u_L)/(S_L-S_star));
      
      //Calculate U_state_L_star
      U_state_L_star.rho = star_coef_left;
      U_state_L_star.momentum = star_coef_left*S_star;
      U_state_L_star.energy = star_coef_left*(U_state_L.energy/U_state_L.rho + (S_star - u_L)*(S_star + P_L/(rho_L*(S_L-u_L))));
      
      Euler::U_state F_L;
      F_L = m.ptr_euler->flux(U_state_L);
     //  Consider overloading the + operator to write this in one line 
      (*itflux).rho = F_L.rho + S_L*(U_state_L_star.rho -U_state_L.rho);
      (*itflux).momentum = F_L.momentum + S_L*(U_state_L_star.momentum -U_state_L.momentum);
      (*itflux).energy = F_L.energy + S_L*(U_state_L_star.energy -U_state_L.energy);
   
							
    }

    if( (S_star <= 0.0) && (S_R >= 0.0)){


      star_coef_right = rho_R*((S_R-u_R)/(S_R-S_star));
      
      //Calculate U_state_R_star
      U_state_R_star.rho = star_coef_right;
      U_state_R_star.momentum = star_coef_right*S_star;
      U_state_R_star.energy = star_coef_right*(U_state_R.energy/U_state_R.rho + (S_star - u_R)*(S_star + P_R/(rho_R*(S_R-u_R))));

      Euler::U_state F_R;    
      F_R = m.ptr_euler->flux(U_state_R);
//      Consider overloading the + operator to write this in one line 
      (*itflux).rho = F_R.rho + S_R*(U_state_R_star.rho -U_state_R.rho);
      (*itflux).momentum = F_R.momentum + S_R*(U_state_R_star.momentum -U_state_R.momentum);
      (*itflux).energy = F_R.energy + S_R*(U_state_R_star.energy -U_state_R.energy);


    }

    if( 0.0 >= S_R){
      
      Euler::U_state F_R;
      F_R = m.ptr_euler->flux(U_state_R);
      *itflux = F_R;
      
    }
					       

    itflux++;
  }

  return flux;
}

void Mesh_update(Mesh &m, std::vector<Euler::U_state> &flux, double dt){
 
  double dt_dx = dt/m.dx;
  std::vector<Euler::U_state>::iterator itflux = flux.begin();
  
  for(int i = m.nGhost; i <m.ncells+ m.nGhost;i++){

    m.data[i].rho = m.data[i].rho - dt_dx*((*(itflux+1)).rho-((*itflux).rho));
    m.data[i].momentum = m.data[i].momentum - dt_dx*((*(itflux+1)).momentum - (*itflux).momentum);
    m.data[i].energy = m.data[i].energy - dt_dx*((*(itflux+1)).energy-((*itflux).energy));
    
    itflux++;

  }

}
*/
//Need to modify this 1D  WAF.Unsure whether to create WAF_x_dir and WAF_y_dir or just permute v and u and use one fcn. 
//This function is not strictly applicable to the 1D case. It includes the existance of a tangential velocity v, which gives
//rise to a sheer wave, thus a new wave is included in the WAF calculation (plus a new limiter -- see page 553 Toro ed.2009)
blitz::Array<Euler::U_state,1> WAF_1D(blitz::Array<Euler::U_state,1> input_data, double dt, double ds, double ncells,double nGhost,std::string limiter,std::string sweep){
 
  //Total vector of fluxes
  blitz::Array<Euler::U_state,1> flux(ncells+1);
  
  //Logic to switch u and v depending on the sweep
  //Need to rever this inversion when the result is calculated at the end. 
  //Otherwise, I'll be calculting only the x-sweep twice.
  if (sweep == std::string("x-sweep")){

  }else if(sweep == std::string("y-sweep")){
    for(int i= nGhost; i < ncells+nGhost; i++){
    input_data(i).moment_u = input_data(i).moment_v;
    input_data(i).moment_v = input_data(i).moment_u;
    }
  }else{
  std::cout <<"Sweep not specified. Fail to compute 1D WAF" << "\n";
  }

  //Given that the mesh is no longer an input argument to access Euler class function, an empty euler object is needed
  // to access the different data members and function members
  Euler e;
  double gamma = e.gamma;
  //Variables used

  double P_star,P_L,P_R; //Presures
  double rho_star,rho_L,rho_R;//densities
  double u_L,u_R;//particle/gas speed in cell

  double a_L, a_R;// sound speed in cell
  double S_L, S_R,S_star;// wave speed in cell left wave, right wave, contact wave

  Euler::U_state U_star_L, U_star_R; // Star states of conserved var
  
  double P_pvrs;
  double rho_bar, a_bar;//average density and average sound speed of left and right cells

  Euler::W_state w_temp_left;
  Euler::W_state w_temp_right;

  double q_R,q_L;

  Euler::U_state U_state_L;
  Euler::U_state U_state_R;
  Euler::W_state W_L;
  Euler::W_state W_R;
  Euler::U_state U_state_L_star;
  Euler::U_state U_state_R_star;

  double star_coef_left; // The coeficient in eq. 10.73 from Toro(ed.2009); (S_k-u_k)/(S_k-u_star_k);
  double star_coef_right; // The coeficient in eq. 10.73 from Toro(ed.2009);
  
  
  //Courant number for each wave speed
  double c_L,c_star,c_R,c_shear;
 

  //Ratio for each wave speed
  double r_L,r_star,r_R,r_shear;
  double dq_l,dq_l_right_interface, dq_l_left_interface;
  double dq_star, dq_star_right_interface, dq_star_left_interface;
  double dq_r, dq_r_right_interface, dq_r_left_interface;
  double dq_shear, dq_shear_right_interface, dq_shear_left_interface;

  //Limiter functions
  double minmod_l,minmod_star,minmod_r,minmod_shear;
  double superbee_l,superbee_star,superbee_r,superbee_shear;

  blitz::Array<Euler::U_state,1> left_interface;
  blitz::Array<Euler::U_state,1> right_interface;
  left_interface.resize(4);
  right_interface.resize(4);

  //Dummy variable used to calculate shear wave max speed
  double shear_speed;
  
  //Middle state used in WAF calculation to account for the presence of a shear wave.
  Euler::W_state w_middle_state;
  Euler::U_state F_ML;
  Euler::U_state F_MR;
  int f_idx = 0;
  std::cout <<"Inside WAF_1D flux function " << "\n";
  //Loop over whole domain
  for(int i = nGhost-1; i < ncells + nGhost; i++){

    std::cout << "Inside the main loop. Iteration: " << i << "\n";
    //Select U_state and initialise W_state

    U_state_L = input_data(i);
    U_state_R = input_data(i+1);

    W_L = e.PfromC(U_state_L);
    W_R = e.PfromC(U_state_R);
    
    //---------Pressure estimate-------------------------

    P_L = W_L.P;
    P_R = W_R.P;

    u_L = W_L.u;
    u_R = W_R.u;

    rho_L = W_L.rho;
    rho_R = W_R.rho;
   
    a_L = e.a(W_L);
    a_R = e.a(W_R);

    // std::cout <<"Inside the HLLC function" << "\n";
        
    rho_bar = 0.5*(rho_L + rho_R);
    a_bar = 0.5*(a_L + a_R);

    P_pvrs = 0.5*(P_L + P_R)-0.5*(u_R-u_L)*rho_bar*a_bar;

    P_star = std::max(0.0,P_pvrs);
      
    //---------------------------------------------------
   
    //----------Wave speed estimate----------------------
    
    //Calculate q_R

    if(P_star <= P_R){
      q_R = 1.0;
    }
    else{
      q_R = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_R)-1));
    }
        
    //Calculate q_L
    if(P_star <= P_L){
      q_L = 1.0;
    }
    else{
      q_L = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_L)-1));
    }
 
    //Calculate S_R and S_L
    S_L = u_L -a_L*q_L;
    S_R = u_R + a_R*q_R;
    
    double numerator = P_R - P_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R); 
    double denominator = rho_L*(S_L-u_L)-rho_R*(S_R-u_R);

    S_star = numerator/denominator; 

    //---------------------------------------------------

    //---------------HLLC fluxes -------------------------------------------------//

    //Calculate F_L     
      Euler::U_state F_L;
      F_L = e.flux(U_state_L);
       

     //Calculate F_L_star
      Euler::U_state F_L_star;
      star_coef_left = rho_L*((S_L-u_L)/(S_L-S_star));
      
      //---Calculate U_state_L_star
      U_state_L_star.rho = star_coef_left;
      U_state_L_star.moment_u = star_coef_left*S_star;
      U_state_L_star.moment_v = star_coef_left*W_L.v;
      U_state_L_star.energy = star_coef_left*(U_state_L.energy/U_state_L.rho + (S_star - u_L)*(S_star + P_L/(rho_L*(S_L-u_L))));
//      /* Consider overloading the + operator to write this in one line 
      F_L_star.rho = F_L.rho + S_L*(U_state_L_star.rho -U_state_L.rho);
      F_L_star.moment_u = F_L.moment_u + S_L*(U_state_L_star.moment_u - U_state_L.moment_u);
      F_L_star.moment_v = F_L.moment_v + S_L*(U_state_L_star.moment_v - U_state_L.moment_v);
      F_L_star.energy = F_L.energy + S_L*(U_state_L_star.energy -U_state_L.energy);
   
      //Calculate F_R
      Euler::U_state F_R;
      F_R = e.flux(U_state_R);
     					

      //Calculate F_R_star
      Euler::U_state F_R_star;
      star_coef_right = rho_R*((S_R-u_R)/(S_R-S_star));
      
      //Calculate U_state_R_star
      U_state_R_star.rho = star_coef_right;
      U_state_R_star.moment_u = star_coef_right*S_star;
      U_state_R_star.moment_v = star_coef_right*W_R.v;
      U_state_R_star.energy = star_coef_right*(U_state_R.energy/U_state_R.rho + (S_star - u_R)*(S_star + P_R/(rho_R*(S_R-u_R))));


//      /* Consider overloading the + operator to write this in one line 
      F_R_star.rho = F_R.rho + S_R*(U_state_R_star.rho -U_state_R.rho);
      F_R_star.moment_u = F_R.moment_u + S_R*(U_state_R_star.moment_u -U_state_R.moment_u);
      F_R_star.moment_v = F_R.moment_v + S_R*(U_state_R_star.moment_v -U_state_R.moment_v);
      F_R_star.energy = F_R.energy + S_R*(U_state_R_star.energy -U_state_R.energy);

      //---------------END of HLLC Fluxes, F_L, F_L_star, F_R, F_R_star------------------//


      //Compute the courant number
      c_L = dt*S_L/ds;
      c_star = dt*S_star/ds;
      c_R = dt*S_R/ds;

      //To calculate the Cfl number associated with the shear wave I need a shear wave speed.From p.553 from Toro
      // I can deduce that the contact wave and the shear wave have the same wave speed. 

      c_shear = dt*S_R/ds ;


      //Compute the ratios
      //-----compute dq, dq_up, dq_down for each wave at the interface. Up for 
      
      //dq for the interface we are trying to solve for i+1/2
      dq_l = U_state_L.rho-U_state_L_star.rho;
      dq_star = U_state_L_star.rho - U_state_R_star.rho;
      dq_shear = (e.PfromC(U_state_L_star)).v-(e.PfromC(U_state_R_star)).v;
      dq_r = U_state_R_star.rho - U_state_R.rho;

      //std::cout << "About to calculate the HLLC values at the right and left interface to our current interface" << "\n";

      left_interface = HLLC_U_state(input_data(i-1), input_data(i));
      right_interface = HLLC_U_state(input_data(i+1),input_data(i+2));
      
      //std::cout << " calculate left and right interfaces " << "\n";

      dq_l_left_interface = left_interface(0).rho-left_interface(1).rho;
      dq_star_left_interface = left_interface(1).rho - left_interface(2).rho;
      dq_r_left_interface = left_interface(2).rho - left_interface(3).rho;
      dq_shear_left_interface = (e.PfromC(left_interface(1))).v-(e.PfromC(left_interface(2))).v;


      dq_l_right_interface = right_interface(0).rho-right_interface(1).rho;
      dq_star_right_interface = right_interface(1).rho - right_interface(2).rho;
      dq_r_right_interface = right_interface(2).rho - right_interface(3).rho;
      dq_shear_right_interface =  (e.PfromC(right_interface(1))).v-(e.PfromC(right_interface(2))).v;
      
      //std::cout << "We finished calculating " << "\n";

      /*
      std::cout << "This is the " << i << " element" <<"\n";

      std::cout << "The 4 states at the interface are: " << "\n";
      U_state_L.print();
      U_state_L_star.print();
      U_state_R_star.print();
      U_state_R.print();
      std::cout << "The jumps are dq_l,dq_star ,dq_r " << dq_l <<"\t" << dq_star << "\t" << dq_r << "\t" << dq_shear << "\n";
      std::cout << " The jumps of the right interfacare are " << dq_l_right_interface <<"\t" << dq_star_right_interface << "\t" << dq_r_right_interface << "\t" << dq_shear_right_interface << "\n"; 
      std::cout << " The jumps of the left interfacare are " << dq_l_left_interface <<"\t" << dq_star_left_interface << "\t" << dq_r_left_interface << "\t" << dq_shear_left_interface << "\n"; 
      */

      /*
      std::cout << "The nearby states are " << "\n";
      std::cout << "The interface to the left of i+1/2 " << "\n";
      std::cout << "U_L : " << "\n";
      left_interface(0).print();
      std::cout << "U_L_star : " << "\n";
      left_interface(1).print();
      std::cout << "U_R_star : " << "\n";
      left_interface(2).print();
      std::cout << "U_R : " << "\n";
      left_interface(3).print();


      std::cout << "The interface to the right of i+1/2 " << "\n";
      std::cout << "U_L : " << "\n";
      right_interface(0).print();
      std::cout << "U_L_star : " << "\n";
      right_interface(1).print();
      std::cout << "U_R_star : " << "\n";
      right_interface(2).print();
      std::cout << "U_R : " << "\n";
      right_interface(3).print();
      */
      std::cout << "\n";
      
   
     
      //---------Computing ratios from d_q jumps--------------------//      
      if(dq_l > 0){
	if (c_L > 0){
	  r_L = dq_l_left_interface/dq_l;
	}
	else{
	  r_L = dq_l_right_interface/dq_l;
	}
      }
      else{
	minmod_l = 0.0;
	superbee_l = 0.0;
	r_L = 0.0;
      }

	if(dq_star > 0){

	  if (c_star > 0){
	    r_star = dq_star_left_interface/dq_star;
	  }
	  else{
	    r_star = dq_star_right_interface/dq_star;
	  }
	}
      else{
	minmod_star = 0.0;
	superbee_star = 0.0;
	r_star = 0.0;
      }

	  if(dq_r > 0){

	    if (c_R > 0){
	      r_R = dq_r_left_interface/dq_r;
	    }
	    else{
	      r_R = dq_r_right_interface/dq_r;
	    }
	  }
	  else{
	    minmod_r = 0.0;
	    superbee_r = 0.0;
	    r_R = 0.0;
	  } 

//----------Shear jump and ratio ---------------------------
	  if(dq_shear > 0){

	    if (c_shear > 0){
	      r_shear = dq_shear_left_interface/dq_shear;
	    }
	    else{
	      r_shear = dq_shear_right_interface/dq_shear;
	    }
	  }
	  else{
	    minmod_shear = 0.0;
	    superbee_shear = 0.0;
	    r_shear = 0.0;
	  } 
   
//----------------------------------------------------------



      //Compute limiter functions
      minmod_l = minmod(r_L,fabs(c_L));
      minmod_star = minmod(r_star, fabs(c_star));
      minmod_r = minmod(r_R,fabs(c_R));
      minmod_shear = minmod(r_shear,fabs(c_shear));

      superbee_l = superbee(r_L, fabs(c_L));
      superbee_star = superbee(r_star, fabs(c_star));
      superbee_r = superbee(r_R, fabs(c_R));
      superbee_shear = superbee(r_shear,fabs(c_shear));
      /*
      std::cout << "The ratios are r_L, r_star, r_R,r_shear " << r_L<<"\t" << r_star << "\t" << r_R << "\t"<<r_shear<<"\n"; 

      std::cout << "The minmod lim m_l, m_star, m_r, m_shear   " << minmod_l <<"\t" <<minmod_star << "\t" << minmod_r << "\t" << minmod_shear << "\n";  
      */
       if(minmod_shear > minmod_star){
	 
	 
	 w_middle_state.rho = U_state_R_star.rho;
	 w_middle_state.u = (e.PfromC(U_state_R_star)).u;
	 w_middle_state.v = W_L.v;
	 w_middle_state.P = (e.PfromC(U_state_R)).P;


	 F_ML = e.flux((e.CfromP(w_middle_state)));

	 flux(f_idx).rho =  0.5*(F_L.rho+F_R.rho)-0.5*(sign(c_L)*minmod_l*(F_L_star.rho-F_L.rho) +\
						   sign(c_shear)*minmod_shear*(F_ML.rho-F_L_star.rho) +\
						   sign(c_star)*minmod_star*(F_R_star.rho-F_ML.rho)+\
						   sign(c_R)*minmod_r*(F_R.rho-F_R_star.rho));

	 flux(f_idx).moment_u =  0.5*(F_L.moment_u+F_R.moment_u)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_u-F_L.moment_u) + \
						   sign(c_shear)*minmod_shear*(F_ML.moment_u-F_L_star.moment_u) +\
						   sign(c_star)*minmod_star*(F_R_star.moment_u-F_ML.moment_u)+\
						   sign(c_R)*minmod_r*(F_R.moment_u-F_R_star.moment_u));

	 flux(f_idx).moment_v =  0.5*(F_L.moment_v+F_R.moment_v)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_v-F_L.moment_v) + \
						   sign(c_shear)*minmod_shear*(F_ML.moment_v-F_L_star.moment_v) +\
						   sign(c_star)*minmod_star*(F_R_star.moment_v-F_ML.moment_v)+\
						   sign(c_R)*minmod_r*(F_R.moment_v-F_R_star.moment_v));
	 

	 flux(f_idx).energy =  0.5*(F_L.energy+F_R.energy)-0.5*(sign(c_L)*minmod_l*(F_L_star.energy-F_L.energy) + \
						   sign(c_shear)*minmod_shear*(F_ML.energy-F_L_star.energy) +\
						   sign(c_star)*minmod_star*(F_R_star.energy-F_ML.energy)+\
						   sign(c_R)*minmod_r*(F_R.energy-F_R_star.energy));

	 std::cout << "The shear wave is before the contact wave" << "\n"; 

       }else if(minmod_shear < minmod_star){
	 std::cout << "The shear wave is behind the contact wave." << "\n";
	 std::cout << "The resulting middle state is rho*l and v_R"<< "\n";

	 w_middle_state.rho = U_state_L_star.rho;
	 w_middle_state.u = (e.PfromC(U_state_R_star)).u;
	 w_middle_state.v = W_R.v;
	 w_middle_state.P = (e.PfromC(U_state_R)).P;



	 F_MR = e.flux((e.CfromP(w_middle_state)));

	 flux(f_idx).rho =  0.5*(F_L.rho+F_R.rho)-0.5*(sign(c_L)*minmod_l*(F_L_star.rho-F_L.rho) + \
						   sign(c_shear)*minmod_shear*(F_MR.rho-F_L_star.rho) +\
						   sign(c_star)*minmod_star*(F_R_star.rho-F_MR.rho)+\
						   sign(c_R)*minmod_r*(F_R.rho-F_R_star.rho));

	 flux(f_idx).moment_u =  0.5*(F_L.moment_u+F_R.moment_u)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_u-F_L.moment_u) + \
								  sign(c_shear)*minmod_shear*(F_MR.moment_u-F_L_star.moment_u) + \
								  sign(c_star)*minmod_star*(F_R_star.moment_u-F_MR.moment_u)+ \
								  sign(c_R)*minmod_r*(F_R.moment_u-F_R_star.moment_u));

	 flux(f_idx).moment_v =  0.5*(F_L.moment_v+F_R.moment_v)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_v-F_L.moment_v) + \
								  sign(c_shear)*minmod_shear*(F_MR.moment_v-F_L_star.moment_v) + \
								  sign(c_star)*minmod_star*(F_R_star.moment_v-F_MR.moment_v)+ \
								  sign(c_R)*minmod_r*(F_R.moment_v-F_R_star.moment_v));

	 	 flux(f_idx).energy =  0.5*(F_L.energy+F_R.energy)-0.5*(sign(c_L)*minmod_l*(F_L_star.energy-F_L.energy) + \
								    sign(c_shear)*minmod_shear*(F_MR.energy-F_L_star.energy) + \
								    sign(c_star)*minmod_star*(F_R_star.energy-F_MR.energy)+	\
								    sign(c_R)*minmod_r*(F_R.energy-F_R_star.energy));
		 

       }else if(minmod_shear == minmod_star){
	 //No middle state exist. Normal WAF calculation.
	 std::cout << "Standard WAF being calculated.No middle state" << "\n";
	 flux(f_idx).rho = 0.5*(F_L.rho+F_R.rho)-0.5*(sign(c_L)*minmod_l*(F_L_star.rho-F_L.rho) + \
						  sign(c_star)*minmod_star*(F_R_star.rho-F_L_star.rho) + \
						  sign(c_R)*minmod_r*(F_R.rho-F_R_star.rho));

	 flux(f_idx).moment_u=0.5*(F_L.moment_u+F_R.moment_u)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_u-F_L.moment_u) + \
							       sign(c_star)*minmod_star*(F_R_star.moment_u-F_L_star.moment_u) + \
							       sign(c_R)*minmod_r*(F_R.moment_u-F_R_star.moment_u));

	 flux(f_idx).moment_v=0.5*(F_L.moment_v+F_R.moment_v)-0.5*(sign(c_L)*minmod_l*(F_L_star.moment_v-F_L.moment_v) + \
							       sign(c_star)*minmod_star*(F_R_star.moment_v-F_L_star.moment_v) + \
							       sign(c_R)*minmod_r*(F_R.moment_v-F_R_star.moment_v));
						     

	 flux(f_idx).energy = 0.5*(F_L.energy + F_R.energy) - 0.5* (sign(c_L)*minmod_l*(F_L_star.energy-F_L.energy) + \
						     sign(c_star)*minmod_star*(F_R_star.energy-F_L_star.energy) + \
						     sign(c_R)*minmod_r*(F_R.energy-F_R_star.energy));


       }else{
	 std::cout <<"Something went wrong in the limiter calcualtion" << "\n";       
       }

       std::cout << "The flux calculated is : " << "\n";
       flux(f_idx).print();

       /*
       flux(i).rho = 0.0;
       flux(i).moment_u = 0.0;
       flux(i).moment_v = 0.0;
       flux(i).energy = 0.0;
       */
       // 
      //---------------------//      
  f_idx++;
  }
  //Logic in case of y-sweep;
  /*
  if (sweep == std::string("x-sweep")){

  }else if(sweep == std::string("y-sweep")){
    for(int i = 0; i < ncells+1; i++){
      flux(i).moment_u = flux(i).moment_v;
      flux(i).moment_v = flux(i).moment_u;
    }
  }else{
    std::cout <<"Sweep not specified. Fail to compute 1D WAF" << "\n";
  }
  */

  for(int i = 0; i < ncells+1; i++){
    std::cout<< "About to return flux. And the value of the 1D ith flux is : " << "\n";
    flux(i).print();
  }
  return flux;
      
      
}

int sign(double c){
  
  if(c >= 0){
    return 1;
  }
  else{
    return -1;
  }

}
double minmod(double r, double c){
  if (r <= 0){
    return 1.0;
  }
  if( r >= 0 && r <= 1){
    return (1-(1-c)*r);
  }
  if (r > 1 ){
    return c;
  }
  else{
    std::cout  << "Error in minmod limiter calculation " << r << "\n";
    return 0;
  }
}

double superbee(double r, double c){

  if(r <= 0){
    return 1;
  }
  else{
    double top = (1-c)*2*r;
    return (1 - top/(1+r));
  }

}



//Modified to account for the exatra speed v in 2D. 
blitz::Array<Euler::U_state,1> HLLC_U_state(Euler::U_state U_state_L, Euler::U_state U_state_R){

  blitz::Array<Euler::U_state,1> U_interface; //vector of size 4 with U_L,U_L_star,U_R_star,U_star;
  
  U_interface.resize(4);
  U_interface(0)=U_state_L;
  U_interface(3)=U_state_R;

  Euler e;
  const double gamma = e.gamma;//This is a massive fudge!!!!! 
  
  //Variables used
  double P_star,P_L,P_R; //Presures
  double rho_star,rho_L,rho_R;//densities
  double u_L,u_R;//particle/gas speed in cell

  double a_L, a_R;// sound speed in cell
  double S_L, S_R,S_star;// wave speed in cell left wave, right wave, contact wave
    
  double P_pvrs;
  double rho_bar, a_bar;//average density and average sound speed of left and right cells

  Euler::W_state w_temp_left;
  Euler::W_state w_temp_right;

  double q_R,q_L;//Part of HLLC Calculation

 
  Euler::W_state W_L;
  Euler::W_state W_R;
  Euler::U_state U_state_L_star;
  Euler::U_state U_state_R_star;

  double star_coef_left; // The coeficient in eq. 10.73 from Toro(ed.2009); (S_k-u_k)/(S_k-u_star_k);
  double star_coef_right; // The coeficient in eq. 10.73 from Toro(ed.2009);
  
  
  W_L = e.PfromC(U_state_L);
  W_R = e.PfromC(U_state_R);


    //---------Pressure estimate-------------------------

    P_L = W_L.P;
    P_R = W_R.P;

    u_L = W_L.u;
    u_R = W_R.u;

    rho_L = W_L.rho;
    rho_R = W_R.rho;
   
    a_L = e.a(W_L);
    a_R = e.a(W_R);

    // std::cout <<"Inside the HLLC function" << "\n";
        
    rho_bar = 0.5*(rho_L + rho_R);
    a_bar = 0.5*(a_L + a_R);

    P_pvrs = 0.5*(P_L + P_R)-0.5*(u_R-u_L)*rho_bar*a_bar;

    P_star = std::max(0.0,P_pvrs);
      
    //---------------------------------------------------
   
    //----------Wave speed estimate----------------------
    
    //Calculate q_R

    if(P_star <= P_R){
      q_R = 1.0;
    }
    else{
      q_R = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_R)-1));
    }
        
    //Calculate q_L
    if(P_star <= P_L){
      q_L = 1.0;
    }
    else{
      q_L = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_L)-1));
    }
 
    //Calculate S_R and S_L
    S_L = u_L -a_L*q_L;
    S_R = u_R + a_R*q_R;
    
    double numerator = P_R - P_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R); 
    double denominator = rho_L*(S_L-u_L)-rho_R*(S_R-u_R);

    S_star = numerator/denominator; 

    //---------------------------------------------------

    //---------------HLLC fluxes -------------------------------------------------//

    //Calculate F_L     
    // Euler::U_state F_L;
    // F_L = m.ptr_euler->flux(U_state_L);
       

     //Calculate F_L_star
     // Euler::U_state F_L_star;
      star_coef_left = rho_L*((S_L-u_L)/(S_L-S_star));
      
      //---Calculate U_state_L_star
      U_state_L_star.rho = star_coef_left;
      U_state_L_star.moment_u = star_coef_left*S_star;
      U_state_L_star.moment_v = star_coef_left*W_L.v;
      U_state_L_star.energy = star_coef_left*(U_state_L.energy/U_state_L.rho + (S_star - u_L)*(S_star + P_L/(rho_L*(S_L-u_L))));
      U_interface(1)=U_state_L_star;

      // /* Consider overloading the + operator to write this in one line 
	// F_L_star.rho = F_L.rho + S_L*(U_state_L_star.rho -U_state_L.rho);
      //F_L_star.momentum = F_L.momentum + S_L*(U_state_L_star.momentum -U_state_L.momentum);
      //F_L_star.energy = F_L.energy + S_L*(U_state_L_star.energy -U_state_L.energy);
   
      //Calculate F_R
      // Euler::U_state F_R;
      // F_R = m.ptr_euler->flux(U_state_R);
     					

      //Calculate F_R_star
      // Euler::U_state F_R_star;
      star_coef_right = rho_R*((S_R-u_R)/(S_R-S_star));
      
      //Calculate U_state_R_star
      U_state_R_star.rho = star_coef_right;
      U_state_R_star.moment_u = star_coef_right*S_star;
      U_state_R_star.moment_v = star_coef_right*W_R.v;
      U_state_R_star.energy = star_coef_right*(U_state_R.energy/U_state_R.rho + (S_star - u_R)*(S_star + P_R/(rho_R*(S_R-u_R))));
      U_interface(2)=U_state_R_star;

//      /* Consider overloading the + operator to write this in one line 
      // F_R_star.rho = F_R.rho + S_R*(U_state_R_star.rho -U_state_R.rho);
      // F_R_star.momentum = F_R.momentum + S_R*(U_state_R_star.momentum -U_state_R.momentum);
      // F_R_star.energy = F_R.energy + S_R*(U_state_R_star.energy -U_state_R.energy);

      return U_interface;

}

