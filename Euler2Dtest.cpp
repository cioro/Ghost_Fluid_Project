#include "Euler2D.hpp"
#include <iostream>



int main(){

  Euler e;
  Euler::W_state w2D(1.0,0.0,0.0,1.0);
  Euler::U_state u2D(1.0,0.0,0.0,2.5); 

  w2D.print();
  u2D.print();
  
  double speed;
  speed = e.a(w2D);
  std::cout << "The speed of sound is " << speed << "\n";

  double int_energy;
  int_energy = e.int_energy(w2D);
  std::cout << "The internal energy is " << int_energy << "\n";

  Euler::U_state F;
  Euler::U_state G;

  F = e.flux_x(u2D);
  G = e.flux_y(u2D);

  std::cout << "These are the fluxes " << "\n";
  F.print();
  G.print();

  Euler::U_state U_new;
  Euler::W_state W_new;

  U_new = e.CfromP(w2D);
  std::cout << "Check transform. This sould be u(1,0,0,2.5) " << "\n";
  U_new.print();

  W_new = e.PfromC(u2D);
  std::cout << "Check transform. This should be w(1,0,0,1)" << "\n";
  W_new.print();

  return 0;
}
