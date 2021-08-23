#include<iostream>
#include"Engine/array.hpp"
#include"cuMesh.hpp"




int main(){

   Mesh2d<Array<double,2>,double> mesh(1.0, 1.0, 110, 110);
   Mesh2d<Array<double,2>,double> mesh2(0.5, 0.5, 110, 110);
   Mesh2d<Array<double,2>,double> mesh3(0.5, 0.5, 110, 110);

   mesh.child = &mesh2;


   std::cout<< mesh.x.getDim(2) <<" esta bien papa?" << std::endl;

   mesh.solve(0.0, 0.0);


   std::cout<< mesh.getNx() <<" "<< mesh.getLx() << std::endl;
   std::cout<< mesh.getNy() <<" "<< mesh.getLy() << std::endl;
   std::cout<< mesh.getDeltaX() <<" "<< mesh.getDeltaY() << std::endl;

   std::cout<< mesh.child->getNx() <<" "<< mesh.child->getLx() 
            << std::endl;
   std::cout<< mesh.child->getNy() <<" "<< mesh.child->getLy() 
            << std::endl;
   std::cout<< mesh.child->getDeltaX() <<" "<< mesh.child->getDeltaY() 
            << std::endl;
   return 0;
}
