#include<string.h>
#include<iostream>
#include<limits>
#include<math.h>

#include"stlIO.hpp"

template<typename Tprec>
StlIO<Tprec>::StlIO(const char _name[])  
{
   
  file_name = new char[strlen(_name)];
  strcpy(file_name, _name);
  std::cout << "Open file: " << file_name << std::endl;
  fileSTL.open(file_name, ios::binary);
  if(fileSTL.good()){std::cout << "    SUCCED" << std::endl;}
  fileSTL.read(header, 80);
  fileSTL.read((char*)(&num_facets), 4);
  std::cout << "Number of triangles: " << num_facets << std::endl;
  gmtry = new triangleSTL<Tprec>[num_facets]; 

}

template<typename Tprec>
void StlIO<Tprec>::read(){
   unsigned short attr;

   triangleSTL<float> trash;

   for(unsigned int i = 0; i < num_facets; i++){

       fileSTL.read((char*)(trash.n), 12);

       gmtry[i].n[0] = static_cast<Tprec>(trash.n[0]);
       gmtry[i].n[1] = static_cast<Tprec>(trash.n[1]);
       gmtry[i].n[2] = static_cast<Tprec>(trash.n[2]);

       fileSTL.read((char*)(trash.vrtx1), 12);

       gmtry[i].vrtx1[0] = static_cast<Tprec>(trash.vrtx1[0]);
       gmtry[i].vrtx1[1] = static_cast<Tprec>(trash.vrtx1[1]);
       gmtry[i].vrtx1[2] = static_cast<Tprec>(trash.vrtx1[2]);



       fileSTL.read((char*)(trash.vrtx2), 12);

       gmtry[i].vrtx2[0] = static_cast<Tprec>(trash.vrtx2[0]);
       gmtry[i].vrtx2[1] = static_cast<Tprec>(trash.vrtx2[1]);
       gmtry[i].vrtx2[2] = static_cast<Tprec>(trash.vrtx2[2]);



       fileSTL.read((char*)(trash.vrtx3), 12);

       gmtry[i].vrtx3[0] = static_cast<Tprec>(trash.vrtx3[0]);
       gmtry[i].vrtx3[1] = static_cast<Tprec>(trash.vrtx3[1]);
       gmtry[i].vrtx3[2] = static_cast<Tprec>(trash.vrtx3[2]);


       fileSTL.read((char*)(&attr), 2);
/*
       printf("%f o %f o %f\n", gmtry[i].n[0], gmtry[i].n[1], gmtry[i].n[2]);
       printf("%f o %f o %f\n", gmtry[i].vrtx1[0], gmtry[i].vrtx1[1], gmtry[i].vrtx1[2]);
       printf("%f o %f o %f\n", gmtry[i].vrtx2[0], gmtry[i].vrtx2[1], gmtry[i].vrtx2[2]);
       printf("%f o %f o %f\n\n", gmtry[i].vrtx3[0], gmtry[i].vrtx3[1], gmtry[i].vrtx3[2]);

*/
   }



/*
   for(unsigned int i = 0; i < num_facets; i++){

       fileSTL.read((char*)(normal1), 12);
       gmtry[i].n[0] = static_cast<Tprec>(normal[0]);
       gmtry[i].n[1] = static_cast<Tprec>(normal[1]);
       gmtry[i].n[2] = static_cast<Tprec>(normal[2]);
       printf("%f o %f o %f\n", normal1[0], normal1[1], normal1[2]);

       fileSTL.read((char*)(normal), 12);
       gmtry[i].vrtx1[0] = static_cast<Tprec>(normal[0]);
       gmtry[i].vrtx1[1] = static_cast<Tprec>(normal[1]);
       gmtry[i].vrtx1[2] = static_cast<Tprec>(normal[2]);
       printf("%f o %f o %f\n", normal[0], normal[1], normal[2]);

       fileSTL.read((char*)(normal), 12);

       gmtry[i].vrtx2[0] = static_cast<Tprec>(normal[0]);
       gmtry[i].vrtx2[1] = static_cast<Tprec>(normal[1]);
       gmtry[i].vrtx2[2] = static_cast<Tprec>(normal[2]);

       printf("%f o %f o %f\n", normal[0], normal[1], normal[2]);

       fileSTL.read((char*)(normal), 12);

       gmtry[i].vrtx3[0] = static_cast<Tprec>(normal[0]);
       gmtry[i].vrtx3[1] = static_cast<Tprec>(normal[1]);
       gmtry[i].vrtx3[2] = static_cast<Tprec>(normal[2]);

       printf("%f o %f o %f\n", normal[0], normal[1], normal[2]);

       fileSTL.read((char*)(&attr), 2);

       printf("%d\n", attr);

   }*/


  fileSTL.close();

}

template<typename Tprec>
void StlIO<Tprec>::getVrtx(triangleSTL<Tprec> *_vrtx){
     //_vrtx = new triangleSTL[num_facets]; 
     //_vrtx = new triangleSTL[1]; 
     memcpy(_vrtx, gmtry, num_facets*(sizeof(triangleSTL<Tprec>)));
}

template<typename Tprec>
void StlIO<Tprec>::getNameSTL() const {
  std::cout << "    "<<std::endl << file_name << std::endl;
}

template<typename Tprec>
unsigned int StlIO<Tprec>::getNumTrngl() const {
  return num_facets;
}


template<typename Tprec>
Tprec StlIO<Tprec>::getMaxAxisX() const {
   Tprec maxVal = -1.0*numeric_limits<Tprec>::max();

   for(int unsigned id = 0; id < num_facets; id++){
      maxVal = max(gmtry[id].vrtx1[0], maxVal);
      maxVal = max(gmtry[id].vrtx2[0], maxVal);
      maxVal = max(gmtry[id].vrtx3[0], maxVal);
   }

   return maxVal;
}

template<typename Tprec>
void StlIO<Tprec>::rotateAxisZ(float alpha){

    for(int id = 0; id < getNumTrngl(); id++){
       Tprec x1 = gmtry[id].vrtx1[0],
             y1 = gmtry[id].vrtx1[1];
              
        gmtry[id].vrtx1[0] = cos(alpha)*x1 - sin(alpha)*y1;
        gmtry[id].vrtx1[1] = sin(alpha)*x1 + cos(alpha)*y1;
        
        x1 = gmtry[id].vrtx2[0],
        y1 = gmtry[id].vrtx2[1];
        
        gmtry[id].vrtx2[0] = cos(alpha)*x1 - sin(alpha)*y1;
        gmtry[id].vrtx2[1] = sin(alpha)*x1 + cos(alpha)*y1;
        
        x1 = gmtry[id].vrtx3[0],
        y1 = gmtry[id].vrtx3[1];
        
        gmtry[id].vrtx3[0] = cos(alpha)*x1 - sin(alpha)*y1;
        gmtry[id].vrtx3[1] = sin(alpha)*x1 + cos(alpha)*y1;
        
        
        x1 = gmtry[id].n[0],
        y1 = gmtry[id].n[1];
        
        gmtry[id].n[0] = cos(alpha)*x1 - sin(alpha)*y1;
        gmtry[id].n[1] = sin(alpha)*x1 + cos(alpha)*y1;
        
              
    }



}





