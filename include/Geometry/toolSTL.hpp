#ifndef __TOOLSTL_HPP
#define __TOOLSTL_HPP

#include<cmath>

#include "stlIO.hpp"
#include "meshCreator.hpp"
#include "vector3d.hpp"

struct At{
    int longID;
    int shortID;
    int newID;
    int i;
};


template<typename Tprec>
struct Point
{
    Tprec x;
    Tprec y;
    Tprec z;
    int longID;
    int shortID;
    int newID;
    int i;
};

template<typename Tprec>
struct PointList
{
    Tprec      x;
    Tprec      y;
    Tprec      z;
    int        clusterDim;
    int        filled;
    int        posId;
    int        id;
    PointList* point;
};

template<typename Tprec>
int getNumOfPoints(PointList<Tprec>* pointList){
    int counter;
    PointList<Tprec>* temp = pointList;
    while(temp != 0){
        counter++;
        temp = temp->point;
    }
    return counter;

}

template<typename Tprec>
void getPosId(PointList<Tprec>* pointList){
    int counter;
    PointList<Tprec>* temp   = pointList;
    
    do{
        temp->point->posId =  temp->posId + temp->clusterDim;
        temp = temp->point;
    }while(temp->point != 0);

}

template<typename Tprec>
void printList(PointList<Tprec>* pointList){
         
    if(pointList->point != 0){
        printList(pointList->point);
    }
}

template<typename Tprec>
void deletePointList(PointList<Tprec>* pointList){
     
    PointList<Tprec>* temporal; 
    if(pointList->point != 0){
        temporal = pointList->point;
    }

    delete [] pointList;

    if(temporal != 0){
      deletePointList(temporal);
    }
}

template<typename Tprec>
void searchForPoints(PointList<Tprec>* pointList, Point<Tprec>* point){
    if(pointList->x == point->x && pointList->y == point->y && pointList->z == point->z){
        pointList->clusterDim ++;
    }else if(pointList->point != 0){
        searchForPoints(pointList->point, point);
    }else{
        pointList->point    = new PointList<Tprec>;
        pointList->point->x = point->x; 
        pointList->point->y = point->y; 
        pointList->point->z = point->z; 
        pointList->point->clusterDim = 0;
        pointList->point->id     = pointList->id + 1;
        pointList->point->point  = 0; 
        pointList->point->filled = 0;  
        searchForPoints(pointList->point, point);    
    }

}

template<typename Tprec>
void makeAnArrayOfPoints(PointList<Tprec>** arrayOfPoints, PointList<Tprec>* pointList){

    arrayOfPoints[pointList->id] = pointList;
    if(pointList->point != 0){
       makeAnArrayOfPoints(arrayOfPoints, pointList->point);
    }
}

template<typename Tprec>
void sortClustersForAp(PointList<Tprec>* pointList, Point<Tprec>* aL, Point<Tprec>* aP){

     if(pointList->x == aL->x && pointList->y == aL->y && pointList->z == aL->z){
         int id = pointList->posId + pointList->filled; 
         
         aP[id].x      = aL->x;
         aP[id].y      = aL->y;
         aP[id].z      = aL->z;
         aP[id].longID = aL->longID;
         aP[id].newID  = pointList->id;
         aP[id].i      = id;
         
         if(pointList->filled == 0){
            aP[id].shortID = aL->longID;          
         } else{
            aP[id].shortID = aP[pointList->posId].shortID;   
         }
         pointList->filled++;
     } else if (pointList->point != 0){
          sortClustersForAp(pointList->point, aL, aP);
     }

}

///////////////////////////////////////////////////////////////7

template<typename Tprec>
void getNeighborList(int             *neighbor, 
                     StlIO<Tprec>    *stl
                  ){

 
   //point list data 
   int               numOfVertex = 3*stl->getNumTrngl();
   PointList<Tprec>* pointList;
   Point<Tprec>*     point;
   Point<Tprec>*     aP;
   int               numOfDiferentPoints = 0;
   
   pointList = new PointList<Tprec>;
   point     = new Point<Tprec>[numOfVertex];
   aP        = new Point<Tprec>[numOfVertex];
   
   for(int triangl = 0; triangl < stl->getNumTrngl(); triangl++){
   
       int id = 3*triangl;
       
       point[id].x      = stl->gmtry[triangl].vrtx1[0];
       point[id].y      = stl->gmtry[triangl].vrtx1[1];
       point[id].z      = stl->gmtry[triangl].vrtx1[2];
       point[id].longID = id;
       
       point[id + 1].x      = stl->gmtry[triangl].vrtx2[0];
       point[id + 1].y      = stl->gmtry[triangl].vrtx2[1];
       point[id + 1].z      = stl->gmtry[triangl].vrtx2[2];
       point[id + 1].longID = id + 1;
       
       point[id + 2].x      = stl->gmtry[triangl].vrtx3[0];
       point[id + 2].y      = stl->gmtry[triangl].vrtx3[1];
       point[id + 2].z      = stl->gmtry[triangl].vrtx3[2];
       point[id + 2].longID = id + 2;
       
   }

   //initial point
   pointList->x          = point[0].x;
   pointList->y          = point[0].y;
   pointList->z          = point[0].z;
   pointList->id         = 0;
   pointList->point      = 0;
   pointList->posId      = 0;
   pointList->clusterDim = 0;
   pointList->filled     = 0;


   for(int id = 0; id < numOfVertex; id++){
      searchForPoints(pointList, &point[id]);
   }
   
   numOfDiferentPoints = getNumOfPoints(pointList);
     
   PointList<Tprec>** arrayOfPoints = new PointList<Tprec>*[numOfDiferentPoints];
   
   makeAnArrayOfPoints(arrayOfPoints, pointList);
   
   getPosId(pointList);
   
   //printList(pointList);
   
   for(int id = 0; id < numOfVertex; id++){
      sortClustersForAp(pointList, &point[id], aP);
   }
   
   for(int id = 0; id < numOfVertex; id++){ 
   
      point[aP[id].longID].i       = aP[id].i;
      point[aP[id].longID].shortID = aP[id].shortID;
      point[aP[id].longID].newID   = aP[id].newID;
   }
   

   
   for(int triangl = 0; triangl < stl->getNumTrngl(); triangl++){
      for(int i = 0; i < 3; i++){
           int id             = triangl*3 + i;         
           int pivotPoint     = point[point[id].shortID].i;
           int pointSelected  = pivotPoint;     
           
           while(aP[pointSelected].newID == aP[pivotPoint].newID){
               int triangleSelected = aP[pointSelected].longID / 3;
               
               int ident[2];
               int match =  0;
               ident[0]  = -1;
               ident[1]  = -1;
                
               //we got two triangles shearing the same vertex
               //we need to shearch for other vertex to fing a sibiling
               
                if(triangl != triangleSelected){
               
               for(int j = 0; j < 3; j++){
                  int idTOri = triangl*3 + j;
                  for(int k = 0; k < 3; k++){
                      int idTSel = triangleSelected*3 + k;  

                      if(point[idTOri].shortID == point[idTSel].shortID){
                          ident[match] = j;
                          match++;
                      }
                  
                  }
               }
               
               if(match == 2){
                  if((ident[0] == 0 && ident[1] == 1) || (ident[1] == 0 && ident[0] == 1)){
                  
                     neighbor[triangl*3 + 0] = triangleSelected;  
                      
                              
                  }else if ((ident[0] == 1 && ident[1] == 2) || (ident[1] == 1 && ident[0] == 2)){
                  
                     neighbor[triangl*3 + 1] = triangleSelected;
                      
                    
                  }else if ((ident[0] == 2 && ident[1] == 0) || (ident[1] == 2 && ident[0] == 0)){                                     

                     neighbor[triangl*3 + 2] = triangleSelected;
                    
                  }
               }
               }//diferent triangle
               pointSelected++;
           }
           
      }
   
   }


   deletePointList(pointList);
   delete [] point;

}

////////////////////////////////////////////////////////////////////////7

template<typename Tprec>
void stlToLevelSet(Tprec           *phiLS,
                   Tprec           *mask, 
                   int             *neighbor,
                   Vector3d<Tprec> *mesh,
                   int             *index1,
                   int             *index2,
                   int             *index3,
                   int              Nx,
                   int              Ny,
                   int              Nz,
                   int*             np,
                   StlIO<Tprec>    *stl
                  ){

   Tprec distance = 0;
   Tprec t  = 0.0;
   Tprec t1 = 0.0;
   Tprec t2 = 0.0;
   Tprec t3 = 0.0;

   Tprec d  = 0.0;
   Tprec d0 = 0.0;
   Tprec d1 = 0.0;
   Tprec d2 = 0.0;

   Tprec dt0 = 0.0;
   Tprec dt1 = 0.0;
   Tprec dt2 = 0.0;

   Tprec st0 = 0.0;
   Tprec st1 = 0.0;
   Tprec st2 = 0.0;

   Tprec rad0 = 0.0;
   Tprec rad1 = 0.0;
   Tprec rad2 = 0.0;


   //using 

   vecMath::Vector3d<Tprec> neighbPlane1;
   vecMath::Vector3d<Tprec> neighbPlane2;
   vecMath::Vector3d<Tprec> neighbPlane3;
   
   Vector3d<Tprec> plane;
   Vector3d<Tprec> vector;

   plane.x = new Tprec[3];
   plane.y = new Tprec[3];
   plane.z = new Tprec[3];

   vector.x = new Tprec[3];
   vector.y = new Tprec[3];
   vector.z = new Tprec[3];

   Tprec sgama =sqrt(2.0)*(mesh->x[1] - mesh->x[0]);
   
   //point list data 
   int numOfVertex = 3*stl->getNumTrngl();   
   int idnp = 0;
   
   //computing planes and inner region 

   for(int k = 0; k <  Nz; k++ ){
      for(int j = 0; j <  Ny; j++ ){
         for(int i = 0; i <  Nx; i++ ){

             int id = Nx*Ny*k + Nx*j + i;
             mask[id]  = 2.0;
             phiLS[id] = 1.0;
             
         }
      }
   }



   for(int triangl = 0; triangl < stl->getNumTrngl() ; triangl++){


      d  = - stl->gmtry[triangl].n[0]*stl->gmtry[triangl].vrtx2[0] 
           - stl->gmtry[triangl].n[1]*stl->gmtry[triangl].vrtx2[1]
           - stl->gmtry[triangl].n[2]*stl->gmtry[triangl].vrtx2[2]; 

      Tprec d1, d2, d3;

      vecMath::Vector3d<Tprec> vector1(0.0,0.0,0.0);
      vecMath::Vector3d<Tprec> vector2(0.0,0.0,0.0);
      vecMath::Vector3d<Tprec> vector3(0.0,0.0,0.0);
      
      vecMath::Vector3d<Tprec> normal(stl->gmtry[triangl].n[0],stl->gmtry[triangl].n[1],stl->gmtry[triangl].n[2]);
      
      vecMath::Vector3d<Tprec> normalN1(stl->gmtry[neighbor[triangl*3 + 0]].n[0],
                                        stl->gmtry[neighbor[triangl*3 + 0]].n[1],
                                        stl->gmtry[neighbor[triangl*3 + 0]].n[2]);
                                        
      vecMath::Vector3d<Tprec> normalN2(stl->gmtry[neighbor[triangl*3 + 1]].n[0],
                                        stl->gmtry[neighbor[triangl*3 + 1]].n[1],
                                        stl->gmtry[neighbor[triangl*3 + 1]].n[2]);
                                                                          
      vecMath::Vector3d<Tprec> normalN3(stl->gmtry[neighbor[triangl*3 + 2]].n[0],
                                        stl->gmtry[neighbor[triangl*3 + 2]].n[1],
                                        stl->gmtry[neighbor[triangl*3 + 2]].n[2]);                                   
                                        
      vecMath::Vector3d<Tprec> normal1(0.0,0.0,0.0);
      vecMath::Vector3d<Tprec> normal2(0.0,0.0,0.0);
      vecMath::Vector3d<Tprec> normal3(0.0,0.0,0.0);
      
      vector1.x = stl->gmtry[triangl].vrtx2[0] - stl->gmtry[triangl].vrtx1[0]; 
      vector1.y = stl->gmtry[triangl].vrtx2[1] - stl->gmtry[triangl].vrtx1[1]; 
      vector1.z = stl->gmtry[triangl].vrtx2[2] - stl->gmtry[triangl].vrtx1[2]; 
      
      vector2.x = stl->gmtry[triangl].vrtx3[0] - stl->gmtry[triangl].vrtx2[0]; 
      vector2.y = stl->gmtry[triangl].vrtx3[1] - stl->gmtry[triangl].vrtx2[1]; 
      vector2.z = stl->gmtry[triangl].vrtx3[2] - stl->gmtry[triangl].vrtx2[2]; 
      
      vector3.x = stl->gmtry[triangl].vrtx1[0] - stl->gmtry[triangl].vrtx3[0]; 
      vector3.y = stl->gmtry[triangl].vrtx1[1] - stl->gmtry[triangl].vrtx3[1]; 
      vector3.z = stl->gmtry[triangl].vrtx1[2] - stl->gmtry[triangl].vrtx3[2]; 
      
      normal1 = normal + normalN1;
      normal2 = normal + normalN2;
      normal3 = normal + normalN3;
      
      normal1 = vector1.vectorProduct(normal1);
      normal2 = vector2.vectorProduct(normal2);
      normal3 = vector3.vectorProduct(normal3);
      
      normal1.normalize();
      normal2.normalize();
      normal3.normalize();
      
      d1 = -normal1.x*stl->gmtry[triangl].vrtx2[0] 
           -normal1.y*stl->gmtry[triangl].vrtx2[1]
           -normal1.z*stl->gmtry[triangl].vrtx2[2];     
           
      d2 = -normal2.x*stl->gmtry[triangl].vrtx3[0] 
           -normal2.y*stl->gmtry[triangl].vrtx3[1]
           -normal2.z*stl->gmtry[triangl].vrtx3[2];
      
      d3 = -normal3.x*stl->gmtry[triangl].vrtx1[0] 
           -normal3.y*stl->gmtry[triangl].vrtx1[1]
           -normal3.z*stl->gmtry[triangl].vrtx1[2]; 

   //   for(idnp = 0; idnp <  *np; idnp++ ){
     printf("%d   %f      \n", triangl, d);

   for(int k = 0; k <  Nz; k++ ){
      for(int j = 0; j <  Ny; j++ ){
         for(int i = 0; i <  Nx; i++ ){
             int id = Nx*Ny*k + Nx*j + i;
      
              //  int i = index1[idnp];
              //  int j = index2[idnp];
              //  int k = index3[idnp];

                
                int idip1 = Nx*Ny*k + Nx*j + i + 1,
                    idim1 = Nx*Ny*k + Nx*j + i - 1,
                    idjp1 = Nx*Ny*k + Nx*(j+1) + i,
                    idjm1 = Nx*Ny*k + Nx*(j-1) + i,
                    idkp1 = Nx*Ny*(k+1) + Nx*j + i,
                    idkm1 = Nx*Ny*(k-1) + Nx*j + i;    
                       
                int idip2 = Nx*Ny*k + Nx*(j+1) + i + 1,
                    idim2 = Nx*Ny*k + Nx*(j-1) + i - 1,
                    idjp2 = Nx*Ny*k + Nx*(j+1) + i - 1,
                    idjm2 = Nx*Ny*k + Nx*(j-1) + i + 1;
                        
                 
                //planes distance

                t = -  normal.x*mesh->x[id] 
                    -  normal.y*mesh->y[id]
                    -  normal.z*mesh->z[id]; 

               
                t1 = - normal1.x*mesh->x[id] 
                     - normal1.y*mesh->y[id]
                     - normal1.z*mesh->z[id];
                
                t2 = -  normal2.x*mesh->x[id] 
                     -  normal2.y*mesh->y[id]
                     -  normal2.z*mesh->z[id]; 
                     
                t3 = -  normal3.x*mesh->x[id] 
                     -  normal3.y*mesh->y[id]
                     -  normal3.z*mesh->z[id];

                
                if((d - t) <= 0.0   && (d - t) >= -sgama
                                    &&  (d1 - t1) <= 0.0
                                    &&  (d2 - t2) <= 0.0 
                                    &&  (d3 - t3) <= 0.0                                                                    
                  ){                      


                     
                     idkp1     = 1.0;

                     //phiLS[id] = (d - t) / sgama;
                     //
                     

                  }else if((d - t) >= 0.0   && (d - t) <= 2.0*sgama
                                    &&  (d1 - t1) <= 0.0
                                    &&  (d2 - t2) <= 0.0 
                                    &&  (d3 - t3) <= 0.0                                                                   
                  ){ }

           
      }//end loop k 
      }
   }
   }//end loop triangle id
   
   idnp = 0;
   
   for(int k = 2; k <  Nz-2; k++ ){
      for(int j = 2; j <  Ny-2; j++ ){
         for(int i = 2; i <  Nx-2; i++ ){
             int id = Nx*Ny*k + Nx*j + i;
             if(mask[id] == 1.0){
                  index1[idnp] = i;
                  index2[idnp] = j;
                  index3[idnp] = k;
                  idnp++;
             }
             
         }
      }
   }
   
   *np = idnp;
   
   delete [] plane.x;
   delete [] plane.y;
   delete [] plane.z;

   delete [] vector.x;
   delete [] vector.y;
   delete [] vector.z;


}


#endif
