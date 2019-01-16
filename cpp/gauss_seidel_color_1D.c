// Author: Wes Kendall
// Copyright 2011 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// An intro MPI hello world program that uses MPI_Init, MPI_Comm_size,
// MPI_Comm_rank, MPI_Finalize, and MPI_Get_processor_name.
//
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define SIZE 4
#define N2 36
int main(int argc, char** argv) {






// ofstream myfile ("example.txt");
//   if (myfile.is_open())
//   {
//     myfile << "This is a line.\n";
//     myfile << "This is another line.\n";
//     myfile.close();
//   }
//   
//   
  
  
  
  
  // Initialize the MPI environment. The two arguments to MPI Init are not
  // currently used by MPI implementations, but are there in case future
  // implementations might need the arguments.
  MPI_Init(NULL, NULL);

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  printf("Hello world from processor %s, rank %d out of %d processors\n",
         processor_name, world_rank, world_size);

  int number;
  if (world_rank == 0) {
    number = -1;
    for(int jj=1;jj<world_size;jj++)    
       MPI_Send(&number, 1, MPI_INT, jj, 0, MPI_COMM_WORLD);
  } else if (world_rank >0 ) {
    MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    printf("Process %d, received number %d from process 0\n",world_rank,number);
}



// MPI_Barrier(MPI_COMM_WORLD);
// printf("\n");
// MPI_Barrier(MPI_COMM_WORLD);
// const int MAX_NUMBERS = 20;
// int numbers[MAX_NUMBERS];
// int number_amount;
// if (world_rank == 0) {
//     // Pick a random amount of integers to send to process one
//     srand(time(NULL));
//     number_amount = (rand() / (float)RAND_MAX) * MAX_NUMBERS;
// 
//     // Send the amount of integers to process one
//     MPI_Send(numbers, number_amount, MPI_INT, 1, 0, MPI_COMM_WORLD);
//  //    for(int ii=0;ii<MAX_NUMBERS;ii++)
// //         printf("COSA MANDO: %d , %d  to 1\n", ii,numbers[ii]);
//     printf("0 sent %d numbers to 1\n", number_amount);
// } else if (world_rank == 1) {
//     MPI_Status status;
//     // Receive at most MAX_NUMBERS from process zero
//     MPI_Recv(numbers, MAX_NUMBERS, MPI_INT, 0, 0, MPI_COMM_WORLD,
//              &status);
// 
//     // After receiving the message, check the status to determine
//     // how many numbers were actually received
//     MPI_Get_count(&status, MPI_INT, &number_amount);
// 
//     // Print off the amount of numbers, and also print additional
//     // information in the status object
//     printf("1 received %d numbers from 0. Message source = %d, "
//            "tag = %d\n",
//            number_amount, status.MPI_SOURCE, status.MPI_TAG);
//      //  MPI_Send(numbers, number_amount, MPI_INT, 1, 0, MPI_COMM_WORLD);
// }
// MPI_Barrier(MPI_COMM_WORLD);
// if (world_rank == 1) {
//     for(int ii=0;ii<MAX_NUMBERS;ii++)
//       printf("COSA RICEVO: %d , %d  to 1\n", ii,numbers[ii]);
//     }
//     int a=1;
//     if(world_rank==0)
//     a=99;
// MPI_Bcast(&a, 1, MPI_INT, 0, MPI_COMM_WORLD);    
//  printf("BBCAST RICEVO: %d\n", a);   
//  MPI_Barrier(MPI_COMM_WORLD);
// printf("\n");   
// printf("\n");   
// printf("\n");   
// 
// MPI_Barrier(MPI_COMM_WORLD);








int send_count=3;
int b[world_size*send_count];
int c[send_count];
int d[world_size*send_count];
int g[world_size*send_count];
 
 
 
if(world_rank==0)
 {for(int ii=0;ii<world_size*send_count;ii++)
     b[ii]=ii*ii;
     }
     
 MPI_Scatter(
    &b,
    send_count,
    MPI_INT,
    &c,
    send_count,
    MPI_INT,
    0,
    MPI_COMM_WORLD );
  // Finalize the MPI environment. No more MPI calls can be made after this
  for(int ii=0;ii<send_count;ii++)
     {c[ii]=c[ii]+3;
     printf(" PRIMA processor %d: %d\n", world_rank,c[ii]);}
     
     
     
//      MPI_Gather(
//     &c,
//     send_count,
//     MPI_INT,
//     &d,
//     send_count,
//     MPI_INT,
//     0,
//     MPI_COMM_WORLD );
    
//  if(world_rank==0)   
//       {for(int ii=0;ii<world_size*send_count;ii++)
//      printf("DOPO processor %d: %d\n", world_rank,d[ii]);
//      }
     
  MPI_Allgather(
    &c,
    send_count,
    MPI_INT,
    g,
    send_count,
    MPI_INT,
    MPI_COMM_WORLD );
  
  printf("vettore c, processore %d %d\n",world_rank,c[0]);
  for(int ii=0;ii<world_size*send_count;ii++)
     printf("MPI_Allgather DOPO processor %d: %d\n", world_rank,g[ii]);
     
  
  
  
  
  
  
  
  
MPI_Barrier(MPI_COMM_WORLD);
printf("\n");   
printf("\n");   
printf("\n");   

MPI_Barrier(MPI_COMM_WORLD);

  
  
  
int N=world_size*send_count;
    int rem = (N)%world_size; // elements remaining after division among processes
    int sum = 0;                // Sum of counts. Used to calculate displacements

 int *sendcounts = malloc(sizeof(int)*world_size);
    int *displs = malloc(sizeof(int)*world_size);

    // calculate send counts and displacements
    for (int i = 0; i < world_size; i++) {
        sendcounts[i] = N/world_size;
        if (rem > 0) {
            sendcounts[i]++;
            rem--;
        }

        displs[i] = sum;
        sum += sendcounts[i];
    }
    
    
    
MPI_Scatterv(&g, sendcounts,displs,
                 MPI_INT, &c,displs[world_rank],
                 MPI_INT,0, MPI_COMM_WORLD );

for(int jj=0;jj< world_size;jj++)  
{           
if(jj==world_rank){
for(int ii=0;ii< sendcounts[world_rank];ii++)
printf("DOPOOOOOO subdivide, id %d number and %d \n",world_rank, c[ii]); 
}  
}

int summa;
int* tmp=malloc(sizeof(int)*1);

MPI_Allreduce(&world_rank,
    &summa,
    SIZE,
    MPI_INT,
    MPI_SUM,
    MPI_COMM_WORLD);
    

    printf("somma: %d \n",summa);

 int summa2[SIZE];
 int tmp2[SIZE];  
 for(int ii=0;ii< SIZE;ii++)
 {tmp2[ii]=1+ii*ii;
printf("tmp2, id %d number and %d \n",world_rank, tmp2[ii]); }

MPI_Allreduce(&tmp2,
    &summa2,
    SIZE,
    MPI_INT,
    MPI_SUM,
    MPI_COMM_WORLD);   

for(int ii=0;ii< SIZE;ii++)
printf("summa2, id %d number and %d \n",world_rank, summa2[ii]); 

    
    

// Split the communicator based on the color and use the
// original rank for ordering
// int color = world_rank / 4; 
// MPI_Comm row_comm;
// MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &row_comm);
// 
// int row_rank, row_size;
// MPI_Comm_rank(row_comm, &row_rank);
// MPI_Comm_size(row_comm, &row_size);
// 
// printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n",
// 	world_rank, world_size, row_rank, row_size);
// 
// MPI_Comm_free(&row_comm);


N=sqrt(N2);
double Amat[N2];

double Adiag[N];
double xvec[N];
double bvec[N];


//Amat=malloc(sizeof(double)*(N*N));
Amat[0]=1;
Adiag[0]=1;
bvec[0]=0;
xvec[0]=0;

for(int jj=1;jj<N;jj++)
    {Amat[jj]=0;
     Adiag[jj]=2;
     bvec[jj]=1;
     xvec[jj]=0;
     }
for(int ii=1;ii<N-1;ii++)
    {
     for(int jj=0;jj<ii;jj++)
         Amat[ii*N+jj]=0;
     Amat[ii*N+ii-1]=1;
     Amat[ii*N+ii]=2;
     Amat[ii*N+ii+1]=1;
     for(int jj=ii+2;jj<N;jj++)
         Amat[ii*N+jj]=0;
    }
for(int jj=0;jj<N;jj++)
    Amat[N*(N-1)+jj]=0;
Amat[N*N-1]=1;    
Adiag[N-1]=1;
bvec[N-1]=0;    
xvec[N-1]=0; 

for(int ii=0;ii<N;ii++)
    {for(int jj=0;jj<N;jj++)
        printf("%f ",Amat[ii*N+jj]);
     printf("\n");
     }        


rem = (N2)%world_size; 
sum = 0;
free(sendcounts);
free(displs);
sendcounts = malloc(sizeof(int)*world_size);
displs = malloc(sizeof(int)*world_size);

int rem_vec = (N)%world_size; 
int sum_vec = 0;
int *sendcounts_vec = malloc(sizeof(int)*world_size);
int *displs_vec = malloc(sizeof(int)*world_size);
// calculate send counts and displacements
    for (int i = 0; i < world_size; i++) {
        
        sendcounts_vec[i] = (N)/world_size;
        
//         if (rem > 0) {
//             sendcounts[i]++;
//             rem--;
//         }
        if (rem_vec > 0) {
            sendcounts_vec[i]++;
            rem_vec--;
        }
        
        sendcounts[i] = (sendcounts_vec[i])*N ;
        sum += sendcounts[i];
        
        displs_vec[i] = sum_vec;
        displs[i] = displs_vec[i]*N;
        sum_vec += sendcounts_vec[i];
        printf("sendcounts_vec[%d]==%d, sendcounts[%d]==%d, displs_vec[%d]=%d,  displs[%d]=%d\n",i,sendcounts_vec[i],i, sendcounts[i],i, displs_vec[i], i, displs[i]) ;  
    } 
    
       
MPI_Barrier(MPI_COMM_WORLD);
printf("\n");   
printf("\n");   
printf("\n");   
printf("\n");   
printf("\n");   
printf("\n"); 
printf("\n");   
printf("\n");   
printf("\n"); 
MPI_Barrier(MPI_COMM_WORLD);
   
         
if(0==world_rank){
for(int ii=0;ii< world_size;ii++)
printf("N=%d , id=%d, sendcounts=%d, displs=%d, sendcounts_vec=%d, displs_vec=%d \n",N,world_rank, sendcounts[ii],displs[ii],sendcounts_vec[ii],displs_vec[ii]); 
}  


double Aloc[sendcounts[world_rank]];
double Adiagloc[sendcounts_vec[world_rank]];
double bloc[sendcounts_vec[world_rank]];
double xloc[sendcounts_vec[world_rank]];
// the coloring in 1D is given by only 2 colors: odd and even numbers
int number_of_colors=2;
int color_vec[N];
int colorloc[sendcounts_vec[world_rank]];
double resloc[sendcounts_vec[world_rank]];
double correctionloc[sendcounts_vec[world_rank]];

if(world_rank==0)
{
for(int ii=0;ii<N;ii++)
   {color_vec[ii]=ii%2;
   printf("color %d %d \n",ii,color_vec[ii]);}

}


// double* Aloc;
// Aloc=malloc(sizeof(double)*sendcounts[world_rank]); 
 MPI_Barrier(MPI_COMM_WORLD); 
// for(int ii=0;ii< N2;ii++)
// printf("Amat %d and %d %f \n",world_rank,ii, Amat[ii]); 
    MPI_Barrier(MPI_COMM_WORLD); 
    
MPI_Scatterv(&Amat, sendcounts,displs,
             MPI_DOUBLE, &Aloc,sendcounts[world_rank],
             MPI_DOUBLE,0, MPI_COMM_WORLD );
 
MPI_Scatterv(&bvec, sendcounts_vec,displs_vec,
             MPI_DOUBLE, &bloc,sendcounts_vec[world_rank],
             MPI_DOUBLE,0, MPI_COMM_WORLD );

MPI_Scatterv(&xvec, sendcounts_vec,displs_vec,
             MPI_DOUBLE, &xloc,sendcounts_vec[world_rank],
             MPI_DOUBLE,0, MPI_COMM_WORLD );    
                 
MPI_Scatterv(&Adiag, sendcounts_vec,displs_vec,
             MPI_DOUBLE, &Adiagloc,sendcounts_vec[world_rank],
             MPI_DOUBLE,0, MPI_COMM_WORLD );  

MPI_Scatterv(&color_vec, sendcounts_vec,displs_vec,
             MPI_INT, &colorloc,sendcounts_vec[world_rank],
             MPI_INT,0, MPI_COMM_WORLD ); 
                    
   
   
int localcolors1[sendcounts_vec[world_rank]];
int localcolors2[sendcounts_vec[world_rank]];
int mat[number_of_colors][sendcounts_vec[world_rank]];
int cont[number_of_colors];  
int contjj[2];  
contjj[0]=0;
contjj[1]=0;         

 for(int ii=0;ii< number_of_colors;ii++)
  cont[ii]=0;      
 for(int ii=0;ii<sendcounts_vec[world_rank] ;ii++)  
    {//printf("COLORLOC: world_rank==%d, colorloc[%d]==%d\n",world_rank,ii,colorloc[ii] );      
 
 
  
  
 for(int jj=0;jj<number_of_colors;jj++)
  {
   if(colorloc[ii]==jj)
   {
   mat[jj][cont[jj]]=ii;
   cont[jj]=cont[jj]+1;
  // printf("world_rank==%d, mat[%d][%d]==%d, cont[%d]==%d\n",world_rank,jj,cont[jj],mat[jj][cont[jj]-1],jj,cont[jj] ); 
   }
  }           

    }                         

 
 MPI_Barrier(MPI_COMM_WORLD);                
//printf("VEDIAMO, id %d, %f \n",world_rank, Aloc);
for(int jj=0;jj<world_size;jj++)
{   
if(jj==world_rank)
{
for(int ii=0;ii< sendcounts[world_rank];ii++)
printf("Aloc id %d, %d, %f \n",world_rank,ii, Aloc[ii]); 
for(int ii=0;ii< sendcounts_vec[world_rank];ii++)
printf("vettore b id %d, %d, %f \n",world_rank,ii, bloc[ii]); 

}
MPI_Barrier(MPI_COMM_WORLD);
}

     
     MPI_Barrier(MPI_COMM_WORLD);
// smoothing steps
int smoothing_steps=100;





for(int ss=0;ss<smoothing_steps;ss++)
{
// loop on the colors: here it is easy because we know we just have two
//for(int cc=0;cc<number_of_colors ;cc++)
//printf("ss==%d,cc==%d, world_rank==%d cont[%d]==%d \n", ss,cc,world_rank,cc, cont[cc]);   
for(int cc=0;cc<number_of_colors;cc++)
 { 
    // printf("SMOOTH==%d, COLORE==%d, CONTC==%d\n", ss,cc, cont[cc]);
     //compute residual
     
        for(int ii=0;ii<sendcounts_vec[world_rank] ;ii++)   
          xloc[ii]=0;
          
//      MPI_Barrier(MPI_COMM_WORLD);  
//      printf("\n");
//      MPI_Barrier(MPI_COMM_WORLD); 
//      printf("\n");
//      MPI_Barrier(MPI_COMM_WORLD);      
//      MPI_Barrier(MPI_COMM_WORLD);  
//      printf("\n");
//      MPI_Barrier(MPI_COMM_WORLD); 
//      printf("\n");
//      MPI_Barrier(MPI_COMM_WORLD);  
    // printf("1 sendcounts_vec[%d]==%d \n",world_rank,sendcounts_vec[world_rank]);
    
     for(int ii=0;ii<sendcounts_vec[world_rank] ;ii++)   
          {
          
           double tmp=0.0;
           for(int jj=0;jj<N;jj++)
               {tmp+=Aloc[ii*N+jj]*xvec[jj];
                }
            //printf(" 2 sendcounts_vec[%d]==%d \n",world_rank,sendcounts_vec[world_rank]);
           resloc[ii]=bloc[ii]-tmp;
           //printf("ss=%d, world_rank==%d, color==%d, bloc[ii]==%f,tmp==%f, resloc[%d]==%f \n",ss,world_rank,cc,bloc[ii],tmp,ii,resloc[ii]);
           }
     
     MPI_Barrier(MPI_COMM_WORLD);      
     MPI_Scatterv(&xvec, sendcounts_vec,displs_vec,
             MPI_DOUBLE, &xloc,sendcounts_vec[world_rank],
             MPI_DOUBLE,0, MPI_COMM_WORLD ); 

        //gauss-seidel on cc color
        
       for(int ii=0;ii<cont[cc] ;ii++)
          {
          int index=mat[cc][ii];
          xloc[index]+=(resloc[index]/Adiagloc[index]);
         //printf("ss==%d,cc==%d, world_rank==%d, resloc[%d]==%f, xloc[%d]== %f \n", ss,cc,world_rank, index,resloc[index],index,xloc[index]);
          }
      
      
      MPI_Barrier(MPI_COMM_WORLD);  
      //printf(" sendcounts_vec[world_rank], "     
   
      
  MPI_Allgatherv(
    &xloc,
    sendcounts_vec[world_rank],
    MPI_DOUBLE,
    &xvec,
    sendcounts_vec,
    displs_vec,
    MPI_DOUBLE,
    MPI_COMM_WORLD );
    
    
  
       // if(world_rank==0)

            
          
 }
 
 
 
 
  //    compute residual
//      for(int ii=0;ii<sendcounts_vec[world_rank] ;ii++)   
//           {
//            double tmp=0.0;
//            for(int jj=0;jj<N;jj++)
//                {tmp+=Aloc[ii*N+jj]*xvec[jj];
//                printf("computing res on rank==%d, with xvec[%d]==%f and tmp==%f  and Aloc[ii*N+jj]==%f \n",world_rank,jj,xvec[jj],tmp,Aloc[ii*N+jj]);
//                }
//             
//            resloc[ii]=bloc[ii]-tmp;
//            printf(" smooth ==%d, on rank==%d, with resloc[%d]==%f and tmp==%f \n",ss,world_rank,ii,resloc[ii],tmp);
//            }
//          MPI_Scatterv(&xvec, sendcounts_vec,displs_vec,
//              MPI_DOUBLE, &xloc,sendcounts_vec[world_rank],
//              MPI_DOUBLE,0, MPI_COMM_WORLD ); 
//      
//  
//          
//      MPI_Barrier(MPI_COMM_WORLD);
//      for(int ii=0;ii<sendcounts_vec[world_rank];ii++)
//         {};//printf("smoothing %d resloc[ii] rank %d  component %d %f \n",ss, world_rank, ii,resloc[ii]);
// 
//      MPI_Barrier(MPI_COMM_WORLD);
//      printf("\n");
//      printf("\n");
//      printf("\n");
//      printf("\n");
//      printf("\n");
//      MPI_Barrier(MPI_COMM_WORLD);
//                 
//      gauss-seidel on first color
//        for(int ii=0;ii<sendcounts_vec[world_rank] ;ii=ii+2)
//           {xloc[ii]+=(resloc[ii]/Adiagloc[ii]);
//           printf("world_rank==%d, correctionloc[%d]== %f \n", world_rank, ii,correctionloc[ii]);
//         
//           }
//     
//     
//     
//     MPI_Barrier(MPI_COMM_WORLD);      
//     for(int ii=0;ii<sendcounts_vec[world_rank];ii++)
//         printf("smoothing %d xloc[ii] rank %d  component %d %f \n",ss, world_rank, ii,xloc[ii]);
//       
//         MPI_Allgather(
//     &xloc,
//     sendcounts_vec[world_rank],
//     MPI_DOUBLE,
//     &xvec,
//     sendcounts_vec[world_rank],
//     MPI_DOUBLE,
//     MPI_COMM_WORLD );
//       if(world_rank==0)
//       {printf(" \n");
//        for(int ii=0;ii<N;ii++)
//             printf("1 color x[ii]  %d %f \n", ii,xvec[ii]);
//         }
//         
//      recompute residual 
//     for(int ii=0;ii<sendcounts_vec[world_rank] ;ii++)   
//           {
//            double tmp=0.0;
//            for(int jj=0;jj<N;jj++)
//                tmp+=Aloc[ii*N+jj]*xvec[jj];
//            resloc[ii]=bloc[ii]-tmp;
//            }
//         MPI_Scatterv(&xvec, sendcounts_vec,displs_vec,
//              MPI_DOUBLE, &xloc,sendcounts_vec[world_rank],
//              MPI_DOUBLE,0, MPI_COMM_WORLD );      
//      update the correction and then the residual
// 
//         
//      MPI_Barrier(MPI_COMM_WORLD);
//      
//      gauss-seidel on second color
//      printf("2 color world_rank==%d, sendcounts_vec[world_rank] ==%d\n",world_rank,sendcounts_vec[world_rank]);
//        for(int ii=1;ii<sendcounts_vec[world_rank] ;ii=ii+2)
//           {xloc[ii]+=(resloc[ii]/Adiagloc[ii]);
//           printf("2 color world_rank==%d, resloc[%d]==%f,Adiagloc[%d]==%f, correctionloc[%d]== %f \n", world_rank, ii,resloc[ii],ii,Adiagloc[ii],ii,xloc[ii]);
//           }
//           
//           
//     MPI_Allgather(
//     &xloc,
//     sendcounts_vec[world_rank],
//     MPI_DOUBLE,
//     &xvec,
//     sendcounts_vec[world_rank],
//     MPI_DOUBLE,
//     MPI_COMM_WORLD );


         
}


// int vecloc[sendcounts_vec[world_rank]];
// int vecglob[N];
// 
// 
// for(int ii=0;ii<sendcounts_vec[world_rank];ii++)
// vecloc[ii]=world_rank*world_rank+1;
// 
// 
//   MPI_Allgatherv(
//     &vecloc,
//     sendcounts_vec[world_rank],
//     MPI_INT,
//     &vecglob,
//     sendcounts_vec,
//     displs_vec,
//     MPI_INT,
//     MPI_COMM_WORLD );
//  
 //    for(int ii=0;ii<N;ii++)
//             printf("world_rank==%d, vecglob[%d] = %d \n", world_rank,ii,vecglob[ii]);
//     
    for(int ii=0;ii<N;ii++)
            printf("world_rank==%d, color x[ii]  %d %f sendcounts_vec[world_rank]==%d\n", world_rank,ii,xvec[ii],sendcounts_vec[world_rank]);
    
 MPI_Barrier(MPI_COMM_WORLD);         
//  if(world_rank==0)
// {printf(" \n");
// for(int ii=0;ii<N;ii++)
// printf("x[ii]  %d %f \n", ii,xvec[ii]);
// }
 
MPI_Finalize();   
}



