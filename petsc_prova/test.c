#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include <stdio.h>


typedef struct node {
    struct node * next;
    PetscInt *val;
    int size;
} patch_list;



void patch_initialize(patch_list* head, int values) {
    head = malloc(sizeof(patch_list));
    head->next=NULL;
    head->size=values; 
    head->val =malloc(sizeof values * values);
    for(int nn=0;nn<values;nn++)
    {
        head->val[nn] = nn;
        printf("ooo: %d \n",head->val[nn]);
        }
     
}

void patch_push(patch_list* head, int values) {
    printf("qui3,  \n");
    patch_list* current = head; 
    printf("qui333,  \n");
    if(current->next==NULL)
     printf("222 e' null");
    while (current->next != NULL) {
        current = current->next;
    }
    printf("qui4, \n");
    /* now we can add a new variable */
    current->next = malloc(sizeof(patch_list));
    current->next->val =malloc(sizeof values * values);
    current->next->size=values;
    printf("VALUE: %d \n",values);
    for(int nn=0;nn<values;nn++)
    {current->next->val[nn] = nn;
     printf("%d \n",current->next->val[nn]);};
    
//     printf(" \n");
    current->next->next = NULL;
}

void print_list(patch_list * head) {
    patch_list * current = head;

    while (current != NULL) {
        printf("%d\n", current->val);
        current = current->next;
    }
}

int main(int argc, char **argv)
{ 


  PetscInitialize(&argc, &argv, NULL, NULL);
  PetscScalar minusone=-1;
  PetscInt i,j,N=4;
  double* a;
  Vec x;
  Vec Ax;
  Vec b;
  Vec res;
  Mat Amat;
  
// initialize vectors 
VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE,N, &b);
VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE,N, &x);
VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE,N, &Ax);
VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE,N, &res);
VecSet(x, 1.);
VecSet(b, 6.);


// initialize matrix
MatCreateAIJ(MPI_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,N,N,N,NULL,N,NULL,&Amat);
//MatSetFromOptions(Amat);
MatSetUp(Amat);
MatSetValue(Amat,0,0,1,INSERT_VALUES);
for (i=1; i<N-1; i++) {
     MatSetValue(Amat,i,i-1,1,INSERT_VALUES);
     MatSetValue(Amat,i,i,2,INSERT_VALUES);
     MatSetValue(Amat,i,i+1,1,INSERT_VALUES);
      }
MatSetValue(Amat,N-1,N-1,1,INSERT_VALUES);


MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);




// Compute the residual res=b-Ax
MatMult(Amat,x,Ax);  
VecWAXPY(res,minusone,Ax,b);


// View vectors and matrix
MatView(Amat,PETSC_VIEWER_STDOUT_WORLD);
VecView(b, PETSC_VIEWER_STDOUT_WORLD); 
VecView(x, PETSC_VIEWER_STDOUT_WORLD);
VecView(Ax, PETSC_VIEWER_STDOUT_WORLD);
VecView(res, PETSC_VIEWER_STDOUT_WORLD);




KSP ksp;
KSPCreate(MPI_COMM_WORLD,&ksp);
KSPSetOperators(ksp,Amat,Amat);
KSPSolve(ksp,b,x);
KSPDestroy(&ksp);

VecView(x, PETSC_VIEWER_STDOUT_WORLD);





//   
//   
//   
//   MatCreate(PETSC_COMM_WORLD,&Amat); 
//   //MatSetType(Amat,MATSEQMAIJ);
//   MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,N,N);  
//   
// 
// 
//   PetscMalloc1(N*N,&a);
  // MatDiagonalSet(Amat,x,INSERT_VALUES);
  //   for(int nn1=0;nn1<N;nn1++)
  //      for(int nn2=0;nn2<N;nn2++)
  //         MatSetValues(Amat,nn1,&nn1,nn2,&nn2, nn1+nn2,INSERT_VALUES) ;
//   for ( PetscInt i=0; i<N; i++) 
// 		for ( PetscInt j=0; j<N; j++) 
// 			a[i*N+j ] = i*N+j;
//   for ( PetscInt  i=0; i<N; i++) 
// 		for ( PetscInt  j=0; j<N; j++) 
//               {printf("Sto inserendo valori nella matrice %d, %d, %f \n",i, j, a[i*N+j]);
//                MatGetSize(Amat,&m,&n);
//                printf("matrix size: %d, %d \n",m,n);
//                MatSetValue(Amat,i,j,1.0,INSERT_VALUES);
//               }
              
// MatSeqDenseSetPreallocation(Amat,a);
// MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY); 
// MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY);
// MatView(Amat, PETSC_VIEWER_STDOUT_WORLD);
// 
// 
// VecSet(x, 1.);
// VecSet(Ax, 2.);
// VecView(x, PETSC_VIEWER_STDOUT_WORLD);
// VecView(Ax, PETSC_VIEWER_STDOUT_WORLD);
// 
// MatGetSize(Amat,&m,&n);
// printf("matrix: %d, %d \n",m,n);
// 
// VecGetSize(x,&m);
// printf("vector: %d \n",m);
// 
// VecGetSize(Ax,&m);
// printf("Ax vector: %d \n",m);
//    


//   
//    v
//   VecSet(v, 2.);
//   VecView(x, PETSC_VIEWER_STDOUT_WORLD);
//   VecView(v, PETSC_VIEWER_STDOUT_WORLD);
//   VecView(Ax, PETSC_VIEWER_STDOUT_WORLD);
//   MatView(Amat, PETSC_VIEWER_STDOUT_WORLD);
//   
//   PetscInt m,n;
//    MatGetSize(Amat,&m,&n);

struct node* patch_mesh[N];
//IS        *is;
for(int n=0;n<N;n++)
    {
    patch_mesh[n]=malloc(sizeof(patch_list));
    patch_mesh[n]->val=malloc(sizeof(n)*(n+1));
    patch_mesh[n]->size=n+1;
     for(int m=0; m<=n;m++)
         {
         patch_mesh[n]->val[m]=m;
         printf("Value: %d \n",patch_mesh[n]->val[m]);
         }
         printf(" \n");
    }
    
    
// loop on the nodes
PetscInt n_loc;
PetscInt *idx;
PetscScalar *v;

for(int n=0;n<N;n++)
{  



  n_loc=patch_mesh[n]->size;
  idx=malloc(sizeof(PetscInt)*n_loc); 
  idx=patch_mesh[n]->val;
  v=malloc(sizeof(PetscScalar)*n_loc*n_loc);
  printf("idx: %d\n",idx[0]);
  printf("n_loc: %d\n",n_loc);
  
 MatGetValues(Amat,n_loc,idx,n_loc,idx, v);
   for(int ii=0;ii<n_loc*n_loc;ii++)
      printf("vvv: %f\n", v[ii]);
  printf("\n");
  printf("\n");
  
  
  Vec x_loc;
  Vec b_loc;

  VecCreate(PETSC_COMM_SELF,&x_loc);
  VecSetSizes(x_loc,PETSC_DECIDE, n_loc);
  
  Mat Aloc;
  MatCreateSeqDense(PETSC_COMM_SELF,n_loc,n_loc,NULL,&Aloc);
  MatSetUp(Aloc);
  
  for(int ii=0;ii<n_loc;ii++)
      for(int jj=0;jj<n_loc;jj++)
          MatSetValue(Aloc,*(idx+ii),*(idx+jj),*(v+ii*n_loc+jj),INSERT_VALUES);
  MatAssemblyBegin(Aloc, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(Aloc, MAT_FINAL_ASSEMBLY);  
  MatView(Aloc, PETSC_VIEWER_STDOUT_WORLD);
  
  
  IS is;
  const PetscInt *indices;
  ISCreateGeneral(PETSC_COMM_SELF,n_loc,idx,PETSC_COPY_VALUES,&is);
  ISGetIndices(is,&indices);
  ISView(is,PETSC_VIEWER_STDOUT_SELF);
  VecGetSubVector(res,is,&b_loc);
  VecView(b_loc, PETSC_VIEWER_STDOUT_WORLD);
  ISRestoreIndices(is,&indices);
  ISDestroy(&is);
  
  KSP ksp_loc;
  KSPCreate(PETSC_COMM_SELF,&ksp_loc);
  KSPSetOperators(ksp_loc,Aloc,Aloc);
  KSPSolve(ksp_loc,b_loc,x_loc);
  KSPDestroy(&ksp_loc);


  
  MatDestroy(&Aloc);
  free(v);
  free(idx);
}





// Free memory
VecDestroy(&b);
VecDestroy(&x);
VecDestroy(&Ax);
VecDestroy(&res);
MatDestroy(&Amat);

// patch_list* head = NULL;
// patch_initialize(head,10);
// printf("qui2 %d, \n",head->val[0]);
// //head = malloc(sizeof(patch_list));
// //head->next = NULL;
//     if(head->next==NULL)
//      printf("3333 e' null");
//     else
//      printf("3333 non e' null");
//  for(int n=1;n<=10;n++)
//       patch_push(head, n );
// printf("qui, \n");
// int cont=0;
// while (head->next != NULL) {
//         printf("%d, \n",head->val[cont]); 
//         cont++;
//         head = head->next;
//     }
//     
// head = malloc(sizeof(patch_list));
// head->val=malloc(23 * sizeof *head->val);
// head->size=1;
// head->val[0]=1000;
// head->next = NULL;
// size_t sz = sizeof head->val ;/// sizeof head->val[0];
// printf("size = %d, %d\n",  m,n);
// //printf("size = %zu\n",  head->val[1]);
// // head->val = 1;
// // head->next = NULL;
// // 
// // 
// for(int n=1;n<=10;n++)
//     push(head, n );
// // if (head == NULL) {
// //     return 1;
// // }
// // print_list(head);
// 
// 
// fprintf(stderr, "Hello, please enter your age\n");
  PetscFinalize();
return 0; 
}