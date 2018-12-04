#include "petscvec.h"
#include "petscmat.h"
#include <stdio.h>


typedef struct node {
    struct node * next;
    int *val;
    int size;
} patch_list;



void push(patch_list* head, int values) {
    patch_list* current = head;
    while (current->next != NULL) {
        current = current->next;
    }
    /* now we can add a new variable */
    current->next = malloc(sizeof(patch_list));
    current->next->val =malloc(sizeof values * values);
    current->next->size=values;
    for(int nn=0;nn<values;nn++)
    {current->next->val[nn] = nn;};
     printf(" \n");
    printf("%d \n",current->next->size);
    printf(" \n");
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



  Vec x;
  Vec v;
  Mat Amat;
  PetscInitialize(&argc, &argv, NULL, NULL);
  VecCreateSeq(PETSC_COMM_SELF, 10, &x);
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,10, &v);
  
  
  int N=9;
  MatCreate(PETSC_COMM_WORLD,&Amat); 
  MatSetType(Amat,MATSEQAIJ);
  MatSetSizes(Amat,PETSC_DECIDE,PETSC_DECIDE,N,N);  
  
  MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY);
  
  
  VecSet(x, 1.);
  VecSet(v, 2.);
  VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  VecView(v, PETSC_VIEWER_STDOUT_WORLD);
  
  
patch_list* head = NULL;
head = malloc(sizeof(patch_list));
head->val=malloc(23 * sizeof *head->val);
head->size=1;
head->val[0]=1000;
head->next = NULL;
size_t sz = sizeof head->val ;/// sizeof head->val[0];
printf("size = %zu\n",  head->val[0]);
//printf("size = %zu\n",  head->val[1]);
// head->val = 1;
// head->next = NULL;
// 
// 
for(int n=1;n<=10;n++)
    push(head, n );
// if (head == NULL) {
//     return 1;
// }
// print_list(head);


fprintf(stderr, "Hello, please enter your age\n");
  PetscFinalize();
return 0; 
}