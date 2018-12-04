#include "petscvec.h"
#include <iostream>


typedef struct node {
    int val;
    struct node * next;
} node_t;

int main(int argc, char **argv)
{ 



  Vec x;
  Vec v;
  PetscInitialize(&argc, &argv, NULL, NULL);
  VecCreateSeq(PETSC_COMM_SELF, 10, &x);
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE,10, &v);
  VecSet(x, 1.);
  VecSet(v, 2.);
  VecView(x, PETSC_VIEWER_STDOUT_WORLD);
  VecView(v, PETSC_VIEWER_STDOUT_WORLD);
  
  
  node_t * head = NULL;
head = malloc(sizeof(node_t));
if (head == NULL) {
    return 1;
}

head->val = 1;
head->next = NULL;

std::cout << "static destructor\n";
  PetscFinalize();
return 0; 
}