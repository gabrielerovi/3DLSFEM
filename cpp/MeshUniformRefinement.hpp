#include <dolfin.h>
#include "Topology.hpp"
#ifndef MESHUNIFORMFENICS_H
#define MESHUNIFORMFENICS_H

class MeshUniformRefinement{

private:
unsigned int _number_of_levels;
std::vector<std::shared_ptr<Mesh> > _mesh;
std::vector<std::shared_ptr<Mesh> > init_mesh(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &number_of_levels);

public:
MeshUniformRefinement(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &number_of_levels);// constructor
//~MeshUniformRefinement(void);

unsigned int number_of_levels(void) const;
std::vector<std::shared_ptr<Mesh> > list(void) const;
};

unsigned int MeshUniformRefinement::number_of_levels(void) const 
{return _number_of_levels;}
std::vector<std::shared_ptr<Mesh> > MeshUniformRefinement::list(void) const 
{return _mesh;}



std::vector<std::shared_ptr<Mesh> > MeshUniformRefinement::init_mesh(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &number_of_levels )
{
std::vector<std::shared_ptr<Mesh> > mesh_list(number_of_levels);
mesh_list[0]=mesh;
for(int lev=1;lev<number_of_levels;lev++)
  mesh_list[lev]= std::make_shared<Mesh>(refine(*(mesh_list[lev-1])));
return mesh_list;
}



MeshUniformRefinement::MeshUniformRefinement(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &number_of_levels)
:
_number_of_levels(number_of_levels),
_mesh(init_mesh(mesh,number_of_levels))
{}



// MeshUniformRefinement::~MeshUniformRefinement(void)
// {
// for(int lev=0;lev<_number_of_levels;lev++)
//    _mesh[lev]->clean();
// }
#endif



