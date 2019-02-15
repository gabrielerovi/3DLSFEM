#include <dolfin.h>
#include "Topology.hpp"

#ifndef UNIFORMFENICS_H
#define UNIFORMFENICS_H

class UniformRefinement{

private:
unsigned int _number_of_levels;
std::vector<std::shared_ptr<Mesh> > _mesh;
std::vector<Topology > _topology;

std::vector<std::shared_ptr<Mesh> > init_mesh(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &number_of_levels);
std::vector<Topology> init_topology(const std::vector<std::shared_ptr<Mesh> > &mesh_list,const unsigned int &number_of_levels,const unsigned int &gdim);

public:
UniformRefinement(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &number_of_levels,const unsigned int gdim);// constructor
~UniformRefinement(void);

unsigned int number_of_levels(void){return _number_of_levels;};
std::vector<std::shared_ptr<Mesh> > mesh(void){return _mesh;};
std::vector<Topology > topology(void){return _topology;};

};

std::vector<std::shared_ptr<Mesh> > UniformRefinement::init_mesh(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &number_of_levels )
{
std::vector<std::shared_ptr<Mesh> > mesh_list(number_of_levels);
mesh_list[0]=mesh;
for(int lev=1;lev<number_of_levels;lev++)
  mesh_list[lev]= std::make_shared<Mesh>(refine(*(mesh_list[lev-1])));
return mesh_list;
}


std::vector<Topology> UniformRefinement::init_topology(const std::vector<std::shared_ptr<Mesh> > &mesh_list,const unsigned int &number_of_levels,const unsigned int &gdim)
{
std::vector<Topology> topology_list;
for(int lev=0;lev<number_of_levels;lev++)
 topology_list.push_back(Topology(mesh_list[lev],gdim));
return topology_list;
}


UniformRefinement::UniformRefinement(const std::shared_ptr<dolfin::Mesh> &mesh,const unsigned int &number_of_levels,const unsigned int gdim)
:
_number_of_levels(number_of_levels),
_mesh(init_mesh(mesh,number_of_levels)),
_topology(init_topology(_mesh,number_of_levels,gdim))
{}



UniformRefinement::~UniformRefinement(void)
{
for(int lev=0;lev<_number_of_levels;lev++)
   _mesh[lev]->clean();
for(int lev=0;lev<_number_of_levels;lev++)
   _topology[lev].clear();
}
#endif



