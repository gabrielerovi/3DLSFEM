#include <dolfin.h>
#include "Topology.hpp"
#include "MeshUniformRefinement.hpp"

#ifndef UNIFORMFENICS_H
#define UNIFORMFENICS_H

class TopologyUniformRefinement{

public:
TopologyUniformRefinement(const MeshUniformRefinement &meshL,const unsigned int &number_of_levels,const unsigned int gdim);// constructor
~TopologyUniformRefinement(void);
unsigned int number_of_levels(void) const ; 
std::vector<Topology > list(void) const;


private:
unsigned int _number_of_levels;
std::vector<Topology > _topology;
std::vector<Topology> init_topology(const MeshUniformRefinement &meshL,const unsigned int &number_of_levels,const unsigned int &gdim);


};

unsigned int TopologyUniformRefinement::number_of_levels(void) const {return _number_of_levels;};
std::vector<Topology > TopologyUniformRefinement::list(void) const {return _topology;}

std::vector<Topology> TopologyUniformRefinement::init_topology(const MeshUniformRefinement &meshL,const unsigned int &number_of_levels,const unsigned int &gdim)
{
std::vector<Topology> topology_list;
for(int lev=0;lev<number_of_levels;lev++)
 topology_list.push_back(Topology(meshL.list().at(lev),gdim));
return topology_list;
}


TopologyUniformRefinement::TopologyUniformRefinement(const MeshUniformRefinement &meshL,const unsigned int &number_of_levels,const unsigned int gdim)
:
_number_of_levels(number_of_levels),
_topology(init_topology(meshL,number_of_levels,gdim))
{}


TopologyUniformRefinement::~TopologyUniformRefinement(void)
{
for(int lev=0;lev<_number_of_levels;lev++)
   _topology[lev].clear();
}
#endif



