from dolfin import *
from scipy.io import savemat
from petsc4py import PETSc
from scipy.sparse import lil_matrix
import numpy as np
import GaussPointsWeights as GPW 
import RaviartThomasBasisFunctions as RT 
import topology_lsfem as top 
#import input_lsfem as input



def boundary(x, on_boundary):
    return on_boundary
def A(sigma):
	return beta*sigma+alpha*tr(sigma)*Identity(gdim)
def epsilon(u):
    return sym(grad(u))


    
class Sample:
  name = ''
    
    
    
number_of_levels=2
mesh_list=list()
Patch_list=list()
W_list=list()
dm_list=list()
mesh_topology=list()
WC2WFnumpy_list=list()
coarse_to_fine_list=list()
##########################################################
################ COARSE MESH #############################  
##########################################################  
aL = 0.0
aR = 1.0
bB = 0.0
bT = 1.0
Na = 2
Nb = 2
Nc = 2
mesh = UnitCubeMesh(Na,Nb,Nc) #,"crossed"

#mesh=Mesh("tetrahedron.xml")
#mesh=refine(mesh)

gdim = mesh.geometry().dim()
mesh_points=mesh.coordinates()
N=mesh.num_vertices()


mesh_list.append(mesh)
for ll in range(1,number_of_levels):
 mesh_list.append(refine(mesh_list[ll-1])) 


# Functioncal constants: F=C_equilibrium||div sigma +f ||^2  + C_constitutive ||A(sigma)-epsilon(u) ||^2
C_equilibrium=1
C_constitutive=1.0
# Lame parameters and other related constants
Lambda = Constant(1.0)
Mu = Constant(1.0)
beta=Constant(1.0/(2*Mu))
alpha = Constant(-beta * Lambda /(2*Lambda +2*Mu))
gamma=Constant(Lambda+2*Mu);
# used tollerance
DOLFIN_TOLL=0.000000000001

# external force f and dirichlet boundary displacement g
f = Expression(("5*pi*pi*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]) - 2*pi*pi*cos(pi*x[0])*cos(pi*x[2])*sin(pi*x[1]) - 2*pi*pi*cos(pi*x[0])*cos(pi*x[1])*sin(pi*x[2])","5*pi*pi*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]) - 2*pi*pi*cos(pi*x[1])*cos(pi*x[2])*sin(pi*x[0]) - 2*pi*pi*cos(pi*x[0])*cos(pi*x[1])*sin(pi*x[2])","5*pi*pi*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]) - 2*pi*pi*cos(pi*x[1])*cos(pi*x[2])*sin(pi*x[0]) - 2*pi*pi*cos(pi*x[0])*cos(pi*x[2])*sin(pi*x[1])"),degree=4)
g = Expression(("0.0","0.0","0.0"),degree=4)

#f =  Expression(("1.0","1.0","1.0"),degree=4)



coordinates_node_="coordinates_node_"
coordinates_edge_="coordinates_edge_"
dofmap_edge_="dofmap_edge_"
dofmap_node_="dofmap_node_"
for ll in range(0,number_of_levels):
    mesh_topology.append(Sample())
    # node to edge        
    mesh_list[ll].init(0,gdim-1)
    mesh_topology[ll].N_2_F=mesh_list[ll].topology()(0,gdim-1)
    # edge to node     
    mesh_list[ll].init(gdim-1,0)
    mesh_topology[ll].F_2_N=mesh_list[ll].topology()(gdim-1,0)
    # element to node
    mesh_list[ll].init(gdim,0)
    mesh_topology[ll].K_2_N=mesh_list[ll].topology()(gdim,0)
    # element to edge
    mesh_list[ll].init(gdim,gdim-1)
    mesh_topology[ll].K_2_F=mesh_list[ll].topology()(gdim,gdim-1)
    # edge to element
    mesh_list[ll].init(gdim-1,gdim)
    mesh_topology[ll].F_2_K=mesh_list[ll].topology()(gdim-1,gdim)
    # node to element
    mesh_list[ll].init(0,gdim)
    mesh_topology[ll].N_2_K=mesh_list[ll].topology()(0,gdim)
            
    lev="{num:{width}}".format(num=ll, width=1)
    savemat(coordinates_node_+lev+".mat",{"coordinates_node":mesh_list[ll].coordinates()})
    Wtmp,dmtmp=top.create_mixed_function_space_and_dofmap(mesh_list[ll])
    W_list.append(Wtmp)
    dm_list.append(dmtmp)
    
    
    
    
    Patch_list.append( top.create_patch(mesh_list[ll],dm_list[ll],mesh_topology[ll].N_2_F)  )   

     
    coordinates_face=np.zeros((mesh_list[ll].num_facets(),gdim))
    dofmap_face=np.zeros((mesh_list[ll].num_facets(),gdim))
    dofmap_node=np.zeros((mesh_list[ll].num_vertices(),gdim))
    for ff in range(0,mesh_list[ll].num_facets()):
        coordinates_face[[ff],[0]]=Facet(mesh_list[ll],ff).midpoint().x()
        coordinates_face[[ff],[1]]=Facet(mesh_list[ll],ff).midpoint().y()
        coordinates_face[[ff],[2]]=Facet(mesh_list[ll],ff).midpoint().z()
        dofmap_face[[ff],:]=dm_list[ll].entity_dofs(mesh_list[ll],gdim-1,np.array([ff],dtype=np.uintp))   
    for nn in range(0,mesh_list[ll].num_vertices()):
        dofmap_node[[nn],:]=dm_list[ll].entity_dofs(mesh_list[ll],0,np.array([nn],dtype=np.uintp))
    savemat(coordinates_edge_+lev+".mat",{"coordinates_edge":coordinates_face})
    savemat(dofmap_edge_+lev+".mat",{"dofmap_edge":dofmap_face})
    savemat(dofmap_node_+lev+".mat",{"dofmap_node":dofmap_node})
    

for ll in range(1,number_of_levels):
     WC2WFnumpy_list.append(RT.mixed_interpolation_matrix(mesh_list[ll-1],mesh_list[ll], W_list[ll-1],W_list[ll],dm_list[ll-1],dm_list[ll],mesh_topology[ll-1],mesh_topology[ll]))
     #coarse_to_fine_list.append(from_coarse_to_fine_faces(mesh_list[ll-1],mesh_list[ll],mesh_topology[ll-1],mesh_topology[ll],  dm_list[ll-1],dm_list[ll],W_list[ll-1].dim(),W_list[ll].dim()))









savemat("WC2WFnumpy_list.mat",{"WC2WFnumpy_list":WC2WFnumpy_list[0]})

bcs=DirichletBC(W_list[0].sub(1), g, DomainBoundary())  

w = Function(W_list[0])
sigma,u = split(w)

(sigma, u) = TrialFunctions(W_list[0])
(tau, v) = TestFunctions(W_list[0])

aaa =  C_constitutive * inner(A(sigma)-epsilon(u),A(tau)-epsilon(v))*dx + C_equilibrium * inner(div(sigma),div(tau))*dx
bbb = C_equilibrium * inner(f,div(tau))*dx
solve(aaa==bbb,w,bcs)




# sigma,u = split(w)
# J =C_constitutive * inner(A(sigma)-epsilon(u),A(sigma)-epsilon(u))*dx + C_equilibrium * inner(div(sigma)+f,div(sigma)+f)*dx
# G = derivative(J,w)
#    
# solve(G==0,w,bcs)



# constant_vector = Constant(("1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0"))
# w = interpolate(constant_vector, W_list[0])
sigmacoarse,ucoarse = w.split()



ufile_pvd = File("lsfem_elasticity_u_coarse.pvd")
ufile_pvd << ucoarse
ufile_pvd = File("lsfem_elasticity_sigma_coarse.pvd")
ufile_pvd << sigmacoarse



solution_list=list()
prova=w.vector().array()
prova=prova[:,None]
solution_list.append(prova)
for ll in range(1,number_of_levels):
   tmp=WC2WFnumpy_list[ll-1]*solution_list[ll-1]
   solution_list.append(tmp)
w3 = Function(W_list[number_of_levels-1])
w3.vector()[:] = solution_list[number_of_levels-1]



# constant_vector = Constant(("1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0","1.0"))
# w3 = interpolate(constant_vector, W_list[1])
# w3.vector()[:] = 1 #
#print(solution_list)
sigma,u = w3.split()
# Gmat = assemble(G)
ufile_pvd << u


ufile_pvd = File("lsfem_elasticity_u.pvd")
ufile_pvd << u
ufile_pvd = File("lsfem_elasticity_sigma.pvd")
ufile_pvd << sigma

print( norm(ucoarse, 'H1', mesh))
print( norm(u, 'H1', mesh))
print( norm(sigmacoarse, 'Hdiv', mesh))
print( norm(sigma, 'Hdiv', mesh))


AAA = PETScMatrix() 
BBB = PETScVector()  
WWW = PETScVector()  

# M_petsc = PETSc.Mat() 
# x_petsc = PETSc.Vec()
# bdr_dofs = PETSc.IS()
# bdr_dofs.createGeneral(Patch_list[0][0])





assemble(aaa,tensor=AAA)
assemble(bbb,tensor=BBB)
solve(AAA,WWW,BBB)

DDD=assemble(aaa)
savemat("A.mat",{"A":AAA.array()})
savemat("B.mat",{"B":BBB.array()})
facce=np.arange(mesh_list[0].num_facets())
faccefine=np.arange(mesh_list[1].num_facets())

# print("coarse dofs")
# print(dm_list[0].entity_dofs(mesh_list[0],gdim-1,np.array(facce,dtype=np.uintp)))
# print("fine dofs")
# print(dm_list[1].entity_dofs(mesh_list[1],gdim-1,np.array(faccefine,dtype=np.uintp)))

# print(type(AAA))
# print(AAA.mat().getSubMatrix(bdr_dofs,bdr_dofs, M_petsc))
# print(M_petsc.assemble())
# print(bdr_dofs)
# AAA.mat().getSubMatrix(bdr_dofs,bdr_dofs, M_petsc)
# AAA.mat().getSubMatrix(bdr_dofs,bdr_dofs, M_petsc)

bdr_dofs = PETSc.IS()
bdr_dofs.createGeneral(Patch_list[0][0])
bdr_dofs_tot = PETSc.IS()
bdr_dofs_tot.createGeneral(np.arange(0,W_list[0].dim(),dtype=np.int32) )
M_petsc = PETSc.Mat()
b_petsc = PETSc.Vec()
x_petsc = PETSc.Vec()



res_petsc = PETSc.Vec()
b_petsc_tot = PETSc.Vec()
BBB.vec().getSubVector(bdr_dofs_tot, b_petsc_tot)

#print(type(BBB))
# BBB=assemble(bbb)
# solve(AAA,
# bcs.apply(AAA, BBB)
# print(type(AAA))
# ss=Function(W_list[0])
# SS=ss.vector()


# print(mesh_list[0].color)
# colore=mesh_list[0].color
# cell_colors(mesh_list[0])
# print(type(colore))
w.vector()[:]=0
for ee in range(0,1):
 print(ee)
 for nn in range(0,mesh_list[0].num_vertices()):
      bdr_dofs.createGeneral(Patch_list[0][nn])
      A_petsc = PETSc.Mat()
      A_mxn_petsc = PETSc.Mat()
      
      AAA.mat().getSubMatrix(bdr_dofs, bdr_dofs, A_petsc)
      AAA.mat().getSubMatrix(bdr_dofs, bdr_dofs_tot, A_mxn_petsc)
      
      BBB.vec().getSubVector(bdr_dofs, b_petsc)
      BBB.vec().getSubVector(bdr_dofs, x_petsc)
      
      A_loc = PETScMatrix(A_petsc)
      b_loc = PETScVector(b_petsc)
      x_loc = PETScVector(x_petsc)
      
      
      #A = as_backend_type(AAA).mat()
      #print(A)
      #prova=A_mxn_petsc*b_petsc_tot
      #prova=prova.mat()

      #res_loc=b_loc-prova
      
      solve(A_loc,x_loc,b_loc)
#       AAA.mat().getSubMatrix(bdr_dofs, bdr_dofs_tot, M2_petsc)
#       MatMult(M2_petsc,x_petsc)
#       A_loc = PETScMatrix(A_petsc)
#       b_loc = PETScVector(b_petsc)
#       x_loc = PETScVector(x_petsc)
#       solve(A_loc,x_loc,b_loc)
#       w.vector()[Patch_list[0][nn]]=x_loc
#       
      
      
    
      
#     print()
#     print(Patch_list[0][nn])
#     a_loc=AAA.array()[np.ix_(Patch_list[0][nn],Patch_list[0][nn])]
#     print("a_loc")
#     print(type(a_loc))
#     b_loc=bbb.array()[Patch_list[0][nn]]
#     #b_loc=b_loc[:,None]
#     print("b_loc")
#     print(b_loc)
#     print("SS[Patch_list[0][nn]]")
#     print(SS[Patch_list[0][nn]])
#     print(type(a_loc))
#     print(type(b_loc))
#     print(type(SS[Patch_list[0][nn]]))
#     print()
#     solve(a_loc, SS[Patch_list[0][nn]], b_loc)
# 


# systemmatrix = derivative(G,w)
# systemmatrix=assemble(systemmatrix)
#print("systemmatrix.array()")
#print(systemmatrix.array())

# savemat("solutioncoarse.mat",{"solc":solutioncoarse})
# savemat("solutionfine.mat",{"solf":solutionfine})
# savemat("WC2WFnumpy.mat",{"WC2WFnumpy":WC2WFnumpy})
# savemat("coordinates_coarse.mat",{"coordinates_coarse":mesh_points})
# savemat("mesh_points2.mat",{"mesh_points2":mesh_points2})
# coordinates_edge_coarse=np.zeros((mesh.num_facets(),2))
# coordinates_edge_fine=np.zeros((mesh2.num_facets(),2))
# dofmap_edge_coarse=np.zeros((mesh.num_facets(),2))
# dofmap_edge_fine=np.zeros((mesh2.num_facets(),2))
# dofmap_node_coarse=np.zeros((mesh.num_vertices(),2))
# dofmap_node_fine=np.zeros((mesh2.num_vertices(),2))
# for ee in range(0,mesh.num_facets()):
#     coordinates_edge_coarse[[ee],[0]]=Facet(mesh,ee).midpoint().x()
#     coordinates_edge_coarse[[ee],[1]]=Facet(mesh,ee).midpoint().y()
#     dofmap_edge_coarse[[ee],[0,1]]=dm.entity_dofs(mesh,1,np.array([ee],dtype=np.uintp)) 
# for ee in range(0,mesh2.num_facets()):
#     coordinates_edge_fine[[ee],[0]]=Facet(mesh2,ee).midpoint().x()
#     coordinates_edge_fine[[ee],[1]]=Facet(mesh2,ee).midpoint().y()
#     dofmap_edge_fine[[ee],[0,1]]=dm2.entity_dofs(mesh2,1,np.array([ee],dtype=np.uintp)) 
#     
# for nn in range(0,mesh.num_vertices()):
#     dofmap_node_coarse[[nn],[0,1]]=dm.entity_dofs(mesh,0,np.array([nn],dtype=np.uintp))
# 
# for nn in range(0,mesh2.num_vertices()):
#     dofmap_node_fine[[nn],[0,1]]=dm2.entity_dofs(mesh2,0,np.array([nn],dtype=np.uintp))
# savemat("coordinates_edge_coarse.mat",{"coordinates_edge_coarse":coordinates_edge_coarse})
# savemat("coordinates_edge_fine.mat",{"coordinates_edge_fine":coordinates_edge_fine})
# savemat("dofmap_edge_coarse.mat",{"dofmap_edge_coarse":dofmap_edge_coarse})
# savemat("dofmap_edge_fine.mat",{"dofmap_edge_fine":dofmap_edge_fine})
# savemat("dofmap_node_coarse.mat",{"dofmap_node_coarse":dofmap_node_coarse})
# savemat("dofmap_node_fine.mat",{"dofmap_node_fine":dofmap_node_fine})






 
# v =Function(W)
# bmesh = BoundaryMesh(mesh2, 'exterior')
# coordinates_bmesh=bmesh.coordinates()


#from scipy.io import savemat
#savemat("mio.mat",Gmat.array())

#sigma1 = interpolate(sigma_raw1,RT)
#sigma2 = interpolate(sigma_raw2,RT)
#sigma_tensor= interpolate(sigma_tot,RTtensor)
#err_tensor=sigma_tensor - sigma
#sigmaxx = interpolate(sigmaxx_ex, feRT)
#sigmaxy = project(sigmaxy_ex, feRT)
#sigmayy = project(sigmayy_ex, feRT)

# Save solution in VTK format
#ufile_pvd = File("lsfem_elasticity_error_sigma.pvd")
#ufile_pvd << err_tensor




# Define Dirichlet boundary (x = 0 or x = 1)
#def boundary_L(x):
#    return x[0] < aL + DOLFIN_EPS
#def boundary_B(x):
#    return x[1] < bB + DOLFIN_EPS
#def boundary_R(x):
#    return x[0] > aR - DOLFIN_EPS
#def boundary_T(x):
#    return x[1] > bB - DOLFIN_EPS

#class BOUNDARY_TOP(SubDomain):
#	def inside(self,x,on_boundary):return ((x[1] > bB - DOLFIN_EPS) and on_boundary)
#class BOUNDARY_RIGHTBOTTOM(SubDomain):
#	def inside(self,x,on_boundary):return (((x[0] > aR - DOLFIN_EPS)or (x[1] < bB + DOLFIN_EPS)) and on_boundary)

# sigmaxx_ex=Expression("G *  pi * cos(pi*x[0])*sin(pi*x[1])+ L * pi * sin(pi*x[0])*cos(pi*x[1]) ", G=gamma,L=Lambda,degree=10)
# sigmaxy_ex=Expression("M *pi*sin(pi*x[0])*cos(pi*x[1])+ M *pi*sin(pi*x[1])*cos(pi*x[0])",M=Mu,degree=4)
# sigmayy_ex=Expression("L * pi * cos(pi*x[0])*sin(pi*x[1])  + G  * pi * sin(pi*x[0])*cos(pi*x[1])",L=Lambda,G=gamma,degree=4)
# sigma_raw1=Expression(("G *  pi * cos(pi*x[0])*sin(pi*x[1])+ L * pi * sin(pi*x[0])*cos(pi*x[1]) ","M *pi*sin(pi*x[0])*cos(pi*x[1])+ M *pi*sin(pi*x[1])*cos(pi*x[0])"), G=gamma,L=Lambda,M=Mu,degree=10)
# sigma_raw2=Expression(("M *pi*sin(pi*x[0])*cos(pi*x[1])+ M *pi*sin(pi*x[1])*cos(pi*x[0])","L * pi * cos(pi*x[0])*sin(pi*x[1])  + G  * pi * sin(pi*x[0])*cos(pi*x[1])"), G=gamma,L=Lambda,M=Mu,degree=10)
# sigma_tot=Expression((("G *  pi * cos(pi*x[0])*sin(pi*x[1])+ L * pi * sin(pi*x[0])*cos(pi*x[1]) ","M *pi*sin(pi*x[0])*cos(pi*x[1])+ M *pi*sin(pi*x[1])*cos(pi*x[0])"),("M *pi*sin(pi*x[0])*cos(pi*x[1])+ M *pi*sin(pi*x[1])*cos(pi*x[0])","L * pi * cos(pi*x[0])*sin(pi*x[1])  + G  * pi * sin(pi*x[0])*cos(pi*x[1])")), G=gamma,L=Lambda,M=Mu,degree=10)
# h = Expression((("0.0","0.0"),("0.0","0.00")),degree=2)
#j = Expression(("1.0","0.00"),degree=2)
# j = Expression((("0.0","-0.005"),("0.0","0.0")),degree=2)
# kk = Expression(("0.0","-0.005"),degree=2)

#pc = mesh2.data().array("parent_cell", mesh2.topology().dim())
#print(pc)
# V1=FunctionSpace(mesh,"P",1)
# V2=FunctionSpace(mesh2,"P",1)
# mat=PETScDMCollection.create_transfer_matrix(V1,V2)
# print(mat.array())



#WC2WF.setValues([W_num_dofs-1],[W_num_dofs-1],1)
#WC2WF.assemble()
#print(WC2WF.getValues(range(15,W_num_dofs), range(15,W_num_dofs)))
#print(type(WC2WF))

#viewer=PETSc.Viewer().createBinary('test.dat', 'w')
#viewer(WC2WF)
#BB=WC2WF.convert("dense")
#BB.getDenseArray()
#print(BB)


#element = W.element()




#dofs_x = W.tabulate_dof_coordinates().reshape(W_num_dofs,gdim)


 
# for dofs, dofs_x in zip(dofs, dofs_x):
#      print(dofs)
#      print( ':') 
#      print(dofs_x)



# for c in cells(mesh):
#  print(W.dofmap().cell_dofs(c.index()))
  #    for v in vertices(c):
   #    print(v)
     #  Prova2[v]=Prova2[v]+1 





#marker = MeshFunction("bool", mesh, mesh.topology().dim(), True)




#print(dm2.entity_dofs(mesh2,0,np.array(nodes,dtype=np.uintp)))
#print(mesh_points)
#print(mesh_points2)
#print(dm.entity_dofs(mesh,0,np.array(nodes,dtype=np.uintp)))
#print(WC2WFnumpy)
#print()
#print("fine edges")

#for ee in range(0,30):
   #provami1=[Facet(mesh2,ee).midpoint().x(),Facet(mesh2,ee).midpoint().y()] 
   #provami2=dm2.entity_dofs(mesh2,1,np.array([ee],dtype=np.uintp))
   #print(ee)
   #print(provami1)
   #print(provami2)
    
#print()
#print("coarse edges")
#for ee in range(0,9):
   #provami1=[Facet(mesh,ee).midpoint().x(),Facet(mesh,ee).midpoint().y()] 
   #provami2=dm.entity_dofs(mesh,1,np.array([ee],dtype=np.uintp))
   #print(provami1)
   #print(provami2)
    
#print()
#print("fine nodes")
#for nn in range(0,15):
   #provami1=mesh_points2[nn]
   #provami2=dm2.entity_dofs(mesh2,0,np.array([nn],dtype=np.uintp))
   #print(provami1)
   #print(provami2)
#print()
#print("coarse nodes")
#for nn in range(0,6):
   #provami1=mesh_points[nn]
   #provami2=dm.entity_dofs(mesh,0,np.array([nn],dtype=np.uintp))
   #print(provami1)
   #print(provami2)
   
   
   
   ############## PETSC MATRIX
# WC2WF= PETSc.Mat().create()
# WC2WF.setSizes([W2_num_dofs, W_num_dofs])
# WC2WF.setType("aij")
# WC2WF.setUp()
   
# WC2WF= PETSc.Mat().create()
# WC2WF.setSizes([W2_num_dofs, W_num_dofs])
# WC2WF.setType("aij")
# WC2WF.setUp()
         #                WC2WF.setValues(tmp[0],[tmpC1[0], tmpC2[0]],[0.5, 0.5])
#                         WC2WF.setValues(tmp[1],[tmpC1[1], tmpC2[1]],[0.5, 0.5])

   #WC2WF.assemble()
   
   
   
   
   
   
   
   
   
   
   
   # boundaries = FacetFunction("size_t", mesh)
# boundaries.set_all(0)
# BOUNDARY_TOP().mark(boundaries, 1)
# BOUNDARY_RIGHTBOTTOM().mark(boundaries, 2)

# ds = Measure('ds', domain = mesh, subdomain_data = boundaries)

#J =C_constitutive * inner(A(sigma)-epsilon(u),A(sigma)-epsilon(u))*dx + C_equilibrium * inner(div(sigma)+f,div(sigma)+f)*dx+ C_boundary* inner(dot(sigma,FacetNormal(mesh))-mm,dot(sigma,FacetNormal(mesh))-mm)*ds(1)+C_boundary* inner(dot(sigma,FacetNormal(mesh))-nn,dot(sigma,FacetNormal(mesh))-nn)*ds(2)

#bc1 = DirichletBC(W.sub(1), g, boundary_L)
#bc2 = DirichletBC(W.sub(1), g, boundary_B)
#bc4 = DirichletBC(W.sub(1), g, boundary_T)
#bc3 = DirichletBC(W.sub(1), g, boundary_R)


#bc1 = DirichletBC(W.sub(1), g, boundary_L)
#bc2 = DirichletBC(W.sub(0), h, boundary_B)
#bc4 = DirichletBC(W.sub(0).sub(1), kk, boundary_T)
#bc3 = DirichletBC(W.sub(0), h, boundary_R)
#bcs=[bc1,bc2,bc4,bc3]



# w.vector()[:] = 1 #0
# w.vector()[1]=0
# w.vector()[2]=0
# w.vector()[12]=0
# w.vector()[14]=0
# w.vector()[13]=0
# w.vector()[15]=0
#constant_vector = Constant(("1.0","1.0","1.0","1.0","1.0","1.0"))
#w = interpolate(constant_vector, W)

