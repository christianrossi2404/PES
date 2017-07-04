% This program solves a convection-diffusion problem 
% in a square domain [0,1]x[0,1] using bilinear elements.
% 

clear, close all, home

global diffusion  h 

disp(' ')
disp('No source term is considered');
diffusion = input('Diffusion coefficient = ');
disp(' ')

% GEOMETRY
% Matrix of nodel coordinates and conectivities

disp(' Number of dimensions')
disp('2: 2D')
disp('3: 3D')
ndim=input('nº Dimensions?   ');
if ndim ~= [2 3] 
    disp('ERROR: This program only works in 2D & 3D');
    clear all
    return
end
disp(' ')
disp(' Select the element type');
if ndim==2
    disp('1: Linear Triangulars');
    disp('2:   ||   Quadrilaterals');
    disp('3: Cuadratic Triangulars');
    disp('4:    ||     Quadrilaterals');
elseif ndim==3
    disp('1: Linear Tetrahedrals');
    disp('2:   ||   Hexahedrals');
    disp('3: Cuadratic Tetrahedrals');
    disp('4:    ||     Hexahedrals');
end
element=input('Element type   ');
if element ~= [1 2 3 4]
    disp ('ERROR: Select an existing element type');
    clear all
    return
end

tic

addpath('Nodes')
addpath('Elements')
addpath('BoundaryConditions')
if ndim==2
    if element == 1
        X=load('nodes_2D_tri_linear');
        T=load('elem_2D_tri_linear');
    elseif element == 2
        X=load('nodes_2D_quad_linear');
        T=load('elem_2D_quad_linear');
    elseif element == 3
        X=load('nodes_2D_tri_quad');
        T=load('elem_2D_tri_quad');
    elseif element == 4
        X=load('nodes_2D_quad_quad');
        T=load('elem_2D_quad_quad');
    end
elseif ndim==3
    if element == 1
        X=load('nodes_3D_tri_lin.txt');
        T=load('elem_3D_tri_lin.txt');
    elseif element == 2
        X=load('nodes_3D_quad_lin.txt');
        T=load('elem_3D_quad_lin.txt');
    elseif element == 3
        X=load('nodes_3D_tri_quad.txt');
        T=load('elem_3D_tri_quad.inp');
    elseif element == 4
        X=load('nodes_3D_quad_quad.txt');
        T=load('elem_3D_quad_quad.txt');
    end
end
        X=X(:,2:end);
        T=T(:,2:end);

[nnode,ncoord] = size(X);

%number of nodes in each element

if ndim==2
    if element == 1
        nelnodes = 3;
        npoints = 1;
    elseif element == 2
        nelnodes = 4;
        npoints = 4;
    elseif element == 3
        nelnodes = 6;
        npoints = 3;
    elseif element == 4
        nelnodes = 8;
        npoints = 9;  
    end
elseif ndim==3
    if element == 1
        nelnodes = 4;
        npoints = 1;
    elseif element == 2
        nelnodes = 8;
        npoints = 8;
    elseif element == 3
        nelnodes = 10;
        npoints = 4;
    elseif element == 4
        nelnodes = 20;
        npoints = 27;  
    end
end

%triangular 1; quadrilateral 2
if element ==1 && 3
    elident = 1;
else
    elident = 2;
end

% NUMERICAL INTEGRATION
% Quadrature,Shape Functions

n = numberofintegrationpoints(ncoord,nelnodes,elident);
if ndim==2
    w = iWeight(ncoord,nelnodes,npoints,elident);
    xi = iPoints(ncoord,nelnodes,npoints,elident);
    if element == 1 
        N = shapefunctions(nelnodes,ncoord,elident,xi,n);
        N=N';
        dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi,n);
    elseif element == 2
        N = shapefunctions(nelnodes,ncoord,elident,xi,n);
        N=N';
        dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi,n);
    elseif element==3
        N = shapefunctions(nelnodes,ncoord,elident,xi,n);
        dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi,n);
    elseif element==4
        N = shapefunctions(nelnodes,ncoord,elident,xi,n);
        dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi,n);
    end
elseif ndim==3
    w = iWeight_3D(ncoord,nelnodes,npoints,elident);
    xi = iPoints_3D(ncoord,nelnodes,npoints,elident);
    if element==1
        N = shapefunctions_3D(nelnodes,ncoord,elident,xi,n);
        dNdxi = shapefunctionderivs_3D(nelnodes,ncoord,elident,xi,n);
        N=N';
        dNdxi=dNdxi';
    elseif element==2
        N = shapefunctions_3D(nelnodes,ncoord,elident,xi,n);
        dNdxi = shapefunctionderivs_3D(nelnodes,ncoord,elident,xi,n);
    elseif element==3
        N = shapefunctions_3D(nelnodes,ncoord,elident,xi,n);
        dNdxi = shapefunctionderivs_3D(nelnodes,ncoord,elident,xi,n);
    elseif element==4
        N = shapefunctions_3D(nelnodes,ncoord,elident,xi,n);
        dNdxi = shapefunctionderivs_3D(nelnodes,ncoord,elident,xi,n);
    end
end
        pospg = xi;
        wpg = w;


% SYSTEM RESULTING OF DISCRETIZING THE WEAK FORM
if ndim==2
    [K,f] = CreateMatrix(X,T,pospg,wpg,N,dNdxi);
elseif ndim==3
    [K,f] = CreateMatrix_3D(X,T,pospg,wpg,N,dNdxi,n);
end

% BOUNDARY CONDITIONS 
% Boundary conditions are imposed using Lagrange multipliers
  % nodesDir1 = nodes on which solution is u=1
  % nodesDir0 = nodes on which solution is u=0

    if ndim==2
        if element == 1
            nodesDir1=load('BC_inlet_2D_tri_lin.txt');
            nodesDir0=load('BC_outlet_2D_tri_lin.txt');
        elseif element == 2
            nodesDir1=load('BC_inlet_2D_quad_lin.txt');
            nodesDir0=load('BC_outlet_2D_quad_lin.txt');
        elseif element == 3
            nodesDir1=load('BC_inlet_2D_tri_quad.txt');
            nodesDir0=load('BC_outlet_2D_tri_quad.txt');
        elseif element == 4
            nodesDir1=load('BC_inlet_2D_quad_quad.txt');
            nodesDir0=load('BC_outlet_2D_quad_quad.txt');
        end
    elseif ndim==3
        if element == 1
            nodesDir1=load('BC_inlet_3D_tri_lin.txt');
            nodesDir0=load('BC_outlet_3D_tri_lin.txt');
        elseif element == 2
            nodesDir1=load('BC_inlet_3D_quad_lin.txt');
            nodesDir0=load('BC_outlet_3D_quad_lin.txt');
        elseif element == 3
            nodesDir1=load('BC_inlet_3D_tri_quad.txt');
            nodesDir0=load('BC_outlet_3D_tri_quad.txt');
        elseif element == 4
            nodesDir1=load('BC_inlet_3D_quad_quad.txt');
            nodesDir0=load('BC_outlet_3D_quad_quad.txt');
        end
    end

    %Boundary condition matrix
    C = [nodesDir1, ones(length(nodesDir1),1);
         nodesDir0, zeros(length(nodesDir0),1)];


ndir = size(C,1);
neq  = size(f,1);
A = zeros(ndir,neq);
A(:,C(:,1)) = eye(ndir);
b = C(:,2);


% SOLUTION OF THE LINEAR SYSTEM
% Entire matrix
Ktot = [K A';A zeros(ndir,ndir)];
ftot = [f;b];

sol = Ktot\ftot;
Temp = sol(1:neq);
multip = sol(neq+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     POSTPROCESS     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% WRITING HEADING FOR VTK FILE %%%%%%%%%%%%%%%%%%%%%%%%%
addpath('ParaviewFiles')
 nnode=length(X);
 nelem=length(T);
	% printing heading to file
    %f=fopen('MyParaviewFile.vtk','w');
    if ndim==2
        if element == 1
            f=fopen('MyParaviewFile_2D_tri_lin.vtk','w');
        elseif element == 2
            f=fopen('MyParaviewFile_2D_cuad_lin.vtk','w');
        elseif element == 3
            f=fopen('MyParaviewFile_2D_tri_cuad.vtk','w');
        elseif element == 4
            f=fopen('MyParaviewFile_2D_cuad_cuad.vtk','w');
        end
    elseif ndim==3  
        if element == 1
            f=fopen('MyParaviewFile_3D_tri_lin.vtk','w');
        elseif element == 2
            f=fopen('MyParaviewFile_3D_cuad_lin.vtk','w');
        elseif element == 3
            f=fopen('MyParaviewFile_3D_tri_cuad.vtk','w');
        elseif element == 4
            f=fopen('MyParaviewFile_3D_cuad_cuad.vtk','w');
        end
    end
fprintf(f,'# vtk DataFile Version 1.0\n');
fprintf(f,'ECM-CELL DIFFUSION-MECHANICS\n');
fprintf(f,'ASCII\n');
fprintf(f,'\n');
fprintf(f,'DATASET UNSTRUCTURED_GRID\n');
fprintf(f,'%s %8i %s\n','POINTS', nnode ,'float');

%%%%%%%%%%%%%%%%%%%%%% WRITING COORDINATES OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%
	% printing coordinates
if ndim==2
    for i1=1:nnode
        fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:2),0.00000000E+000);
    end
    fprintf(f,'\n');
elseif ndim==3
    for i1=1:nnode
        fprintf(f,'%14.8E %14.8E %14.8E\n',X(i1,1:end));
    end
    fprintf(f,'\n');
end

%%%%%%%%%%%%%%%%%%%%%% WRITING CONNECTIVITY OF NODES %%%%%%%%%%%%%%%%%%%%%%%%%
if ndim==2
    if element==1        
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*4);
        for i1=1:nelem
            fprintf(f,'%4i  %10i  %10i  %10i\n',3,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 5);
        end
        fprintf(f,'\n'); 
    elseif element==2
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*5);
        for i1=1:nelem
            fprintf(f,'%4i  %10i  %10i  %10i  %10i\n',4,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 9);        
        end
        fprintf(f,'\n');  
    elseif element==3
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*7);
        for i1=1:nelem
            fprintf(f,'%4i %10i  %10i  %10i  %10i  %10i  %10i\n',6,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 22);
        end
        fprintf(f,'\n');            

    elseif element==4
        fprintf(f,'%s %8i %8i\n','CELLS', nelem , nelem*9);
        for i1=1:nelem
            fprintf(f,'%4i %10i  %10i %10i %10i %10i  %10i %10i %10i\n',8,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 23);
        end
        fprintf(f,'\n');   
    end
elseif ndim==3
    if element==1        
        fprintf(f,'%s %8i %8i %8i\n','CELLS', nelem , nelem*5);
        for i1=1:nelem
            fprintf(f,'%4i  %10i  %10i  %10i  %10i\n',4,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 10); % 10 de diapo
        end
        fprintf(f,'\n'); 
    elseif element==2
        fprintf(f,'%s %8i %8i %8i\n','CELLS', nelem , nelem*9);
        for i1=1:nelem
            fprintf(f,'%4i %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i\n',8,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
        fprintf(f,' %4i ', 12); % 12 de diapo
        end
        fprintf(f,'\n');            
    elseif element==3
        fprintf(f,'%s %8i %8i %8i\n','CELLS', nelem , nelem*11);
        for i1=1:nelem
            fprintf(f,'%4i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i\n',10,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1,T(i1,9)-1,T(i1,10)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 24); % 24 magico
        end
        fprintf(f,'\n');   
    elseif element==4
        fprintf(f,'%s %8i %8i %8i\n','CELLS', nelem , nelem*21);
        for i1=1:nelem
            fprintf(f,'%4i %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i  %10i\n',20,T(i1,1)-1,T(i1,2)-1,T(i1,3)-1,T(i1,4)-1,T(i1,5)-1,T(i1,6)-1,T(i1,7)-1,T(i1,8)-1,T(i1,9)-1,T(i1,10)-1,T(i1,11)-1,T(i1,12)-1,T(i1,13)-1,T(i1,14)-1,T(i1,15)-1,T(i1,16)-1,T(i1,17)-1,T(i1,18)-1,T(i1,19)-1,T(i1,20)-1);
        end
        fprintf(f,'\n');
        fprintf(f,'%s %8i\n','CELL_TYPES', nelem);
        for i1=1:nelem
            fprintf(f,' %4i ', 25); % 25 magico
        end
        fprintf(f,'\n');   
    end
end

%%%%%%%%%%%%%%%%%%%%%% WRITING VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(f,'%s %8i\n','POINT_DATA', nnode);
fprintf(f,'SCALARS Ux float\n');
fprintf(f,'LOOKUP_TABLE default\n');
for i1=1:nnode
    fprintf(f,'%14.8E\n', Temp(i1) );
end

toc

fclose(f)
