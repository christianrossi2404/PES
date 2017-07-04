%================= SHAPE FUNCTIONS ==================================
%
%        Calculates shape functions for various element types
%
function N = shapefunctions(nelnodes,ncoord,elident,xi,n)

N = zeros(n,nelnodes);
%  2D elements
 
xi = xi';
% for i1=1:n
%    Triangular element
if ( nelnodes == 3 ) 
    for i1=1:n
       N(1,1) = xi(1,1);
       N(1,2) = xi(1,2);
       N(1,3) = 1.-xi(1,1)-xi(1,2);               
    end
elseif ( nelnodes == 6 ) 
    for i1=1:n
         xi3 = 1.-xi(i1,1)-xi(i1,2);
       N(i1,1) = (2.*xi(i1,1)-1.)*xi(i1,1);
       N(i1,2) = (2.*xi(i1,2)-1.)*xi(i1,2);
       N(i1,3) = (2.*xi3-1.)*xi3;
       N(i1,4) = 4.*xi(i1,1)*xi(i1,2);
       N(i1,5) = 4.*xi(i1,2)*xi3;
       N(i1,6) = 4.*xi3*xi(i1,1);
    end
%    Rectangular element                
elseif ( nelnodes == 4 ) 
    for i1=1:n
       N(i1,1) = 0.25*(1.-xi(i1,1))*(1.-xi(i1,2));
       N(i1,2) = 0.25*(1.+xi(i1,1))*(1.-xi(i1,2));
       N(i1,3) = 0.25*(1.+xi(i1,1))*(1.+xi(i1,2));
       N(i1,4) = 0.25*(1.-xi(i1,1))*(1.+xi(i1,2));
    end
elseif (nelnodes == 8)
    for i1=1:n
       N(i1,1) = -0.25*(1.-xi(i1,1))*(1.-xi(i1,2))*(1.+xi(i1,1)+xi(i1,2));
       N(i1,2) = 0.25*(1.+xi(i1,1))*(1.-xi(i1,2))*(xi(i1,1)-xi(i1,2)-1.);
       N(i1,3) = 0.25*(1.+xi(i1,1))*(1.+xi(i1,2))*(xi(i1,1)+xi(i1,2)-1.);
       N(i1,4) = 0.25*(1.-xi(i1,1))*(1.+xi(i1,2))*(xi(i1,2)-xi(i1,1)-1.);
        
       N(i1,5) = 0.5*(1.-xi(i1,1)*xi(i1,1))*(1.-xi(i1,2));
       N(i1,6) = 0.5*(1.+xi(i1,1))*(1.-xi(i1,2)*xi(i1,2));
       N(i1,7) = 0.5*(1.-xi(i1,1)*xi(i1,1))*(1.+xi(i1,2));
       N(i1,8) = 0.5*(1.-xi(i1,1))*(1.-xi(i1,2)*xi(i1,2));
    end
end
