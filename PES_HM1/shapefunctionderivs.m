

%================= SHAPE FUNCTION DERIVATIVES ======================
%
function dNdxi = shapefunctionderivs(nelnodes,ncoord,elident,xi,n)
xi=xi';
  dNdxi = zeros(ncoord*n,nelnodes);
%
% 1D elements
%
%   if (ncoord == 1) 
%     if (nelnodes==2) 
%       dNdxi(1,1) = 0.5;
%       dNdxi(2,1) = -0.5;
%     elseif (nelnodes == 3) 
%       dNdxi(1,1) = -0.5+xi(1);
%       dNdxi(2,1) =  0.5+xi(1);
%       dNdxi(3,1) = -2.*xi(1);
%     end
% %
% %  2D elements
% %

% %    Triangular element

if ( nelnodes == 3 )
     for i1=1:n
       dNdxi(i1*2-1,1) = 1.;
       dNdxi(i1*2,2) = 1.;
       dNdxi(i1*2-1,3) = -1.;
       dNdxi(i1*2,3) = -1.;  
     end
elseif ( nelnodes == 6 ) 
     for i1=1:n
         xi3 = 1.-xi(i1,1)-xi(i1,2); 
       dNdxi(i1*2-1,1) = 4.*xi(i1,1)-1.;
       dNdxi(i1*2,2) = 4.*xi(i1,2)-1.;
       dNdxi(i1*2-1,3) = -(4.*xi3-1.);
       dNdxi(i1*2,3) = -(4.*xi3-1.);
       dNdxi(i1*2-1,4) = 4.*xi(i1,2) ;
       dNdxi(i1*2,4) = 4.*xi(i1,1) ;
       dNdxi(i1*2-1,5) =  - 4.*xi(i1,2);
       dNdxi(i1*2,5) =  4.*xi3 - 4.*xi(i1,2);
       dNdxi(i1*2-1,6) = 4.*xi3 - 4.*xi(i1,1);
       dNdxi(i1*2,6) = -4.*xi(i1,1);        
     end
 
%    Rectangular element

elseif ( nelnodes == 4 ) 
     for i1=1:n
       dNdxi(i1*2-1,1) = -0.25*(1.-xi(i1,2));
       dNdxi(i1*2,1) = -0.25*(1.-xi(i1,1));
       dNdxi(i1*2-1,2) = 0.25*(1.-xi(i1,2));
       dNdxi(i1*2,2) = -0.25*(1.+xi(i1,1));
       dNdxi(i1*2-1,3) = 0.25*(1.+xi(i1,2));
       dNdxi(i1*2,3) = 0.25*(1.+xi(i1,1));
       dNdxi(i1*2-1,4) = -0.25*(1.+xi(i1,2));
       dNdxi(i1*2,4) = 0.25*(1.-xi(i1,1));
     end
elseif (nelnodes == 8) 
     for i1=1:n
       dNdxi(i1*2-1,1) = 0.25*(1.-xi(i1,2))*(2.*xi(i1,1)+xi(i1,2));
       dNdxi(i1*2,1) = 0.25*(1.-xi(i1,1))*(xi(i1,1)+2.*xi(i1,2));
       dNdxi(i1*2-1,2) = 0.25*(1.-xi(i1,2))*(2.*xi(i1,1)-xi(i1,2));
       dNdxi(i1*2,2) = 0.25*(1.+xi(i1,1))*(2.*xi(i1,2)-xi(i1,1));
       dNdxi(i1*2-1,3) = 0.25*(1.+xi(i1,2))*(2.*xi(i1,1)+xi(i1,2));
       dNdxi(i1*2,3) = 0.25*(1.+xi(i1,1))*(2.*xi(i1,2)+xi(i1,1));
       dNdxi(i1*2-1,4) = 0.25*(1.+xi(i1,2))*(2.*xi(i1,1)-xi(i1,2));
       dNdxi(i1*2,4) =  0.25*(1.-xi(i1,1))*(2.*xi(i1,2)-xi(i1,1));
       dNdxi(i1*2-1,5) = -xi(i1,1)*(1.-xi(i1,2));
       dNdxi(i1*2,5) = -0.5*(1.-xi(i1,1)*xi(i1,1));
       dNdxi(i1*2-1,6) = 0.5*(1.-xi(i1,2)*xi(i1,2));
       dNdxi(i1*2,6) = -(1.+xi(i1,1))*xi(i1,2);
       dNdxi(i1*2-1,7) = -xi(i1,1)*(1.+xi(i1,2));
       dNdxi(i1*2,7) = 0.5*(1.-xi(i1,1)*xi(i1,1));
       dNdxi(i1*2-1,8) = -0.5*(1.-xi(i1,2)*xi(i1,2));
       dNdxi(i1*2,8) = -(1.-xi(i1,1))*xi(i1,2);
     end
end



