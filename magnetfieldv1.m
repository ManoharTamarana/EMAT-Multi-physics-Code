% version1 only magnetic field

% quasi 3d problem trail for simulation of receiver simulation of only
% magnetic field
% 2 ppm magnet with coil case with time variation input
% fine mesh


clc                                               % for clearing display of command window
clear all                                         % for clearing all variables in workspace 
close all                                         % for closing all other windows in matlab
            
tic                                               % for finding the run time of matlab code

xp          =       0:0.05e-3:8e-3             ;   % creating coordinates in x direction
yp          =       0:0.05e-3:8e-3             ;   % creating coordinates in y direction
deltax      =       0.05e-3                    ;
deltay      =       0.05e-3                    ;
ab          =       size(xp)                  ;   
ab1         =       size(yp)                  ;
nn          =       ab(1,2)*ab1(1,2)          ;   % number of nodes
p           =       zeros(2,nn)               ;   % gives the size of p matrix       
[X,Y]       =       meshgrid(xp,yp)           ;   

X1          =       transpose(X)              ;
Y1          =       transpose(Y)              ;

p           =       [X1(:),Y1(:)]             ;

p1          =       transpose(p)              ;

x           =       p1(1,:)                   ;
y           =       p1(2,:)                   ;


ne          =       2*(ab(1,2)-1)*(ab1(1,2)-1);  % number of elements
t           =       zeros(3,ne)               ;  % triangle matrix for storing all the node numbers of each triangle element 

% TRI         =      delaunay(x,y)              ;   % meshing using delaunay function 
% t           =      transpose(TRI)             ;
% 


t(1,1)      =       1                         ; 
t(2,1)      =       2                         ;
t(3,1)      =       ab(1,2)+1                 ;
t(1,2)      =       2                         ;
t(2,2)      =       ab(1,2)+2                 ;
t(3,2)      =       ab(1,2)+1                 ;
for i       =       1:(ne/2)-1
    j       =       ab(1,2)-1                 ;
     if rem(i,j)   ==      0  
        t(:,2*i+1)  =       t(:,2*i-1)+2      ;
        t(:,2*i+2)  =       t(:,2*i)+2        ;
     else
        t(:,2*i+1)  =     t(:,2*i-1)+1        ;
        t(:,2*i+2)  =     t(:,2*i)+1          ;
    end 
end



endof       =       6*3                       ;  % endof is degrees of freedom in each element
tndof       =       6*nn                      ;  % total number of DOF in entire assembly 
      

ek           =       zeros(endof,endof)       ;    % stiffness matrix for each element
gk           =       sparse(tndof,tndof)       ;    % total global stiffness matrix for total assembly

% element damping matrix is not needed
ec           =       zeros(endof,endof)        ;
gc           =       sparse(tndof,tndof)        ; % global damping matrix for entire assembly


% element force matrix is not needed

gf           =       zeros(tndof,1)      ;    % global force column matrix for total assembly

% element potential matrix is not needed

ga           =       zeros(tndof,1)      ;    % global output potential column matrix for total assembly


% 1 - number corresponds to permenant magnet
% 2 - number corresponds to copper coil
% 3 - number corresponds to air medium
% 4 - number corresponds to aluminium conductor

% properties of different media   

% pr = permeability 

pr1         =       0.00000132                ;  % permeability of neodymium magnet
pr2         =       0.00000125                ;  % permeability of copper coil
pr3         =       0.00000125                ;  % permeability of air
pr4         =       0.00000125                ;  % permeability of aluminium conductor


% ec = electrical conductivity

ec1         =       1000000                   ;  % electrical conductivity of neodymium magnet
ec2         =       59600000                  ;  % electrical conductivity of copper coil 
ec3         =       10.^(-15)                 ;  % electrical conductivity of air
ec4         =       35000000                  ;  % electrical conductivity of aluminium conductor


cnst1       =       ec1/12                ;  % used for evaluation of damping matrix in permenant magnet region 
cnst2       =       ec2/12                ;  % used for evaluation of damping matrix in copper coil region
cnst3       =       ec3/12                ;  % used for evaluation of damping matrix in air region
cnst4       =       ec4/12                ;  % used for evaluation of damping matrix in aluminium conductor
  
ths=1e-6;
thnss=ths;
 B1 = [2 0 0 0          0          0          1 0 0 0           0           0           1 0 0 0          0          0           ;
       0 2 0 0          0          0          0 1 0 0           0           0           0 1 0 0          0          0           ;
       0 0 2 0          0          0          0 0 1 0           0           0           0 0 1 0          0          0           ;
       0 0 0 ths*ths/6  0          0          0 0 0 ths*ths/12  0           0           0 0 0 ths*ths/12 0          0           ;
       0 0 0 0          ths*ths/6  0          0 0 0 0           ths*ths/12  0           0 0 0 0          ths*ths/12 0           ;
       0 0 0 0          0          ths*ths/6  0 0 0 0           0           ths*ths/12  0 0 0 0          0          ths*ths/12  ;
       1 0 0 0          0          0          2 0 0 0           0           0           1 0 0 0          0          0           ;
       0 1 0 0          0          0          0 2 0 0           0           0           0 1 0 0          0          0           ;
       0 0 1 0          0          0          0 0 2 0           0           0           0 0 1 0          0          0           ;
       0 0 0 ths*ths/12 0          0          0 0 0 ths*ths/6   0           0           0 0 0 ths*ths/12 0          0           ;
       0 0 0 0          ths*ths/12 0          0 0 0 0           ths*ths/6   0           0 0 0 0          ths*ths/12 0           ;
       0 0 0 0          0          ths*ths/12 0 0 0 0           0           ths*ths/6   0 0 0 0          0          ths*ths/12  ;
       1 0 0 0          0          0          1 0 0 0           0           0           2 0 0 0          0          0           ;
       0 1 0 0          0          0          0 1 0 0           0           0           0 2 0 0          0          0           ;
       0 0 1 0          0          0          0 0 1 0           0           0           0 0 2 0          0          0           ;
       0 0 0 ths*ths/12 0          0          0 0 0 ths*ths/12  0           0           0 0 0 ths*ths/6  0          0           ;
       0 0 0 0          ths*ths/12 0          0 0 0 0           ths*ths/12  0           0 0 0 0          ths*ths/6  0           ;
       0 0 0 0          0          ths*ths/12 0 0 0 0           0           ths*ths/12  0 0 0 0          0          ths*ths/6  ];

   % evaluation of stiffness matrix for each element and evaluation of


disp('assembly started')

for i       =       1:ne

       ek   =       zeros(endof,endof)        ;
       ec   =       zeros(endof,endof)        ;
       B    =       zeros(endof,endof)        ;       
       u    =       t(1,i)                    ;
       v    =       t(2,i)                    ;
       w    =       t(3,i)                    ;
       xi   =       x(1,u)                    ;
       yi   =       y(1,u)                    ;
       xj   =       x(1,v)                    ;
       yj   =       y(1,v)                    ;
       xk   =       x(1,w)                    ;
       yk   =       y(1,w)                    ;
       xc   =      (xi+xj+xk)/3               ;
       yc   =      (yi+yj+yk)/3               ;
       a    =      [xi yi 1; xj yj 1; xk yk 1];  % area formula 
       A1   =      det(a)/2                   ;
       A    =      abs(A1)                    ;
       ai   =      (xi*yk-xk*yi)/(2*A)        ;
       aj   =      (xj*yi-xi*yj)/(2*A)        ;
       ak   =      (xk*yj-xj*yk)/(2*A)        ;
       bi   =      (yj-yk)/(2*A)              ;
       bj   =      (yk-yi)/(2*A)              ;
       bk   =      (yi-yj)/(2*A)              ;
       ci   =      (xj-xk)/(2*A)              ;
       cj   =      (xk-xi)/(2*A)              ;
       ck   =      (xi-xj)/(2*A)              ;

       
       sfn1 =    ai*bi*xc+ci*yc    ;  % shape function N1
       sfn2 =    aj+bj*xc+cj*yc    ;  % shape function N2
       sfn3 =    ak+bk*xc+ck*yc    ;  % shape function N3
         
         %1         2       3                  4 %     5        6                 7          8         9               10                              11                            12                                               13                             14                               15                                               16                              17                                  %18%%%
         
    B =[ ci*ci  -bi*ci      0              cj*ci   -bj*ci      0               ck*ci     -bk*ci       0                0                               0                              0                                                0                              0                               0                                                 0                              0                                     0                                                    ;
        -bi*ci   bi*bi      0             -cj*bi    bj*bi      0              -ck*bi      bk*bi       0                0                               0                              0                                                0                              0                               0                                                 0                              0                                     0                                                    ;
         0       0          bi*bi+ci*ci    0        0          bj*bi+cj*ci     0          0           bk*bi+ck*ci     -bi/3                           -ci/3                           0                                              -bi/3                           -ci/3                            0                                                -bi/3                          -ci/3                                  0                                                    ;
         ci*cj  -bi*cj      0               cj*cj  -bj*cj      0               ck*cj     -bk*cj       0                0                               0                              0                                                 0                             0                               0                                                 0                              0                                     0                                                    ;
        -ci*bj   bi*bj      0             -cj*bj    bj*bj      0              -ck*bj      bk*bj       0                0                               0                              0                                                 0                             0                               0                                                 0                              0                                     0                                                    ;
         0        0         bi*bj+ci*cj     0       0          bj*bj+cj*cj      0         0           bk*bj+ck*cj     -bj/3                           -cj/3                           0                                              -bj/3                           -cj/3                            0                                                -bj/3                           -cj/3                                 0                                                    ;
         ci*ck  -bi*ck       0             cj*ck    -bj*ck     0               ck*ck     -bk*ck       0                0                               0                              0                                                 0                             0                               0                                                  0                             0                                     0                                                    ;
        -ci*bk   bi*bk       0            -cj*bk     bj*bk     0              -ck*bk      bk*bk        0               0                               0                              0                                                 0                             0                               0                                                 0                               0                                    0                                                    ;
         0       0          bi*bk+ci*ck     0         0        bj*bk+cj*ck      0         0           bk*bk+ck*ck     -bk/3                           -ck/3                           0                                               -bk/3                          -ck/3                            0                                                -bk/3                           -ck/3                                 0                                                    ;
         0       0         -bi/3            0         0       -bj/3             0         0          -bk/3             1/6                            -bi*ci*thnss*thnss/12           0                                               (cj*ci*thnss*thnss/12)+(1/12)  -bj*ci*thnss*thnss/12            0                                                (ck*ci*thnss*thnss/12)+(1/12)   -bk*ci*thnss*thnss/12                 0                                                    ; 
         0       0         -ci/3            0         0       -cj/3             0         0          -ck/3            -ci*bi*thnss*thnss/12            (1/6)+(bi*bi*thnss*thnss/12)   0                                              -cj*bi*thnss*thnss/12            (1/12)+(bj*bi*thnss*thnss/12)   0                                                -ck*bi*thnss*thnss/12           (1/12)+(bk*bi*thnss*thnss/12)         0                                                    ;
         0       0           0              0         0        0                0         0           0                0                                0                             (bi*bi*thnss*thnss/12)+(ci*ci*thnss*thnss/12)    0                               0                              (bj*bi*thnss*thnss/12)+(cj*ci*thnss*thnss/12)     0                               0                                   (ck*ci*thnss*thnss/12)+(bk*bi*thnss*thnss/12)         ;
         0       0         -bi/3            0         0       -bj/3             0         0           -bk/3            (ci*cj*thnss*thnss/12)+(1/12)   -bi*cj*thnss*thnss/12          0                                               (cj*cj*thnss*thnss/12)+(1/6)   -bj*cj*thnss*thnss/12            0                                                 (ck*cj*thnss*thnss/12)+(1/12)  -bk*cj*thnss*thnss/12                 0                                                    ;
         0       0         -ci/3            0         0       -cj/3             0         0           -ck/3            -ci*bj*thnss*thnss/12           (bi*bj*thnss*thnss/12)+(1/12)  0                                               -cj*bj*thnss*thnss/12           (bj*bj*thnss*thnss/12)+(1/6)    0                                                 -ck*bj*thnss*thnss/12           (bk*bj*thnss*thnss/12)+(1/12)        0                                                    ;
         0       0          0               0         0        0                0         0            0                0                               0                             (bi*bj*thnss*thnss/12)+(ci*cj*thnss*thnss/12)    0                               0                              (bj*bj*thnss*thnss/12)+(cj*cj*thnss*thnss/12)      0                               0                                   (bk*bj*thnss*thnss/12)+(cj*ck*thnss*thnss/12)        ;
         0       0         -bi/3            0         0       -bj/3             0         0            -bk/3            (ci*ck*thnss*thnss/12)+(1/12)  -bi*ck*thnss*thnss/12          0                                                (cj*ck*thnss*thnss/12)+(1/12)  -bj*ck*thnss*thnss/12           0                                                  (ck*ck*thnss*thnss/12)+(1/6)   -bk*ck*thnss*thnss/12                 0                                                   ;
         0       0         -ci/3            0         0       -cj/3             0         0            -ck/3           -ci*bk*thnss*thnss/12          (bi*bk*thnss*thnss/12)+(1/12)   0                                               -cj*bk*thnss*thnss/12            (bj*bk*thnss*thnss/12)+(1/12)  0                                                 -ck*bk*thnss*thnss/12            (bk*bk*thnss*thnss/12)+(1/6)         0                                                   ;
         0       0          0               0         0        0                0         0             0                0                                 0                         (bi*bk*thnss*thnss/12)+(ci*ck*thnss*thnss/12)     0                               0                              (bj*bk*thnss*thnss/12)+(cj*ck*thnss*thnss/12)       0                              0                                    (bk*bk*thnss*thnss/12)+(ck*ck*thnss*thnss/12)]      ;

if (xc<=6e-3 && xc>=2e-3)&&(yc<=2e-3 && yc>=1e-3)      % for region of aluminium conductor 
              ek   =      (A*B*thnss)/pr4                      ;
              ec   =       A*B1*cnst4*thnss                    ;

     elseif (xc<=6e-3 && xc>=2e-3)&&(yc<=7e-3 && yc>=4e-3)   % region of copper coil
              ek   =      (A*B*thnss)/pr1                      ;
              ec   =       A*B1*cnst1*thnss                    ;
  
     else
              ek   =      (A*B*thnss)/pr3                      ;
              ec   =       A*B1*cnst3*thnss                    ;

     end
     

 
% globalization of stiffness matrix for entire assmebly domain

v1 = [6*u-5 6*u-4 6*u-3 6*u-2 6*u-1 6*u 6*v-5 6*v-4 6*v-3 6*v-2 6*v-1 6*v 6*w-5 6*w-4 6*w-3 6*w-2 6*w-1 6*w  ];



for iv1 = 1:endof
     ii1 = v1(1,iv1);
   
    for iv2 = 1:endof
        ii2=v1(1,iv2);
        gk(ii1,ii2) = gk(ii1,ii2)+ek(iv1,iv2) ;
        gc(ii1,ii2) = gc(ii1,ii2)+ec(iv1,iv2) ;
    end
end
     

end
disp('assembly ended')
 
disp('imposing of boundary conditions started')  
for i = 1:nn
 % application of boundary conditions
    
         if x(1,i)     ==    0                      % boundary condition A = 0 on left edge of domain
         
         gk(6*i-5,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-5)      =    zeros(tndof,1)       ;
         gk(6*i-5,6*i-5)  =    1                    ;
         
         gk(6*i-4,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-4)      =    zeros(tndof,1)       ;
         gk(6*i-4,6*i-4)  =    1                    ;
         
         gk(6*i-3,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-3)      =    zeros(tndof,1)       ;
         gk(6*i-3,6*i-3)  =    1                    ;
         
         gk(6*i-2,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-2)      =    zeros(tndof,1)       ;
         gk(6*i-2,6*i-2)  =    1                    ;
         
         gk(6*i-1,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-1)      =    zeros(tndof,1)       ;
         gk(6*i-1,6*i-1)  =    1                    ;
         
         
         gk(6*i,:)      =    zeros(1,tndof)         ;
         gk(:,6*i)      =    zeros(tndof,1)         ;
         gk(6*i,6*i)  =    1                        ;
         
         
         
         gc(6*i-5,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-5)      =    zeros(tndof,1)       ;
         
         gc(6*i-4,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-4)      =    zeros(tndof,1)       ;
         
         gc(6*i-3,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-3)      =    zeros(tndof,1)       ;
         
         gc(6*i-2,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-2)      =    zeros(tndof,1)       ;
         
         gc(6*i-1,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-1)      =    zeros(tndof,1)       ;
         
         gc(6*i,:)      =    zeros(1,tndof)         ;
         gc(:,6*i)      =    zeros(tndof,1)         ;
         
         
         
         
       end
     
       if y(1,i)     ==    0                      % boundary condition A = 0 on bottom edge of domain
        
        gk(6*i-5,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-5)      =    zeros(tndof,1)       ;
         gk(6*i-5,6*i-5)  =    1                    ;
         
         gk(6*i-4,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-4)      =    zeros(tndof,1)       ;
         gk(6*i-4,6*i-4)  =    1                    ;
         
         gk(6*i-3,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-3)      =    zeros(tndof,1)       ;
         gk(6*i-3,6*i-3)  =    1                    ;
         
         gk(6*i-2,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-2)      =    zeros(tndof,1)       ;
         gk(6*i-2,6*i-2)  =    1                    ;
         
         gk(6*i-1,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-1)      =    zeros(tndof,1)       ;
         gk(6*i-1,6*i-1)  =    1                    ;
         
         
         gk(6*i,:)      =    zeros(1,tndof)         ;
         gk(:,6*i)      =    zeros(tndof,1)         ;
         gk(6*i,6*i)  =    1                        ;
         
          gc(6*i-5,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-5)      =    zeros(tndof,1)       ;
         
         gc(6*i-4,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-4)      =    zeros(tndof,1)       ;
         
         gc(6*i-3,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-3)      =    zeros(tndof,1)       ;
         
         gc(6*i-2,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-2)      =    zeros(tndof,1)       ;
         
         gc(6*i-1,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-1)      =    zeros(tndof,1)       ;
         
         gc(6*i,:)      =    zeros(1,tndof)         ;
         gc(:,6*i)      =    zeros(tndof,1)         ;
         
         
 
      
       end
       
       if x(1,i)     ==    8e-3                    % boundary condition A = 0 on right edge of domain
        
         gk(6*i-5,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-5)      =    zeros(tndof,1)       ;
         gk(6*i-5,6*i-5)  =    1                    ;
         
         gk(6*i-4,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-4)      =    zeros(tndof,1)       ;
         gk(6*i-4,6*i-4)  =    1                    ;
         
         gk(6*i-3,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-3)      =    zeros(tndof,1)       ;
         gk(6*i-3,6*i-3)  =    1                    ;
         
         gk(6*i-2,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-2)      =    zeros(tndof,1)       ;
         gk(6*i-2,6*i-2)  =    1                    ;
         
         gk(6*i-1,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-1)      =    zeros(tndof,1)       ;
         gk(6*i-1,6*i-1)  =    1                    ;
         
         
         gk(6*i,:)      =    zeros(1,tndof)         ;
         gk(:,6*i)      =    zeros(tndof,1)         ;
         gk(6*i,6*i)  =    1                        ;
         
          gc(6*i-5,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-5)      =    zeros(tndof,1)       ;
         
         gc(6*i-4,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-4)      =    zeros(tndof,1)       ;
         
         gc(6*i-3,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-3)      =    zeros(tndof,1)       ;
         
         gc(6*i-2,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-2)      =    zeros(tndof,1)       ;
         
         gc(6*i-1,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-1)      =    zeros(tndof,1)       ;
         
         gc(6*i,:)      =    zeros(1,tndof)         ;
         gc(:,6*i)      =    zeros(tndof,1)         ;
         
         
               
       
       end 
       
       if y(1,i)     ==    8e-3                    % boundary condition A = 0 on top edge of domain 
        
         gk(6*i-5,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-5)      =    zeros(tndof,1)       ;
         gk(6*i-5,6*i-5)  =    1                    ;
         
         gk(6*i-4,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-4)      =    zeros(tndof,1)       ;
         gk(6*i-4,6*i-4)  =    1                    ;
         
         gk(6*i-3,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-3)      =    zeros(tndof,1)       ;
         gk(6*i-3,6*i-3)  =    1                    ;
         
         gk(6*i-2,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-2)      =    zeros(tndof,1)       ;
         gk(6*i-2,6*i-2)  =    1                    ;
         
         gk(6*i-1,:)      =    zeros(1,tndof)       ;
         gk(:,6*i-1)      =    zeros(tndof,1)       ;
         gk(6*i-1,6*i-1)  =    1                    ;
         
         
         gk(6*i,:)      =    zeros(1,tndof)         ;
         gk(:,6*i)      =    zeros(tndof,1)         ;
         gk(6*i,6*i)  =    1                        ;
         
         gc(6*i-5,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-5)      =    zeros(tndof,1)       ;
         
         gc(6*i-4,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-4)      =    zeros(tndof,1)       ;
         
         gc(6*i-3,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-3)      =    zeros(tndof,1)       ;
         
         gc(6*i-2,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-2)      =    zeros(tndof,1)       ;
         
         gc(6*i-1,:)      =    zeros(1,tndof)       ;
         gc(:,6*i-1)      =    zeros(tndof,1)       ;
         
         gc(6*i,:)      =    zeros(1,tndof)         ;
         gc(:,6*i)      =    zeros(tndof,1)         ;
         
         
         
         
   
       end
       
end
disp('imposing of boundary conditions ended')
% input to magnet

disp('loading started')

tgfv        =     795774.715*(A/3)*ths        ; % 1 TESLA input

% for dividing input among nodes
 m          =     61                          ;

 nnfv       =     -tgfv/m                     ; % nnfv = negative node force value
 pnfv       =     -nnfv                      ; % pnfv = positive node force value
   
 gfm        =     zeros(tndof,1)              ;
   
      % giving current density as input to permenant magnet
      
   for i=1:nn    
          u1          =    x(1,i)                ;
          v1          =    y(1,i)                ;
          
       if  u1==2e-3 && (v1<=7e-3 && v1>=4e-3)  
           
          gfm(6*i-3,1)     =     pnfv              ;
          
       end
      
       if  u1==4e-3 && (v1<=7e-3 && v1>=4e-3) 
           
          gfm(6*i-3,1)     =     2*nnfv            ;
          
       end 
      
       if  u1==6e-3 && (v1<=7e-3 && v1>=4e-3)  
          
           gfm(6*i-3,1)     =     pnfv             ;
           
       end 
   end 
   

    
gf=gfm;

disp('loading ended')

 
 
 disp('LU decomposition started')
 
 [L,U]    =     lu(gk)                   ; % lu decomposition
 
 disp('LU decomposition ended')
 
 disp('solution evaluation started')
 
 ga          =      (U\(L\(gf)))  ;

disp('solution calculation done')

toc



% storing x direction MVP seperately

gaxnf        =        zeros(nn,1)               ;
   
      for i=1:nn
            gaxnf(i,1)    =     ga(6*i-5,1)       ;
      end

% storing y direction MVP seperately

gaynf        =        zeros(nn,1)               ;

   
      for i=1:nn
            gaynf(i,1)    =     ga(6*i-4,1)       ;
      end
      
      
% storing z direction MVP seperately

gaznf        =        zeros(nn,1)               ;

   
      for i=1:nn
            gaznf(i,1)    =     ga(6*i-3,1)       ;
      end
      
     
      
 
% storing total MVP seperately

gatnf       =         sqrt(gaxnf.^2+gaynf.^2+gaznf.^2)          ;


% storing derivative of x  direction MVP w r to z direction  seperately

 dgaxwznf        =        zeros(nn,1)               ;


      for i=1:nn
            dgaxwznf(i,1)    =     ga(6*i-2,1)       ;
      end

      
      
 % storing derivative of y  direction MVP w r to z direction  seperately

 dgaywznf        =        zeros(nn,1)               ;

   
      for i=1:nn
            dgaywznf(i,1)    =     ga(6*i-1,1)       ;
      end  
   
      
  % storing derivative of z  direction MVP w r to z direction  seperately

 dgazwznf        =        zeros(nn,1)               ;

  
      for i=1:nn
            dgazwznf(i,1)    =     ga(6*i,1)       ;
      end
  

% storing x direction MVP in grid form

gaxgf       =         zeros(ab1(1,2),ab(1,2))   ;
 

  
     for i=1:ab1(1,2)
     
            gaxgf(i,:)   =       gaxnf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1)   ;
   
     end



% storing y direction MVP in grid form 
gaygf       =         zeros(ab1(1,2),ab(1,2))   ;
 

     
     for i=1:ab1(1,2)
     
            gaygf(i,:)   =       gaynf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1)   ;
   
     end



% storing z direction MVP in grid form 
gazgf       =         zeros(ab1(1,2),ab(1,2))   ;
 

    
     for i=1:ab1(1,2)
     
            gazgf(i,:)   =       gaznf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1)   ;
   
     end

% storing total direction MVP in grid form 

gatgf       =         zeros(ab1(1,2),ab(1,2))   ;
 

     
     for i=1:ab1(1,2)
     
            gatgf(i,:)   =       gatnf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1)   ;
   
     end
     
   
     
% storing x direction MVP derivative w r t z in grid form 
dgaxwzgf       =         zeros(ab1(1,2),ab(1,2))   ;
 

    
     for i=1:ab1(1,2)
     
            dgaxwzgf(i,:)   =       dgaxwznf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1)   ;
   
     end
  
     
     

 % storing y direction MVP derivative w r t z in grid form 
dgaywzgf       =         zeros(ab1(1,2),ab(1,2))   ;
 

    
     for i=1:ab1(1,2)
     
            dgaywzgf(i,:)   =       dgaywznf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1)   ;
   
     end

         
 
 % storing z direction MVP derivative w r t z in grid form
dgazwzgf       =         zeros(ab1(1,2),ab(1,2))   ;
 

     
     for i=1:ab1(1,2)
     
            dgazwzgf(i,:)   =       dgazwznf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1)   ;
   
     end
   

% finding magnetic flux density of permenant magnet ( which does not change with time )
 
 
mx          =     [-1 -1 -1; 0 0 0; 1 1 1]    ;
my          =     [1 0 -1; 1 0 -1; 1 0 -1]    ;

% grid form of variable is achieved only in fine mesh region 


    
dgazwygf      =     conv2(gazgf,my,'same')/deltax  ;
dgazwxgf      =     conv2(gazgf,mx,'same')/deltax  ;
dgaywxgf      =     conv2(gaygf,mx,'same')/deltax ;
dgaxwygf      =     conv2(gaxgf,my,'same')/deltax  ;

    
%cnvsfactor=1e6;% amphere/meeter  to tesla of B(magnetic flux density)
cnvsfactor=1;

bxgf       =    (dgazwygf-dgaywzgf)/cnvsfactor               ;
bygf       =   (-dgazwxgf+dgaxwzgf)/cnvsfactor               ;
bzgf       =    (dgaywxgf-dgaxwygf)/cnvsfactor               ;
%bzgf      =     conv2(frmvpgf,my,'same')/w  ;

rbgf      =     sqrt(bxgf.^2+bygf.^2+bzgf.^2)   ; % node wise resultant magnetic field
thgf      =     atan2(bygf,bxgf)                ; % node wise angle


% process of getting node wise data of magnetic flux density from grid data
% of magnetic flux density



    
tw1         =     transpose(bxgf)           ;
tw2         =     transpose(bygf)           ;
tw3         =     transpose(bzgf)           ;
tw4         =     transpose(rbgf)           ;
tw5         =     transpose(thgf)           ;


    
bxnf        =     tw1(:)                           ;
bynf        =     tw2(:)                           ;
bznf        =     tw3(:)                           ;
rbnf        =     tw4(:)                           ;
thnf        =     tw5(:)                          ;




%   
% plotting results




%      figure(1)
%        poy1        =     polyshape([2e-3; 6e-3; 6e-3; 2e-3],[1e-3; 1e-3; 2e-3; 2e-3])     ;                  % permenant magnet coordinates
%        poy2        =     polyshape([0e-3; 8e-3; 8e-3; 0e-3],[0e-3 ;  0e-3;  8e-3; 8e-3 ])                           ;                 % aluminum conductor coordinates
%       %poy3        =     polyshape([2e-3; 6e-3; 6e-3; 2e-3],[4e-3 ;  4e-3;  7e-3; 7e-3 ])                           ;                 % aluminum conductor coordinates
%       plot(poy1)
%       hold on
%       plot(poy2)
%       hold on
%       r1=bzgf;
% % %     % imagesc((gdygf(:,:,i)));
% %       %imshow(mat2gray(gatgf(:,:,i)));
% %     %  colormap(gray);
%       
%        contour(X,Y,r1,10)
% %       caxis([min(r1(:)) max(r1(:))])
% %       %caxis([-1e4 1e4 ])
% %      disp(i)
% %       %surf(X,Y,gdtgf(:,:,i))
% %       %colormap(gray)
% % %       mat2gray(gdtgf(:,:,i));
% %      % contour(X,Y,gdygf(:,:,i),'100')
% %       %colormap(gray)
% %       %caxis([0 0.004 ])
% %       pause()
% %       clear figure(61)
% %       close all
% 
% %  
% % % variation of various parameters along the edge of the plate
% % 
% % pt2=1e-3:0.1e-3:6e-3;
% % k=1631;
% % %n=1500;
% % n=5;
% % %r2=gaznf;
% % r2=dgazwznf;
% % for i=1:ab1(1,2)-30
% %     
% % plf3(i,1)=r2(k,1,n);
% % k=k+1;
% % end
% % figure(15)
% % plot(pt2,plf3);
% % grid on
% % 
% % 
% % 
% % 
% %   tttt=0:ss:ft;
% % sz1=size(tttt);
% % %r1=gaznf;
% % r1=dgaxwznf;
% % j=1660; %  node number required variation with time
% % for i=1:sz1(1,2)
% %     oooo(1,i)=r1(j,1,i);
% % end
% % figure()
% % plot(tttt,oooo)
% 

