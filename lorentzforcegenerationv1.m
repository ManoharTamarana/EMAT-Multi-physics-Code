%  final lorentz force generation code version 1

% quasi 3d problem
% 2 ppm magnet with coil case with time variation input
% permenant magnet field strength used - 0.3 Tesla - 1 tesla
% input excitation of coil current - 50-100 ampheres


clc                                               % for clearing display of command window
clear all                                         % for clearing all variables in workspace 
close all                                         % for closing all other windows in matlab
            
tic                                               % for finding the run time of matlab code

xp          =       0:0.05e-3:8e-3            ;   % creating coordinates in x direction
yp          =       0:0.05e-3:8e-3            ;   % creating coordinates in y direction
deltax      =       0.05e-3                   ;
deltay      =       0.05e-3                   ;
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

%  initialization for applying  time varying input

st          =      0                          ;  % st = start time
ft          =      5e-6                       ;  % ft = final time
ss          =      0.01e-6                    ;  % ss = step size for transient analysis
nts         =     (ft-st)/ss                  ;  % nts = number of time steps
nts         =      round(nts)                 ;


endof       =       6*3                       ;  % endof is degrees of freedom in each element
tndof       =       6*nn                      ;  % total number of DOF in entire assembly 
      

ek           =       zeros(endof,endof)       ;    % stiffness matrix for each element
gk           =       sparse(tndof,tndof)      ;    % total global stiffness matrix for total assembly

% element damping matrix is not needed
ec           =       zeros(endof,endof)        ;
gc           =       sparse(tndof,tndof)        ; % global damping matrix for entire assembly


% element force matrix is not needed

gf           =       zeros(tndof,1,nts+1)      ;    % global force column matrix for total assembly

% element potential matrix is not needed

ga           =       zeros(tndof,1,nts+1)      ;    % global output potential column matrix for total assembly
dga          =       zeros(tndof,1,nts+1)      ; % dga = derivative of global potential vector


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

     elseif (xc<=6e-3 && xc>=2e-3)&&(yc<=7e-3 && yc>=4e-3)   % region of permenant magnet
              ek   =      (A*B*thnss)/pr1                      ;
              ec   =       A*B1*cnst1*thnss                    ;
  
      else                                                    % region of air
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

   disp('assembly completed')
   
   disp('started imposing boundary conditions')
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
disp('end of boundary conditions imposing')




disp('loading started')

% input to the magnet

tgfv        =     795774.715*(A/3)*ths        ; % 1 TESLA input

% for dividing input among nodes
 m          =     61                          ;
nnfv        =     -tgfv/m                     ; % nnfv = negative node force value
pnfv        =     -nnfv                       ; % pnfv = positive node force value
  

 
 
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

   
   
% input to the coil

gfc     =      zeros(tndof,1,nts+1)  ;
thou    =      ss*251                ;


  for i=1:nn    
          u1          =    x(1,i)                ;
          v1          =    y(1,i)                ;
          
if  v1==3e-3 && (u1<=6e-3 && u1>=2e-3)  
             
     for n=1:nts+1

          it            =  ss*(n-1)                                                                        ; % it = instant time
        
          tbc           =      50*exp(-5e11*(it-thou)^2)*cos(2*pi*1e6*(it-thou))*ths;

         gfc(6*i-5,1,n) = + tbc    ;

     end
 
       end
       
  end

      
  
for n=1:nts+1
    
   gf(:,:,n)   =   gfm+gfc(:,:,n)   ;
   
end


% % % implimentation of newmark scheme for time varying case
% % 
% % % finding magnetic vector potential at every node point at each time step
% % 
% % initialization of newmark constants 
% 
% % gmma=0.5; % gmma = gamma in newmark method
% % bta=0.25; % bta = beta used in newmark method
% % ct1=gmma/(bta*ss); % ct1=constant 1
% % ct2=gmma/bta; % ct2 = constant 2
% % ct3=ss*((gmma/(2*bta))-1); % ct3 = constant 3
% 
% alpa     =     0.38                      ;
% delt     =     0.75                      ;
% aa0      =     1/(alpa*ss*ss)            ;
% aa1      =     delt/(alpa*ss)            ;
% aa2      =     1/(alpa*ss)               ;
% aa3      =     (1/2*alpa)-1              ;
% aa4      =     (delt/alpa)-1             ;
% aa5      =     (ss*0.5)*((delt/alpa)-2)  ;
% aa6      =     ss*(1-delt)               ;
% aa7      =     delt*ss                   ;
% 
% 
% 
% gkk      =     gk+aa1*gc                 ; % updated newmark gk matrix
% %invofgkk =     inv(gkk)                  ;

% disp('LU decomposition started')

% [L,U]    =     lu(gkk)                   ; % lu decomposition

% disp('LU decomposition ended')

% %U=chol(gkk);
% 
% % gf already has force input whereas gff used for updatation of force
% 
% gff      =     zeros(tndof,1,nts+1)        ;

%  disp('solution evaluation started')
% 
% for n=1:nts+1
%            if n==1
% %                  ga(:,:,n)       =      zeros(tndof,1,n)       ;
%                  ga(:,:,n)         =      (gk)\gf(:,:,n)         ; % inv(gk)
%                  dga(:,:,n)        =      zeros(tndof,1,n)       ;
%                  ddga(:,:,n)       =      zeros(tndof,1,n)       ;
%            else
%                 
%                 gff(:,:,n)         =      gf(:,:,n)+gc*(aa1*ga(:,:,n-1)+aa4*dga(:,:,n-1)+aa5*ddga(:,:,n-1))   ; % updated newmark gf matrix
%                 ga(:,:,n)          =      (U\(L\(gff(:,:,n))))                                                 ; 
%                 ddga(:,:,n)        =      aa0*(ga(:,:,n)-ga(:,:,n-1))-aa2*dga(:,:,n-1)-aa3*ddga(:,:,n-1)      ;
%                 dga(:,:,n)         =      dga(:,:,n-1)+aa6*ddga(:,:,n-1)+aa7*ddga(:,:,n)                       ; 
%            end
%            n
%            
% end
%  disp('solution evaluation ended')
% toc


%% wilson theta method

theta    =     3                         ;
aa0      =     6/(theta*theta*ss*ss)     ;
aa1      =     3/(theta*ss)              ;
aa2      =     2*aa1                     ;
aa3      =     theta*ss*0.5              ;
aa4      =     aa0/theta                 ;
aa5      =     -aa2/theta                ;
aa6      =     1-(3/theta)               ;
aa7      =     0.5*ss                    ;
aa8      =     ss*ss/6                   ;


gkk      =     gk+aa1*gc                 ;    % updated newmark gk matrix
% invofgkk =     inv(gkk)                  ;

disp('LU decomposition started')

 [L,U]   =     lu(gkk)                   ;    % lu decomposition
 
 disp('LU decomposition ended')
%U=chol(gkk); 
 
 disp('solution evaluation started')
 
 for n=1:nts+1
      
        if n==1
             % initial conditions 
               % ga(:,:,n)     =    zeros(tndof,1,n)   ;
                ga(:,:,n)     =    (gk\gf(:,:,n))  ; 
                dga(:,:,n)    =    zeros(tndof,1,n)   ;
                ddga(:,:,n)   =    zeros(tndof,1,n)   ;
         else
  
                ngfff(:,:,n)  =    gf(:,:,n-1)+theta*(gf(:,:,n)-gf(:,:,n-1))+gc*(aa1*ga(:,:,n-1)+2*dga(:,:,n-1)+aa3*ddga(:,:,n-1))   ;  % updated newmark gf matrix
                nga(:,:,n)    =    (U\(L\(ngfff(:,:,n))))                                                                             ;
                 
                ddga(:,:,n)   =    aa4*(nga(:,:,n)-ga(:,:,n-1))+aa5*dga(:,:,n-1)+aa6*ddga(:,:,n-1)                                   ;
                dga(:,:,n)    =    dga(:,:,n-1)+aa7*(ddga(:,:,n)+ddga(:,:,n-1))                                                      ;
                ga(:,:,n)     =    ga(:,:,n-1)+ss*dga(:,:,n-1)+aa8*(ddga(:,:,n)+2*ddga(:,:,n-1))                                     ;
                
%                 ddgd(:,:,n)   =    aa0*(gd(:,:,n)-gd(:,:,n-1))-aa2*dgd(:,:,n-1)-aa3*ddgd(:,:,n-1)       ;  % finding velocity at each and every step
%                 dgd(:,:,n)    =    dgd(:,:,n-1)+aa6*ddgd(:,:,n-1)+aa7*ddgd(:,:,n)                       ;  % finding acceleration at each and every step
%  
         end
         n
         

end
 
 
disp('solution calculation done')

toc



% storing x direction MVP seperately

gaxnf        =        zeros(nn,1,nts+1)               ;

for n=1:nts+1    
      for i=1:nn
            gaxnf(i,1,n)    =     ga(6*i-5,1,n)       ;
      end
end
% storing y direction MVP seperately

gaynf        =        zeros(nn,1,nts+1)               ;

for n=1:nts+1    
      for i=1:nn
            gaynf(i,1,n)    =     ga(6*i-4,1,n)       ;
      end
end      
      
% storing z direction MVP seperately

gaznf        =        zeros(nn,1,nts+1)               ;

for n=1:nts+1    
      for i=1:nn
            gaznf(i,1,n)    =     ga(6*i-3,1,n)       ;
      end
      
end     
      
 
% storing total MVP seperately

gatnf       =         sqrt(gaxnf.^2+gaynf.^2+gaznf.^2)          ;


% storing derivative of x  direction MVP w r to z direction  seperately

 dgaxwznf        =        zeros(nn,1,nts+1)               ;

 for n=1:nts+1
      for i=1:nn
            dgaxwznf(i,1,n)    =     ga(6*i-2,1,n)       ;
      end
 end
      
      
 % storing derivative of y  direction MVP w r to z direction  seperately

 dgaywznf        =        zeros(nn,1,nts+1)               ;

 for n=1:nts+1   
      for i=1:nn
            dgaywznf(i,1,n)    =     ga(6*i-1,1,n)       ;
      end  
 end     
      
  % storing derivative of z  direction MVP w r to z direction  seperately

 dgazwznf        =        zeros(nn,1,nts+1)               ;

 for n=1:nts+1   
      for i=1:nn
            dgazwznf(i,1,n)    =     ga(6*i,1,n)       ;
      end
 end     

% storing x direction MVP in grid form by each layer for each
% time step

gaxgf       =         zeros(ab1(1,2),ab(1,2),nts+1)   ;
 

 for n=1:nts+1     
     for i=1:ab1(1,2)
     
            gaxgf(i,:,n)   =       gaxnf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)   ;
   
     end
 end 


% storing y direction MVP in grid form by each layer for each
% time step

gaygf       =         zeros(ab1(1,2),ab(1,2),nts+1)   ;
 

 for n=1:nts+1     
     for i=1:ab1(1,2)
     
            gaygf(i,:,n)   =       gaynf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)   ;
   
     end
 end


% storing z direction MVP in grid form by each layer for each
% time step

gazgf       =         zeros(ab1(1,2),ab(1,2),nts+1)   ;
 

 for n=1:nts+1     
     for i=1:ab1(1,2)
     
            gazgf(i,:,n)   =       gaznf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)   ;
   
     end
 end    
     
% storing total direction MVP in grid form by each layer for each
% time step

gatgf       =         zeros(ab1(1,2),ab(1,2),nts+1)   ;
 

for n=1:nts+1      
     for i=1:ab1(1,2)
     
            gatgf(i,:,n)   =       gatnf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)   ;
   
     end
     
end     
     
% storing x direction MVP derivative w r t z in grid form by each layer for each
% time step

dgaxwzgf       =         zeros(ab1(1,2),ab(1,2),nts+1)   ;
 

 for n=1:nts+1     
     for i=1:ab1(1,2)
     
            dgaxwzgf(i,:,n)   =       dgaxwznf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)   ;
   
     end
 end   
     
     

 % storing y direction MVP derivative w r t z in grid form by each layer for each
% time step

dgaywzgf       =         zeros(ab1(1,2),ab(1,2),nts+1)   ;
 

 for n=1:nts+1     
     for i=1:ab1(1,2)
     
            dgaywzgf(i,:,n)   =       dgaywznf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)   ;
   
     end
 end  
         
 
 % storing z direction MVP derivative w r t z in grid form by each layer for each
% time step

dgazwzgf       =         zeros(ab1(1,2),ab(1,2),nts+1)   ;
 

 for n=1:nts+1     
     for i=1:ab1(1,2)
     
            dgazwzgf(i,:,n)   =       dgazwznf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)   ;
   
     end
       
 end    

% finding magnetic flux density of permenant magnet ( which does not change with time )
 
 
mx          =     [-1 -1 -1; 0 0 0; 1 1 1]    ;
my          =     [1 0 -1; 1 0 -1; 1 0 -1]    ;

% grid form of variable is achieved only in fine mesh region 

for n=1:nts+1
    
dgazwygf(:,:,n)      =     conv2(gazgf(:,:,n),my,'same')/deltax  ;
dgazwxgf(:,:,n)      =     conv2(gazgf(:,:,n),mx,'same')/deltax  ;
dgaywxgf(:,:,n)      =     conv2(gaygf(:,:,n),mx,'same')/deltax ;
dgaxwygf(:,:,n)      =     conv2(gaxgf(:,:,n),my,'same')/deltax  ;
end
    
%cnvsfactor=1e6;% amphere/meter  to tesla of B(magnetic flux density)
cnvsfactor=1;

bxgf       =    (dgazwygf-dgaywzgf)/cnvsfactor               ;
bygf       =   (-dgazwxgf+dgaxwzgf)/cnvsfactor               ;
bzgf       =    (dgaywxgf-dgaxwygf)/cnvsfactor               ;
%bzgf      =     conv2(frmvpgf,my,'same')/w  ;

rbgf      =     sqrt(bxgf.^2+bygf.^2+bzgf.^2)   ; % node wise resultant magnetic field
thgf      =     atan2(bygf,bxgf)                ; % node wise angle


% process of getting node wise data of magnetic flux density from grid data
% of magnetic flux density

% 
for n=1:nts+1
    
tw1(:,:,n)         =     transpose(bxgf(:,:,n))           ;
tw2(:,:,n)         =     transpose(bygf(:,:,n))           ;
tw3(:,:,n)         =     transpose(bzgf(:,:,n))           ;
tw4(:,:,n)         =     transpose(rbgf(:,:,n))           ;
tw5(:,:,n)         =     transpose(thgf(:,:,n))           ;
end
for n=1:nts+1
      tw6          =     tw1(:,:,n)                       ;
      tw7          =     tw2(:,:,n)                       ;
      tw8          =     tw3(:,:,n)                       ;
      tw9          =     tw4(:,:,n)                       ;
      tw10         =     tw5(:,:,n)                       ;
        
bxnf(:,:,n)        =     tw6(:)                           ;
bynf(:,:,n)        =     tw7(:)                           ;
bznf(:,:,n)        =     tw8(:)                           ;
rbnf(:,:,n)        =     tw9(:)                           ;
thnf(:,:,n)        =     tw10(:)                          ;

end



%cnvfct2=1e6; % considering thickness for quasi static thickness
cnvfct2=1;
% eddy current calculation

ko        =       0       ;
ki        =       0       ;
kj        =       0       ;
kk        =       0       ;
for n=1:nts+1
    
     for i=1:nn    
          u1          =    x(1,i)             ;
          v1          =    y(1,i)             ;
         if (u1>=2e-3 &&  u1<=6e-3) && (v1<=2e-3 && v1>=1e-3)
                  ec   =     ec4              ;
                  ko   =     ko+1             ;

         elseif (v1==3e-3) && (u1>=2e-3 &&  u1<=6e-3)
                  ec   =     ec2              ;
                  kk   =     kk+1             ;
         elseif   (u1>=2e-3 &&  u1<=6e-3) && (v1<=7e-3 && v1>=4e-3)
                  ec   =     ec1              ;
                  kj   =     kj+1             ;              
         else
                  ec   =     ec3              ;
                  ki   =     ki+1             ;
         end
               
     % calculation of eddy current
     
         jexnf(i,1,n)   =    -ec*dga(6*i-5,1,n)/cnvfct2    ;
         jeynf(i,1,n)   =    -ec*dga(6*i-4,1,n)/cnvfct2    ;
         jeznf(i,1,n)   =    -ec*dga(6*i-3,1,n)/cnvfct2    ;
     end
     
end
% 
% %  jexgf = eddy current x component

   for n=1:nts+1

     for i=1:ab1(1,2)
         
            jexgf(i,:,n)  =   jexnf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n) ;
     end

   end
  
%    %  jeygf = eddy current y component

   for n=1:nts+1

     for i=1:ab1(1,2)
         
            jeygf(i,:,n) =   jeynf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)  ;
     end

   end
   
% %  jezgf = eddy current z component

   for n=1:nts+1

     for i=1:ab1(1,2)
         
            jezgf(i,:,n) =   jeznf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)  ;
     end

   end
   

lfxnf        =   zeros(nn,1,nts+1)                       ;
lfynf        =   zeros(nn,1,nts+1)                       ;
lfznf        =   zeros(nn,1,nts+1)                       ;
 
% %  for n=1:nts+1
% %        for i=1:nn
% %            lfznf(i,1,n)= jenf(3*i-2,1,n)*bynf(i,1);
% %            lfxnf(i,1,n)= jenf(3*i-1,1,n)*bznf(i,1);
% %            lfynf(i,1,n)= jenf(3*i,1,n)*bxnf(i,1);
% %        end
% %  end
% %  
%  
%  
for n=1:nts+1
    
       for i=1:nn
           lfznf(i,1,n)   =   jexnf(i,1,n)*bynf(i,1,n)-jeynf(i,1,n)*bxnf(i,1,n)   ;
           lfxnf(i,1,n)   =   jeynf(i,1,n)*bznf(i,1,n)-jeznf(i,1,n)*bynf(i,1,n)   ;
           lfynf(i,1,n)   =   jeznf(i,1,n)*bxnf(i,1,n)-jexnf(i,1,n)*bznf(i,1,n)   ;
       end
       
end

lftnf=sqrt(lfxnf.*lfxnf+lfynf.*lfynf+lfznf.*lfznf);

 
%  lfx = lorentz force x component

   for n=1:nts+1

     for i=1:ab1(1,2)
         
            lfxgf(i,:,n)  =   lfxnf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n) ;
     end

  end


%  lfy = lorentz force y component

   for n=1:nts+1

     for i=1:ab1(1,2)
         
            lfygf(i,:,n) =   lfynf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)  ;
     end

   end
   
   %  lfz = lorentz force y component

   for n=1:nts+1

     for i=1:ab1(1,2)
         
            lfzgf(i,:,n) =   lfznf(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)  ;
     end

   end
   
   
   toc
%    pause()

  



% extraction of lorentz force of this stage to keep in next stage
j=1;
stg1r=zeros(1701,1,nts+1);
%for n=1:nts+1
    
    for i=1:nn
       u1=x(1,i);
       v1=y(1,i);
       if ( (u1<=6e-3 && u1>=2e-3) && (v1<=2e-3 && v1>=1e-3) )
           stg1r(j,1,:) = lfznf(i,1,:);
           j=j+1;
       end
    end
%end




xp1          =       0:0.05e-3:4e-3            ;   % creating coordinates in x direction
yp1          =       0:0.05e-3:1e-3            ;   % creating coordinates in y direction
deltax      =       0.05e-3                   ;
deltay      =       0.05e-3                   ;
ab2          =       size(xp1)                  ;   
ab3         =       size(yp1)                  ;
%nn          =       ab(1,2)*ab1(1,2)          ;   % number of nodes
%p           =       zeros(2,nn)               ;   % gives the size of p matrix       
[X2,Y2]       =       meshgrid(xp1,yp1)           ;   



 for n=1:nts+1

     for i=1:ab3(1,2)
         
            stg1rgf(i,:,n) =   stg1r(ab2(1,2)*i-(ab2(1,2)-1):ab2(1,2)*i,1,n)  ;
     end

   end

% 
% 
% 
% % plotting results
% 
% 
%  figure(1)
% 
%  %r1=gazgf;
%  r1=dgazwzgf;
%  for i=1:10:nts+1
%      figure(1)
%        poy1        =     polyshape([2e-3; 6e-3; 6e-3; 2e-3],[1e-3; 1e-3; 2e-3; 2e-3])     ;                  % permenant magnet coordinates
% %  %poy2        =     polyshape([0e-3; 8e-3; 8e-3; 0e-3],[0e-3 ;  0e-3;  8e-3; 8e-3 ])                           ;                 % aluminum conductor coordinates
%       poy3        =     polyshape([2e-3; 6e-3; 6e-3; 2e-3],[4e-3 ;  4e-3;  7e-3; 7e-3 ])                           ;                 % aluminum conductor coordinates
%       plot(poy1)
%       hold on
%       plot(poy3)
%       hold on
% %     % imagesc((gdygf(:,:,i)));
%       %imshow(mat2gray(gatgf(:,:,i)));
%     %  colormap(gray);
%     %figure()
%       contour(X,Y,r1(:,:,i),10)
%       caxis([min(r1(:)) max(r1(:))])
%       %caxis([-1e4 1e4 ])
%      disp(i)
%       %surf(X,Y,gdtgf(:,:,i))
%       %colormap(gray)
% %       mat2gray(gdtgf(:,:,i));
%      % contour(X,Y,gdygf(:,:,i),'100')
%       %colormap(gray)
%       %caxis([0 0.004 ])
%       pause(0.5)
%       clear figure(61)
%       close all
%       
%       
%  end
%  
% % variation of various parameters along the edge of the plate
% 
% pt2=2e-3:0.05e-3:6e-3;
% k=6481;
% %n=1500;
% n=5;
% %r2=gaznf;
% r2=dgazwznf;
% for i=1:ab1(1,2)-80
%     
% plf3(i,1)=r2(k,1,n);
% k=k+1;
% end
% figure(20)
% plot(pt2,plf3);
% grid on
% 
% 
% tttt=0:ss:ft;
% sz1=size(tttt);
% r1=gaznf;
% %r1=dgaxwznf;
% j=6521; %  node number required variation with time
% for i=1:sz1(1,2)
%     oooo(1,i)=r1(j,1,i);
% end
% figure()
% plot(tttt,oooo)
% 
% 
% %    % plotting of eddy current
% %    
% %    % pt=0:ss:ft;
% %   pt=0:ss:ft;
% %   r1=gaxnf;
% %    for n=1:nts+1
% %    plf2(1,n)=r1(4091,1,n);
% %  % plf(1,n)=lfynf(1661,1,n);
% %    end
% %   figure(127)
% % plot(pt,plf2);
% % grid on
% %    
% % 
%  
% % % pt=0:ss:ft;
% %   pt=0:ss:ft;
% %   r1=lfznf;
% %    for n=1:nts+1
% %    %plf(1,n)=jexnf(1661,1,n);
% %   plf(1,n)=r1(1661,1,n);
% %    end
% %    
% %   figure(114)
% %   plot(pt,plf);
% %   grid on
% 
% 
% figure(101)
% 
%  %poy1        =     polyshape([2e-3; 6e-3; 6e-3; 2e-3],[1e-3; 1e-3; 2e-3; 2e-3])     ;                  % permenant magnet coordinates
%  %poy2        =     polyshape([0e-3; 8e-3; 8e-3; 0e-3],[0e-3 ;  0e-3;  8e-3; 8e-3 ])                           ;                 % aluminum conductor coordinates
%  %poy3        =     polyshape([2.5e-3; 5.5e-3; 5.5e-3; 2.5e-3],[4e-3 ;  4e-3;  7e-3; 7e-3 ])                           ;                 % aluminum conductor coordinates
%  %plot(poy1)
%  %hold on
%  %plot(poy2)
% % hold on
%  %plot(poy3)
%  %hold on
%  
%  %r1=gaxgf;
%  %r1=dgazwzgf;
%  r1=lfzgf;
%  for i=1:50:nts+1
%      figure(1)
%        poy1        =     polyshape([2e-3; 6e-3; 6e-3; 2e-3],[1e-3; 1e-3; 2e-3; 2e-3])     ;                  % permenant magnet coordinates
% %  %poy2        =     polyshape([0e-3; 8e-3; 8e-3; 0e-3],[0e-3 ;  0e-3;  8e-3; 8e-3 ])                           ;                 % aluminum conductor coordinates
%       %poy3        =     polyshape([2e-3; 6e-3; 6e-3; 2e-3],[4e-3 ;  4e-3;  7e-3; 7e-3 ])                           ;                 % aluminum conductor coordinates
%       plot(poy1)
%       hold on
% %       plot(poy3)
% %       hold on
% %     % imagesc((gdygf(:,:,i)));
%       %imshow(mat2gray(gatgf(:,:,i)));
%     %  colormap(gray);
%     %figure()
%       contour(X,Y,r1(:,:,i),10)
%       caxis([min(r1(:)) max(r1(:))])
%       %caxis([-1e4 1e4 ])
%      disp(i)
%      hold off
%       %surf(X,Y,gdtgf(:,:,i))
%       %colormap(gray)
% %       mat2gray(gdtgf(:,:,i));
%      % contour(X,Y,gdygf(:,:,i),'100')
%       %colormap(gray)
%       %caxis([0 0.004 ])
%       pause()
%       clear figure(61)
%       %close all
%       
%       
%  end
%  


% variation of various parameters along the edge of the plate

pt2=2e-3:0.05e-3:6e-3;
%k=6481;
k=6521
%n=1500;
n=350;
r2=lfznf;
%r2=dgazwznf;
for i=1:ab1(1,2)-80
    
plf3(i,1)=r2(k,1,n);
k=k+1;
end
figure(15)
plot(pt2,plf3);
grid on




tttt=0:ss:ft;
sz1=size(tttt);
r3=lfznf;
%r3=gfc;
oooo1=zeros(1,501);
%r3=dgaxwznf;
%j=6*9741-5;
j=6521; %  node number required variation with time
j=6521-2*161;
%j=6541;
%j=3301;
for i=1:sz1(1,2)
    oooo(1,i)=r3(j,1,i);
end
figure()
plot(tttt,oooo)
hold on
%plot(tttt,oooo1)





r1=stg1gf;
 for i=1:50:nts+1
     figure(1)
%        poy1        =     polyshape([0e-3; 4e-3; 4e-3; 0e-3],[0e-3; 0e-3; 1e-3; 1e-3])     ;                  % permenant magnet coordinates
% %  %poy2        =     polyshape([0e-3; 8e-3; 8e-3; 0e-3],[0e-3 ;  0e-3;  8e-3; 8e-3 ])                           ;                 % aluminum conductor coordinates
%       %poy3        =     polyshape([2e-3; 6e-3; 6e-3; 2e-3],[4e-3 ;  4e-3;  7e-3; 7e-3 ])                           ;                 % aluminum conductor coordinates
%       plot(poy1)
%       hold on
%       plot(poy3)
%       hold on
%     % imagesc((gdygf(:,:,i)));
      %imshow(mat2gray(gatgf(:,:,i)));
    %  colormap(gray);
    %figure()
   % caxis([-6.9e14 4.9e13 ])
     %caxis([-1e11 7e11 ])
    %caxis([min(r1(:)) max(r1(:))])
      contour(X2,Y2,r1(:,:,i),10)
      %caxis([min(r1(:)) max(r1(:))])
%        caxis([-8.6e10 6.8e11 ])
%caxis([-6.9e14 4.9e13 ])
%caxis([-1e5 1e10 ])
     disp(i)
     hold off
      %surf(X,Y,gdtgf(:,:,i))
      %colormap(gray)
%       mat2gray(gdtgf(:,:,i));
     % contour(X,Y,gdygf(:,:,i),'100')
      %colormap(gray)
      %caxis([0 0.004 ])
      pause()
     % clear figure(1)
      %close all
      
      
 end




% 
% 
% 
% % 
% % % code for extraction of this stage results to keep in next stage
% % j=1;
% % stg1=zeros(81,1,nts+1);
% % 
% % 
% %     
% %     for i=1:nn
% %        u1=x(1,i);
% %        v1=y(1,i);
% %        if ( (u1<=6e-3 && u1>=2e-3) && (v1<=2e-3 && v1>=1e-3) )
% %            stg1(j,1,:) = lfznf(i,1,:);
% %            j=j+1;
% %            disp('kate')
% %        end
% %     end
% 
% 
%            