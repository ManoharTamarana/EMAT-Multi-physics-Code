% quasi 3d Receiver 

% total time run code of sh wave propagation

% clc                                               % for clearing display of command window
% clear all                                         % for clearing all variables in workspace 
% close all                                         % for closing all other windows in matlab
            
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

%  initialization for applying  time varying input

st          =      0                          ;  % st = start time
ft          =      40e-6                      ;  % ft = final time
ss          =      10e-9                      ;  % ss = step size for transient analysis
nts         =     (ft-st)/ss                  ;  % nts = number of time steps
nts         =      round(nts)                 ;


endof       =       6*3                       ;  % endof is degrees of freedom in each element
tndof       =       6*nn                      ;  % total number of DOF in entire assembly 
      

ek           =       zeros(endof,endof)       ;    % stiffness matrix for each element
gk           =       sparse(tndof,tndof)      ;    % total global stiffness matrix for total assembly

ec           =       zeros(endof,endof)        ;
gc           =      sparse(tndof,tndof)        ; % global damping matrix for entire assembly


% element force matrix is not needed

gf           =       zeros(tndof,1,nts+1)      ;    % global force column matrix for total assembly
gfc          =       zeros(tndof,1,nts+1)      ; 
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

     elseif (xc<=6e-3 && xc>=2e-3)&&(yc<=7e-3 && yc>=4e-3)   % region of magnet
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

disp('assembly completed')

                                         ;
  
 disp('loading started')
  
 
%  % extract displacements from propagation code

% load('dgdznfv1.mat')
% load('bxnfv1')
% load('bynfv1')
% load('bznfv1')

j=1;
   for i=1:nn    
          u1          =    x(1,i)                ;
          v1          =    y(1,i)                ;
          
if  (u1<=6e-3 && u1>=2e-3)  && (v1<=2e-3 && v1>=1e-3)  
                 
                 

   for n=1:nts+1
        it           = ss*(n-1) ; % it = instant time

           gfc(6*i-5,1,n)   =   ec4*(-bynf(i,1)*dgdznfv(j,1,n))  ;
           gfc(6*i-4,1,n)   =   ec4*(bxnf(i,1)*dgdznfv(j,1,n))  ;
        
                   
   end
       j=j+1;
             end
             
   end          
            
for n=1:nts+1
    
   gf(:,:,n)   =   gfc(:,:,n)*ths*(A/3)  ;
   
   
end

disp('loading ended')   



% % % implimentation of newmark scheme for time varying case
% % 
% disp('integration scheme used is Newmark method')

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
% alpa     =     0.25                      ;
% delt     =     0.5                       ;
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
%disp('LU started')
% %invofgkk =     inv(gkk)                  ;
% [L,U]    =     lu(gkk)                   ; % lu decomposition
% %U=chol(gkk);
% disp('LU ended')
% % gf already has force input whereas gff used for updatation of force
% 
% gff      =     zeros(tndof,1,nts+1)        ;
% 

 
%  
%   disp('started imposing boundary conditions')  
% for i = 1:nn
%  % application of boundary conditions
%     
%          if x(1,i)     ==    0                      % boundary condition A = 0 on left edge of domain
%          
%          gk(6*i-5,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-5)      =    zeros(tndof,1)       ;
%          gk(6*i-5,6*i-5)  =    1                    ;
%          
%          gk(6*i-4,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-4)      =    zeros(tndof,1)       ;
%          gk(6*i-4,6*i-4)  =    1                    ;
%          
%          gk(6*i-3,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-3)      =    zeros(tndof,1)       ;
%          gk(6*i-3,6*i-3)  =    1                    ;
%          
%          gk(6*i-2,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-2)      =    zeros(tndof,1)       ;
%          gk(6*i-2,6*i-2)  =    1                    ;
%          
%          gk(6*i-1,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-1)      =    zeros(tndof,1)       ;
%          gk(6*i-1,6*i-1)  =    1                    ;
%          
%          
%          gk(6*i,:)        =    zeros(1,tndof)       ;
%          gk(:,6*i)        =    zeros(tndof,1)       ;
%          gk(6*i,6*i)      =    1                    ;
%          
%          
%          
%          gc(6*i-5,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-5)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-4,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-4)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-3,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-3)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-2,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-2)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-1,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-1)      =    zeros(tndof,1)       ;
%          
%          gc(6*i,:)        =    zeros(1,tndof)       ;
%          gc(:,6*i)        =    zeros(tndof,1)       ;
%          
%          gkk(i,:)         =    zeros(1,tndof)       ;
%          gkk(:,i)         =    zeros(tndof,1)       ;
%          gkk(i,i)         =    1                    ;
%          
%          gff(i,1,:)       =    0                    ;
%          
%          
%          
%          
%        end
%      
%        if y(1,i)     ==    0                      % boundary condition A = 0 on bottom edge of domain
%         
%          gk(6*i-5,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-5)      =    zeros(tndof,1)       ;
%          gk(6*i-5,6*i-5)  =    1                    ;
%          
%          gk(6*i-4,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-4)      =    zeros(tndof,1)       ;
%          gk(6*i-4,6*i-4)  =    1                    ;
%          
%          gk(6*i-3,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-3)      =    zeros(tndof,1)       ;
%          gk(6*i-3,6*i-3)  =    1                    ;
%          
%          gk(6*i-2,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-2)      =    zeros(tndof,1)       ;
%          gk(6*i-2,6*i-2)  =    1                    ;
%          
%          gk(6*i-1,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-1)      =    zeros(tndof,1)       ;
%          gk(6*i-1,6*i-1)  =    1                    ;
%          
%          
%          gk(6*i,:)        =    zeros(1,tndof)       ;
%          gk(:,6*i)        =    zeros(tndof,1)       ;
%          gk(6*i,6*i)      =    1                    ;
%          
%          gc(6*i-5,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-5)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-4,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-4)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-3,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-3)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-2,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-2)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-1,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-1)      =    zeros(tndof,1)       ;
%          
%          gc(6*i,:)        =    zeros(1,tndof)       ;
%          gc(:,6*i)        =    zeros(tndof,1)       ;
%          
%          gkk(i,:)         =    zeros(1,tndof)       ;
%          gkk(:,i)         =    zeros(tndof,1)       ;
%          gkk(i,i)         =    1                    ;
%          
%          gff(i,1,:)       =    0                    ;
%          
%          
%  
%       
%        end
%        
%        if x(1,i)     ==    8e-3                    % boundary condition A = 0 on right edge of domain
%         
%          gk(6*i-5,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-5)      =    zeros(tndof,1)       ;
%          gk(6*i-5,6*i-5)  =    1                    ;
%          
%          gk(6*i-4,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-4)      =    zeros(tndof,1)       ;
%          gk(6*i-4,6*i-4)  =    1                    ;
%          
%          gk(6*i-3,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-3)      =    zeros(tndof,1)       ;
%          gk(6*i-3,6*i-3)  =    1                    ;
%          
%          gk(6*i-2,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-2)      =    zeros(tndof,1)       ;
%          gk(6*i-2,6*i-2)  =    1                    ;
%          
%          gk(6*i-1,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-1)      =    zeros(tndof,1)       ;
%          gk(6*i-1,6*i-1)  =    1                    ;
%          
%          
%          gk(6*i,:)        =    zeros(1,tndof)       ;
%          gk(:,6*i)        =    zeros(tndof,1)       ;
%          gk(6*i,6*i)      =    1                    ;
%          
%          gc(6*i-5,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-5)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-4,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-4)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-3,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-3)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-2,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-2)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-1,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-1)      =    zeros(tndof,1)       ;
%          
%          gc(6*i,:)        =    zeros(1,tndof)       ;
%          gc(:,6*i)        =    zeros(tndof,1)       ;
%          
%          gkk(i,:)         =    zeros(1,tndof)       ;
%          gkk(:,i)         =    zeros(tndof,1)       ;
%          gkk(i,i)         =    1                    ;
%          
%          gff(i,1,:)       =    0                    ;
%                
%        
%        end 
%        
%        if y(1,i)     ==    8e-3                    % boundary condition A = 0 on top edge of domain 
%         
%          gk(6*i-5,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-5)      =    zeros(tndof,1)       ;
%          gk(6*i-5,6*i-5)  =    1                    ;
%          
%          gk(6*i-4,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-4)      =    zeros(tndof,1)       ;
%          gk(6*i-4,6*i-4)  =    1                    ;
%          
%          gk(6*i-3,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-3)      =    zeros(tndof,1)       ;
%          gk(6*i-3,6*i-3)  =    1                    ;
%          
%          gk(6*i-2,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-2)      =    zeros(tndof,1)       ;
%          gk(6*i-2,6*i-2)  =    1                    ;
%          
%          gk(6*i-1,:)      =    zeros(1,tndof)       ;
%          gk(:,6*i-1)      =    zeros(tndof,1)       ;
%          gk(6*i-1,6*i-1)  =    1                    ;
%          
%          
%          gk(6*i,:)        =    zeros(1,tndof)       ;
%          gk(:,6*i)        =    zeros(tndof,1)       ;
%          gk(6*i,6*i)      =    1                    ;
%          
%          gc(6*i-5,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-5)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-4,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-4)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-3,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-3)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-2,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-2)      =    zeros(tndof,1)       ;
%          
%          gc(6*i-1,:)      =    zeros(1,tndof)       ;
%          gc(:,6*i-1)      =    zeros(tndof,1)       ;
%          
%          gc(6*i,:)        =    zeros(1,tndof)       ;
%          gc(:,6*i)        =    zeros(tndof,1)       ;
%          
%          gkk(i,:)         =    zeros(1,tndof)       ;
%          gkk(:,i)         =    zeros(tndof,1)       ;
%          gkk(i,i)         =    1                    ;
%            
%          gff(i,1,:)       =    0                    ;
%          
%          
%    
%        end
%        
% end
% disp('end of boundary conditions imposing')
%  

 %disp('solution evaluation started')
 
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
%                  
% 
% 
% %                   gff(:,:,n)               =      gf(:,:,n)+gc*(aa1*ga(:,:,n-1)+aa4*dga(:,:,n-1))                     ; % updated newmark gf matrix
% %                   ga(:,:,n)          =      invofgkk*gff(:,:,n)                                                 ; 
% %                  %  ga(:,:,n)         =      U\(transpose(U)\gff)                                                           ; % finding gk inverse using lu decomposition
% %                  dga(:,:,n)        =      dga(:,:,n-1)+(aa7*aa0*(ga(:,:,n)-ga(:,:,n-1)))-aa7*aa2*dga(:,:,n-1) ; % finding first derivative of potential at each step
% %                  
% 
% 
% 
%            end
%            n
%            
% end
% disp('solution calculation done')

%% wilson theta method

theta    =     2                         ;
aa0      =     6/(theta*theta*ss*ss)     ;
aa1      =     3/(theta*ss)              ;
aa2      =     2*aa1                     ;
aa3      =     theta*ss*0.5              ;
aa4      =     aa0/theta                 ;
aa5      =     -aa2/theta                ;
aa6      =     1-(3/theta)               ;
aa7      =     0.5*ss                    ;
aa8      =     ss*ss/6                   ;



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
         
         
         gk(6*i,:)        =    zeros(1,tndof)       ;
         gk(:,6*i)        =    zeros(tndof,1)       ;
         gk(6*i,6*i)      =    1                    ;
         
         
         
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
         
         gc(6*i,:)        =    zeros(1,tndof)       ;
         gc(:,6*i)        =    zeros(tndof,1)       ;
         
         %gkk(i,:)         =    zeros(1,tndof)       ;
        % gkk(:,i)         =    zeros(tndof,1)       ;
        % gkk(i,i)         =    1                    ;
         
         gf(i,1,:)        =    0                    ;
          
         
         
         
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
         
         
         gk(6*i,:)        =    zeros(1,tndof)       ;
         gk(:,6*i)        =    zeros(tndof,1)       ;
         gk(6*i,6*i)      =    1                    ;
         
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
         
         gc(6*i,:)        =    zeros(1,tndof)       ;
         gc(:,6*i)        =    zeros(tndof,1)       ;
         
        % gkk(i,:)         =    zeros(1,tndof)       ;
         %gkk(:,i)         =    zeros(tndof,1)       ;
         %gkk(i,i)         =    1                    ;
         
         gf(i,1,:)        =    0                    ;
         
         
 
      
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
         
         
         gk(6*i,:)        =    zeros(1,tndof)       ;
         gk(:,6*i)        =    zeros(tndof,1)       ;
         gk(6*i,6*i)      =    1                    ;
         
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
         
         gc(6*i,:)        =    zeros(1,tndof)       ;
         gc(:,6*i)        =    zeros(tndof,1)       ;
         
         %gkk(i,:)         =    zeros(1,tndof)       ;
         %gkk(:,i)         =    zeros(tndof,1)       ;
        % gkk(i,i)         =    1                    ;
         
         gf(i,1,:)        =    0                    ;
               
       
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
         
         
         gk(6*i,:)        =    zeros(1,tndof)       ;
         gk(:,6*i)        =    zeros(tndof,1)       ;
         gk(6*i,6*i)      =    1                    ;
         
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
         
         gc(6*i,:)        =    zeros(1,tndof)       ;
         gc(:,6*i)        =    zeros(tndof,1)       ;
         
         %gkk(i,:)         =    zeros(1,tndof)       ;
         %gkk(:,i)         =    zeros(tndof,1)       ;
         %gkk(i,i)         =    1                    ;
           
         gf(i,1,:)        =    0                    ;
         
         
   
       end
       
end
disp('end of boundary conditions imposing')
 
 
 



%gc=zeros(tndof,tndof);
disp('LU started')
gkk      =     gk+aa1*gc                 ;    % updated newmark gk matrix
% invofgkk =     inv(gkk)                  ;
 [L,U]   =     lu(gkk)                   ;    % lu decomposition
%U=chol(gkk); 
 disp('LU ended')
 
 
 disp('solution evaluation started')
 
 for n=1:nts+1
      
        if n==1
             % initial conditions 
                ga(:,:,n)     =    zeros(tndof,1,n)   ;
                %ga(:,:,n)     =    (gk\gf(:,:,n))  ; 
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
dgaxwtnf        =        zeros(nn,1,nts+1)               ;

for n=1:nts+1    
      for i=1:nn
            gaxnf(i,1,n)    =     ga(6*i-5,1,n)       ;
            dgaxwtnf(i,1,n)    =     dga(6*i-5,1,n)       ;
      end
end
% storing y direction MVP seperately

gaynf        =        zeros(nn,1,nts+1)               ;
dgaywtnf        =        zeros(nn,1,nts+1)               ;

for n=1:nts+1    
      for i=1:nn
            gaynf(i,1,n)    =     ga(6*i-4,1,n)       ;
           dgaywtnf(i,1,n)    =     dga(6*i-4,1,n)       ;

      end
end      
      
% storing z direction MVP seperately

gaznf        =        zeros(nn,1,nts+1)               ;
dgazwtnf        =        zeros(nn,1,nts+1)               ;

for n=1:nts+1    
      for i=1:nn
            gaznf(i,1,n)    =     ga(6*i-3,1,n)       ;
            dgazwtnf(i,1,n)    =     dga(6*i-3,1,n)       ;

      end
      
end     
      
 
% storing total MVP seperately

gatnf       =         sqrt(gaxnf.^2+gaynf.^2+gaznf.^2)          ;
gat1nf       =         sqrt(gaxnf.^2+gaynf.^2)          ;
dgat1nf       =       sqrt(dgaxwtnf.^2)*1e-6 ;




% computing voltage for 1mm coil length

 yu=9701;
 zu=9721;
 lnth=1e-3;
 voltg=zeros(21,nts+1);
%  tvoltg=zeros(1,nts+1);
uio=1;
 
     
     for i=1:nn
         if y(1,i)==3e-3 &&( x(1,i)<=3e-3 && x(1,i)>=2e-3)
             for n=1:nts+1
              voltg(uio,n)=dgaxwtnf(i,1,n);
             end
   uio=uio+1;
         end
     end


tvoltg1mm=zeros(1,nts+1);
lloc=2e-3;
uloc=3e-3;
ser=lloc:deltax:uloc;
ser1=size(ser);
ert=9701;

for i=1:ser1(1,2)
    for n=1:nts+1
   
    tvoltg1mm(1,n)=tvoltg1mm(1,n)+voltg(i,n);
end
    ert=ert+1;
   
end



% computing voltage for 4mm coil length
 yu=9701;
 zu=9781;
 lnth=4e-3;
 voltg=zeros(81,nts+1);
%  tvoltg=zeros(1,nts+1);
uio=1;
 
     
     for i=1:nn
         if y(1,i)==3e-3 &&( x(1,i)<=6e-3 && x(1,i)>=2e-3)
             for n=1:nts+1
              voltg(uio,n)=dgaxwtnf(i,1,n);
             end
   uio=uio+1;
         end
     end


tvoltg4mm=zeros(1,nts+1);
lloc=2e-3;
uloc=6e-3;
ser=lloc:deltax:uloc;
ser1=size(ser);
ert=9701;
%ert=9863;
%for n=1:nts+1;
for i=1:ser1(1,2)
    %varb(1,:)=dgat1nf(ert,1,:);
    tvoltg4mm(1,:)=tvoltg4mm(1,:)+voltg(i,:);
    ert=ert+1;
end




% 
% plotting results

% 
% %  tvoltg1=tvoltg*lnth;
%  tttt=0:ss:ft;
%  sz1=size(tttt);
%  tttt1=tttt*3162;
% %r3=tvoltg1;
% %r3=voltg;
% %r3=dgdznf;
% r3=gf;
%  %r3=dgat1nf;
% % r3=gd;
% %r3=gf;
% %r3=dgaywtnf;
% %dgdznf
%  %r3=voltg;
% %r3=tvoltg;
% %j=2455; 
% %j=2485; %  node number required variation with time
% %j=1661;
% %j=3462;
% 
% %j=6*6330-5;
% 
% %j=22850;
% %j=9691;
% j=4911;
%  for i=1:sz1(1,2)
%      oooo(1,i)=r3(6*j-4,1,i);
%   %oooo(1,i)=r3(j,1,i);
% %oooo(1,i)=r3(j,i);
% %oooo(1,i)=r3(1,i);
% end
%  figure()
%  %plot(tttt,r3)
%  plot(tttt,oooo)
%  grid on
% 

% plotting voltage computated 1 mm length coil
%end
lnth=1e-3;
%voltg1=voltg*lnth*1e-6;
tvoltg1=tvoltg1mm*lnth;
  tttt=0:ss:ft;
 sz1=size(tttt);
 tttt1=tttt*3162*1e3;
%r3=voltg;
r3=tvoltg1;


 %r3=dgaxwtnf;
 %r3=gf;
  %r3=dgat1nf;
%r3=dgd;
%j=4911;
%r3=tvoltg1;
%j=24061;
%j=6*6320-5;
%j=6*6561-5;
j=12021;
j=9741;
j=11;
j=9711;
%j=4911;
%j=6*9691-5;
 for i=1:sz1(1,2)
%     oooo(1,i)=r3(6*j-5,1,i);
  %oooo(1,i)=r3(j,1,i);
%oooo(1,i)=r3(j,i);
oooo(1,i)=r3(1,i);
end
 figure(9)
 %plot(tttt,r3)
 plot(tttt,oooo)
grid on




% plotting voltage computated 4 mm length coil
%end
lnth=4e-3;
%voltg1=voltg*lnth*1e-6;
tvoltg2=tvoltg4mm*lnth;
  tttt=0:ss:ft;
 sz1=size(tttt);
 tttt1=tttt*3162*1e3;
%r3=tvoltg1;
r3=tvoltg2;


 %r3=dgaxwtnf;
 %r3=gf;
  %r3=dgat1nf;
%r3=gd;
%j=4911;
%j=24061;
%j=6*6320-5;
% j=6*6550-5;
%j=6*6491-5;
%j=6*5072-5;
%j=4881;
%j=9651;
%j=6*9691-5;
 for i=1:sz1(1,2)
%     oooo(1,i)=r3(6*j-5,1,i);
  %oooo(1,i)=r3(j,1,i);
%oooo(1,i)=r3(j,i);
oooo(1,i)=r3(1,i);
end
 figure(6)
 %plot(tttt,r3)
 plot(tttt,oooo)
grid on

% 
%  tvoltg=zeros(1,nts+1);
%  j=1;
%   n=101;
%   spc=2e-3:deltax:6e-3;
%   rest=zeros(1,81);
%  for i=9701:1:9781
%     rest(1,j)=dgaxwtnf(i,1,n);
%     j=j+1;
%  end
%  figure()
%  plot(spc,rest)
%  
%  
%  
%  
%  
%  
%  
%  
%  
% tvoltg=zeros(1,nts+1);
% lloc=2e-3;
% uloc=3e-3;
% ser=lloc:deltax:uloc;
% ser1=size(ser);
% ert=9701;
% %ert=9863;
% %for n=1:nts+1;
% for i=1:ser1(1,2)
%     %varb(1,:)=dgat1nf(ert,1,:);
%     tvoltg(1,:)=tvoltg(1,:)+voltg(i,:);
%     ert=ert+1;
% end
% 
% %end
% lnth=1e-3;
% %voltg1=voltg*lnth*1e-6;
% tvoltg1=tvoltg1mm*lnth;
%   tttt=0:ss:ft;
%  sz1=size(tttt);
%  tttt1=tttt*3162*1e3;
% %r3=tvoltg1;
% r3=tvoltg1;
% 
% 
%  %r3=dgaxwtnf;
%  %r3=gf;
%   %r3=dgat1nf;
% %r3=gd;
% %j=4911;
% %j=24061;
% %j=6*6320-5;
% % j=6*6550-5;
% %j=6*6491-5;
% %j=6*5072-5;
% %j=4881;
% %j=9651;
% %j=6*9691-5;
%  for i=1:sz1(1,2)
% %     oooo(1,i)=r3(6*j-5,1,i);
%   %oooo(1,i)=r3(j,1,i);
% %oooo(1,i)=r3(j,i);
% oooo(1,i)=r3(1,i);
% end
%  figure(8)
%  %plot(tttt,r3)
%  plot(tttt,oooo)
% grid on
% 
%  
% figure(2)
% % hold on
% % plot(tttt,iluo2)
% % hold on
% % plot(tttt,iluo3)
% hold on
% plot(tttt,iuo5)
%  
% figure(4)
% % hold on
% % plot(tttt,iluo)
% hold on
% plot(tttt,iuo5)











% 
% 
% 
% transpvvoltg1mm=transpose(tvoltg1mm);
% j=3000;
% for i=1:501
% %     dst1(1,i)=transpvvoltg1mm(1,j);
% dst1(1,i)=tvoltg1mm(1,j);
%     j=j+1;
%     
% end
% 
% dst2=fft(dst1)
% %dst3=fft(transpvvoltg1mm);
% 
% 
% % plotting fft
% 
% %voltg1=voltg*lnth*1e-6;
% tvoltg2=tvoltg4mm*lnth;
%    tttt=30e-6:ss:35e-6;
%  tttt=0:ss:ft;
%  sz1=size(tttt);
%  tttt1=tttt*3162*1e3;
% %r3=tvoltg1;
% r41=dst2;
% %r41=transpose(dst3);
% 
%  for i=1:sz1(1,2)
% %     oooo(1,i)=r3(6*j-5,1,i);
%   %oooo(1,i)=r3(j,1,i);
% %oooo(1,i)=r3(j,i);
% oooo51(1,i)=r41(1,i);
% end
%  figure(113)
%  %plot(tttt,r3)
%  %plot(tttt,oooo51)
% grid on
% 
% % 
% %  tvoltg=zeros(1,nts+1);
% semilogy(tttt,oooo51)
% 
% 
% 
% 
tttt3=0:10e-9:40e-6;
hh3=fft(tvoltg1mm);
hh5=abs(hh3);
hh6=fftshift(hh5);
figure(106)
plot(tttt3,hh6)
hold on


hh4=fft(tvoltg4mm);
hh7=abs(hh4);
hh9=fftshift(hh7);
figure(106)
plot(tttt3,hh9)
hold on