%with out defect code for propagation 

% sh wave propagation case
% wave propagation analysis of specimen with lorentz force as input 

% Domain - 60 mm * 1 mm
% element size - 0.05 mm
% time step - 10 nano sec

% clc                                               % for clearing command window
% clear all                                         % for clearing all variables in workspace 
% close all                                         % for closing all other windows in matlab

tic

xp          =       0:0.05e-3:60e-3           ;   % creating coordinates in x direction
yp          =       0:0.05e-3:1e-3            ;   % creating coordinates in y direction
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

ne          =       2*(ab(1,2)-1)*(ab1(1,2)-1);   % number of elements
t           =       zeros(3,ne)               ;

% TRI         =      delaunay(x,y)              ;   % meshing using delaunay function 
% t           =      transpose(TRI)             ;


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

 %figure(1)
% trimesh(TRI,x,y)                                  % plotting mesh
% 
%for i=1:nn
%    str=string(i);
%text(x(i),y(i),  str)
%end

% initialization  for applying  time varying input

st          =      0                          ;   % st = start time
ft          =      40e-6                      ;   % ft = final time
ss          =      10e-9                      ;   % ss = step size for transient analysis
nts         =      (ft-st)/ss                 ;   % nts = number of time steps
nts         =      round(nts)                 ;


endof       =      3                          ;   % endof is degrees of freedom in each element
tndof       =      nn                         ;   % total number of DOF in entire assembly 


% initilization of required variables 

% 3 dimension array initialization is needed because wave propagation displacements changes with time

ek          =       zeros(endof,endof)        ;   % stiffness matrix for each element
gk          =       zeros(tndof,tndof)        ;   % global stiffness matrix for entire assembly

em          =       zeros(endof,endof)        ;   % element mass matrix for each element
gm          =       zeros(tndof,tndof)        ;   % global mass matrix for entire assembly

gd          =       zeros(tndof,1,nts+1)      ;   % global displacement column matrix for total assembly

gf          =       zeros(tndof,1)            ;   % global force  matrix (here it is global flux matrix)
gff         =       zeros(tndof,1,nts+1)      ;   % global force matrix for each time step
   

dgd         =       zeros(tndof,1,nts+1)      ;   % derivative displacement
ddgd        =       zeros(tndof,1,nts+1)      ;   % double derivative displacement

gdxnf       =      zeros(nn,1,nts+1)          ;   % gdxnf = global displacement x direction node form
gdynf       =      zeros(nn,1,nts+1)          ;   % gdxnf = global displacement y direction node form
gdtnf       =      zeros(nn,1,nts+1)          ;   % gdxnf = global displacement resultant node form

gdxgf       =      zeros(ab1(1,2),ab(1,2),nts+1); % gdxgf = global displacement x direction grid form
gdygf       =      zeros(ab1(1,2),ab(1,2),nts+1); % gdygf = global displacement y direction grid form
gdtgf       =      zeros(ab1(1,2),ab(1,2),nts+1); % gdygf = global displacement resultant grid form



E           =      69e9                       ;   % E  = youngs modulus of aluminium
nu          =      0.32                       ;   % nu = poisson ratio of aluminium
 
D           =     (E/(1-(nu*nu)))*[1 nu 0; nu 1 0; 0 0 (1-nu)*0.5]; % D =matrix used in calculation of stress 
den         =     2700                        ;   % den= desnity of aluminium

% consistent mass matrix 
% em          =     M*[2 0 1 0 1 0; 0 2 0 1 0 1; 1 0 2 0 1 0; 0 1 0 2 0 1; 1 0 1 0 2 0; 0 1 0 1 0 2] ;

% evaluation of stiffness matrix for each element
% evaluation of mass matrix for each element
thnss       =       1e-6                      ;  % thickness of specimen in z direction ( quasi 3D )
poi=0;
disp(' assembly started')
for i       =       1:ne
   




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
       
       
%       if (xc<=30.25e-3 && xc>=30e-3)&&(yc<=0.75e-3 && yc>=0.25e-3)      % for region of void
%             
%            ek   =       zeros(endof,endof)        ;
%            em   =       zeros(endof,endof)        ;
%            poi=poi+1;
%       else

         
       ek   =       zeros(endof,endof)        ;
       B    =       zeros(endof,endof)        ;
       
       a    =      [xi yi 1; xj yj 1; xk yk 1];
       A1   =      det(a)/2                   ;   % area of element 
       A    =      abs(A1)                    ;
       bi   =      (yj-yk)/(2*A)              ;
       bj   =      (yk-yi)/(2*A)              ;
       bk   =      (yi-yj)/(2*A)              ;
       ci   =      (xj-xk)/(2*A)              ;
       cj   =      (xk-xi)/(2*A)              ;
       ck   =      (xi-xj)/(2*A)              ;
       
        B   =      [bi*bi+ci*ci  bi*bj+ci*cj bi*bk+ci*ck ; bj*bi+cj*ci bj*bj+cj*cj bj*bk+cj*ck ; bk*bi+ck*ci bk*bj+ck*cj bk*bk+ck*ck ] ;

   

       ek   =      B*A*thnss                  ;   % ek = element stiffness matrix
       
       % lumped mass matrix
       
       em   =    (den*A/(27e9))*(1/3)*[1 0 0;0 1 0;0 0 1]*thnss ;  % em = element mass matrix
       
     
% globalization of stiffness matrix

     
       gk(u,u)    =      gk(u,u)+ek(1,1)           ;
       gk(u,v)    =      gk(u,v)+ek(1,2)           ;
       gk(u,w)    =      gk(u,w)+ek(1,3)           ;
       
       gk(v,u)    =      gk(v,u)+ek(2,1)           ;
       gk(v,v)    =      gk(v,v)+ek(2,2)           ;
       gk(v,w)    =      gk(v,w)+ek(2,3)           ;
       
       gk(w,u)    =      gk(w,u)+ek(3,1)           ;
       gk(w,v)    =      gk(w,v)+ek(3,2)           ;
       gk(w,w)    =      gk(w,w)+ek(3,3)           ;
       
       gm(u,u)    =      gm(u,u)+em(1,1)           ;
       gm(u,v)    =      gm(u,v)+em(1,2)           ;
       gm(u,w)    =      gm(u,w)+em(1,3)           ;
       
       gm(v,u)    =      gm(v,u)+em(2,1)           ;
       gm(v,v)    =      gm(v,v)+em(2,2)           ;
       gm(v,w)    =      gm(v,w)+em(2,3)           ;
       
       gm(w,u)    =      gm(w,u)+em(3,1)           ;
       gm(w,v)    =      gm(w,v)+em(3,2)           ;
       gm(w,w)    =      gm(w,w)+em(3,3)           ;
       
       
     
end
disp('assembly ended')

%load('stg1rv3.mat')

% 
j1     =       1     ;  

disp('loading forcing term started')

    for i=1:nn
        
        u1=x(1,i);
        v1=y(1,i);
        
         if (u1<=4e-3 && u1>=0e-3)&&(v1<=1e-3 && v1>=0e-3)      % for region of aluminium conductor 
              
             for n=1:501

                gff(i,1,n)=stg1r(j1,1,n)*A*thnss*((thnss*A)/3);  % JA/3 term ( ist factor for lorentz force density conversion next factor is for force term )

              end
               j1=j1+1;
         end
      
    end

disp('loading forcing term finished')




% % newmark scheme
%  
% disp('integration scheme used is Newmark method')
% 
% initialization of newmark constants 
% 
% gmma=0.5; % gmma = gamma in newmark method
% bta=0.25; % bta = beta used in newmark method
% ct1=gmma/(bta*ss); % ct1=constant 1
% ct2=gmma/bta; % ct2 = constant 2
% ct3=ss*((gmma/(2*bta))-1); % ct3 = constant 3
% 
% 
% 
% alpa     =     0.75                      ;
% delt     =     0.75                       ;
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
% disp('evaluation of LU decomposition started ')
% 
%  gkk     =     gk+aa0*gm                 ;    % updated newmark gk matrix
%  [L,U]   =     lu(gkk)                   ;    % lu decomposition
%  
% disp('LU decomposition done') 
%  
% 
% 
%     
% 
% 
% disp('solution evaluation started')
% 
%  for n=1:nts+1
%          if n==1
%             % initial conditions 
%                 gd(:,:,n)     =    zeros(tndof,1,n)   ; 
%                 dgd(:,:,n)    =    zeros(tndof,1,n)   ;
%                 ddgd(:,:,n)   =    zeros(tndof,1,n)   ;
%          else
%    
%                 gfff          =    gff(:,:,n)+gm*(aa0*gd(:,:,n-1)+aa2*dgd(:,:,n-1)+aa3*ddgd(:,:,n-1))   ;  % updated newmark gf matrix
%                 gd(:,:,n)     =    U\(L\gfff)                                                           ;  % finding gk inverse using lu decomposition
%                 ddgd(:,:,n)   =    aa0*(gd(:,:,n)-gd(:,:,n-1))-aa2*dgd(:,:,n-1)-aa3*ddgd(:,:,n-1)       ;  % finding velocity at each and every step
%                 dgd(:,:,n)    =    dgd(:,:,n-1)+aa6*ddgd(:,:,n-1)+aa7*ddgd(:,:,n)                       ;  % finding acceleration at each and every step
%  
%          end
% n
%  end
% disp('solution evaluation done')
% toc
% 



% % wilson thetha scheme
% 
% disp('integration scheme used is wilson theta method')
% 
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

disp(' evaluation of LU decompostion started')

 gkk     =     gk+aa0*gm                 ;    % updated newmark gk matrix
 [L,U]   =     lu(gkk)                   ;    % lu decomposition
 
 disp('LU decomposition done')
 
 
 
 disp('started imposing boundary conditions')
for i = 1:nn
%  application of boundary conditions
    
         if x(1,i)     ==   0e-3 && y(1,i) == 0e-3                     % boundary condition A = 0 on left edge of domain
         
         gk(i,:)      =    zeros(1,tndof)       ;
         gk(:,i)      =    zeros(tndof,1)       ;
         gk(i,i)      =    1                    ;
         
         gm(i,:)      =    zeros(1,tndof)       ;
         gm(:,i)      =    zeros(tndof,1)       ;
         
         gkk(i,:)     =    zeros(1,tndof)       ;
         gkk(:,i)     =    zeros(tndof,1)       ;
         gkk(i,i)     =    1                    ;
         
         gff(i,1,:)    =    0                    ;
         
         end
end

disp('imposing boundary conditions done')
 
 disp('solution evaluation started')
 
 for n=1:nts+1
         if n==1
             % initial conditions 
                gd(:,:,n)     =    zeros(tndof,1,n)   ; 
                dgd(:,:,n)    =    zeros(tndof,1,n)   ;
                ddgd(:,:,n)   =    zeros(tndof,1,n)   ;
         else
   
                ngfff(:,:,n)  =    gff(:,:,n-1)+theta*(gff(:,:,n)-gff(:,:,n-1))+gm*(aa0*gd(:,:,n-1)+aa2*dgd(:,:,n-1)+2*ddgd(:,:,n-1))   ;  % updated newmark gf matrix
                ngd(:,:,n)    =    U\(L\ngfff(:,:,n))                                                   ;  % finding gk inverse using lu decomposition
                ddgd(:,:,n)   =    aa4*(ngd(:,:,n)-gd(:,:,n-1))+aa5*dgd(:,:,n-1)+aa6*ddgd(:,:,n-1)      ;
                dgd(:,:,n)    =    dgd(:,:,n-1)+aa7*(ddgd(:,:,n)+ddgd(:,:,n-1))                         ;
                gd(:,:,n)     =    gd(:,:,n-1)+ss*dgd(:,:,n-1)+aa8*(ddgd(:,:,n)+2*ddgd(:,:,n-1))        ;
                
%                 ddgd(:,:,n)   =    aa0*(gd(:,:,n)-gd(:,:,n-1))-aa2*dgd(:,:,n-1)-aa3*ddgd(:,:,n-1)       ;  % finding velocity at each and every step
%                 dgd(:,:,n)    =    dgd(:,:,n-1)+aa6*ddgd(:,:,n-1)+aa7*ddgd(:,:,n)                       ;  % finding acceleration at each and every step
%  
         end
         n

 end
disp('solution evaluation done')

% 


% conversion of x displacement in node form into x displacement in grid
% form
 
for  n=1:nts+1

      for i=1:ab1(1,2)
            gdtgf(i,:,n) =     gd(ab(1,2)*i-(ab(1,2)-1):ab(1,2)*i,1,n)  ;
   
      end
   
end
%pause()



% extract displacements from propagation code
j=1;
 for i=1:nn
          u1          =    x(1,i)                ;
          v1          =    y(1,i)                ;
          
if  (u1<=4e-3 && u1>=0e-3)  && (v1<=1e-3 && v1>=0e-3)  
            
        dgdznfv(j,1,:)=dgd(i,1,:);
        j=j+1;
        
 end

 end
toc

% % plotting of results
% 
% 
% figure(17)
%  for i=1:100:nts+1
%     % imagesc((gdygf(:,:,i)));
%     % imshow(mat2gray(gdtgf(:,:,i)));
%     % colormap(gray);
%     % contourf(X,Y,gdtgf(:,:,i),100)
%     % surf(X,Y,gdtgf(:,:,i))
%     % colormap(gray)
%     % mat2gray(gdtgf(:,:,i));
%     contour(X,Y,gdtgf(:,:,i),10)
%      caxis([min(gdtgf(:)) max(gdtgf(:)) ])
%     pause()
%     disp(i)
%  
%  end
%  
%  
%  
% % plot for variation of displacement with time
% % at particular location
% 
% tttt   =    0:ss:ft      ;
% sz1    =    size(tttt)   ;
% 
% 
% r3     =     gd          ;
% %r3=rst3;
% %r4     =     gff         ;
% j      =     24261-11*1201       ; %  node number required variation with time
% for i=1:sz1(1,2)
%     
%         oooo(1,i)  =   r3(j,1,i)   ;
%         %oooo1(1,i) =   r4(j,1,i)   ;
%         
% end
% figure(1)
% grid on
% plot(tttt,oooo)
% grid on
% hold on
% 
% 
% 
% % plot for variation of displacement along the edge of the plate
% % at particular instant of time
% 
% pt2    =    0e-3:deltax:60e-3  ;
% k      =    24021-11*1201              ;
% n      =   2000               ;  % time step at which the variation along the length is needed 
% 
% %r2     =    gazgf;
% %r2     =    dgazwzgf;
% %r2     =    jexgf;
% r2     =    gd                 ;
% %r2=rst2;
% plf3   =    zeros(1,1201)      ;
% 
% for i=1:1201
%     
%       plf3(1,i)   =     r2(k,1,n)   ;
%       k           =     k+1         ;
%       
% end
% figure(2)
% plot(pt2,plf3);
% grid on
% hold on
% 
% 
% % plot for 
% % 
% % 
% % % Extraction of stage 2 results
% % 
% % j5    =       1                       ;
% % stg2  =       zeros(81*21,1,nts+1)    ;
% % 
% % % n2=n+3500;
% % 
% %   for i=1:nn
% %       
% %           u1          =    x(1,i)                ;
% %           v1          =    y(1,i)                ;
% %            
% %           if (u1<=4e-3 && u1>=0e-3)&&(v1<=1e-3 && v1>=0e-3)
% %     
% %               stg2(j5,1,:)   =    gd(i,1,:)     ; 
% %               j5             =    j5+1          ;
% %  
% %               i 
% %           end
% % 
% %   end
% % 
% %    
% %  % matching of the time of return travel of wave front
% %  
% %  j6    =       1                       ;
% %  stg2m =       zeros(81*21,1,501)      ; % stage 2 modified results matching the time
% % % time step size see
% % 
% % for n=1:501
% %     
% %              n3  = n+500 ;
% % 
% %             stg2m(:,1,n)=    stg2(:,1,n3)           ;
% % 
% % end
% 
% 
% % toc
%  tttt=0:ss:ft;
%  sz1=size(tttt);
%  
% %r3=tvoltg1;
% %r3=voltg;
% %r3=dgdznf;
% r3=gd;
%  %r3=dgat1nf;
% % r3=dgazwtnf;
% %r3=dgdznf;
% %r3=dgaxwtnf;
% %dgdznf
% %r3=voltg;
% %j=2455; 
% %j=2485; %  node number required variation with time
% %j=1661;
% j=12051;
% %j=12111; % 5mm from center of left edge of plate
% % j=12211; % 10mm from center of left edge of plate
% %j=12311; % 15mm from center of left edge of plate
%  %j=12411; % 20mm from center of left edge of plate
%  %j=12511; % 25mm from center of left edge of plate
%  j=12611; % 30mm from center of left edge of plate
%  for i=1:sz1(1,2)
% %     oooo(1,i)=r3(6*j-5,1,i);
%   oooo(1,i)=r3(j,1,i);
% %oooo(1,i)=r3(j,i);
% end
%  figure()
%  plot(tttt,oooo)
%  grid on
%  %plot(tttt,tvoltg)
% % plot(tttt,voltg(2451,:))
%  %plot(tttt,voltg(1661,:))
% %plot(tttt,dgat1nf(2451,1,:))
%   
