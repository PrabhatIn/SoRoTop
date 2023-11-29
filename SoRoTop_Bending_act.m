function SoRoTop_Bending_act(nelx,nely,volfrac,sefrac,penal,rmin,ks,etaf,betaf,lst,betamax,delrbst,maxit)
%% ___PART 1.____________________________MATERIAL AND FLOW PARAMETERS
E1 = 1;
Emin = E1*1e-6;
nu = 0.30;
[Kv,epsf,r,Dels] = deal(1,1e-7,0.1,2);                       % flow parameters
[Ds, Kvs]= deal((log(r)/Dels)^2*epsf,Kv*(1 - epsf));  % flow parameters
%% ____PART 2._______________FINITE ELEMENT and NON-DESIGN DOMAIN PREPARATION
[nel,nno] = deal(nelx*nely, (nelx+1)*(nely+1));
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
Udofs = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
[Lnode,Rnode]= deal(1:nely+1, (nno-nely):nno);
[Bnode,Tnode]= deal((nely+1):(nely+1):nno, 1:(nely+1):(nno-nely));
[Pdofs,allPdofs,allUdofs] = deal(Udofs(:,2:2:end)/2,1:nno,1:2*nno);
Kp = 1/6*[4 -1 -2 -1;-1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4]; % flow matrix: Darcy Law
KDp = 1/36*[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]; % Drainage matrix
Te = 1/12*[-2 2 1 -1;-2 -1 1 2;-2 2 1 -1;-1 -2 2 1;-1 1 2 -2; -1 -2 2 1; -1 1 2 -2; -2 -1 1 2]; % transformation matrix
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
Ke = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);             %stiffness matrix
iP = reshape(kron(Pdofs,ones(4,1))',16*nel,1);
jP = reshape(kron(Pdofs,ones(1,4))',16*nel,1);
iT = reshape(kron(Udofs,ones(4,1))',32*nel,1);
jT = reshape(kron(Pdofs,ones(1,8))',32*nel,1);
iK = reshape(kron(Udofs,ones(8,1))',64*nel,1);
jK = reshape(kron(Udofs,ones(1,8))',64*nel,1);
IFprj=@(xv,etaf,betaf)((tanh(betaf*etaf) + tanh(betaf*(xv-etaf)))/...            %projection function
    (tanh(betaf*etaf) + tanh(betaf*(1 - etaf))));
dIFprj=@(xv,etaf,betaf) betaf*(1-tanh(betaf*(xv-etaf)).^2)...
    /(tanh(betaf*etaf)+tanh(betaf*(1-etaf)));                    % derivative of the projection function
elNrs = reshape(1:nel,nely,nelx);                                 % element grid
s1 = elNrs(19*nely/20:nely,1:nelx/10);                            % Solid element or element with rho =1
v1 = elNrs(nely/10:nely,9*nelx/10:nelx);
[NDS, NDV ] = deal( s1, v1 );
act = setdiff((1 : nel)', union( NDS, NDV ));
opdof = 2*Tnode(end); % output degree of freedom
%% ____PART 3.__PRESSURE & STRUCTURE B.C's, LOADs, DOFs, Lag. Multi.,sL
[PF, Pin] =deal(0.00001*ones(nno,1),1); %pressure-field preparation
PF([Lnode, Tnode]) = 0; PF([Rnode, Bnode]) = Pin; % applying pressure load
fixedPdofs = allPdofs(PF~=0.00001);
freePdofs  = setdiff(allPdofs,fixedPdofs);
pfixeddofsv = [fixedPdofs' PF(fixedPdofs)]; % p-fixed and its value
fixedUdofs = [2*Lnode(end:-1:19*nely/20+1)-1  2*Lnode(end:-1:19*nely/20+1) 2*Rnode-1]; %fixed displ.
freeUdofs = setdiff(allUdofs,fixedUdofs);
[L,U,lam2] = deal(zeros(2*nno,1));
[lam1,mu1] = deal(zeros(nno,1)); %initialize lambda
[L(opdof)] = 1 ;                  % dummy load and spring constant
%% ___PART 4._________________________________________FILTER PREPARATION
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = H./sum(H,2);                             % matrix of weights (filter)
%% ___PART 5.__________________________MMA OPTIMIZATION PREPARATION & INITIALIZATION
x = zeros(nel,1); % design variable
x(act) = (volfrac*(nel-length(NDV))-length(NDS) )/length(act); x(NDS) = 1;
[nMMA,pMMA,qMMA] = deal(length(act),2,2);
[mMMA,xMMA,xTilde,mvLt] = deal(pMMA+qMMA,x(act),x,0.1);
[xminvec,xmaxvec] = deal(zeros(nMMA,1),ones(nMMA,1)); %Min. & Max
[low, upp] = deal(xminvec,xmaxvec); % Low and Upp limits MMA
[cMMA,dMMA, a0] = deal(1000*ones(mMMA,1),zeros(mMMA,1),1);
aMMA = [ones(pMMA,1); zeros(qMMA,1)];
[xold1,xold2] = deal(xMMA);
[betap,loop, change] = deal(1,0,1);
[etab,etae] =deal(0.5,0.5+delrbst);
costadd = 10000;
%% ____PART 6_____________________________________MMA OPTIMIZATION LOOP
while(loop<maxit && change>0.0001)
    loop = loop + 1;  % Updating the opt. iteration
    %___PART 6.1___________Compute blueprint and eroded physical desing variables
    xTilde(act) = xMMA;  xTilde =Hs'*xTilde; xTilde(NDS)=1;  xTilde(NDV)=0;
    xphysb  =  IFprj(xTilde,etab,betap); xphysb(NDS)=1; xphysb(NDV)=0;
    xphyse = IFprj(xTilde,etae,betap); xphyse(NDS)=1; xphyse(NDV)=0;
    %___PART 6.2__________Performing blueprint design analysis using ObjObjSens function
    [objb, objsensb,volb,volsensb, ~, ~,Ub,Fb,PFb] = ObjObjSens(xphysb,nel,E1,Emin,penal,Kv,Kvs,epsf,Ds,etaf,betaf,Udofs,freeUdofs, ...
        Pdofs,pfixeddofsv,fixedPdofs,freePdofs,iP,jP,iT,jT,iK,jK,Kp,KDp,Te,Ke,opdof,ks,loop,IFprj,dIFprj,L,U,lam1,lam2,mu1,lst,volfrac,sefrac);
    %___PART 6.3___________Performing eroded design analysis using ObjObjSens function
    [obje, objsense,~,~, SEe, SEsense,~,~,~] = ObjObjSens(xphyse,nel,E1, Emin,penal,Kv,Kvs,epsf,Ds,etaf,betaf,Udofs,freeUdofs, ...
        Pdofs,pfixeddofsv,fixedPdofs,freePdofs,iP,jP,iT,jT,iK,jK,Kp,KDp,Te,Ke,opdof,ks,loop,IFprj,dIFprj,L,U,lam1,lam2,mu1,lst,volfrac,sefrac);
    %___PART 6.4_____________________Filtering and projecting objective and constraints sensitivities
    objsensb = Hs'*(objsensb.*dIFprj(xphysb,etab,betap)); % blueprint sensitiivty
    objsense= Hs'*(objsense.*dIFprj(xphyse,etae,betap));
    volsensb = Hs'*(volsensb.*dIFprj(xphysb,etab,betap));
    SEsense = Hs'*(SEsense.*dIFprj(xphyse,etae,betap));
    %___PART 6.5_________________Stacking constraints and their sensitivities
    constr =[volb SEe];
    constrsens = [volsensb SEsense];
    normf = 1;
    %___PART 6.6______________________SETTING and CALLING MMA OPTIMIZATION
    fval = [costadd + objb*normf,costadd + obje*normf,constr]';
    dfdx  = [objsensb(act,:)*normf, objsense(act,:)*normf, constrsens(act,:)]';
    [xminvec, xmaxvec]= deal(max(0, xMMA - mvLt),min(1, xMMA + mvLt));
    [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(mMMA,nMMA,loop,xMMA,xminvec,xmaxvec,xold1,xold2, ...
        0,0,0,fval,dfdx,0*dfdx,low,upp,a0,aMMA,cMMA,dMMA);
    %___PART 6.7____________Updating__________
    [xold2,xold1, xnew]= deal(xold1, xMMA,xmma);
    change = max(abs(xnew-xMMA)); % Calculating change
    xMMA = xnew;
    if(mod(loop,50)==0 && betap<=betamax), betap= betap*2;end % beta updation
    %___PART 6.8____________________________Printing and plotting results
    fprintf(' It.:%5i Obji.:%11.4f  Obje.:%11.4f Voli.:%7.3f ch.:%7.3f\n',loop,fval(1),fval(2),mean(xphysb),change);
    colormap(gray); imagesc(1-reshape(xphysb, nely, nelx));caxis([0 1]);axis equal off;drawnow;
end
%% ______PART 7____plotting results with final pressure field___________
PFP = figure(2); set(PFP,'color','w'); axis equal off, hold on; colormap('gray');
node = [ (1:nno)' reshape(repmat(0:nelx,nely+1,1),nno,1) repmat(0:-1:-nely,1, nelx+1)']; % nodal coordinates
elem(:,1) = (1:nel)'; elem(:,2:5) = Pdofs; % element and connectivity information
X = reshape(node(elem(:,2:5)',2),4,nel); Y = reshape(node(elem(:,2:5)',3),4,nel);
X1 = 2*max(node(:,2))-reshape(node(elem(:,2:5)',2),4,nel); % for x-symmetry
for i = 1:nel,elemP(i) = sum(PFb(elem(i,2:5)))/4/Pin;end
patch(X, Y, 1-xphysb,'EdgeColor','none');caxis([0 1]);
patch(X1, Y, 1-xphysb,'EdgeColor','none');caxis([0 1]);
 for i = 1:nel
     if (xphysb(i)<0.2 && elemP(i)>0.70)
          patch(X(:,i), Y(:,i), 1-elemP(i),'FaceColor',[0 0.8078 0.90],'EdgeColor','none')
          patch(X1(:,i), Y(:,i), 1-elemP(i),'FaceColor',[0 0.8078 0.90],'EdgeColor','none')
     elseif (xphysb(i)<0.2 && elemP(i)<0.70)
         patch(X(:,i), Y(:,i), 1-elemP(i),'FaceColor','w','EdgeColor','none')
          patch(X1(:,i), Y(:,i), 1-elemP(i),'FaceColor','w','EdgeColor','none')
     end
 end
%% ____________ Plotting deformed profile____________________________
DFP = figure(3); set(DFP,'color','w'); axis equal off, hold on;colormap(gray);
xn= node;                             % defomed nodal position
xn(:,2) = node(:,2) + 0.00025*Ub(1:2:end); xn(:,3) = node(:,3) + 0.00025*Ub(2:2:end);
Xn = reshape(xn(elem(:,2:5)',2),4,nel); Yn = reshape(xn(elem(:,[2:5])',3),4,nel);
Xn1 = 2*max(node(:,2))-reshape(xn(elem(:,2:5)',2),4,nel); % for symmetry about x-axis
patch(Xn, Yn,  1-xphysb,'EdgeColor','none');caxis([0 1]);
patch(Xn1, Yn, 1-xphysb,'EdgeColor','none');caxis([0 1]);
for i = 1:nel
     if (xphysb(i)<0.2 && elemP(i)>0.70)
          patch(Xn(:,i), Yn(:,i), 1-elemP(i),'FaceColor',[0 0.8078 0.90],'EdgeColor','none')
          patch(Xn1(:,i), Yn(:,i), 1-elemP(i),'FaceColor',[0 0.8078 0.90],'EdgeColor','none')
     elseif (xphysb(i)<0.2 && elemP(i)<0.70)
         patch(Xn(:,i), Yn(:,i), 1-elemP(i),'FaceColor','w','EdgeColor','none')
          patch(Xn1(:,i), Yn(:,i), 1-elemP(i),'FaceColor','w','EdgeColor','none')
     end
 end
%% ___________PART 8__________________
function[obj,objsens,vol,volsens, SEc, SEsens,U,F,PF] = ObjObjSens(xphys,nel,E1, Emin,penal,Kv,kvs,epsf,Ds,etaf,betaf,Udofs,freeUdofs,Pdofs,pfixeddofsv, ...
    fixedPdofs,freePdofs,iP,jP,iT,jT,iK,jK,Kp,KDp,Te,Ke,opdof,ks,loop,IFprj,dIFprj,L,U,lam1,lam2,mu1,lst,volfrac,sefrac)
% ___PATT 8.1_______SOLVING FLOW BALANCE EQUATION
Kc = Kv*(1-(1-epsf)*IFprj(xphys,etaf,betaf));         %Flow coefficient
Dc = Ds*IFprj(xphys,etaf,betaf);                            %Drainage coefficient
Ae = reshape(Kp(:)*Kc' + KDp(:)*Dc',16*nel,1);    %Elemental flow matrix in vector form
AG = (sparse(iP,jP,Ae)+ sparse(iP,jP,Ae)')/2;         %Global flow matrix
Aff = AG(freePdofs,freePdofs);     %AG for free pressure dofs
dAff_ldl = decomposition(Aff,'ldl'); % Decomposing Aff matrix
PF(freePdofs,1) = dAff_ldl\(-AG(freePdofs,fixedPdofs)*pfixeddofsv(:,2));
PF(pfixeddofsv(:,1),1) = pfixeddofsv(:,2);              % Final P-field
%__PART 8.2_DETERMINING CONSISTENT NODAL LOADS and GLOBAL Disp. Vector
Ts = reshape(Te(:)*ones(1,nel), 32*nel, 1);        %Elemental transformation matrix in vector form
TG = sparse(iT, jT, Ts);                                       %Global transformation matrix
F = -TG*PF;                                                       % Dertmining nodal forces
E = Emin + xphys.^penal*(E1 - Emin);                %Material interpolation
Ks = reshape(Ke(:)*E',64*nel,1);                         %Elemental stiffness matrix in vector form
KG = (sparse(iK,jK,Ks) + sparse(iK,jK,Ks)')/2;    %Global stiffnes matrix
KG(opdof,opdof) = KG(opdof,opdof) + ks;          % adding the workpiece stiffness
dKG_chol = decomposition(KG(freeUdofs,freeUdofs),'chol','lower'); % decomposed freedofs stiffness
U(freeUdofs) = dKG_chol\F(freeUdofs); %Global Disp. Vect.
%__PART 8.3______objective evaluation
obj = L'*U; % maximizing the output deformation
%__PART 8.4__________sensitivity analysis
lam2(freeUdofs) = -dKG_chol\L(freeUdofs);
lam1(freePdofs) = -(lam2(freeUdofs)'*TG(freeUdofs,freePdofs))/dAff_ldl;
objsT1 = (E1 - Emin)*penal*xphys.^(penal - 1).*sum((lam2(Udofs)*Ke).*U(Udofs),2);
dC1k = -dIFprj(xphys,etaf,betaf).* sum((lam1(Pdofs)*(kvs*Kp)) .* PF(Pdofs),2);
dC1d =  dIFprj(xphys,etaf,betaf).* sum((lam1(Pdofs)*(Ds*KDp)) .* PF(Pdofs),2);
objsT2 = dC1k + dC1d;
objsens = (objsT1 + lst*objsT2); % final sensitivities
%__PART 8.5____volume sensitivities
vol = sum(xphys)/(nel*volfrac)-1;
volsens = 1/(volfrac*nel)*ones(nel,1);
%___PART 8.6____Strain energy sensitivities
if(loop==1) , SE_perm = sefrac*(0.5*U'*KG*U); save SE_perm SE_perm; end; load SE_perm
SEc =   0.5*U'*KG*U/SE_perm -1;
SET1 = -0.5*(E1 - Emin)*penal*xphys.^(penal - 1).*sum(([U(Udofs)]*Ke).*[U(Udofs)],2);
mu1(freePdofs) =   (U(freeUdofs)'*TG(freeUdofs,freePdofs))/dAff_ldl;
dSEk = -dIFprj(xphys,etaf,betaf).* sum((mu1(Pdofs)*(kvs*Kp)) .* PF(Pdofs),2);
dSEd =  dIFprj(xphys,etaf,betaf).* sum((mu1(Pdofs)*(Ds*KDp)) .* PF(Pdofs),2);
SET2 = dSEk + dSEd;
SEsens = (SET1 + SET2)/SE_perm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    SoRoTop is written for  pedagogical purposes. A  detailed description can be  %
%    found in the paper:"SoRoTop: a hitchhiker's guide to topology optimization    % 
%    MATLAB code for design-dependent pneumatic-driven soft robots" Optimization   % 
%    and Engineering, 2023.                                                        %
%                                                                                  %
%    Code and its extensions are available  online as supplementary material       %  
%    of the paper and also available at:                                           %
%                                      https://github.com/PrabhatIn/SoRoTop        %
%                                                                                  %
%    Please send your comment to: pkumar@mae.iith.ac.in                            %
%                                                                                  %
%    One may also refer to the following two papers for more detail:               % 
%                                                                                  %
%    1. Kumar P, Frouws JS, Langelaar M (2020) Topology optimization of fluidic    %
%    pressure-loaded structures and compliant mechanisms using the Darcy method.   %
%    Structural and Multidisciplinary Optimization 61(4):1637-1655                 %
%    2. Kumar P, Langelaar M (2021) On topology optimization of design-dependent   % 
%    pressure-loaded three-dimensional structures and compliant mechanisms.        %
%    International Journal for Numerical Methods in Engineering 122(9):2205-2220   %
%    3. P. Kumar (2023) TOPress: a MATLAB implementation for topology optimization %
%    of structures subjected to design‑dependent pressure loads, Structural and    % 
%    Multidisciplinary Optimization, 66(4), 2023                                   %
%                                                                                  %   
%                                                                                  %
%                                                                                  %
%    Disclaimer:                                                                   %
%    The author does not guarantee that the code is free from erros but reserves   %
%    all rights. Further, the author shall not be liable in any event caused by    % 
%    use of the above code and its extensions                                      %
%                                                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%