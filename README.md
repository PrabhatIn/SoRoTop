# SoRoTop
## SoRoTop: a hitchhiker’s guide to topology optimization MATLAB code for design-dependent pneumatic-driven soft robots
`SoRoTop.m` provides a MATLAB implementation for topology optimization for design-dependent pneumatic-driven soft robots. Demands for pneumatic-driven soft robots are constantly rising for various applications. However, they are often designed manually due to the lack of systematic methods. Moreover, design-dependent characteristics of pneumatic actuation pose distinctive challenges. SoRoTop is developed for designing such pneumatic-driven soft robots using topology optimization.
## About author
 Prabhat Kumar, Department of Mechanical and Aerospace Engineering, Indian Institute of Technology Hyderabad, India. Please send your comments and suggestions to  pkumar@mae.iith.ac.in or prabhatkumar.rns@gmail.com
## How to use
1. Please see the paper: P. Kumar (2023) [SoRoTop: a hitchhiker’s guide to topology optimization MATLAB code for design-dependent pneumatic-driven soft robots,  Optimization and Engineering, 2023](https://link.springer.com/article/10.1007/s11081-023-09865-1)
2. MMA setting:
(i) SoRoTop uses the MMA written in 1999 and updated in the 2002 version. The mmasub function has the following form
## [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);
 SoRoTop code calls the MMA on line 110 as
## [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(mMMA,nMMA,loop,xMMA,xminvec,xmaxvec,xold1,xold2,0,0,0,fval,dfdx,0*dfdx,low,upp,a0,aMMA,cMMA,dMMA);
(ii) 
With the 2006 version of MMA, one can modify MMA call (line 110) to:
## [xmma,~,~,~,~,~,~,~,~,low,upp]=mmasub(mMMA,nMMA,loop,xMMA,xminvec,xmaxvec,xold1,xold2,0,0,fval,dfdx,low,upp,a0,aMMA,cMMA,dMMA);
## Citation
For citing the paper, please use the following bibtex format:
```
@article{kumar2023sorotop,
  title={{SoRoTop}: a hitchhiker’s guide to topology optimization {MATLAB} code for design-dependent pneumatic-driven soft robots},
  author={Kumar, Prabhat},
  journal={Optimization and Engineering},
  pages={1--35},
  year={2023},
  publisher={Springer}
}
```
