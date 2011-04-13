function [x,u,d,l] = BackwardTridSolve(tu,tc,tl,b)
% BackwardTridSolve - Solves Tridiagonal system Tx=b by UDL factorization
% The UDL factorization is obtained by Backward Gaussian elimination
% withouth pivoting
%
% Syntax:  
%   [x,u,d,l] = BackwardTridSolve(tu,tc,tl,b);
%
% Inputs:
%   tu - Upper diagonal of the Tridiagonal matrix T   
%   tc - Main  diagonal of the Tridiagonal matrix T   
%   tl - Lower diagonal of the Tridiagonal matrix T   
%   b  - Right-side of the equation to solve
% 
% Outputs:
%   u  - Upper diagonal of U in T=UDL
%   d  - Main  diagonal of D in T=UDL
%   l  - Lower diagonal of L in T=UDL
%   x  - Solution of the tridiagonal system.
%
% Example:
%   N=100;
%   tc=randn(1,N);
%   tu=randn(1,N-1); 
%   x=randn(1,N);
%   b=TriSymGaxpy(tc,tu,x);
%   [Fx,Fu,Fd,Fl] = MexForwardTridSolve(tu,tc,tu,b);
%   [Bx,Bu,Bd,Bl] = BackwardTridSolve(tu,tc,tu,b);
%   Fe=sum(abs(Fx-x))
%   Be=sum(abs(Bx-x))
%   FB=sum(abs(Fx-Bx))
%   FBe=sum(abs(0.5*Fx+0.5*Bx-x))
%   
% Implemented in mex built-in function MexBackwardTridSolve.c
%  
% See also: ForwardTridSolve TridSolve TriSymGaxpy

% Author: Roger Pique-Regi
% Ming Hsieh Department of Electrical Engineering, University of Southern California 
% Department of Hemathology/Oncology, Childrens Hospital Los Angeles 
% email: rpique@gmail.com
% Website: http://sipi.usc.edu/~piquereg
% March 2008; 

%   GADA -- Genome Alteration Detection Algorithm 
%   Copyright (C) 2008,2009  Childrens Hospital of Los Angeles
%   Author: Roger Pique-Regi    rpique@gmail.com


N=length(b);
assert(length(tc)==N,'Dimensions do not match');
assert(length(tu)==N-1,'Dimensions do not match');
assert(length(tl)==N-1,'Dimensions do not match');

[x,u,d,l] = MexBackwardTridSolve(tu,tc,tl,b);
