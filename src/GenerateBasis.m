#!/Applications/Mathematica . app/Contents/MacOS/MathKernel -script

(* ::Package:: *)
args = $CommandLine[[5;;]];
l=ToExpression[args[[1]]]
m=ToExpression[args[[2]]]
Nb=ToExpression[args[[3]]]

Print["l=",l]
Print["m=",m]
Print["Nb=",Nb]



 
 p = 1;
 h = 1;
 
 I1[p1_, p2_] := Beta[p1 + 1, p2 - p1 - 1]*2^(p1 - p2 + 1);
 I2[p1_, p2_] := I1[p1, p2 + 1];

 
 C0[n_] := 8*n^2;
 C1[n_] := 2 (1 + 2*n*(n + 1) - 3 (2 n + 1));
 C2[n_] := 2 (1 - 2 (n + 1));

 
 K[k_, n_] := -Which[n == 0 && m == 0,
     C1[n]*I1[k + n, k + n + 2 p + 2] +
      C2[n]*I1[k + n + 1, k + n + 2 p + 2] -
      l (l + 1) I1[k + n, k + n + 2 p], n > 0 && m == 0,
     C0[n]*I1[k + n - 1, k + n + 2 p + 2] +
      C1[n]*I1[k + n, k + n + 2 p + 2] +
      C2[n]*I1[k + n + 1, k + n + 2 p + 2] -
      l (l + 1) I1[k + n, k + n + 2 p], n > 0 && m > 0,
     C0[n]*I1[k + n - 1, k + n + 2 p + 2] +
      C1[n]*I1[k + n, k + n + 2 p + 2] +
      C2[n]*I1[k + n + 1, k + n + 2 p + 2] -
      l (l + 1) I1[k + n, k + n + 2 p] -
      m^2*I2[k + n - 1, k + n + 2 p]];
 
 DotC[k_, n_] :=
   Which[k == 0 && n == 0,
    K[0, 0] + p K[0, 1] + p K[1, 0] + p^2 K[1, 1], k == 0 && n > 0,
    K[0, n + 1] + p K[1, n + 1], k > 0 && n == 0,
    K[k + 1, 0] + p K[k + 1, 1], k > 0 && n > 0, K[k + 1, n + 1]];
 
 
 (*https://reference.wolfram.com/language/ParallelTools/tutorial/\
 ParallelEvaluation.html*)
  
  Print["Compute Gram Matrix"]
  
 n0 = If[m == 0, 0, 1];
 MatC = Table[DotC[k, n], {k, n0, Nb}, {n, n0, Nb}];

  
 invdet = N[(1/Det[MatC])];
  
  Print[invdet];
  
  Print["Compute Cholesky Decomposition"]
  
  R1 = CholeskyDecomposition[MatC];
  
  Print["Compute Orthogonal coefficients"]
  
  Y1 = Inverse[R1];
  
  fbare[x_, n_] := If[n > 0, (x)^n, 1];
  ft[x_, n_] := If[n == 0, fbare[x, 0] + p fbare[x, 1], fbare[x, n + 1]];

  
  fGS[x_, n_] :=
    Module[{y, z}, Clear[y]; Clear[z];
     y = Sum[Y1[[k + 1 - n0, n + 1 - n0]]*ft[x, k], {k, n0, n}];
      z = Simplify[y]; Return[z]];
  
  Print["Sampling the orthogonal basis function"]
  
Xi[x_] := (1 + h x)/(1 - x);
           
  TabSingle[Nr_] :=
    Module[{fg, nbx, tab}, Clear[fg]; Clear[nbx];
     fg = fGS[x, Nr]/(Xi[x] + h); nbx = Max[Nr*100, 100];
     tab = N[Table[{i/nbx, If[i < nbx, fg /. x -> i/nbx, 0]}, {i,
         0, nbx}]];Print["n="<>ToString[Nr]]; Return[tab]];
 
           
  TabAll[Nmax_] := Table[TabSingle[Nr], {Nr, n0, Nmax}];
  
  nmax = Nb;
  tabAll = TabAll[nmax];
 
  
  Print["Exporting the interpolated orthogonal basis functions"]
  
  Do[Export["test/test_l_" <> ToString[l] <> "_m_" <> ToString[m] <>
     "_n_" <> ToString[n] <> ".hdf5",
    tabAll[[n + 1 - n0]], {"Datasets", "InterpolationTable"}], {n, n0,
     nmax}]

 
  
 
 Exit[]
 
 
