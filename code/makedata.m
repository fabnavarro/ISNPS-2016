function [X, Y, sig] = makedata(n,signal_name)
% Generates a data sample of size n corresponding to experiment.
%
%%%
% This function is based on MakeSinalNewb written by
% Anestis Antoniadis and Jeremie Bigot
% which is based on a code provided to them by Buckheit, Chen, Donoho,
% Johnstone & Scargle.
%%%
%
% 	Input
%		n	        Desired sampling length.
%		signal_name	String: 'Angles','Parabolas','Corner'.
%
%	Output
%		X   sorted X_i.
%		Y   observed data Y_i.
%       sig s_star.
%
%
%   Copyright (c) 2017 Fabien Navarro and Adrien Saumard

X = sort(rand(1,n));
epsilon = randn(1,n);
if strcmp(signal_name,'Angles'),
  sig = ((2*X + 0.5).*(X <= 0.15)) + ...
        ((-12*(X-0.15) + 0.8).*(X > 0.15 & X <= 0.2)) + ...
        0.2*(X > 0.2 & X <= 0.5) + ...
        ((6*(X - 0.5) + 0.2).*(X > 0.5 & X <= 0.6)) + ...
        ((-10*(X - 0.6) + 0.8).*(X > 0.6 & X <= 0.65)) + ...
        ((-0.5*(X - 0.65) + 0.3).*(X > 0.65 & X <= 0.85)) + ...
        ((2*(X - 0.85) + 0.2).*(X > 0.85));
elseif strcmp(signal_name,'Parabolas'),
  pos = [0.1 0.2 0.3 0.35 0.37 0.41 0.43 0.5 0.7 0.9];
  hgt = [(-30) 60 (-30) 500 (-1000) 1000 (-500) 7.5 (-15) 7.5];
  sig = zeros(size(X));
  for j =1:length(pos)
    sig = sig + hgt(j).*((X-pos(j)).^2).*(X > pos(j));
  end
  sig = sig + 0.8;
elseif strcmp(signal_name,'Corner'),
  sig = X;
  sig(X <= 0.5) = 62.387.*10.*X(X <= 0.5).^3.*(1-4.*X(X <= 0.5).^2);
  sig(0.5 < X & X <= 0.8) = 62.387.*3.*(0.125-X(0.5 < X & X <= 0.8).^3).*X(0.5 < X & X <= 0.8).^4;
  sig(X > 0.8) = 62.387.*59.443.*(X(X > 0.8)-1).^3;
  sig=(0.6/range(sig)).*sig+0.6;
else
  disp('Unknown signal names');
  disp('Allowable Names are:')
  disp('Angles'),
  disp('Parabolas'),
  disp('Corner'),
end

sigma_v = abs(cos(10*X))/10;
Y = sig+sigma_v.*epsilon;
end
