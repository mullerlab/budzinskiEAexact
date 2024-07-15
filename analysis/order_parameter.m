function r = order_parameter( ang, NN )
% INPUT
% ang - TxN matrix of angular values
% NN - number of nodes in the cv-NN
%
% OUTPUT
% r - order parameter \in [0,1]
%

r = (1/NN) * sum( exp(1i*ang), 2 );
r = abs(r);
