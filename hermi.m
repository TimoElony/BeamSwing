function Y = hermi(X,w_k,h)
%Hermi nimmt einen Vektor mit 2n Lösungskoeffizienten aus FEM
% und bestimmt mit deren Hilfe alle Werte Y an den Stellen X
i_Stuetz = floor(X./h+1); % Bestimmt Nummer der Stützstellen unter gesuchtem x
x_i = (i_Stuetz-1)*h; % x-Position dieser Stützstellen
Ksi = (X-x_i)/h; % Abstand von diesen Stützstellen als relatives X
% Formfunktionen(Polynome)
 phi_1 = [2 -3 0 1];
 phi_2 = [1 -2 1 0];
 phi_3 = [-2 3 0 0];
 phi_4 = [1 -1 0 0];
Y = zeros(size(X));
for k = 1:size(X,2)
    Y(k) = w_k(2*i_Stuetz(k)-1)*polyval(phi_1,Ksi(k))...
    +      w_k(2*i_Stuetz(k))*h*polyval(phi_2,Ksi(k));
    if i_Stuetz(k)<size(w_k,1)/2
    Y(k) = Y(k) +w_k(2*i_Stuetz(k)+1)*polyval(phi_3,Ksi(k))...
                +w_k(2*i_Stuetz(k)+2)*h*polyval(phi_4,Ksi(k));
    end
end
end

