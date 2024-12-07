function D = daempfTimo(d,L,n)
% Skript zum erstellen einer Steifigkeitsmatrix
% Achtung das ganze gilt nur f?r diese Formpolynome
% Formfunktionen(Polynome)
phi_1 = [2 -3 0 1];
phi_2 = [1 -2 1 0];
phi_3 = [-2 3 0 0];
phi_4 = [1 -1 0 0];

% Quadrierte Stammf.
st_int_phi{1,1} = polyint(conv(phi_1,phi_1));
st_int_phi{2,2} = polyint(conv(phi_2,phi_2));
st_int_phi{3,3} = polyint(conv(phi_3,phi_3));
st_int_phi{4,4} = polyint(conv(phi_4,phi_4));

% Relevante Produkte Stammf.
st_int_phi{1,2} = polyint(conv(phi_1,phi_2));
st_int_phi{1,3} = polyint(conv(phi_1,phi_3));
st_int_phi{2,3} = polyint(conv(phi_2,phi_3));
st_int_phi{2,4} = polyint(conv(phi_2,phi_4));
st_int_phi{3,4} = polyint(conv(phi_3,phi_4));
st_int_phi{1,4} = polyint(conv(phi_1,phi_4));

% Von 0 bis 1 integrieren
int_phi = zeros(4,4);
for i = 1:4
    for j = i:4
        int_phi(i,j) = polyval(st_int_phi{i,j},1)-polyval(st_int_phi{i,j},0);
    end
end

h = (L)/(n-1); % Schrittweite
x = zeros(1,n);
for i = 1:n
    x(i) = h*(i-1);
end
D = zeros(2*n);

for i = 1:n
   D(2*i-1,2*i-1) = h*((i<n)*int_phi(1,1)+(i>1)*int_phi(3,3)); 
   D(2*i-1,2*i) = h^2*((i<n)*int_phi(1,2)+(i>1)*int_phi(3,4));
   D(2*i,2*i) = h^3*((i<n)*int_phi(2,2)+(i>1)*int_phi(4,4));
   if (i<n)
   D(2*i,2*i+1) = h^2*int_phi(2,3);
   D(2*i,2*i+2) = h^3*int_phi(2,4);
   D(2*i-1,2*i+1) = h*int_phi(1,3);
   D(2*i-1,2*i+2) = h^2*int_phi(1,4);
   end
   
end

D = d*D;
D = D+D'-diag(diag(D)); % Macht S symmetrisch voll

end
