function S = stiffTimo(EI,L,n)
% Skript zum erstellen einer Steifigkeitsmatrix
% Achtung das ganze gilt nur für diese Formpolynome
% Formfunktionen(Polynome)
phi_1 = [2 -3 0 1];
phi_2 = [1 -2 1 0];
phi_3 = [-2 3 0 0];
phi_4 = [1 -1 0 0];

% Zweifache Ableitungen davon(Strichstrich)
phi_1_Ss = polyder(polyder(phi_1));
phi_2_Ss = polyder(polyder(phi_2));
phi_3_Ss = polyder(polyder(phi_3));
phi_4_Ss = polyder(polyder(phi_4));

% Quadrierte zweiteAbl Stammf.
st_int_phi{1,1} = polyint(conv(phi_1_Ss,phi_1_Ss));
st_int_phi{2,2} = polyint(conv(phi_2_Ss,phi_2_Ss));
st_int_phi{3,3} = polyint(conv(phi_3_Ss,phi_3_Ss));
st_int_phi{4,4} = polyint(conv(phi_4_Ss,phi_4_Ss));

% Relevante Produkte zweiter Ableitungen Stammf.
st_int_phi{1,2} = polyint(conv(phi_1_Ss,phi_2_Ss));
st_int_phi{1,3} = polyint(conv(phi_1_Ss,phi_3_Ss));
st_int_phi{2,3} = polyint(conv(phi_2_Ss,phi_3_Ss));
st_int_phi{2,4} = polyint(conv(phi_2_Ss,phi_4_Ss));
st_int_phi{3,4} = polyint(conv(phi_3_Ss,phi_4_Ss));
st_int_phi{1,4} = polyint(conv(phi_1_Ss,phi_4_Ss));

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
S = zeros(2*n);

for i = 1:n
   S(2*i-1,2*i-1) = ((i<n)*int_phi(1,1)+(i>1)*int_phi(3,3)); 
   S(2*i-1,2*i) = h*((i<n)*int_phi(1,2)+(i>1)*int_phi(3,4));
   S(2*i,2*i) = h^2*((i<n)*int_phi(2,2)+(i>1)*int_phi(4,4));
   if (i<n)
   S(2*i,2*i+1) = h*int_phi(2,3);
   S(2*i,2*i+2) = h^2*int_phi(2,4);
   S(2*i-1,2*i+1) = int_phi(1,3);
    S(2*i-1,2*i+2) = h*int_phi(1,4);
    end
    
end


S = (EI/h^3)*S;
S = S+S'-diag(diag(S)); % Macht S symmetrisch voll

end

