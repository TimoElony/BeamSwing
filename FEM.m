clear
clf
tic
%% Parameter
% Eigenschaften des Bauteils
L = 0.5; % L?nge von 50cm
EI = 210*10^9*(8*10^-3)^4/12; % E-Modul mal Fl?chentr?gheitsmoment Rechteckprofil Stahl 3x3mm
mu = 1; % Masse pro Laenge
d = 10; % Daempfungsfaktor
% ZBen
% ZBen fuer Ort(Position im Balken ueber Wert der Auslenkung)
x_ZB_Ort = [0]';
ZB_Ort = [0]'; 
% ZBen f?r Neigung
x_ZB_Neig = [0]';
ZB_Neig = [0]';

% Lasten
% Kr?fte dort wo keine Lager sind
x_Kraft = [L]';
Kraft = [20]';
% Momente dort wo keine Lager sind
x_Mom = [L]';
Mom = [0]';
% Streckenlast als Polynom
q = [1];
% Anzahl St?tzstellen also doppelte Anzahl davon Basispolynome
n = 100; 
% Schrittweite
h = (L)/(n-1); 
% Intervall auf dem die Funktion sp?ter geplottet werden soll
X = 0:h/4:L;
% Anzahl Zeitschritte fuer dynamischen Fall
n_NM = 150;
% Schrittweite fuer dynamischen Fall in Sekunden
h_NM = 0.003;
% Faktoren fuer DGLVerfahren nach NEWMARK
beta_NM = 1/4;
gamma_NM = 1/2;
%% Formfunktionen(Polynome)
phi_1 = [2 -3 0 1];
phi_2 = [1 -2 1 0];
phi_3 = [-2 3 0 0];
phi_4 = [1 -1 0 0];
 %% Streckenlast
% Streckenlast f?r uns nur an den St?tzstellen relevant
% Hier als Polynom vorrausgesetzt
q_vec = zeros(2*n,1);
for i = 1:n;
    if i>1
        q_vec(2*i-1) = polyval(polyint(conv(q,phi_3)),h*(i-1))-polyval(polyint(conv(q,phi_3)),h*(i-2));
        q_vec(2*i) = polyval(polyint(conv(q,phi_4)),h*(i-1))-polyval(polyint(conv(q,phi_4)),h*(i-2));
    end
    if i<n
        q_vec(2*i-1) = q_vec(2*i-1) + polyval(polyint(conv(q,phi_1)),h*(i))-polyval(polyint(conv(q,phi_1)),h*(i-1));
        q_vec(2*i) = q_vec(2*i-1) + polyval(polyint(conv(q,phi_2)),h*(i))-polyval(polyint(conv(q,phi_2)),h*(i-1));
    end
end
 %% Erweiterte Steifheitsmatrix und rechte Seite
 % Anzahl der ZBen
 n_ZB = size(x_ZB_Ort,1)+size(x_ZB_Neig,1);
 C = zeros(2*n,n_ZB);
 rS = zeros(2*n+n_ZB,1);
 % L?sst die Lager/Kr?fte in die n?chstgelegene St?tzstelle fallen
 k_Ort = 2*round(x_ZB_Ort/h+1)-1;
 k_Neig = 2*round(x_ZB_Neig/h+1);
 k_Kraft = 2*round(x_Kraft/h+1)-1;
 k_Mom = 2*round(x_Mom/h+1);
 % Hier werden diese ganzen Bedingungen an ihre Pl?tze verwiesen
 for j = 1:size(x_ZB_Ort,1)
     C(k_Ort(j),j) = 1;
 end
 for j = 1:size(x_ZB_Neig,1)
     C(k_Neig(j),j+size(x_ZB_Ort,1)) = 1;
 end
 for j = 1:size(x_Kraft,1)
     rS(k_Kraft(j),1) = rS(k_Kraft(j),1)+Kraft(j);
 end
 for j = 1:size(x_Mom,1)
     rS(k_Mom(j),1) = rS(k_Mom(j),1)+Mom(j);
 end
 rS(2*n+1:(2*n+n_ZB),1) = [ZB_Ort;ZB_Neig];
 %% Steifheitsmatrix
S = stiffTimo(EI,L,n);
S_e = [S C;C' zeros(size(C,2))];
%% Loesung des LGS
Loesung_STAT = S_e\rS;
%% Plotten statisches System(ABD)
%Koeffizienten fuer Funktionen
w_k = Loesung_STAT(1:2*n); % So werden Zwangskr?fte ausgeblendet
% Einsetzen der Koeffizienten in Ritz-Ansatz fuer Plot
Y = hermi(X,w_k,h);
% Plot des Balken
subplot(3,1,1) , plot(X,Y),title('statische Auslenkung');
hold on
% Lager plotten
xlim([-0.0001,L])
i_Ort = round(x_ZB_Ort/h+1);
plot((i_Ort-1)*h,ZB_Ort,'or')

%% Eliminieren ZK Balken
P_dach = C*((C'*C)\C');
P = eye(2*n)-P_dach;
a_vec = [ZB_Ort;ZB_Neig]; % Vektor enth?lt die Konstanten der Zwangbedingung
f_vec = rS(1:2*n);
a_dach = C*((C'*C)\a_vec);
x = (P_dach+P*S*P)\(P*(f_vec-S*a_dach)+a_dach);

%% dynamischer Fall

%% Massenmatrix
M = massTimo(mu,L,n);
M_e = [M,zeros(size(C));zeros(size(C))',zeros(size(C,2))];

%% Dämpfungsmatrix
D = daempfTimo(d,L,n);
D_e = [D,zeros(size(C));zeros(size(C))',zeros(size(C,2))];

%% Loesen der DGL mit Newmarkmethode
u = zeros(2*n+n_ZB,n_NM);
u_p = zeros(2*n+n_ZB,n_NM);
u_pp = zeros(2*n+n_ZB,n_NM);
u_s = zeros(2*n+n_ZB,n_NM);
u_ps = zeros(2*n+n_ZB,n_NM);
% AnfangsBedingungen
% Auslenkung aus statischem Fall als aben
u(:,1) = Loesung_STAT;
 % ABen f?r einfachen fall
 u_p(:,1) = zeros(2*n+n_ZB,1);
 %u_p(:,1) = zeros(2*n+2,1);
 u_pp(:,1) = zeros(2*n+n_ZB,1);
 % Matrix der linken Seite(?ndert sich nicht)
 NM_Matrix = M_e+gamma_NM*h_NM*D_e+beta_NM*h_NM^2*S_e;
% Loesung
 for j = 1:(n_NM-1)
     % Algorithmus von NEWMARK:
     % (Wir sollten vielleicht lieber eine LR Zerlegung anstellen
     % dann m?ssen wir in der Schleife nur noch vorrueck machen,
     % nicht mehr jedesmal ein echtes LGS l?sen)
     % a)
        u_s(:,j) = u(:,j)+u_p(:,j)*h_NM+(1/2-beta_NM)*u_pp(:,j)*h_NM^2;
        u_ps(:,j) = u_p(:,j)+(1-gamma_NM)*u_pp(:,j)*h_NM;
     % b) (dieser Teil kann optimiert werden)
        %u_pp(:,j+1) = NM_Matrix\(rechte_Seite-S_e*u_s(:,j));
        u_pp(:,j+1) = NM_Matrix\([zeros(2*n,1);ZB_Ort;ZB_Neig]-D_e*u_ps(:,j)-S_e*u_s(:,j));

     % c)
        u(:,j+1) = u_s(:,j)+beta_NM*u_pp(:,j+1)*h_NM^2;
        u_p(:,j+1) = u_ps(:,j)+gamma_NM*u_pp(:,j+1)*h_NM;
 end
 
%% Dynamische Energien und Hermite-Interpolation
 E_el_dyn = zeros(n_NM,1);
 E_kin_dyn = zeros(n_NM,1);
 x_range = linspace(0,L,30);
 w_hermi = zeros(size(x_range,2),n_NM);
 for i = 1:n_NM
     % Zu jedem Zeitschritt wird Energie berechnet
     E_el_dyn(i) = EI/2*u(1:2*n,i)'*S*u(1:2*n,i);
     E_kin_dyn(i) = EI/2*u_p(1:2*n,i)'*M*u_p(1:2*n,i);
     % Und auch die Auslenkung des Balken in z-Richtung
     w_hermi(:,i) = hermi(x_range,u(1:2*n,i),h)';
 end
 E_ges_dyn = E_el_dyn+E_kin_dyn;
 toc
%% Plot
subplot(3,1,2)
 hold on
 xlim([-0.0001,L])
 ylim([min(min(w_hermi)),max(max(w_hermi))])
 for i = 1:n_NM
     cla
     fig = plot(x_range,w_hermi(:,i));
     plot((i_Ort-1)*h,ZB_Ort,'or')
     titelsek = ['Auslenkung nach: ' num2str((i-1)*h_NM)];
     titelsek = [titelsek ' sekunden'];
     title(titelsek);
     %legend(['E_ges :' num2str(E_ges_dyn)])
     pause(0.03)
     drawnow
     if(ishandle(fig)~=1)
         break;
     end
 end
 cla
 subplot(3,1,1)
 hold on
xlim([-0.0001,L])
plot(X,Y);
plot((i_Ort-1)*h,ZB_Ort,'or')
title('statische Auslenkung');
 
 subplot(3,1,2)
hold on
 plot(x_range,w_hermi(:,n_NM));
 plot((i_Ort-1)*h,ZB_Ort,'or')
 title(titelsek);
 subplot(3,1,3)
hold on
 plot(h_NM*2:h_NM:h_NM*(n_NM),E_kin_dyn(2:n_NM),'r');
 plot(h_NM*2:h_NM:h_NM*(n_NM),E_el_dyn(2:n_NM));
 plot(h_NM*2:h_NM:h_NM*(n_NM),E_ges_dyn(2:n_NM),'k');
 xlabel ('Zeit in s')
 ylabel ('Energie in J ')
 
 legend('Kinetische Energie','Elastische Energie','Gesamt Energie')