function [Y] = form_admittance(phase,Zbase,Sbase)

Nnode = 36;
Ybase = 1/Zbase;
%25 loads
P_l = [0 0 0;
           140 140 350; 
           0 0 0; 
           0 0 0; 
           0 0 85; 
           8 85 0; 
           0 0 85; 
           0 0 0;
           17 21 0;
           85 0 0; 
           0 0 85; 
           0 0 0; 
           0 42 0; 
           0 140 21; 
           0 0 0; 
           0 42 0;
           0 0 0;
           0 0 42;
           42 0 0; 
           42 0 0; 
           42 42 42; 
           0 0 85; 
           0 0 0; 
           0 85 0;
           0 0 0;
           0 0 42; 
           85 0 0;
           0 0 42; 
           140 0 0;
           126 0 0;
           0 0 0;
           0 0 42;
           0 0 85; 
           0 0 0;
           0 42 0;
           0 0 85].';

Q_l = [0 0 0; 
       70 70 175; 
       0 0 0; 
       0 0 0;
       0 0 40; 
       4 40 0; 
       0 0 40;
       0 0 0;
       8 10 0;
       40 0 0;
       0 0 40;
       0 0 0;
       0 21 0;
       0 70 10;
       0 0 0;
       0 21 0;
       0 0 0;  
       0 0 21;
       21 0 0;
       21 0 0;
       21 21 21 ;
       0 0 40;
       0 0 0;
       0 40 0;
       0 0 0;
       0 0 21;
       40 0 0;
       0 0 21;
       70 0 0;
       62 0 0;
       0 0 0;
       0 0 21; 
       0 0 40;
       0 0 0;
       0 21 0;
       0 0 40].';
   
   
   
%------------------------------------------------------------------------
% impedance matrix
%------------------------------------------------------------------------

% Configuration 721
Zs1 = [0.2926+0.1973i 0.0673-0.0368i 0.0337-0.0417i;
       0.0673-0.0368i 0.2646+0.1900i 0.0673-0.0368i;
       0.0337-0.0417i 0.0673-0.0368i 0.2926+0.1973i]./Zbase;
Ys1 = sqrt(-1)*159.7919*(10^-6).*eye(3)./Ybase;

% Configuration 722
Zs2 = [0.4751+0.2973i 0.1629-0.0326i 0.1234-0.0607i;
       0.1629-0.0326i 0.4488+0.2678i 0.1629-0.0326i;
       0.1234-0.0607i 0.1629-0.0326i 0.4751+0.2973i]./Zbase;
Ys2 = sqrt(-1)*127.8306*(10^-6).*eye(3)./Ybase;

% Configuration 723
Zs3 = [1.2936+0.6713i 0.4871+0.2111i 0.4585+0.1521i;
       0.4871+0.2111i 1.3022+0.6326i 0.4871+0.2111i;
       0.4585+0.1521i 1.2936+0.6713i 1.2936+0.6713i]./Zbase;
Ys3 = sqrt(-1)*74.8405*(10^-6).*eye(3)./Ybase;

% Configuration 724
Zs4 = [2.0952+0.7758i 0.5204+0.2738i 0.4926+0.2123i;
       0.5204+0.2738i 2.1068+0.7398i 0.5204+0.2738i;
       0.4926+0.2123i 0.5204+0.2738i 2.0952+0.7758i]./Zbase;
Ys4 = sqrt(-1)*60.2483*(10^-6).*eye(3)./Ybase;


%--------------------
% line matrices 
%--------------------

% mile = 5280 feet 
convfm = (1/5280);

Z12 = Zs1*(1850)*convfm;
Z12i = pinv(Z12);
Y12 = .5.*Ys1*(1850)*convfm;

Z23 = Zs2*(960)*convfm;
Z23i = pinv(Z23);
Y23 = .5.*Ys2*(960)*convfm;

Z34 = Zs4*(400)*convfm;
Z34i = pinv(Z34);
Y34 = .5.*Ys4*(400)*convfm;

Z45 = Zs4*(240)*convfm;
Z45i = pinv(Z45);
Y45 = .5.*Ys4*(240)*convfm;

Z46 = Zs4*(320)*convfm;
Z46i = pinv(Z46);
Y46 = .5.*Ys4*(320)*convfm;

Z37 = Zs3*(360)*convfm;
Z37i = pinv(Z37);
Y37 = .5.*Ys3*(360)*convfm;

Z78 = Zs3*(520)*convfm;
Z78i = pinv(Z78);
Y78 = .5.*Ys3*(520)*convfm;

Z89 = Zs4*(80)*convfm;
Z89i = pinv(Z89);
Y89 = .5.*Ys4*(80)*convfm;

Z910 = Zs4*(520)*convfm;
Z910i = pinv(Z910);
Y910 = .5.*Ys4*(520)*convfm;

Z811 = Zs3*(800)*convfm;
Z811i = pinv(Z811);
Y811 = .5.*Ys3*(800)*convfm;

Z1112 = Zs4*(920)*convfm;
Z1112i = pinv(Z1112);
Y1112 = .5.*Ys4*(920)*convfm;

Z1213 = Zs4*(760)*convfm;
Z1213i = pinv(Z1213);
Y1213 = .5.*Ys4*(760)*convfm;

Z1214 = Zs4*(120)*convfm;
Z1214i = pinv(Z1214);
Y1214 = .5.*Ys4*(120)*convfm;

Z1115 = Zs3*(600)*convfm;
Z1115i = pinv(Z1115);
Y1115 = .5.*Ys3*(600)*convfm;

Z1516 = Zs4*(280)*convfm;
Z1516i = pinv(Z1516);
Y1516 = .5.*Ys4*(280)*convfm;

Z317 = Zs2*(1320)*convfm;
Z317i = pinv(Z317);
Y317 = .5.*Ys2*(1320)*convfm;

Z1722 = Zs3*(600)*convfm;
Z1722i = pinv(Z1722);
Y1722 = .5.*Ys3*(600)*convfm;

Z2223 = Zs3*(200)*convfm;
Z2223i = pinv(Z2223);
Y2223 = .5.*Ys3*(200)*convfm;

Z1718 = Zs4*(240)*convfm;
Z1718i = pinv(Z1718);
Y1718 = .5.*Ys4*(240)*convfm;

Z1819 = Zs3*(280)*convfm;
Z1819i = pinv(Z1819);
Y1819 = .5.*Ys3*(280)*convfm;

Z1920 = Zs4*(280)*convfm;
Z1920i = pinv(Z1920);
Y1920 = .5.*Ys4*(280)*convfm;

Z1921 = Zs4*(200)*convfm;
Z1921i = pinv(Z1921);
Y1921 = .5.*Ys4*(200)*convfm;

Z2324 = Zs3*(600)*convfm;
Z2324i = pinv(Z2324);
Y2324 = .5.*Ys3*(600)*convfm;

Z2325 = Zs3*(320)*convfm;
Z2325i = pinv(Z2325);
Y2325 = .5.*Ys3*(320)*convfm;

Z2526 = Zs4*(320)*convfm;
Z2526i = pinv(Z2526);
Y2526 = .5.*Ys4*(320)*convfm;

Z2527 = Zs3*(320)*convfm;
Z2527i = pinv(Z2527);
Y2527 = .5.*Ys3*(320)*convfm;

Z2728 = Zs3*(560)*convfm;
Z2728i = pinv(Z2728);
Y2728 = .5.*Ys3*(560)*convfm;

Z2829 = Zs3*(640)*convfm;
Z2829i = pinv(Z2829);
Y2829 = .5.*Ys3*(640)*convfm;

Z2930 = Zs3*(400)*convfm;
Z2930i = pinv(Z2930);
Y2930 = .5.*Ys3*(400)*convfm;

Z3031 = Zs3*(400)*convfm;
Z3031i = pinv(Z3031);
Y3031 = .5.*Ys3*(400)*convfm;

Z3132 = Zs3*(400)*convfm;
Z3132i = pinv(Z3132);
Y3132 = .5.*Ys3*(400)*convfm;

Z3133 = Zs4*(200)*convfm;
Z3133i = pinv(Z3133);
Y3133 = .5.*Ys4*(200)*convfm;

Z2834 = Zs4*(520)*convfm;
Z2834i = pinv(Z2834);
Y2834 = .5.*Ys4*(520)*convfm;

Z3435 = Zs4*(1280)*convfm;
Z3435i = pinv(Z3435);
Y3435 = .5.*Ys4*(1280)*convfm;

Z3436 = Zs4*(200)*convfm;
Z3436i = pinv(Z3436);
Y3436 = .5.*Ys4*(200)*convfm;


% network admittance matrix
oo = zeros(3);

Y_net = [Z12i+Y12 -Z12i zeros(3,3*(Nnode-2));
         -Z12i Z12i+Z23i+Y12+Y23 -Z23i zeros(3,3*(Nnode-3));
         oo -Z23i Z23i+Z34i+Z37i+Z317i+Y23+Y34+Y37+Y317 -Z34i oo oo -Z37i oo oo oo oo oo oo oo oo oo -Z317i zeros(3,3*(Nnode-17));
         oo oo -Z34i Z34i+Z45i+Z46i+Y34+Y45+Y46 -Z45i -Z46i zeros(3,3*(Nnode-6));
         oo oo oo -Z45i Z45i+Y45 zeros(3,3*(Nnode-5));
         oo oo oo -Z46i oo Z46i+Y46 zeros(3,3*(Nnode-6));
         oo oo -Z37i oo oo oo Z37i+Z78i+Y37+Y78 -Z78i zeros(3,3*(Nnode-8));
         oo oo oo oo oo oo -Z78i Z78i+Y78+Z89i+Y89+Z811i+Y811 -Z89i oo -Z811i zeros(3,3*(Nnode-11));
         zeros(3,3*7) -Z89i Z89i+Y89+Z910i+Y910 -Z910i zeros(3,3*(Nnode-10));
         zeros(3,3*8) -Z910i Z910i+Y910 zeros(3,3*(Nnode-10));
         zeros(3,3*7) -Z811i oo oo Z811i+Y811+Z1112i+Y1112+Z1115i+Y1115 -Z1112i oo oo -Z1115i zeros(3,3*(Nnode-15));
         zeros(3,3*10) -Z1112i Z1112i+Z1213i+Z1214i+Y1112+Y1213+Y1214 -Z1213i -Z1214i zeros(3,3*(Nnode-14));
         zeros(3,3*11) -Z1213i Z1213i+Y1213 zeros(3,3*(Nnode-13));
         zeros(3,3*11) -Z1214i oo Z1214i+Y1214 zeros(3,3*(Nnode-14))
         zeros(3,3*10) -Z1115i oo oo oo Z1115i+Z1516i+Y1115+Y1516 -Z1516i zeros(3,3*(Nnode-16));
         zeros(3,3*14) -Z1516i Z1516i+Y1516 zeros(3,3*(Nnode-16))
         oo oo -Z317i zeros(3,3*13) Z317i+Y317+Z1718i+Y1718+Z1722i+Y1722 -Z1718i oo oo oo -Z1722i zeros(3,3*(Nnode-22));
         zeros(3,3*16) -Z1718i Z1718i+Y1718+Z1819i+Y1819 -Z1819i zeros(3,3*(Nnode-19));
         zeros(3,3*17) -Z1819i Z1819i+Y1819+Z1920i+Y1920+Z1921i+Y1921 -Z1920i -Z1921i zeros(3,3*(Nnode-21));
         zeros(3,3*18) -Z1920i Z1920i+Y1920 zeros(3,3*(Nnode-20)); 
         zeros(3,3*18) -Z1921i oo Z1921i+Y1921 zeros(3,3*(Nnode-21));
         zeros(3,3*16) -Z1722i oo oo oo oo Z1722i+Y1722+Z2223i+Y2223 -Z2223i zeros(3,3*(Nnode-23));
         zeros(3,3*21) -Z2223i Z2223i+Y2223+Z2324i+Y2324+Z2325i+Y2325 -Z2324i -Z2325i zeros(3,3*(Nnode-25));
         zeros(3,3*22) -Z2324i Z2324i+Y2324 zeros(3,3*(Nnode-24));
         zeros(3,3*22) -Z2325i oo Z2325i+Y2325+Z2526i+Y2526+Z2527i+Y2527 -Z2526i -Z2527i zeros(3,3*(Nnode-27));
         zeros(3,3*24) -Z2526i Z2526i+Y2526 zeros(3,3*(Nnode-26));
         zeros(3,3*24) -Z2527i oo Z2527i+Y2527+Z2728i+Y2728 -Z2728i zeros(3,3*(Nnode-28));
         zeros(3,3*26) -Z2728i Z2728i+Y2728+Z2829i+Y2829+Z2834i+Y2834 -Z2829i oo oo oo oo -Z2834i zeros(3,3*(Nnode-34));
         zeros(3,3*27) -Z2829i Z2829i+Y2829+Z2930i+Y2930 -Z2930i zeros(3,3*(Nnode-30));
         zeros(3,3*28) -Z2930i Z2930i+Y2930+Z3031i+Y3031 -Z3031i zeros(3,3*(Nnode-31));
         zeros(3,3*29) -Z3031i Z3031i+Y3031+Z3132i+Y3132+Z3133i+Y3133 -Z3132i -Z3133i zeros(3,3*(Nnode-33))
         zeros(3,3*30) -Z3132i Z3132i+Y3132 zeros(3,3*(Nnode-32));
         zeros(3,3*30) -Z3133i oo Z3133i+Y3133 zeros(3,3*(Nnode-33));
         zeros(3,3*27) -Z2834i zeros(3,3*5) Z2834i+Y2834+Z3435i+Y3435+Z3436i+Y3436 -Z3435i -Z3436i;
         zeros(3,3*33) -Z3435i Z3435i+Y3435 oo;
         zeros(3,3*33) -Z3436i oo Z3436i+Y3436];


%P = P_l(phase,:);
%Q = Q_l(phase,:);
Y = Y_net(phase:3:end,phase:3:end);
   
fac = Sbase/1000;

%P = P./fac;
%Q = Q./fac;


end
   