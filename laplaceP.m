function [pBar ] = laplaceP ( Gamma, alpha, beta, gamma, lmf, lvf, lmv, wf, wv, ...
    wm, kD, Sw, Slv, u, rlvD, reD1, reD2, CwD, ClvD, hD, key, WR, LT )

Y = -Gamma*beta/( (gamma + u)^alpha );

m1 =  -( lmf^2/(lmv + lmf + u*wm) - lmf - lvf - wf*u );
m2 =  -( lmv^2/(lmv + lmf + u*wm) - lmv - lvf - wv*u );
m =  -( lmf*lmv/(lmv + lmf + u*wm) + lvf  );

nu1Squared = ( kD*m2+(1-kD)*m1 ) / ( 2*kD*(1-kD) ) + ...
    sqrt( ( kD*m2 + (1-kD)*m1 )^2 + 4*(m^2 - m1*m2)*kD*(1-kD) ) / ( 2*kD*(1-kD) );
nu2Squared = ( kD*m2+(1-kD)*m1 ) / ( 2*kD*(1-kD) ) - ...
    sqrt( ( kD*m2 + (1-kD)*m1 )^2 + 4*(m^2 - m1*m2)*kD*(1-kD) ) / ( 2*kD*(1-kD) );

% nu1Squared = ( kD*m2+(1-kD)*m1 ) / ( 2*kD*(1-kD) ) + ...
%     sqrt( ( kD*m2 + (-1+kD)*m1 )^2 + 4*(m^2)*kD*(1-kD) ) / ( 2*kD*(1-kD) );
% nu2Squared = ( kD*m2+(1-kD)*m1 ) / ( 2*kD*(1-kD) ) - ...
%     sqrt( ( kD*m2 + (-1+kD)*m1 )^2 + 4*(m^2)*kD*(1-kD) ) / ( 2*kD*(1-kD) );

nu1 = sqrt(nu1Squared);
nu2 = sqrt(nu2Squared);

A1 = -( (1-kD)*nu1Squared - m2 ) / m;
A2 = -( (1-kD)*nu2Squared - m2 ) / m;
F1 = A1;
F2 = A2;

M1 = A1*besseli(0, nu1);
M2 = A2*besseli(0, nu2);
M3 = F1*besselk(0, nu1);
M4 = F2*besselk(0, nu2);
M5 = A1*nu1*besseli(1, nu1);
M6 = A2*nu2*besseli(1, nu2);
M7 = -F1*nu1*besselk(1, nu1);
M8 = -F2*nu2*besselk(1, nu2);

N1 = besseli(0, nu1*rlvD);
N2 = besseli(0, nu2*rlvD);
N3 = besselk(0, nu1*rlvD);
N4 = besselk(0, nu2*rlvD);
N5 = nu1*rlvD*besseli(1, nu1*rlvD);
N6 = nu2*rlvD*besseli(1, nu2*rlvD);
N7 = -nu1*rlvD*besselk(1, nu1*rlvD);
N8 = -nu2*rlvD*besselk(1, nu2*rlvD);  %%%%%%%%% WARNING! nu2 or nu1?

a1 = M1 - Sw*M5 + Slv*N5 - N1;
a2 = M2 - Sw*M6 + Slv*N6 - N2;
a3 = M3 - Sw*M7 + Slv*N7 - N3;
a4 = M4 - Sw*M8 + Slv*N8 - N4;
b1 = u*CwD*(M1 - Sw*M5) + u*ClvD*(N1 - Slv*N5) - hD*Sw*M5 - (1-hD)*Slv*N5;
b2 = u*CwD*(M2 - Sw*M6) + u*ClvD*(N2 - Slv*N6) - hD*Sw*M6 - (1-hD)*Slv*N6;
b3 = u*CwD*(M3 - Sw*M7) + u*ClvD*(N3 - Slv*N7) - hD*Sw*M7 - (1-hD)*Slv*N7;
b4 = u*CwD*(M4 - Sw*M8) + u*ClvD*(N4 - Slv*N8) - hD*Sw*M8 - (1-hD)*Slv*N8;
d1 = nu1*A1*besseli(1, nu1*reD1);
d2 = nu2*A2*besseli(1, nu2*reD1);
d3 = -nu1*F1*besselk(1, nu1*reD1);
d4 = -nu2*F2*besselk(1, nu2*reD2);
d5 = nu1*besseli(1, nu1*reD2);
d6 = nu2*besseli(1, nu2*reD2);
d7 = -nu1*besselk(1, nu1*reD2);
d8 = -nu2*besselk(1, nu2*reD2);

switch key
    case 2 % Finite reservoir with closed outer boundary
        bVector = [1/u Y 0 0]';
        aMatrix = [a1 a2 a3 a4; b1 b2 b3 b4; d1 d2 d3 d4; d5 d6 d7 d8];
        Dvector = aMatrix\bVector;
        D1 = Dvector(1);
        D2 = Dvector(2);
        D3 = Dvector(3);
        D4 = Dvector(4);
        
        pBar = D1*(M1-Sw*M5) + D2*(M2-Sw*M6) + D3*(M3-Sw*M7) + D4*(M4-Sw*M8);
    case 1 % Infinite reservoir
        D3 = ( a2/u - a4*Y ) / ( a2*a3 - a1*a4 );
        D4 = ( a1/u - a3*Y ) / ( a1*a4 - a2*a3 );
            pBar = D3*(M3 - M7) + D4*(M4 - M8);
    case 3 % Camacho-Velazquez (2005)
        if kD == 1
            g = m1 - m^2/m2;
            noStorageP = besselk(0, sqrt(g))/ ( u*sqrt(g)*besselk(1, sqrt(g)) );
            pBar = 1/ ( u^2*CwD + 1/noStorageP );
        else
            C1 = 1/u * ( nu1*besselk(1, nu1)*(kD*A1+1-kD)+ (nu2*besselk(1, nu2)*(kD*A2+1-kD))*...
                ( ( besselk(0, nu1)*(A1-1)+besselk(1, nu1)*nu1*(Sw*A1-Slv) ) /...
                ( besselk(0,nu2)*(1-A2)+besselk(1,nu2)*nu2*(Slv-Sw*A2) ) ) )^(-1);
            
            C2 = 1/u * ( nu2*besselk(1, nu2)*(kD*A2+1-kD) + (nu1*besselk(1, nu1)*(kD*A1+1-kD))*...
                ( ( besselk(0, nu2)*(1-A2)+besselk(1, nu2)*nu2*(Slv-Sw*A2) ) /...
                ( besselk(0,nu1)*(A1-1)+besselk(1,nu1)*nu1*(Sw*A1-Slv) ) ) )^(-1);
%             C2 = 1/u* (  ( ((1-A2)*besselk(0, nu2))/((A1-1)*besselk(0, nu1)) )*nu1*besselk(1, nu1)...
%                 *(kD*A1+1-kD) + nu2*besselk(1, nu2)*(kD*A2+1-kD))^(-1);
%             
%             C1 = 1/u* (  nu1*besselk(1, nu1)*(kD*A1+1-kD) + ...
%                 nu2*besselk(1, nu2)*(kD*A2+1-kD)*( ((A1-1)*besselk(0, nu1))/((1-A2)*besselk(0,nu2)) ))^(-1);
            noStorageP = C1*besselk(0, nu1)+C2*besselk(0,nu2)+ ...
                Slv* ( C1*nu1*besselk(1, nu1)+C2*nu2*besselk(1, nu2) );
            pBar = 1/ ( u^2*CwD + 1/noStorageP );
%             pBar = C1*besselk(0,nu1) + C2*besselk(0,nu2);
        end
    case 4 
        b1 = lvf*(lmv+lmf)+lmf*lmv;
        b2 = lvf*wm;
        b3 = lmv*(lvf+lmf)+lmf*lvf;
        b4 = wv*(lmv+lmf)+wm*(lmv+lvf);
        b5 = wm*wv;
        denom = b3*(lmv+lmf)+u*((lmv+lmf)*b4+...
            wm*b3)+u^2*((lmv+lmf)*b5+wm*b4)+u^3*wm*b5;
        g = lmf*(1-lmv*b1*lmf*b3 + u*(lmv*b2+lmf*b4) + u^2*lmf*b5/denom)...
            + lvf*(1-(b1+b2*u)/(b3+b4*u+b5*u^2)) + wf*u;
        if strcmp(WR, 'on') == 1
            g = u*(lmf+u*wf*(1-wf))/(lmf+u*(1-wf)); % warren and root
        elseif strcmp(LT, 'on') == 1
            t1 = wm*lvf*(lmv+lmf)+(1-wf)*lmv*(lmv+lmf);
            t2 = lvf*(lmv+lmf)^2+lmv*lmf*(lmv+lmf);
            t3 = wm*2*(lmv+lmf)*lvf+wv*lmf^2+(1-wf)*(lmv^2+2*lmf*lmv);
            t4 = wv*lmf+lmv*(1-wf);
            t5 = lmv*(lvf+lmf)+lmf*lvf;
            t6 = wv*(lmf+lvf)+(1-wf)*(lmv+lvf);
            g = u*(lmf*(t1/(t2+u*t3))+lvf*(t4/(t5+u*t6)) + wf);
        end
        rD = 1;
        pBar = besselk(0, sqrt(g)*rD)/ ( u*sqrt(g)*besselk(1, sqrt(g)) );
end
end