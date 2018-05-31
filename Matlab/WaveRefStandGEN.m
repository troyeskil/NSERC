%
%  Movie of Normal Incident plane Wave (lossless to lossy dielectric)
%  G. Bridges, modifed Feb, 2014 
%  Course ECE3590
%  Modified verion of Matlab exercise 8.12 by 
%           [B.M. Notaros, MATLAB-Based Electromagnetics, Pearson, 2014]
%
%          WaveRefStandGEN(Er2,Sig2)  with     Er1=1.0, Sig1=0.0
%
%  Example input values: 9.0     0.0
%                        36.0    0.0
%                        1.0     0.0
%                        9.0     0.02
%                        9.0     0.2
%  Output:
%  E1inc = BLUE
%  E1ref = RED
%  E2trans = GREEN
%  E1tot = MAGENTA
%  Envelope = BLACK
%
%  to make a move  set [M]=WaveMovie(..)
%  then to play movie(M,3) repeating 3 times
%
function WaveRefStandGEN(epsr2,sigma2)
    %
    PlotType=1;
    Ntime=2;
    epsr2=9.0;
    sigma2=0.2;
    %
    f=750*10^6;
    w=2*pi*f;
    %
    eps0=8.854e-12;
    u0=4*pi*1.0e-7;
    c0=1.0/sqrt(u0*eps0);
    eta0=sqrt(u0/eps0);
    %
    Ei0=1.0/sqrt(2.0);
    theta0=0.0;
    epsr1=1;
    mur1=1;
    sigma1=0.0;
    %epsr2=4;
    mur2=1;
    %sigma2=0.01;
    deg2rad=180.0/pi;
    theta0rad=theta0*deg2rad;
    %
    [gamma1,alpha1,beta1,eta1,abseta1,phi1]=basicPropParam(epsr1,mur1,sigma1,f);
    eta1=abseta1*exp(i*phi1);
    lambda1=2*pi/beta1;
    [gamma2,alpha2,beta2,eta2,abseta2,phi2]=basicPropParam(epsr2,mur2,sigma2,f);
    eta2=abseta2*exp(i*phi2);
    lambda2=2*pi/beta2;
    GAMMA=gammaReflCoef(eta1,eta2,'r');
    TAU=tauTranCoef(eta1,eta2,'r');
    GAMMA, TAU
    er0=GAMMA*Ei0*exp(i*theta0rad);
    et0=TAU*Ei0*exp(i*theta0rad);
    zlim=max(2*lambda1,2*lambda1);
    dz=zlim/100;
    z=[-zlim:dz:0];
    T1=zlim*sqrt(eps0*epsr1*u0*mur1);
    T2=zlim*sqrt(eps0*epsr2*u0*mur2);
    dt=max(T1,T2)/400;
    t1=[0:dt:T1];
    t2=[0:dt:Ntime*T2];
    eim=Ei0*sqrt(2);
    erm=abs(er0)*sqrt(2);
    etm=abs(et0)*sqrt(2);
    emax=max([eim,erm,etm]);
    elim=emax+emax/2;
    L=[-elim:elim];
    boundary=zeros(1,length(L));
    %
    for k=1:length(t1);
        Ei=eim.*exp(-alpha1*z).*cos(w*t1(k)-beta1.*z+theta0rad).*heaviside(w*t1(k)-beta1.*(z+zlim));
        [errori,ni]=min(heaviside(w*t1(k)-beta1.*(z-zlim)));
        plot(z(1:ni-2),Ei(1:ni-2),'b','linewidth',2); hold on;
        plot(boundary,L,'k'); hold off;
        M(k) = getframe;
    end;
    %
    E1env=abs(eim.*exp(-alpha1*z - i*beta1.*z) + erm.*exp(alpha1*z + i*beta1.*z + i*angle(er0)));
    E2env=abs(etm.*exp(-alpha2*(-z) - i*beta2.*z));
    %
    for k=1:length(t2);
        Ei=eim.*exp(-alpha1*z).*cos(w*(t2(k)+T1)-beta1.*z+theta0rad);
        Er=erm.*exp(alpha1*z).*cos(w*(t2(k)+T1)+beta1.*z+angle(er0)).*heaviside(w*t2(k)+beta1.*(z));
        Etot=Ei+Er;
        Et=etm.*exp(-alpha2*(-z)).*cos(w*(t2(k)+T1)-beta2.*(-z)+angle(et0)).*heaviside(w*t2(k)-beta2.*(-z));
        [errorr,nr]=max(heaviside(w*t2(k)+beta1.*(z)));
        [errort,nt]=max(heaviside(w*t2(k)-beta2.*(-z)));
        if PlotType == 1; plot(z,Ei,'b','linewidth',2); hold on; end;
        if PlotType == 1; plot(z(nr:length(z)),Er(nr:length(z)),'r','linewidth',2); end;
        plot(-z(nt:length(z)),Et(nt:length(z)),'g-','linewidth',2);
        if PlotType == 2; hold on; end;
        plot(z,Etot,'c','linewidth',4);
        plot(z,E1env,'k','linewidth',2);
        plot(z,-E1env,'k','linewidth',2);
        plot(-z,E2env,'k','linewidth',2);
        plot(-z,-E2env,'k','linewidth',2);
        plot(boundary,L,'k'); hold off;
        M(k) = getframe;
    end;
%
%
%
%
%  calculate medium propagation parameters  
%
function [gamma,alpha,beta,eta,etaabs,etaphi] = basicPropParam(EPSR,MUR,SIGMA,f)
    %
    eps0=8.854e-12;
    u0=4*pi*1.0e-7;
    eta0=sqrt(u0/eps0);
    wrad=2*pi*f;
    %
    var1=1i*wrad*u0*MUR;
    var2=SIGMA+1i*wrad*eps0*EPSR;
    gamma=sqrt(var1*var2);
    eta=sqrt(var1/var2);
    alpha=real(gamma);
    beta=imag(gamma);
    etaabs=abs(eta);
    etaphi=angle(eta);
%
%
%
%
%
%  calculate GAMMA for plane interface
%
function [result1,result2] = gammaReflCoef(eta1,eta2,k)
    %
    if k == 'r';
        result1=(eta2-eta1)/(eta2+eta1);
        result2=' ';
    end;
    if k == 'p';
        result1=abs((eta2-eta1)/(eta2+eta1));
        result2=angle((eta2-eta1)/(eta2+eta1));
    end;
    if k == 'd';
        result1=20*log10(abs((eta2-eta1)/(eta2+eta1)));
        result2=' ';
    end;
%
%
%
%  calculate TAU for plane interface
%
function [result1,result2] = tauTranCoef(eta1,eta2,k)
    %
    if k == 'r';
        result1=(2*eta2)/(eta2+eta1);
        result2=' ';
    end;
    if k == 'p';
        result1=abs((2*eta2)/(eta2+eta1));
        result2=angle((2*eta2)/(eta2+eta1));
    end;
    if k == 'd';
        result1=20*log10(abs((2*eta2)/(eta2+eta1)));
        result2=' ';
    end;
%
%
%