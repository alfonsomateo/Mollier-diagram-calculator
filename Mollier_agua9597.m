format long;

%   Etiquetas:
%   Densidad                   'rho'   (si hay dos fases es a media)
%   Temperatura                'T'
%   Entropía                   's'
%   Entalpía                   'h'
%   Presión de saturación      'psat'
%   Densidad del líquido       'rhow'
%   Densidad del vapor         'rhov'
%   Entropía del líquido       'sw'
%   Entropía del vapor         'sv'
%   Entalpía del líquido       'hw'
%   Entalpía del vapor         'hv'
%   Fracción másica del vapor  'frac'

% LAS PRESIONES ESTÁN EN kPa Y LAS ENTALPÍAS Y ENTROPÍAS USAN kJ
% EL RESTO DE MAGNITUDES ESTÁN EN UNIDADES DEL S.I.

%Variable 1

var1='s';

val1=7;

%Variable 2

var2='h';

val2=2500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vars={var1,var2};
vals=[val1,val2];

[sysflag,rho,T,p,s,h,psat,rhow,rhov,sw,sv,hw,hv,frac]=Magia(vars,vals);

%Imprimir resultados
disp('----Resultados----');
disp(sysflag);
fprintf('---------------------------------------------------------\n');
fprintf('Densidad                     %.6f kg/m³\n',rho);
fprintf('---------------------------------------------------------\n');
fprintf('Temperatura                  %.6f K\n',T);
fprintf('---------------------------------------------------------\n');
fprintf('Presión                      %.6f kPa\n',p);
fprintf('---------------------------------------------------------\n');
fprintf('Entropía                     %.6f kJ/(kg K)\n',s);
fprintf('---------------------------------------------------------\n');
fprintf('Entalpía                     %.6f kJ/kg\n',h);
fprintf('---------------------------------------------------------\n');
fprintf('Presión de saturación        %.6f kPa\n',psat);
fprintf('---------------------------------------------------------\n');
if frac<1
    fprintf('Densidad del líquido         %.6f kg/m³\n',rhow);
    fprintf('---------------------------------------------------------\n');
    fprintf('Densidad del vapor           %.6f kg/m³\n',rhov);
    fprintf('---------------------------------------------------------\n');
    fprintf('Entropía del líquido         %.6f kJ/(kg K)\n',sw);
    fprintf('---------------------------------------------------------\n');
    fprintf('Entropía del vapor           %.6f kJ/(kg K)\n',sv);
    fprintf('---------------------------------------------------------\n');
    fprintf('Entalpía del líquido         %.6f kJ/kg\n',hw);
    fprintf('---------------------------------------------------------\n');
    fprintf('Entalpía del vapor           %.6f kJ/kg\n',hv);
    fprintf('---------------------------------------------------------\n');
    fprintf('Fracción másica del vapor    %.6f\n',frac);
    fprintf('---------------------------------------------------------\n');
end

function [fphi0,fphir,phi0d,phi0t,phird,phirt,p,s,h,psat,rhow,rhov,sw,sv,hw,hv,frac]=create_Functions()
dh=1e-5; %intervalo en las derivadas
O=4; %orden de las derivadas
N=O+1;
rhoc=322;
Tc=647.096;
R=0.46151805;
n0=[-8.3204464837497,...
    6.6832105275932,...
    3.00632,...
    0.012436,...
    0.97315,...
    1.27950,...
    0.96956,...
    0.24873];
gamma0=[1.28728967,...
    3.53734222,...
    7.74073708,...
    9.24437796,...
    27.5075105];

nr=[0.12533547935523e-1,...
    0.78957634722828e1,...
    -0.87803203303561e1,...
    0.31802509345418,...
    -0.26145533859358,...
    -0.78199751687981e-2,...
    0.88089493102134e-2,...
    -0.66856572307965,...
    0.20433810950965,...
    -0.66212605039687e-4,...
    -0.19232721156002,...
    -0.25709043003438,...
    0.16074868486251,...
    -0.40092828925807e-1,...
    0.39343422603254e-6,...
    -0.75941377088144e-5,...
    0.56250979351888e-3,...
    -0.15608652257135e-4,...
    0.11537996422951e-8,...
    0.36582165144204e-6,...
    -0.13251180074668e-11,...
    -0.62639586912454e-9,...
    -0.10793600908932,...
    0.17611491008752e-1,...
    0.22132295167546,...
    -0.40247669763528,...
    0.58083399985759,...
    0.49969146990806e-2,...
    -0.31358700712549e-1,...
    -0.74315929710341,...
    0.47807329915480,...
    0.20527940895948e-1,...
    -0.13636435110343,...
    0.14180634400617e-1,...
    0.83326504880713e-2,...
    -0.29052336009585e-1,...
    0.38615085574206e-1,...
    -0.20393486513704e-1,...
    -0.16554050063734e-2,...
    0.19955571979541e-2,...
    0.15870308324157e-3,...
    -0.16388568342530e-4,...
    0.43613615723811e-1,...
    0.34994005463765e-1,...
    -0.76788197844621e-1,...
    0.22446277332006e-1,...
    -0.62689710414685e-4,...
    -0.55711118565645e-9,...
    -0.19905718354408,...
    0.31777497330738,...
    -0.11841182425981,...
    -0.31306260323435e2,...
    0.31546140237781e2,...
    -0.25213154341695e4,...
    -0.14874640856724,...
    0.31806110878444];
cr=zeros(1,44);
cr(1:15)=1;
cr(16:35)=2;
cr(36:39)=3;
cr(40)=4;
cr(41:44)=6;
dr=[1,1,1,2,2,3,4,1,1,1,...
    2,2,3,4,4,5,7,9,10,11,...
    13,15,1,2,2,2,3,4,4,4,...
    5,6,6,7,9,9,9,9,9,10,...
    10,12,3,4,4,5,14,3,6,6,...
    6,3,3,3];
tr=[-0.5,0.875,1,0.5,0.75,0.375,1,4,6,12,...
    1,5,4,2,13,9,3,4,11,4,...
    13,1,7,1,9,10,10,3,7,10,...
    10,6,10,10,1,2,3,4,8,6,...
    9,8,16,22,23,23,10,50,44,46,...
    50,0,1,4];
alphar=20;
betar=[150,150,250,0.3,0.3];
gammar=[1.21,1.21,1.25];
epsilonr=1;
ar=3.5;
br=[0.85,0.95];
Br=0.2;
Cr=[28,32];
Dr=[700,800];
Ar=0.32;

nsat=[0.11670521452767e4,...
    -0.72421316703206e6,...
    -0.17073846940092e2,...
    0.12020824702470e5,...
    -0.32325550322333e7,...
    0.14915108613530e2,...
    -0.48232657361591e4,...
    0.40511340542057e6,...
    -0.23855557567849,...
    0.65017534844798e3];

Coef=finitedifferences(N,1,O/2+1);

%Derivadas

fphi0=@(delta,tau) phi0(delta,tau);
fphir=@(delta,tau) phir(delta,tau);
phi0d=@(delta,tau) df(delta,tau,fphi0,1);
phi0t=@(delta,tau) df(delta,tau,fphi0,2);
phird=@(delta,tau) df(delta,tau,fphir,1);
phirt=@(delta,tau) df(delta,tau,fphir,2);

p=@(rho,T) pressure(rho,T);
s=@(rho,T) entropy(rho,T);
h=@(rho,T) enthalpy(rho,T);
psat=@(T) pressuresat(T);
rhow=@(T) densitywater(T);
rhov=@(T) densityvapour(T);
sw=@(T) entropywater(T);
sv=@(T) entropyvapour(T);
hw=@(T) enthalpywater(T);
hv=@(T) enthalpyvapour(T);
frac=@(rho,T) fraction(rho,T);

    function phi=phi0(delta,tau)
    phi=0;
    for i=4:8
        phi=phi+n0(i)*log(1-exp(-gamma0(i-3)*tau));
    end
    phi=phi+log(delta)+n0(1)+n0(2)*tau+n0(3)*log(tau);
    end
    function phi=phir(delta,tau)
    phi=0;
    for i=1:7
       phi=phi+nr(i)*delta^dr(i)*tau^tr(i);
    end
    for i=8:51
       phi=phi+nr(i)*delta^dr(i)*tau^tr(i)*exp(-delta^cr(i-7));
    end
    for i=52:54
        phi=phi+nr(i)*delta^dr(i)*tau^tr(i)*exp(-alphar*(delta-epsilonr)^2-betar(i-51)*(tau-gammar(i-51))^2);
    end
    for i=55:56
        theta=(1-tau)+Ar*((delta-1)^2)^(1/(2*betar(i-51)));
        Delta=theta^2+Br*((delta-1)^2)^ar;
        psi=exp(-Cr(i-54)*(delta-1)^2-Dr(i-54)*(tau-1)^2);

        phi=phi+nr(i)*Delta^br(i-54)*delta*psi;
    end
    end
    function phi=phirdA(delta,tau)
    phi=0;
    for i=1:7
       phi=phi+nr(i)*dr(i)*delta^(dr(i)-1)*tau^tr(i);
    end
    for i=8:51
       phi=phi+nr(i)*exp(-delta^cr(i-7))*(delta^(dr(i)-1)*tau^tr(i)*(dr(i)-cr(i-7)*delta^cr(i-7)));
    end
    for i=52:54
        phi=phi+nr(i)*delta^dr(i)*tau^tr(i)*exp(-alphar*(delta-epsilonr)^2-betar(i-51)*(tau-gammar(i-51))^2)*...
            (dr(i)/delta-2*alphar*(delta-epsilonr));
    end
    for i=55:56
        theta=(1-tau)+Ar*((delta-1)^2)^(1/(2*betar(i-51)));
        Delta=theta^2+Br*((delta-1)^2)^ar;
        psi=exp(-Cr(i-54)*(delta-1)^2-Dr(i-54)*(tau-1)^2);
        psid=-2*Cr(i-54)*psi;
        Deltad=(delta-1)*(Ar*theta*(2/betar(i-51))*((delta-1)^2)^(1/(2*betar(i-51))-1)+2*Br*ar*((delta-1)^2)^(ar-1));

        phi=phi+nr(i)*(Delta^br(i-54)*(psi+delta*psid)+(br(i-54)*Delta^(br(i-54)-1)*Deltad)*delta*psi);
    end
    end
    function res=df(x,y,fun,pos)
        vec=zeros(1,N);
        if pos==1
            for i=1:N
                vec(i)=fun(x+dh*(i-O/2-1),y);
            end
        else
            for i=1:N
                vec(i)=fun(x,y+dh*(i-O/2-1));
            end
        end
        res=(vec*Coef)/dh;
    end
    function res=pressure(rho,T,checkfrac)
        if nargin<3
           checkfrac=true; 
        end
        if checkfrac
            fracf=fraction(rho,T);
            if (fracf<1) && (fracf>0)
                res=pressuresat(T);
                return;
            end
        end
        delta=rho/rhoc;
        tau=Tc/T;
        res=rho*R*T*(1+delta*phird(delta,tau));
    end
    function res=entropy(rho,T,checkfrac)
        if nargin<3
           checkfrac=true; 
        end
        if checkfrac
            fracf=fraction(rho,T);
            if (fracf<1) && (fracf>0)
                res=fracf*entropyvapour(T)+(1-fracf)*entropywater(T);
                return;
            end
        end
        delta=rho/rhoc;
        tau=Tc/T;
        res=R*(tau*(phi0t(delta,tau)+phirt(delta,tau))-phi0(delta,tau)-phir(delta,tau));
    end
    function res=enthalpy(rho,T,checkfrac)
        if nargin<3
           checkfrac=true; 
        end
        if checkfrac
            fracf=fraction(rho,T);
            if (fracf<1) && (fracf>0)
                res=fracf*enthalpyvapour(T)+(1-fracf)*enthalpywater(T);
                return;
            end
        end
        delta=rho/rhoc;
        tau=Tc/T;
        res=R*T*(1+tau*(phi0t(delta,tau)+phirt(delta,tau))+delta*phird(delta,tau));
    end
    function res=pressuresat(T)
        theta=T+nsat(9)/(T-nsat(10));
        A=theta^2+nsat(1)*theta+nsat(2);
        B=nsat(3)*theta^2+nsat(4)*theta+nsat(5);
        C=nsat(6)*theta^2+nsat(7)*theta+nsat(8);
        res=1e3*((2*C)/(-B+sqrt(B^2-4*A*C)))^4;
    end
    function res=densitywater(T)
        fun=@(rho) pressure(rho,T,false)-pressuresat(T);
        rho=1e3;
        options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-14,'StepTolerance',1e-14);
        res=fsolve(fun,rho,options);
    end
    function res=densityvapour(T)
        fun=@(rho) pressure(rho,T,false)-pressuresat(T);
        rho=1e-2;
        options = optimoptions('fsolve','Display','none','FunctionTolerance',1e-14,'StepTolerance',1e-14);
        res=fsolve(fun,rho,options);
    end
    function res=entropywater(T)
        res=entropy(densitywater(T),T,false);
    end
    function res=entropyvapour(T)
        res=entropy(densityvapour(T),T,false);
    end
    function res=enthalpywater(T)
        res=enthalpy(densitywater(T),T,false);
    end
    function res=enthalpyvapour(T)
        res=enthalpy(densityvapour(T),T,false);
    end
    function res=fraction(rho,T)
        rhowf=densitywater(T);
        rhowm=rhowf/rho;
        rhowv=rhowf/densityvapour(T);
        if (1-rhowm)<=(1-rhowv)
            res=1;
        elseif rhowv<=1
            res=0;
        else
            res=(1-rhowm)/(1-rhowv);
        end
    end
end

function [sysflag,rho,T,p,s,h,psat,rhow,rhov,sw,sv,hw,hv,frac]=Magia(vars,vals)
nullflag='indeterminado';
[fphi0,fphir,phi0d,phi0t,phird,phirt,p,s,h,psat,rhow,rhov,sw,sv,hw,hv,frac]=create_Functions();
funcs={p,s,h,psat,rhow,rhov,sw,sv,hw,hv,frac};
variables={'rho','T','p','s','h','psat','rhow','rhov','sw','sv','hw','hv','frac'};
magnitudes=[1e1,4e2,1e2,7e0,3e3,1e2,1e3,1e0,1e0,1e0,1e3,1e3,1e0];
X=zeros(1,13);

if strcmp(vars(1),variables(4))
    aux=vals(2);
    vars(1)=vars(2);
    vals(2)=vals(1);
    vars(2)=variables(4);
    vals(1)=aux;
end

n=1;
code=zeros(1,2);
for j=1:2
    for i=1:13
        if strcmp(vars{j},variables{i})
            code(j)=i;
        end
    end
end
if code(1)==code(2)
    sysflag='Variable repetida';
    X(1:13)=nullflag;
elseif (code(1)==1) && (code(2)==2)
    sysflag='Sistema compatible determinado';
    X(1)=vals(1);
    X(2)=vals(2);
    for i=3:5
        X(i)=funcs{i-2}(vals(1),vals(2));
    end
    for i=6:12
        X(i)=funcs{i-2}(vals(2));
    end
    X(13)=funcs{11}(vals(1),vals(2));
elseif (code(1)==2) && (code(2)==1)
    sysflag='Sistema compatible determinado';
    X(1)=vals(2);
    X(2)=vals(1);
    for i=3:5
        X(i)=funcs{i-2}(vals(2),vals(1));
    end
    for i=6:12
        X(i)=funcs{i-2}(vals(1));
    end
    X(13)=funcs{11}(vals(2),vals(1));
elseif (code(1)==2) && (6<=code(2))&&(code(2)<=12)
    fun=@(T) funcs{code(2)-2}(T)-vals(2);
    T=magnitudes(2);
    options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-14,'StepTolerance',1e-14);
    T=fsolve(fun,T,options);
    if T-vals(1)<=1e-4
        sysflag='Sistema indeterminado';
        X(1)=nullflag;
        X(2)=T;
        X(3:5)=nullflag;
        for i=6:12
            X(i)=funcs{i-2}(T);
        end
        X(13)=nullflag;
    else
        sysflag='Sistema incompatible';
        X(1:13)=nullflag;
    end
elseif (6<=code(1))&&(code(1)<=12) && (code(2)==2)
    fun=@(T) funcs{code(1)-2}(T)-vals(1);
    T=magnitudes(2);
    options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-14,'StepTolerance',1e-14);
    T=fsolve(fun,T,options);
    if T-vals(2)<=1e-4
        sysflag='Sistema indeterminado';
        X(1)=nullflag;
        X(2)=T;
        X(3:5)=nullflag;
        for i=6:12
            X(i)=funcs{i-2}(T);
        end
        X(13)=nullflag;
    else
        sysflag='Sistema incompatible';
        X(1:13)=nullflag;
    end
elseif ((6<=code(1))&&(code(1)<=12)) && ((6<=code(2))&&(code(2)<=12))
    fun=@(T) funcs{code(1)-2}(T)-vals(1);
    T=magnitudes(2);
    options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-14,'StepTolerance',1e-14);
    T=fsolve(fun,T,options);
    if funcs{code(2)-2}(T)-vals(2)<=1e-4
        sysflag='Sistema indeterminado';
        X(1)=nullflag;
        X(2)=T;
        X(3:5)=nullflag;
        for i=6:12
            X(i)=funcs{i-2}(T);
        end
        X(13)=nullflag;
    else
        sysflag='Sistema incompatible';
        X(1:13)=nullflag;
    end
elseif (code(1)==13 && not((0<=vals(1))&&(vals(1)<=1))) || (code(2)==13 && not((0<=vals(2))&&(vals(2)<=1)))
    sysflag='Sistema incompatible o indeterminado al no ser la fracción molar mayor que 0 o menor que 1: es necesario especificar otra variable en su lugar';
    X(1:13)=nullflag;
else
%     sysflag='Sistema compatible determinado';
%     fun=cell(1,2);
%     for i=1:2
%         if (code(i)==1) || (code(i)==2)
%             fun{i}=@(Y) vals(i)-Y(code(i));
%         elseif (6<=code(i))&&(code(i)<=12)
%             fun{i}=@(Y) vals(i)-funcs{code(i)-2}(Y(2));
%         else
%             fun{i}=@(Y) vals(i)-funcs{code(i)-2}(Y(1),Y(2));
%         end
%     end
%     system=@(Y) [fun{1}([Y(1)*magnitudes(code(1)),Y(2)*magnitudes(code(2))]),fun{2}([Y(1)*magnitudes(code(1)),Y(2)*magnitudes(code(2))])];
%     Y=[1,1];
%     options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-14,'StepTolerance',1e-14,'MaxFunctionEvaluations',1000);
%     Y=fsolve(system,Y,options);
%     Y=[Y(1)*magnitudes(code(1)),Y(2)*magnitudes(code(2))];
%     X(1)=Y(1);
%     X(2)=Y(2);
%     for i=3:5
%         X(i)=funcs{i-2}(Y(1),Y(2));
%     end
%     for i=6:12
%         X(i)=funcs{i-2}(Y(2));
%     end
%     X(13)=funcs{11}(Y(1),Y(2));
    sysflag='Sistema compatible determinado';
    for i=1:2
        if ((3<=code(i))&&(code(i)<=5)) || (code(i)==13)
            vfun=@(rho,T) vals(i)-funcs{code(i)-2}(rho,T);
            pick=2/i;
            break;
        end
    end
    
    if (code(1)==1)
        X(1)=vals(1);
        X(2)=tempsolver(vals(1));
    elseif (code(2)==1)
        X(1)=vals(2);
        X(2)=tempsolver(vals(2));
    elseif (code(1)==2)
        X(1)=denssolver(vals(1));
        X(2)=vals(1);
    elseif (code(2)==2)
        X(1)=denssolver(vals(2));
        X(2)=vals(2);
    else
        if (6<=code(pick))&&(code(pick)<=12)
            fun=@(rho) vals(pick)-funcs{code(pick)-2}(tempsolver(rho));
            dens=magnitudes(1);
            options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-14,'StepTolerance',1e-14);
            X(1)=fsolve(fun,dens,options);
            X(2)=tempsolver(X(1));
        else
            fun=@(rho) vals(pick)-funcs{code(pick)-2}(rho,tempsolver(rho));
            dens=magnitudes(1);
            options = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-14,'StepTolerance',1e-14);
            X(1)=fsolve(fun,dens,options);
            X(2)=tempsolver(X(1));
        end
    end
    
    for i=3:5
        X(i)=funcs{i-2}(X(1),X(2));
    end
    for i=6:12
        X(i)=funcs{i-2}(X(2));
    end
    X(13)=funcs{11}(X(1),X(2));
end

rho=X(1);
T=X(2);
p=X(3);
s=X(4);
h=X(5);
psat=X(6);
rhow=X(7);
rhov=X(8);
sw=X(9);
sv=X(10);
hw=X(11);
hv=X(12);
frac=X(13);

    function dens=denssolver(T)
        ofun=@(rho) vfun(rho,T);
        dens=magnitudes(1);
        soptions = optimoptions('fsolve','Display','none','FunctionTolerance',1e-14,'StepTolerance',1e-14);
        dens=fsolve(ofun,dens,soptions);
    end
    function temp=tempsolver(rho)
        ofun=@(T) vfun(rho,T);
        temp=magnitudes(2);
        soptions = optimoptions('fsolve','Display','none','FunctionTolerance',1e-14,'StepTolerance',1e-14);
        temp=fsolve(ofun,temp,soptions);
    end
end
function res=finitedifferences(N,d,e)
for i=1:N
    V(i)=i-e;
end
S=zeros(N);
b=zeros(N,1);
for i=1:N
    S(i,:)=V.^(i-1);
    b(i)=factorial(d)*delta_kronecker(i-1,d);
end
res=linsolve(S,b);
    function d=delta_kronecker(a,b)
        if a==b
            d=1;
        else
            d=0;
        end
    end
end