clear all
% Numero de individuos
Nind = 80;
% Numero de variables
Nvar = 2;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = [0.00001 0.00001];
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = [16 16];
%Numero de generaciones del genetico
Ngen = 80000;
% Inercia
W = 0.3;
% Cognitivo
c1 = 1.8;
% Social
c2 = 1.5;

rng('shuffle');

pob = crearPob(Li, Ls, Nind, Nvar);
vel = zeros(size(pob));

for j = 1:Nvar
%    vel(:, j) = ones(Nind, 1) * (Li(j) + Ls(j)) / 2;
end

FO = zeros(Nind, 1);
S = zeros(Nind, 1);
for i = 1:Nind
    FO(i) = funcionObjetivo(pob(i,:));
    g = restdes(pob(i,:));
    h = restigu(pob(i,:));
    S(i) = SVR(g, h);
end

pbest = pob;
pbestFO = FO;
pbestSVR = S;

gbest = pob(1,:);
gbestFO = FO(1);
gbestSVR = S(1);

for i =1:Nind
    if DEB(pbestFO(i), pbestSVR(i), gbestFO, gbestSVR)
        gbest = pbest(i,:);
        gbestFO = pbestFO(i);
        gbestSVR = pbestSVR(i);
    end
end


for p = 1:Ngen
    newvel = zeros(size(vel));
    newpob = zeros(size(pob));
    newFO = zeros(size(FO));
    newSVR = zeros(size(S));

    for i= 1:Nind
        for var = 1:Nvar
            newvel(i, var) = W*vel(i, var) + c1*rand()*(pbest(i,var) - ...
                pob(i,var)) + c2*rand()*(gbest(1,var) - pob(i,var));
            newpob(i, var) = ajustar(pob(i,var) + newvel(i, var), ...
                Li(var), Ls(var));
        end

        newFO(i) = funcionObjetivo(pob(i,:));
        g = restdes(pob(i,:));
        h = restigu(pob(i,:));
        newSVR(i) = SVR(g, h);

        if DEB(newFO(i), newSVR(i), pbestFO(i), pbestSVR(i))
            pbest(i,:) = newpob(i,:);
            pbestFO(i) = newFO(i);
            pbestSVR(i) = newSVR(i);
        end

    end

    for i =1:Nind
        if DEB(pbestFO(i), pbestSVR(i), gbestFO, gbestSVR)
            gbest(1,:) = pbest(i,:);
            gbestFO = pbestFO(i);
            gbestSVR = pbestSVR(i);
        end
    end

    FO = newFO;
    S = newSVR;
    vel = newvel;
    pob = newpob;

    gbest
    gbestFO
    gbestSVR
end


gbest
gbestFO
gbestSVR


function FO = funcionObjetivo(p)
    k1 = 0.09755988;
    k2 = 0.99*k1;
    k3 = 0.0391908;
    k4 = 0.9*k3;
    numerador = k2*p(1)*k3*p(2) + k1*p(1)+k2*p(2);
    denominador = (1+k3*p(1))*(1+k2*p(2)+k1*p(1))*(1+k4*p(2));
    FO = -numerador/denominador;
end

function g = restdes(p)
    g = zeros(1,1);
    g(1) = p(1).^0.5 + p(2).^0.5 - 4;
end

function h = restigu(p)
    h = 0;
end

function s = SVR(g, h)
    s = 0;
    for i = 1:size(g,2)
        s = s + max([0 g(i)]);
    end
    for i = 1:size(h,2)
        s = s + max([0 abs(h(i))-0.00001]);
    end
end

function mejor1= DEB(FO1, SVR1, FO2, SVR2)
FOtrial = FO1;
FOtarget = FO2;
Strial = SVR1;
Starget = SVR2;
if Starget == 0 && Strial == 0
    if FOtrial < FOtarget
        mejor1 = true;
        return
    else
        mejor1 = false;
        return
    end
elseif Starget ~= 0 && Strial ~= 0
    if Starget > Strial
        mejor1 = true;
        return
    else
        mejor1 = false;
        return
    end
else
    if Strial == 0
        mejor1 = true;
        return
    else
        mejor1 = false;
        return
    end
end
end

function pob = crearPob(li, ls, Nind, Nvar)
    pob = zeros(Nind, Nvar);
    for i=1:Nvar
        pob(:,i) = li(i) + (ls(i)- li(i))*rand(Nind, 1);
    end
end

function ajustado = ajustar(valor, li, ls)
    while true
        if valor < li
            valor = 2*li - valor;
        elseif valor > ls
            valor = 2*ls - valor;
        else
            break
        end
    end
    ajustado = valor;
end

