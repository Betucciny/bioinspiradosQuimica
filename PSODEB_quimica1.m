clear all
format longG
% Numero de individuos
Nind = 30;
% Numero de variables
Nvar = 8;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = [100 1000 1000 10 10 10 10 10];
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = [10000 10000 10000 1000 1000 1000 1000 1000];
%Numero de generaciones del genetico
Ngen = 500000;
% Inercia
W = 0.9;
% Cognitivo
c1 = 0.8;
% Social
c2 = 0.6;
turb = 0.1;

rng('shuffle');
colFo = Nvar + 1;

pob = crearPob(Li, Ls, Nind, Nvar);
vel = zeros(size(pob));

for j = 1:Nvar
   vel(:, j) = ones(Nind, 1) * (Li(j) + Ls(j)) / 2;
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
        gbest(1,:) = pbest(i,:);
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
            newvel(i, var) = W*vel(i, var) + c1*rand()*(pbest(i,var) - pob(i,var)) ...
                + c2*rand()*(gbest(1,var) - pob(i,var)) ...
                + turb*(2*rand()-1);
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
    FO = p(1) + p(2) + p(3);
end

function g = restdes(p)
    g = zeros(1,3);
    g(1) = 100*p(1) - p(1)*p(6) + 833.33252*p(4) - 83333.333;
    g(2) = p(2)*p(4) - p(2)*p(7) - 1250*p(4) + 1250*p(5);
    g(3) = p(3)*p(5) - p(3)*p(8) - 2500*p(5) + 1250000;
end

function h = restigu(p)
    h = zeros(1,3);
    h(1) = 0.0025*(p(4) + p(6)) - 1;
    h(2) = 0.0025*(-p(4) + p(5) + p(7)) - 1;
    h(3) =  0.01*(-p(5) + p(8)) - 1; 
end

function s = SVR(g, h)
    s = 0;
    for i = 1:size(g,2)
        s = s + max([0 g(i)]);
    end
    for i = 1:size(h,2)
        s = s + max([0 abs(h(i))-0.0001]);
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
