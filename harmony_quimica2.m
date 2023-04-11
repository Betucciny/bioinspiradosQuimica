clear all
format longG
% Ireaciones
g_max = 1200000;
% Harmony numbers
NH = 80;
% Numero de variables
Nvar = 6;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = [0 0 0 0 0.00001 0.00001];
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = [1 1 1 1 16 16];
% Acceptance rate
rac = 0.9;
% Pitch adjust rate
rpa = 0.8;
% Intelligent acceptance
ria = 0.5;
% Factor para cambiar bw
a = 1;

harmonies = crearPob(Li, Ls, NH, Nvar);
FO = zeros([1 NH]);
S = zeros([1 NH]);

% Calculo de las primeras funciones objetivo
for i = 1:NH
    FO(i) = funcionObjetivo(harmonies(i,:));
    S(i) = SVR(restdes(harmonies(i,:)), restigu(harmonies(i,:)));
end


% Encontremos el peor de nuestras harmonias
windex = indicePeor(FO, S, NH);
bindex = indiceMejor(FO, S, NH);

for g = 1:g_max
    bw = (Ls - Li)/g.^a;
    newH = zeros([1 Nvar]);
    for v = 1:Nvar
        if rand() < rac
            index = randi(NH);
            if rand() < rpa
                newH(v) = ajustar(harmonies(index, v) + bw(v) * ...
                    (-1 + 2*rand()), Li(v), Ls(v));
            else
                if rand() < ria
                    newH(v) = harmonies(bindex, v);
                else
                    newH(v) = harmonies(index, v);
                end
            end
        else
            newH(v) = crearPob(Li(v), Ls(v), 1, 1);
        end
    end
    newFO = funcionObjetivo(newH);
    newS = SVR(restdes(newH), restigu(newH));
    
    if DEB(newFO, newS, FO(windex), S(windex))
        harmonies(windex, :) = newH;
        FO(windex) = newFO;
        S(windex) = newS;
        windex = indicePeor(FO, S, NH);
        bindex = indiceMejor(FO, S, NH);
    end
    harmonies(bindex, :)
    FO(bindex)
    S(bindex)
end

harmonies(bindex, :)
FO(bindex)
S(bindex)


function FO = funcionObjetivo(p)
    FO = -p(4);
end

function g = restdes(p)
    g = zeros(1,1);
    g(1) = p(5).^0.5 + p(6).^0.5 - 4;
end

function h = restigu(p)
    k1 = 0.09755988;
    k2 = 0.99*k1;
    k3 = 0.0391908;
    k4 = 0.9*k3;
    h = zeros(1,4);
    h(1) = p(1) + k1*p(2)*p(5) - 1;
    h(2) = p(2) - p(1) + k2*p(2)*p(6);
    h(3) = p(3) + p(1) + k3*p(3)*p(5) - 1; 
    h(4) = p(4) - p(3) + p(2) - p(1) + k4*p(4)*p(6);
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


function windex = indicePeor(FO, S, NH)
windex = 1;
for i = 1:NH
    if DEB(FO(windex),S(windex), FO(i), S(i))
        windex = i;
    end
end
end

function bindex = indiceMejor(FO, S, NH)
bindex = 1;
for i = 1:NH
    if DEB(FO(i), S(i), FO(bindex), S(bindex))
        bindex = i;
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