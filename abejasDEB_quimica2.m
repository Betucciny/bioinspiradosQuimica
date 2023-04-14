clear all
format longG
% Numero de fuentes
Nf = 100;
% Numero de variables
Nvar = 2;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = [0.00001 0.00001];
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = [16 16];
%Numero de iteraciones del genetico
Niter = 30000;

limite = round(Niter / (2*Nf));

rng('shuffle');
fuentes = crearPob(Li, Ls, Nf, Nvar);
colFO = Nvar + 1;

FO = zeros(Nf, 1);
S = zeros(Nf, 1);
L = zeros(Nf, 1);
for i = 1:Nf
    FO(i) = funcionObjetivo(fuentes(i,:));
    g = restdes(fuentes(i,:));
    h = restigu(fuentes(i,:));
    S(i) = SVR(g, h);
end

bindex = indiceMejor(FO, S, Nf);

for p=1:Niter
    asignacion = randperm(Nf);
    for i=1:Nf
        posible = zeros(1,Nvar);
        k = asignacion(i);
        for j=1:Nvar
            phi = -1 + 2 * rand();
            posible(j) = ajustar(fuentes(i,j) + phi*(fuentes(i,j)-fuentes(k,j)), Li(j), Ls(j));
        end
        FOposible = funcionObjetivo(posible);
        g = restdes(posible);
        h = restdes(posible);
        Sposible = SVR(g, h);
        if DEB(FOposible, Sposible, FO(i), S(i))
            fuentes(i,:) = posible;
            FO(i) = FOposible;
            S(i) = Sposible;
            L(i) = 0;
        else
            L(i) = L(i) + 1;
        end
    end

    %     Por torneo tomamos la mejor solucion
    bindex = indiceMejor(FO, S, Nf);

    for i=1:Nf
        posible = zeros(1,Nvar);
        k = bindex;
        for j=1:Nvar
            phi = -1 + 2 * rand();
            posible(j) = ajustar(fuentes(i,j) + phi*(fuentes(i,j)-fuentes(k,j)), Li(j), Ls(j));
        end
        FOposible = funcionObjetivo(posible);
        g = restdes(posible);
        h = restdes(posible);
        Sposible = SVR(g, h);
        if DEB(FOposible, Sposible, FO(i), S(i))
            fuentes(i,:) = posible;
            FO(i) = FOposible;
            S(i) = Sposible;
            L(i) = 0;
        else
            L(i) = L(i) + 1;
        end
    end
    bindex = indiceMejor(FO, S, Nf);

%     Vuelo inteligente
    for i=1:Nf
        if L(i) > limite
            L(i) = 0;
            if i == bindex
                continue
            end
            k = randi([1, Nf]);
            temp = zeros([1, Nvar]);
            for j = 1:Nvar
                phi = -1 + 2 * rand();
                temp(j) = fuentes(i,j) + phi * (fuentes(k,j) - fuentes(i,j)) ...
                + (1-phi) * (fuentes(bindex,j) - fuentes(i,j));
            end
            fuentes(i,:) = ajustar(temp, Li, Ls);
            g = restdes(fuentes(i,:));
            h = restigu(fuentes(i,:));
            S(i) = SVR(g, h);
        end
    end
    
    fuentes(bindex,:)
    FO(bindex)
    S(bindex)
end

fuentes(bindex,:)
FO(bindex)
S(bindex)


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