Tema 1 ASC - Optimizari

Descriere generala
==================

In aceasta tema am implementat urmatorul calcul cu matrice:

C = A^T * B  
D = C * C^T  
y = 0  
pentru fiecare linie i din C:  
    y = y + D * linie(C, i)  
la final:  
    y = y + x  

Am facut cele 3 variante cerute in enunt:

- blas - varianta care foloseste functii BLAS;
- neopt - varianta simpla, scrisa manual, fara optimizari;
- opt_m - varianta optimizata manual, fara BLAS.

Toate cele 3 variante tin cont de faptul ca matricea D este simetrica. Nu am folosit algoritmi de inmultire cu alta complexitate, cum ar fi Strassen. Am pastrat inmultirea clasica de matrice, asa cum cere enuntul.

Structura arhivei
=================

Arhiva are urmatoarea structura:

src/
  solver_blas.c
  solver_neopt.c
  solver_opt.c

cache/
  blas.cache
  neopt.cache
  opt_m.cache

memory/
  blas.memory
  neopt.memory
  opt_m.memory

grafice/
  grafic1.png
  grafic2.png
  ...
  grafic13.png

README
LLMprompts

Observatie matematica folosita
==============================

Pentru partea in care se aduna rezultatele pentru fiecare linie din C, am folosit o observatie simpla.

Expresia initiala este:

y = suma_i D * linie(C, i) + x

Aceasta se poate scrie si ca:

y = D * s + x

unde s este vectorul care contine sumele pe coloane ale matricei C.

Cu alte cuvinte, in loc sa inmultesc D cu fiecare linie din C separat, pot sa calculez o singura data vectorul s si apoi sa fac un singur produs D * s. Aceasta idee am folosit-o mai ales in varianta blas.

Pentru varianta opt_m am mai folosit si faptul ca:

D * v = (C * C^T) * v = C * (C^T * v)

Astfel nu mai trebuie sa construiesc explicit matricea D. In loc de asta, folosesc doua produse matrice-vector. Acest lucru ajuta deoarece se reduce memoria folosita si se evita parcurgerea unei matrice suplimentare mari.

Varianta neopt
==============

Fisier: solver_neopt.c

Varianta neopt este implementarea de baza. Am scris calculele direct, cu bucle for, fara BLAS si fara optimizari speciale.

Pasii principali sunt:

1. initializez vectorul y cu 0;
2. calculez C = A^T * B;
3. calculez D = C * C^T;
4. pentru D calculez doar triunghiul superior si copiez valoarea si in pozitia simetrica;
5. pentru fiecare linie din C calculez D * linie(C, i);
6. adun rezultatul in y;
7. la final adaug vectorul x.

Aceasta varianta este mai usor de urmarit, dar nu este eficienta. De exemplu, la calculul lui C apar accesari pe coloane in matricea A, iar acest lucru nu este foarte bun pentru cache. Am pastrat varianta cat mai simpla pentru a avea un punct de comparatie cu opt_m si blas.

Varianta blas
=============

Fisier: solver_blas.c

In varianta blas am folosit functii din OpenBLAS pentru operatiile principale. Ideea a fost sa folosesc rutinele deja optimizate din biblioteca, in loc sa fac manual inmultirile mari.

Functiile folosite au fost:

- cblas_dgemm pentru C = A^T * B;
- cblas_dsyrk pentru D = C * C^T;
- cblas_dgemv pentru calculul vectorului cu sumele pe coloane;
- cblas_dsymv pentru produsul cu matricea simetrica D;
- cblas_daxpy pentru adunarea vectorului x.

Am folosit dsyrk si dsymv pentru ca D este simetrica. Astfel nu tratez D ca pe o matrice generala si folosesc o functie potrivita pentru cazul acesta.

Aceasta varianta este cea mai rapida in testele mele, ceea ce era de asteptat. OpenBLAS foloseste optimizari mult mai bune decat o implementare manuala simpla: blocking, vectorizare si acces mai bun la cache.

Varianta opt_m
==============

Fisier: solver_opt.c

Varianta opt_m porneste de la varianta neopt, dar am modificat modul in care sunt accesate datele. Scopul a fost sa reduc numarul de accesari lente in memorie si sa folosesc mai bine cache-ul.

Optimizari folosite:

- am transpus A intr-o matrice separata AT, ca accesul pentru A^T sa fie liniar;
- am folosit inmultire blocata pentru calculul lui C;
- dimensiunea blocului este 40, deoarece in enunt se spune ca N este multiplu de 40;
- am folosit pointeri in unele bucle, pentru a evita indexari repetate;
- am facut unrolling pe 4 elemente in buclele interioare;
- am folosit pragma-uri GCC pentru unrolling si vectorizare;
- nu mai construiesc explicit matricea D;
- folosesc relatia D * v = C * (C^T * v);
- am transpus C intr-o matrice CT pentru a avea acces mai bun la memorie la produsul C^T * v.

Cea mai importanta schimbare este faptul ca nu mai aloc si nu mai parcurg explicit D. Matricea D ar avea dimensiune N x N, deci pentru valori mari ale lui N devine costisitor sa fie construita si citita de mai multe ori. Prin forma C * (C^T * v), calculele raman corecte, dar accesul la memorie este mai bun.

Compilare si rulare
===================

Compilarea se face cu:

make

Rularea variantelor:

./tema1_blas input  
./tema1_neopt input  
./tema1_opt_m input  

Pentru rularea pe cluster am folosit partitia haswell. Pentru varianta blas am rulat in containerul pus la dispozitie, folosind apptainer.

Exemplu:

apptainer exec optimizari.sif bash test.sh

Verificare cu Valgrind memcheck
===============================

Pentru verificarea acceselor la memorie am rulat:

valgrind --tool=memcheck --leak-check=full ./tema1_neopt ../input/input_valgrind  
valgrind --tool=memcheck --leak-check=full ./tema1_blas ../input/input_valgrind  
valgrind --tool=memcheck --leak-check=full ./tema1_opt_m ../input/input_valgrind  

Rezultatele sunt salvate in directorul memory.

| Varianta | Alocari | Eliberari | Memorie ramasa la final | Erori |
|---|---:|---:|---:|---:|
| neopt | 9 | 9 | 0 bytes | 0 |
| blas | 31 | 31 | 0 bytes | 0 |
| opt_m | 11 | 11 | 0 bytes | 0 |

Pentru toate cele 3 variante, memoria alocata a fost eliberata corect. Valgrind nu a raportat memory leaks sau accesari invalide.

Analiza Cachegrind
==================

Pentru analiza cache am folosit:

valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes ./tema1_neopt ../input/input_valgrind  
valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes ./tema1_blas ../input/input_valgrind  
valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes ./tema1_opt_m ../input/input_valgrind  

Inputul de valgrind foloseste N = 400.

| Varianta | I refs | D refs | D1 misses | LLd misses | Branches | Mispredicts | D1 miss rate | Branch miss rate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| neopt | 4,195,263,931 | 2,388,625,383 | 136,444,557 | 80,871 | 130,102,230 | 405,224 | 5.7% | 0.3% |
| blas | 83,718,145 | 31,995,024 | 1,710,411 | 137,482 | 4,306,380 | 90,317 | 5.3% | 2.1% |
| opt_m | 1,914,039,678 | 502,136,455 | 17,201,137 | 100,928 | 54,945,160 | 1,966,975 | 3.4% | 3.6% |

Interpretare Cachegrind
=======================

La neopt se vede clar ca sunt foarte multe instructiuni si foarte multe accesari de date. Varianta este corecta, dar face multe parcurgeri ale matricelor si nu are o localitate foarte buna. Din cauza asta are si cel mai mare numar de D1 misses.

La blas, numarul de instructiuni si de referinte la date este mult mai mic. Asta se intampla pentru ca OpenBLAS foloseste implementari optimizate. Chiar daca rata de branch miss este mai mare decat la neopt, numarul total de branch-uri este mult mai mic, deci impactul nu este asa mare.

La opt_m se vede imbunatatirea fata de neopt. Numarul de D refs scade de la aproximativ 2.38 miliarde la aproximativ 502 milioane. Si D1 misses scade foarte mult, de la aproximativ 136 milioane la 17 milioane. Pentru mine asta a fost cel mai clar semn ca optimizarile pe accesul la memorie au avut efect. Totusi, opt_m nu ajunge la blas, pentru ca BLAS foloseste kernel-uri mult mai optimizate.

Timpi de executie pe haswell
============================

Pentru comparatia principala am rulat cele 3 variante pe haswell. Timpii sunt in secunde.

| N | neopt | blas | opt_m |
|---:|---:|---:|---:|
| 400 | 0.51 | 0.12 | 0.18 |
| 800 | 4.05 | 0.08 | 1.20 |
| 1200 | 13.80 | 0.28 | 4.34 |

Pentru N = 1200, varianta opt_m ruleaza in aproximativ 4.34 secunde, deci se incadreaza sub limita de 9 secunde si si sub limita de 6 secunde.

Se observa ca blas este mult mai rapida decat celelalte variante. Neopt este cea mai lenta, iar opt_m este intre cele doua. Pentru N = 1200, opt_m este de peste 3 ori mai rapida decat neopt.

Am mai rulat opt_m si pentru dimensiuni mai mari, ca sa vad cum scaleaza:

| N | opt_m haswell |
|---:|---:|
| 400 | 0.18 |
| 800 | 1.20 |
| 1200 | 4.34 |
| 1400 | 7.01 |
| 1600 | 10.10 |
| 1800 | 14.72 |

Timpul creste destul de mult cand creste N, ceea ce este normal pentru inmultiri de matrice. Cresterea nu este perfect cubica, pentru ca apar si efecte de cache, zgomot de rulare pe cluster si diferente date de optimizarile compilatorului.

Comparatie haswell vs ucsx
==========================

Pentru bonus am rulat varianta opt_m si pe ucsx. Valorile de mai jos sunt din rularile folosite pentru graficele comparative.

| N | opt_m haswell | opt_m ucsx | Speedup ucsx fata de haswell |
|---:|---:|---:|---:|
| 400 | 0.16 | 0.15 | 1.07x |
| 800 | 1.25 | 0.94 | 1.33x |
| 1200 | 4.14 | 3.33 | 1.24x |
| 1400 | 6.65 | 6.07 | 1.10x |
| 1600 | 10.45 | 10.64 | 0.98x |
| 1800 | 14.41 | 14.09 | 1.02x |

In majoritatea cazurilor, ucsx a iesit putin mai rapid decat haswell. Diferenta cea mai mare este la N = 800, unde speedup-ul este 1.33x. La N = 1600, ucsx este foarte putin mai lent, dar diferenta este mica. Probabil aici conteaza mai mult memoria si cache-ul, nu doar procesorul.

Per total, rezultatele pe ucsx sunt apropiate de cele de pe haswell, dar putin mai bune pentru mai multe valori ale lui N.

Grafice incluse
===============

Am pus graficele in directorul grafice. Ele sunt facute pentru a compara variantele intre ele si pentru a vedea mai clar diferentele de performanta.

Primele grafice compara timpii de rulare pentru neopt, blas si opt_m. Din ele se vede imediat ca neopt creste cel mai mult odata cu N. Blas ramane foarte jos pe grafic, pentru ca timpii sunt mult mai mici. Opt_m este vizibil mai buna decat neopt, dar tot ramane peste blas.

Un alt set de grafice arata cum se comporta opt_m pentru valori mai mari ale lui N. Aici se vede cresterea timpului pe masura ce N devine mai mare. Graficul este util pentru ca arata ca optimizarea nu functioneaza doar pentru testele oficiale 400, 800 si 1200, ci si pentru valori suplimentare.

Am inclus si grafice pentru comparatia haswell vs ucsx. In acestea se vede ca ucsx este de obicei putin mai rapid pentru opt_m, dar nu cu aceeasi diferenta la toate dimensiunile. De exemplu, pentru N = 800 diferenta este mai vizibila, iar pentru N = 1600 cele doua arhitecturi sunt aproape la egalitate.

Graficele cu D1 misses si LLd misses sunt importante pentru partea de cache. Se vede ca neopt are cele mai multe miss-uri in L1, ceea ce explica timpul mai mare. Opt_m reduce mult aceste miss-uri, iar blas ramane cea mai eficienta varianta. La LLd misses diferentele sunt mai mici, deci problema principala se vede mai ales la nivelul L1 data cache.

Am pus si grafice pentru D refs, rata de miss L1 si branch misprediction. La D refs, diferenta dintre neopt si opt_m este foarte clara: opt_m acceseaza mult mai putine date. La branch misprediction, opt_m are o rata mai mare, dar acest lucru nu strica foarte mult performanta, pentru ca reducerea acceselor la memorie conteaza mai mult in cazul acesta.

Ultimul grafic important este cel cu speedup-ul ucsx fata de haswell. Acolo se vede ca ucsx aduce un castig pentru unele valori ale lui N, mai ales la N = 800 si N = 1200. Pentru valorile mari, speedup-ul ajunge aproape de 1, deci cele doua arhitecturi se comporta asemanator.

Pe scurt, graficele confirma ce se vede si in tabele: blas este cea mai rapida varianta, opt_m imbunatateste clar neopt, iar diferentele dintre haswell si ucsx depind destul de mult de dimensiunea problemei.

Efectul optimizarilor manuale
=============================

Optimizarea cea mai importanta a fost imbunatatirea accesului la memorie. In neopt, unele accesari sunt pe coloane sau repeta multe zone mari din matrice. In opt_m, prin transpunerea lui A si C, mai multe accesari devin liniare. Astfel procesorul foloseste mai bine cache-ul si prefetching-ul.

Blocking-ul ajuta la calculul lui C, pentru ca se lucreaza pe bucati mai mici din matrice. Aceste bucati au sanse mai mari sa ramana in cache in timp ce sunt refolosite.

Unrolling-ul reduce overhead-ul din bucle si poate ajuta compilatorul sa genereze cod mai bun. Directivele GCC au fost folosite tot pentru a incuraja vectorizarea si unrolling-ul.

Evitarea matricei D ajuta mult, deoarece D este o matrice mare. Daca o construiesc explicit, trebuie sa aloc mai multa memorie si sa o citesc de multe ori. In opt_m am folosit forma C * (C^T * v), ceea ce reduce numarul de accesari la date.

Rezultatele confirma acest lucru: opt_m are mult mai putine D refs si D1 misses decat neopt, iar timpul pentru N = 1200 scade de la 13.80 secunde la 4.34 secunde.

Detalii despre rularea pe arhitecturi
=====================================

Rularile principale au fost facute pe partitia haswell, asa cum cere enuntul. Pentru bonus am folosit si coada ucsx.

Informatiile despre procesor, cache si compilator au fost obtinute cu:

lscpu  
lscpu -C  
gcc --version  

Pe ambele arhitecturi am folosit acelasi cod sursa si aceleasi tipuri de input. Astfel, diferentele din rezultate vin in principal din arhitectura si din conditiile de rulare de pe noduri.

Utilizarea uneltelor LLM
========================

Am folosit ChatGPT ca sa clarific unele idei si ca sa verific directia de implementare. Nu am folosit raspunsurile ca solutie finala copiata direct, ci mai mult ca punct de plecare.

Prompturile complete sunt puse in fisierul LLMprompts.

Ordinea in care am folosit prompturile a fost:

1. mai intai am cerut o explicatie pentru expresia matematica din enunt;
2. apoi am intrebat cum se poate folosi faptul ca D este simetrica;
3. dupa aceea am cerut o directie pentru varianta blas si functiile potrivite;
4. apoi am verificat ideea pentru varianta neopt;
5. la final am cerut sugestii pentru optimizari manuale, dar codul a fost ajustat prin teste.

Modificarile facute manual au fost mai ales in opt_m: ordinea buclelor, transpunerea matricelor, blocking-ul, unrolling-ul si eliminarea matricei D. Am testat mai multe variante si am pastrat varianta care a dat timp mai bun pe haswell.

Concluzii
=========

Rezultatele obtinute sunt cele asteptate. Varianta neopt este cea mai lenta, deoarece face multe accesari la memorie si nu foloseste optimizari de cache. Varianta blas este cea mai rapida, pentru ca foloseste OpenBLAS, care are implementari foarte optimizate.

Varianta opt_m este compromisul dintre cele doua. Nu ajunge la viteza BLAS, dar imbunatateste clar varianta neopt. Pentru N = 1200, timpul scade de la 13.80 secunde la 4.34 secunde, ceea ce arata ca optimizarile manuale au avut efect.

Din Cachegrind se vede ca timpul nu depinde doar de numarul de operatii, ci si de modul in care sunt accesate datele. Reducerea D refs si D1 misses a avut cel mai mare impact asupra performantei variantei opt_m.