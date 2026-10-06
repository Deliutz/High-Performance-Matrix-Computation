LLMprompts:


1. Spune-mi exact ce comenzi trebuie sa rulez ca sa compilez si sa testez varianta neopt si blas in container, pe haswell

R: srun -p haswell --pty bash
cd ~/asc-public/assignments/1-optimizari
apptainer shell optimizari.sif
cd src
make clean
make
cd ..
./src/tema1_neopt input/input
./src/compare out1 references/out1
./src/compare out2 references/out2
./src/compare out3 references/out3
./src/tema1_blas input/input
./src/compare out1 references/out1
./src/compare out2 references/out2
./src/compare out3 references/out3

2. Nu vreau sa calculezi direct D. Vreau sa folosesti relatia D * v = C * (C^T * v), astfel incat in bucla mare sa folosesti doua inmultiri matrice-vector.

R:Da, asta este o optimizare buna si pastreaza acelasi rezultat matematic.

Daca ai D = C * C^T, atunci pentru orice vector v poti scrie:

D * v = (C * C^T) * v = C * (C^T * v)

Asta inseamna ca nu mai e nevoie sa construiesti explicit matricea D, care ocupa memorie multa si adauga accesari inutile. In loc de asta, in fiecare pas din bucla mare, pentru vectorul v = linie(C, i), faci doua operatii:
mai intai calculezi un vector intermediar t1 = C^T * v, apoi calculezi t2 = C * t1, iar la final aduni t2 in y.

Practic, partea:

y += D * linie(C, i)

devine:

t1 = C^T * linie(C, i)
t2 = C * t1
y += t2

Avantajul este ca:

elimini complet alocarea si calculul lui D
folosesti doar matrici si vectori deja existenti
accesul la memorie este mai bun, mai ales daca ai si C^T precomputat
codul devine mai prietenos pentru optimizari de tip cache, pointeri, unrolling

Ca idee de implementare, dupa ce calculezi C, poti sa faci si CT o singura data. Apoi, pentru fiecare linie i din C:

iei row = C + i * N
calculezi t1[j] = sum_k CT[j][k] * row[k]
calculezi t2[j] = sum_k C[j][k] * t1[k]
faci y[j] += t2[j]

E important de retinut ca nu schimbi complexitatea asimptotica a partii din bucla mare: tot faci doua produse matrice-vector, deci ramai in zona corecta pentru cerinta, doar ca varianta asta este mult mai eficienta in practica decat sa formezi explicit D si apoi sa inmultesti cu vectorul.

3. Uita-te la valorile din cachegrind si explica-mi pe scurt ce inseamna I refs, D refs, D1 misses, LLd misses, Branches si Mispredicts si cum se vad diferentele intre neopt, blas si opt_m.

R:Pe scurt, metricele din cachegrind spun cat de mult “munceste” programul si cat de bine se intelege cu memoria si cu predictorul de ramificatii.

I refs inseamna instruction references, adica numarul de instructiuni citite/executate. Cu cat valoarea asta e mai mare, cu atat programul face mai multa treaba la nivel de cod. In general, varianta neopt are cele mai multe I refs, pentru ca face multe operatii manuale si multe iteratii, in timp ce BLAS are mult mai putine pentru ca muta mare parte din munca in rutine foarte optimizate.

D refs inseamna data references, adica numarul de accesari de date din memorie, fie citiri, fie scrieri. Asta e foarte important la inmultiri de matrici, pentru ca aici se pierde mult timp daca accesul la memorie nu este cache-friendly. La tine se vede foarte clar: neopt are enorm de multe D refs, opt_m are mult mai putine, iar blas are de departe cele mai putine. Asta sugereaza ca BLAS foloseste mult mai eficient datele deja aduse in cache.

D1 misses reprezinta cate accesari de date nu au fost gasite in cache-ul de nivel 1 pentru date. Cu alte cuvinte, procesorul a cautat rapid in L1 si nu a gasit acolo ce voia. In cazul tau, neopt are foarte multe astfel de ratari, opt_m reduce mult numarul lor, iar blas are cele mai putine. Asta e exact genul de comportament pe care vrei sa-l vezi dupa optimizari.

LLd misses inseamna last-level data cache misses, adica ratari si in nivelurile mai mari de cache pentru date, care de obicei inseamna acces mai scump, spre RAM. Aici diferentele nu mai sunt la fel de mari ca la D1 misses, dar tot sunt relevante. Faptul ca blas si opt_m tin valorile jos inseamna ca localitatea in memorie este mai buna decat in neopt.

Branches inseamna numarul total de ramificatii executate, de tip if, for, while si altele similare. Un numar mare de branches apare de obicei cand ai multe bucle si controale de flux. Varianta neopt are multe branches pentru ca are multe bucle clasice, opt_m are tot multe, dar mai putine, iar blas are foarte putine la nivelul codului tau, fiindca o mare parte din logica este ascunsa in apelurile BLAS.

Mispredicts inseamna de cate ori predictorul de ramificatii al procesorului a ghicit gresit directia unei ramificatii. Asta costa timp, pentru ca procesorul trebuie sa arunce o parte din executia speculativa si sa reia corect. La tine, blas si opt_m au rate de misprediction mai mari procentual decat neopt, dar asta nu inseamna automat ca sunt mai lente. Conteaza mult mai mult volumul total de lucru si accesul la memorie.

Daca ne uitam comparativ la cele trei variante, imaginea generala este destul de clara. Neopt face cel mai mult lucru brut: multe instructiuni, foarte multe accesari de date si multe miss-uri in cache, deci este varianta cea mai slaba din punct de vedere al performantei. Opt_m reduce mult numarul de accesari si de ratari in cache fata de neopt, ceea ce arata ca optimizarile manuale chiar au efect. Blas este cea mai eficienta, pentru ca are cele mai mici valori la majoritatea metricalor importante, mai ales la D refs si D1 misses, ceea ce explica de ce in practica merge cel mai bine.

Pe valorile tale se vede foarte bine asta:

neopt: D refs foarte mari, D1 misses foarte mari
opt_m: scadere mare fata de neopt, deci optimizarile manuale au imbunatatit localitatea
blas: cel mai mic volum de accesari si cele mai putine miss-uri, deci cea mai buna utilizare a cache-ului

Concluzia scurta pentru README ar fi ca diferentele dintre variante se vad mai ales in metricele de memorie. Neopt pierde mult timp din cauza acceselor ineficiente, opt_m imbunatateste semnificativ comportamentul in cache, iar blas obtine cele mai bune rezultate pentru ca foloseste rutine foarte bine optimizate pentru arhitectura.

4. Ce functii BLAS s-ar potrivi pentru fiecare pas din formula?

R: Pentru formula din tema, functiile BLAS potrivite ar fi astea:
Pentru C = A^T * B, cea mai potrivita este dgemm, pentru ca ai o inmultire matrice-matrice generala. Aici setezi A ca transpusa si B normal.
Pentru D = C * C^T, tot o inmultire matrice-matrice este, dar pentru ca rezultatul D este simetric, se potriveste foarte bine dsyrk. Asta e mai potrivita decat un dgemm normal, tocmai fiindca exploateaza simetria si calculeaza doar jumatatea relevanta din matrice.
Pentru partea cu bucla:
y += D * linie(C, i)
aici ai doua variante.
Daca respecti formula literal si lucrezi cu D deja calculata, atunci pentru fiecare iteratie poti folosi dsymv pentru produsul matrice simetrica - vector, fiindca D este simetrica. Apoi pentru acumularea in y folosesti daxpy.
Daca vrei varianta mai eleganta, fara sa construiesti explicit D, folosesti relatia:
D * v = C * (C^T * v)
si atunci, in fiecare pas:


t1 = C^T * v cu dgemv


t2 = C * t1 cu dgemv


y += t2 cu daxpy


Pentru y = y + x, functia potrivita este daxpy, pentru ca exact asta face: aduna un vector scalat la alt vector, iar aici scalarul este 1.0.
Deci, pe scurt:


C = A^T * B → dgemm


D = C * C^T → dsyrk


D * vector → dsymv


C^T * vector sau C * vector → dgemv


y += ceva sau y = y + x → daxpy


O mapare foarte buna pentru tema ar fi:


varianta BLAS clasica: dgemm + dsyrk + dsymv + daxpy


varianta BLAS fara construire explicita a lui D: dgemm + dgemv + dgemv + daxpy


A doua e adesea mai eleganta conceptual, dar prima se potriveste mai direct cu cerinta si foloseste explicit faptul ca D este simetrica.

5. Ce grafice ar fi bune ca sa compar clar neopt, blas si opt_m, astfel incat sa se vada si diferentele de timp, si partea de cache, dar si comparatia intre haswell si ucsx?

R: Cele mai bune ar fi cateva grafice simple, fiecare cu un rol clar.

Primul si cel mai important ar fi un grafic de tip linie pentru timp de executie, cu N pe axa Ox si timpul in secunde pe Oy. Aici poti pune cele trei variante pentru valorile 400, 800 si 1200, iar pentru 1600, 1800 si 2000 doar blas si opt_m. Graficul asta arata imediat diferenta reala de performanta.

Pentru comparatia intre arhitecturi, as face inca un grafic de timp, tot cu N pe Ox, dar separat pe haswell si ucsx. De exemplu, fie doua grafice separate, fie acelasi grafic cu serii de forma blas_haswell, blas_ucsx, opt_m_haswell, opt_m_ucsx. Asa se vede direct daca una dintre arhitecturi scaleaza mai bine.

Pe partea de cache, un grafic bun este unul pentru D refs si I refs, pentru ca arata cat de mult “munceste” fiecare varianta. Aici merg bine grafice de tip bara, pentru ca ai doar cateva variante de comparat. Se vede usor ca neopt are de obicei mult mai multe accesari, iar blas mult mai putine.

Apoi as face un grafic separat pentru D1 misses si LLd misses. Asta e important fiindca spune cat de bine foloseste fiecare varianta cache-ul. Daca vrei sa fie mai clar, poti pune scala logaritmica pe Oy, pentru ca diferentele dintre variante pot fi foarte mari.

Pentru branch prediction, un grafic bun este unul cu Branches si Mispredicts, sau chiar mai bine cu Mispred rate. Acolo se vede daca optimizarile au redus sau nu costul ramificatiilor. E util mai ales cand explici de ce o varianta poate avea mai multe mispredicts procentual, dar totusi sa fie mai rapida.

In README as zice ca minimul bun ar fi:

un grafic cu timpii de executie
un grafic cu comparatia haswell vs ucsx
un grafic cu D refs / I refs
un grafic cu D1 misses / LLd misses
un grafic cu Branches / Mispred rate

Daca vrei sa fie si mai curat, poti organiza asa:

grafice de performanta
grafice de cache
grafice de comparatie intre arhitecturi

Asa acoperi si timpul, si cache-ul, si bonusul cu alta arhitectura, fara sa faci prea multe figuri inutile.