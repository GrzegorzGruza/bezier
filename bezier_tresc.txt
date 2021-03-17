Zadanie polega na zaimplementowaniu w paradygmacie programowania funkcyjnego
biblioteki do rysowania kawałkami sześciennych krzywych (ang. piecewise cubic
curves). Pojedynczy segment takiej krzywej ma być odrębną krzywą Béziera
stopnia 3, czyli funkcją B : [0,1] -> R^{2} postaci:

B(t) = p0*(1 - t)^{3} + 3*p1*t*(1 - t)^{2} + 3*p2*t^{2}*(1 - t) + p3*t^{3},

dla 0 <= t <= 1, gdzie p0, p1, p2 i p3 są punktami w R^{2}. Położenie węzłów
p0, p1, p2 i p3 na płaszczyźnie determinuje kształt krzywej, zatem punkty te
zwane są często punktami kontrolnymi Béziera. Do obliczania punktu B(t) leżącego
na krzywej Béziera, zdefiniowanej węzłami p0, p1, p2 i p3, stosuje się algorytm
de Casteljau, bazujący na liczeniu kolejnych kombinacji liniowych:

b0 = (1 - t)*p0 + t*p1
b1 = (1 - t)*p1 + t*p2
b2 = (1 - t)*p2 + t*p3

b0 = (1 - t)*b0 + t*b1
b1 = (1 - t)*b1 + t*b2

b0 = (1 - t)*b0 + t*b1

aż do otrzymania wyniku B(t) = b0.

Punkt na płaszczyźnie R^{2} powinien być w bibliotece reprezentowany typem
point_2d, przechowującym współrzędne punktu jako dwie niemodyfikowalne liczby
rzeczywiste (typu real_t, będącego aliasem typu zmiennoprzecinkowego,
np. double). Współrzędne punktu mogą być inicjowane jedynie w konstruktorze
i nie może to być konstruktor domyślny. Typ ten powinien ponadto przeciążać
tylko następujące operatory:
 + dodawanie dwóch punktów (suma po współrzędnych),
 + mnożenie skalar razy punkt i punkt razy skalar,
 + sprawdzanie równości dwóch punktów (równość po współrzędnych),
 + wysyłanie do strumienia wyjściowego danych punktu w postaci „(wsp_x, wsp_y)”.

Do reprezentowania punktów kontrolnych kawałkami sześciennej krzywej
(zbudowanej z jednego lub kilku segmentów Béziera) należy zamiast zwykłych
danych (wartości) stosować funkcję, która pobiera w parametrze indeks węzła
(liczba całkowita nieujemna typu node_index_t) i zwraca niemodyfikowalny punkt.
Zakładamy przy tym, że indeksy 0, 1, 2, 3 odpowiadają punktom p0, p1, p2 i p3
pierwszego segmentu Béziera, indeksy 4, 5, 6, 7 odpowiadają punktom p0, p1, p2
i p3 segmentu drugiego itd.

W bibliotece powinny się znaleźć specjalne funkcje tworzące pewne wybrane
układy węzłów pojedynczego (tylko pierwszego) segmentu Béziera:
 + Cup() – zwraca funkcję generującą punkty kontrolne
           p0 = (-1, 1), p1 = (-1, -1), p2 = (1, -1) i p3 = (1, 1),
           definiujące krzywą Béziera w kształcie litery U;
 + Cap() – zwraca funkcję generującą punkty kontrolne
           p0 = (-1, -1), p1 = (-1, 1), p2 = (1, 1) i p3 = (1, -1),
           dające krzywą w kształcie litery U obróconej o 180 stopni;
 + ConvexArc() – zwraca funkcję generującą węzły
                 p0 = (0, 1), p1 = (ARC, 1), p2 = (1, ARC) i p3 = (1, 0),
                 gdzie ARC = 4*(sqrt(2) - 1)/3, tworzące łuk wypukły zbliżony
                 do ćwiartki okręgu jednostkowego;
 + ConcaveArc() - zwraca funkcję generującą węzły
                  p0 = (0, 1), p1 = (0, 1 - ARC), p2 = (1 - ARC, 0) i p3 = (1, 0),
                  gdzie ARC j.w., tworzące łuk wklęsły zbliżony do ćwiartki
                  okręgu;
 + LineSegment(p, q) – zwraca funkcję generującą węzły p0 = p1 = p i p2 = p3 = q,
                       definiujące odcinek prostej łączący punkty p i q.
Funkcje generujące węzły, zwracane przez powyższe funkcje, powinny w sytuacji
podania im nieprawidłowego indeksu węzła zgłaszać wyjątek typu std::out_of_range
z komunikatem „a curve node index is out of range”. Prawidłowy indeks węzła ma
wartość mniejszą od liczby węzłów w sześciennej krzywej Béziera, reprezentowanej
stałą NUM_OF_CUBIC_BEZIER_NODES równą 4.

Biblioteka powinna również dawać możliwość wykonywania podstawowych
przekształceń węzłów całej krzywej kawałkami sześciennej. W tym celu powinna
udostępniać następujące funkcje:
 + MovePoint(f, i, x , y) - przesuwa i-ty węzeł krzywej f o x wzdłuż osi X
                            i o y wzdłuż osi Y;
 + Rotate(f, a) - obraca wszystkie węzły krzywej f wokół punktu (0, 0) o kąt
                  a wyrażony w stopniach, przeciwnie do ruchu wskazówek zegara;
 + Scale(f, x, y) - skaluje węzły krzywej f względem punktu (0, 0) skalą x
                    wzdłuż osi X i skalą y wzdłuż osi Y;
 + Translate(f, x, y) - przesuwa wszystkie węzły krzywej f o x wzdłuż osi X
                        i o y wzdłuż osi Y.

Oprócz tego biblioteka powinna udostępniać funkcję generyczną Concatenate,
pozwalającą na łączenie dwóch lub więcej krzywych, poprzez operację składania
funkcji. W szczególności, Concatenate(f1, f2) ma dołączać punkty kontrolne
generowane funkcją f2 do pierwszego segmentu Béziera, generowanego funkcją f1.
Podczas operacji konkatenacji końce sąsiadujących segmentów Béziera winny
w miejscu łączenia pozostać odrębnymi węzłami (punktami). Wersja funkcji
Concatenate z większą liczbą parametrów powinna korzystać z wersji
dwuparametrowej.

Za obsługę rysowania obrazu krzywej leżącej wewnątrz obszaru [-1, 1] x [-1, 1]
ma odpowiadać klasa P3CurvePlotter. Do publicznego interfejsu tej klasy mają
należeć:
 + konstruktor – pobierający kolejno w parametrach funkcję dostarczającą węzły
                 krzywej, liczbę segmentów tej krzywej (domyślnie 1) oraz
                 rozdzielczość wydruku (domyślnie 80 znaków); powinien on
                 obliczać i przygotowywać dane o wszystkich punktach na krzywej,
                 które mają być zaznaczone na wydruku; tych danych nie można
                 już później modyfikować; należy obliczyć tyle punktów na
                 krzywej, ile może maksymalnie zmieścić się na wydruku, tzn.
                 rozdzielczość podniesiona do kwadratu;
 + Print(s, fb, bg) – funkcja stała, drukująca do strumienia wyjściowego s
                      (domyślnie std::cout) kwadrat o boku rozdzielczość,
                      zbudowany ze znaków tak, że punkty leżące na krzywej
                      drukowane są znakiem fb (domyślnie '*'), zaś punkty tła
                      drukowane są znakiem bg (domyślnie spacja);
 + operator() - funkcja stała, pobierająca kolejno w parametrach funkcję
                dostarczającą węzły krzywej, parametr zmiennoprzecinkowy t oraz
                numer segmentu (licząc od 1); celem tego operatora jest
                zwrócenie punktu B(t) leżącego na wskazanym segmencie Béziera.
Wszelkie pomocnicze, dodatkowe metody w klasie P3CurvePlotter powinny być
prywatne. Jeżeli jakiś punkt krzywej znajdzie się poza obszarem wydruku, tzn.
poza kwadratem o wielkości rozdzielczość razy rozdzielczość znaków, to nie może
być drukowany.

Klasy, funkcje, typy pomocnicze i stałe powinny zostać umieszczone odpowiednio
w przestrzeniach nazw:

bezier
bezier::types
bezier::constants

                                      * * *

Przykład użycia biblioteki można znaleźć w pliku bezier_example.cc, zaś efekt
jego działania w plikach tekstowych bezier_example.out i clover.out.

Rozwiązanie będzie kompilowane z użyciem polecenia:

  g++ -Wall -Wextra -pedantic -O3 -std=c++17

Rozwiązanie powinno zawierać plik bezier.h, który należy umieścić
w repozytorium w katalogu

  grupaN/zadanie7/ab123456

gdzie N jest numerem grupy, a ab123456 jest identyfikatorem autora
umieszczającego to rozwiązanie. Katalog z rozwiązaniem nie powinien zawierać
innych plików, ale może zawierać podkatalog prywatne, gdzie można umieszczać
różne pliki, np. swoje testy. Pliki umieszczone w tym podkatalogu nie będą
oceniane. Nie wolno umieszczać w repozytorium plików dużych, binarnych,
tymczasowych (np. *.o) ani innych zbędnych.
