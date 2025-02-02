# Modelowanie ogrzewania domu - Symulacja temperatury w pokojach

Projekt ten implementuje model numeryczny dla przepływu ciepła w domu z grzejnikami, oknami, drzwiami oraz z warunkami brzegowymi Dirichleta i Neumanna.

## Opis projektu

Projekt symuluje rozprzestrzenianie się ciepła w trzech pokojach, z uwzględnieniem grzejników oraz okien. Modeluje również wymianę ciepła przez drzwi, a także nakłada warunki brzegowe na okna i ściany. Symulacja jest przeprowadzana za pomocą metody różnic skończonych, a wyniki wyświetlane są w formie mapy ciepła oraz wykresu całkowitej zużytej energii.

## Struktura projektu

### Klasy:

- **Heater** - reprezentuje grzejnik z parametrami: pozycja, moc, maksymalna temperatura oraz powierzchnia.
- **Window** - reprezentuje okno, które ma określoną temperaturę.
- **Door** - reprezentuje drzwi w domu, które łączą różne pomieszczenia.
- **Room** - reprezentuje pojedynczy pokój w domu, uwzględnia grzejnik i warunki brzegowe.
- **House** - reprezentuje cały dom, składający się z różnych pokoi.
- **Solver** - główny solver, który przeprowadza symulację, nakłada warunki brzegowe i generuje wyniki.

## Wymagania

- Python 3.x
- Zainstalowane biblioteki: `numpy`, `matplotlib`, `ffmpeg`

Aby zainstalować wymagane biblioteki, uruchom:

```
pip install -r requirements.txt
```

## Jak uruchomić projekt

1. Skopiuj pliki `kod.py`, `run_experiments.py` oraz `requirements.txt` do swojego projektu.
2. Zainstaluj wymagane biblioteki przy pomocy `pip install -r requirements.txt`.
3. Uruchom eksperyment, wywołując plik `run_experiments.py`:
```
python run_experiments.py
```
4. Eksperymenty zostaną uruchomione, a wyniki będą dostępne w postaci mapy ciepła i wykresu energii.

5. Możesz wygenerować mapy ciepła oraz wykresy energii używając metod:

```python
alg.plot_result()  # Mapa ciepła
alg.total_energy_use()  # Całkowita energia używana
```

6. Możesz również uzyskać animację przedstawiającą ewolucję temperatury w czasie:
```python
alg.result_animated()
```


## Parametry wejściowe

- **Dyskretyzacja przestrzeni**: Parametry hx, xmax, ymax definiują rozdzielczość siatki przestrzennej i rozmiar obszaru (np. pokój, dom).
- **Dyskretyzacja czasu**: Parametry ht, Tmax definiują krok czasowy oraz czas trwania symulacji.
- **Grzejniki, okna, drzwi**: Możesz dodać różne urządzenia w pokoju, podając ich współrzędne oraz parametry fizyczne.

## Przykładowe wyniki

Po zakończeniu symulacji, możesz zobaczyć mapę ciepła przedstawiającą rozkład temperatury w pomieszczeniu oraz wykres całkowitej zużytej energii w czasie

## Licencja

Projekt jest udostępniany na zasadach licencji MIT.

