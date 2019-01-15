Testiranje:
    - Na pocetku izvrsavanja programa, korisnik bira da li zeli da unosi fajl putem komandne linije ili iz fajla.
        1) Komandna linija:
            - Primer kako praviti graf i ubacivati grane na slici command_line.png
        2) Iz fajla:
            - Primer kako praviti graf i ubacivati grane na slici from_file.png

    - Program se izvrsava dok se na pitanje Continue testing? (y/n) ne odgovori sa 'n'.

    Formatiranje fajla koji sadrzi graf:
        - Prva linija fajla sadrzi broj cvorova.
        - Sledecih n linija predstavlja grane u formatu:
            - Prvi broj u liniji predstavlja cvor od koga ce polaziti grane do cvorova koji su navedeni nakon njega,
            u istoj toj liniji. Sledeci brojevi dakle predstavljaju cvorove do kojih postoji grana iz tog cvora
            koji je naveden na pocetku linije.
