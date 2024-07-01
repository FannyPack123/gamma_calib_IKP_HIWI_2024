# gamma_calib_IKP_HIWI_2024
Gamma-Kalibrierungsprojekt HIWI-Stelle am IKP 2024


quickParsingForJena:
1. auszuwertende Daten als ein Ordner mit beliebigen Namen in Verzeichnis quickData legen
2. in parseData.py den String dataSet zum Namen des Ordners mit den Daten umbenennen
3. Programm ausführen
4. Ausgewertete Daten liegen im quickOut Verzeichnis:
   - Histogramm über alle 2^14 Kanäle des SIS (nicht die "Channels", an denen die Detektoren angeschlossen sind)
   - .csv file mit allen Datenpunkte, welche auch im Histogramm zu sehen sind
   - extra .txt file für debugging Gründe, enthält Anzahl aller aufgenommenen Ereignisse und Kanal mit den meisten Ereignissen
