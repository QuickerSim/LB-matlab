Repozytorium zawiera implementację metody Lattice-Bolztmann na jednorodnych oraz niejednorodnych siatkach w oprogramowaniu MATLAB.

W celu uruchomienia przykładów należy uruchomić jeden ze skryptów:
- 'bfs_uni.m' - Backward facing step na jednorodnej siatce
- 'cavity_uni.m' - Cavity na jednorodnej siatce
- 'cavity_nonuni.m' - Cavity na niejednorodnej siatce

W przypadku jednorodnych przykładów można sterować zagęszczeniem siatki za pomocą parametrów 'nx' oraz 'ny'.

W przypadku niejednorodnym aby zagęszczać siatkę nalezy skorzystać z metody 'createTree', gdzie jako argumenty podaje się punkty, w których ma następować zagęszczenie oraz poziom zagęszczenia.