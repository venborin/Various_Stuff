Various useful stuff

1. get_atoms

   prints out the list of sidechain/solvent/ligand atoms that have at least one atom within a certain cutoff from a residue/solvent/ligand.
   Usefull for MM, QM, an QM/MM otimization your stuff within the relaxed binding pocket/bulk of solvent.
   
   ./get_atoms -h for help
2. knight_move
   
   Have you ever wondered how many moves it would take for knight to move from square A to square B on a chessboard? Well, this program is for you.

3. convolute.cpp

   performs gaussian broadening of the calculated UV-Vis spectrum.
   Building: g++ convolute.cpp -o convolute -std=c++17
   ./convolute -h for help
