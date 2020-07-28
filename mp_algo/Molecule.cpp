#include "Molecule.h"



using namespace bpp;


Molecule::Molecule()
{
   this->ID = 0;
   this->pe = 0;
   this->ke = 100;
   this->mp = 0;
   this->ml = 0;
   this->num_of_hits = 0;

   
}
Molecule::Molecule(Tree* st)
{
   this->ID = 0;
   this->pe = 0;
   this->ke = 100;
   this->mp = 0;
   this->ml = 0;
   this->num_of_hits = 0;
   this->structure = st;
   this->minStruct = st;
   this->minPE = pe;
   this->minHit = 0;
}

Molecule::~Molecule()
{
    //dtor
    //delete structure;
}

void Molecule::update()
{
    
}

void Molecule::seeMolecule()
{
   
   printf("ID-%d, pe-%.2lf, ke-%.2lf, mp-%.2lf, ml-%.2lf, nHit-%d, mPE-%.2lf, mHit-%d\n",  ID , pe, ke, mp, ml, num_of_hits, minPE, minHit); 
}
