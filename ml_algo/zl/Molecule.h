#ifndef MOLECULE_H
#define MOLECULE_H
#include <iostream>

using namespace std;

#include <Bpp/Phyl/Tree.h> /* this includes classes for tree manipulations */
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>

using namespace bpp;

class Molecule
{
    public:
        Molecule();
        Molecule(Tree* tree);
        virtual ~Molecule();
        void update();
        void seeMolecule();
        double pe,minPE, ke, mp, ml;
        int num_of_hits, ID, minHit;
        //pair<float, float> structure;
        Tree* structure;

        Tree* minStruct;

        bool operator<(const Molecule & b)const{
            return (this->ml < b.ml);// < will give increasing sort
        }
    protected:
    private:

};

#endif // MOLECULE_H
