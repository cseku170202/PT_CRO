    #ifndef CRO_H
    #define CRO_H
    #include "Molecule.h"
    #include <iostream>
    #include <bits/stdc++.h>

    using namespace std;


    #include <Bpp/Seq/Alphabet/AlphabetTools.h>  /* this includes all alphabets in one shot */
    #include <Bpp/Seq/Container/VectorSiteContainer.h>
    #include <Bpp/Seq/Container/SiteContainerTools.h>
    #include <Bpp/Seq/Io/Fasta.h>


    /*
     * From Bpp-Phyl:
     */
    #include <Bpp/Phyl/Tree.h> /* this includes classes for tree manipulations */
    #include <Bpp/Phyl/TreeTools.h>
    #include <Bpp/Phyl/Io/Newick.h>
    #include <Bpp/Phyl/Model/Nucleotide/K80.h>
    #include <Bpp/Phyl/Model/Nucleotide/GTR.h>
    #include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
    #include <Bpp/Phyl/Distance/DistanceEstimation.h>
    #include <Bpp/Phyl/Distance/BioNJ.h>
    #include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
    #include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
    #include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>
    #include <Bpp/Phyl/OptimizationTools.h>
    #include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
    #include <Bpp/Phyl/Node.h>
    #include <Bpp/Phyl/TreeTemplate.h>
    #include <Bpp/Phyl/TreeTemplateTools.h>
    #include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>


    /*
     * All Bio++ functions are also in a namespace, so we'll use it:
     */
    #include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
    using namespace bpp;

    class CRO
    {
        public:
            CRO();
            virtual ~CRO();
            void init_pop_with_bionj();
            void init_pop_with_rand_tree();
            Tree* TBR(Tree* tr);
            Tree* NNI(Tree * tr);
            Tree* SPR(Tree* tr);
            bool NNIValidate(Node * Nodo);
            int  SPRvalide (Node* N1, Node* N2);
            float updatePE(float mp, float ml);
            float random1();
            void run();
            void run1();
            void branch_optimization();
            void update(int m);
            void decomposition(int m);
            void on_wall(int m);
            void interaction(int m1, int m2);
            void synthesis(int m1, int m2);
            void synthesisNew(int m1, int m2);
            Tree *  synthesis_with_contribution(Tree * tr1, Tree * tr2);
            Molecule syn(Molecule & m1, Molecule & m2);
            pair<Molecule, Molecule> inter(Molecule & m1, Molecule & m2);
            Tree* wall( Molecule & m);
            pair<Molecule,Molecule> dec(Molecule & m);
            float MP_value(Molecule & m);
            float ML_value(Molecule & m);
            float ML_value_from_tree(Tree * t);
            float MP_value_from_tree(Tree * tr);
            float ofv(Tree* t);
            void erase_molecule( int m);
            void printPopulation();
            Node * selectNodeToCross(TreeTemplate<Node> * tree_, vector<int> nodosIDs );
            Tree* PDG(Tree * t1, Tree * t2 ) ;
            Tree* random_tree(Tree * tr);
            int random_num(int start, int end);
            void readingTreesFromAFile();
            void branch_optimize_parameter();
            void lastHope();
            void branch(Tree * t);
            Tree*  decomposition_with_contribution(Tree * tr1);

            Molecule optimal;

            //input and objective function calculation

            vector<Molecule> pop;

            SubstitutionModel* model;

            DiscreteDistribution* rateDist;

            Fasta seqReader;
                //correct fasta file link
            SequenceContainer* sequences ;
            SiteContainer* sites ;

            NNIHomogeneousTreeLikelihood* tl;
            RHomogeneousTreeLikelihood* tLK;
        

        protected:
        private:
    };

    #endif // CRO_H
