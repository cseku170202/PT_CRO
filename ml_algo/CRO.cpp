        #include "CRO.h"
        #include "Molecule.h"
        #include <bits/stdc++.h>

        #define PI 3.14

        using namespace std;

        



        int cnt = 0;
        double limit = 1000, KELossRate = 0.2, mole_coll = 0.3;
        double alpha = 20, beta = 20, buffer = 0, init_ke = 100;
        int init_pop = 0;
        int nD = 0, nS = 0, nI = 0, nO = 0;
        string st, choice;
        vector< pair<float, float> > result;
        vector<string> species;
        vector<int> numOfHitsCollection;
        double threshold = 0;
        int nos = 0;
        int is_first_molecule = 0;
        int molecule_count = 0;
        ofstream myfile;
        Newick treeWriter;
        clock_t t_ini, t_fin;
        double mini = 999999;
        int blo_tree = 0;
            
        CRO::CRO()
        {

        time_t tic = clock();

        string report, dataset;

        freopen("input.txt", "r", stdin);
        //freopen("rbclOut.txt", "w", stdout);
        cout << "Enter report file name: ";
        cin >> report;
        myfile.open (report);
        cout << "Enter dataset" << endl;
        cin >> dataset;
        cout << "Hello World :) " << endl;
        if(dataset == "rb")
        {
            model = new GTR(&AlphabetTools::DNA_ALPHABET, 1.6325, 0.1290, 0.1950, 0.3371, 0.3812, 0.3104,0.1588 ,  0.1669 , 0.3639  );
            rateDist = new GammaDiscreteDistribution(4,  0.3650 ,  0.3650 );

        }else if(dataset == "mt")
        {
            model = new GTR(&AlphabetTools::DNA_ALPHABET, 0.676539, 0.0201724, 0.0170779, 0.0270668, 0.0283801, 0.3104,  0.3166,  0.1296, 0.2434 );
            rateDist = new GammaDiscreteDistribution(4,  0.05 ,  0.05 );

        }else if(dataset == "rp")
        {
            model = new GTR(&AlphabetTools::DNA_ALPHABET,1.50316, 0.515473, 0.436472, 0.406966, 0.492165, 0.2405, 0.2204,  0.2865,  0.2526  );
            rateDist = new GammaDiscreteDistribution(4,  0.5450 ,  0.5450 );

        }else if(dataset == "zl")
        {
            model = new GTR(&AlphabetTools::DNA_ALPHABET, 0.8085, 0.0864304, 0.248328, 0.308679, 0.430362,0.28759,  0.17822,  0.14682, 0.38737);
            rateDist = new GammaDiscreteDistribution(4,  0.867 ,  0.867 );

        }
        
        //rateDist = new ConstantRateDistribution();
        
        

        cout << "Hello !! Enter name of the dataset to be used" << endl;
        cin >> st;
        cout << "Enter parameters" << endl;
        cout << "Limit mole_coll alpha beta KELossRate"<<endl;
        cin >> limit >> mole_coll >> alpha >> beta >> KELossRate;
        // cout << "Enter threshold value for termination" << endl;
        // cin >> threshold;


        sequences = seqReader.readSequences(st, &AlphabetTools::DNA_ALPHABET);
        sites = new VectorSiteContainer(*sequences);
        delete sequences;
        SiteContainerTools::removeGapOnlySites(*sites); 
        SiteContainerTools::changeGapsToUnknownCharacters(*sites);

        species = sites->getSequencesNames();

        for (vector<string>::iterator it=species.begin(); it!=species.end(); ++it)
             {
                
                cout << ' ' << *it;
                //printf("%s\n",*it );
             }
             nos = species.size();
             cout << "number of species" << species.size() << endl;

             cout << "Enter initial population" << endl;
             cin >> init_pop;

             

            // cout << "Computing tree..." << endl;
            // DistanceEstimation distanceMethod(model, rateDist, sites);
            // DistanceMatrix* distances = distanceMethod.getMatrix();
 
            // BioNJ bionj(*distances);
            // Tree* tre = bionj.getTree();
            // // NeighborJoining nj(*distances, false, true); // unrooted tree.
            // // Tree* tre = nj.getTree();

            // Molecule bioTree(tre);
            // myfile << "Bionj tree : " << MP_value(bioTree) << "(MP) and " << ML_value(bioTree) << "(ML)" << endl;
            // pop.push_back(bioTree);


            

            // for(int i = 0; i < init_pop ; )
            // {
                
            //    Tree* t_new = tre->clone();
            //    t_new = random_tree(t_new);
            //    //cout << TreeTemplateTools::treeToParenthesis(t_new) << endl;
            //    Molecule mn(t_new);
               
            //     if(mn.structure->getNumberOfLeaves() == nos)
            //     {
            //             i++;
            //             pop.push_back(mn);
                    
            //             update(i+1);    
            //     }
                
            // }

            // for(int i = 0; i < init_pop; i++)
            // {
                
            //     TreeTemplate<Node>* tree = TreeTemplateTools::getRandomTree(species, false);
            //     tree->setBranchLengths(0.05);
            //     Tree* t = (Tree*) tree;
                
            //     myfile << "Random tree : "<< i << " " << ofv(t) << endl;
            //     Molecule mn(t);
            //     pop.push_back(mn);

            // }
             readingTreesFromAFile();
             
            int fi = 0;
            for(vector<Molecule>::iterator it = pop.begin(); it != pop.end(); ++it)
            {

                //it->mp = this->MP_value(*it);
                it->mp = 0;
                it->ml = this->ML_value(*it);
                it->ID = molecule_count++;

                //it->pe = it->mp;

                //for ml signle
                it->pe = updatePE(it->mp,it->ml);

                it->minPE = it->pe;

                it->ke = 100;
                //update(fi++);
                if(!is_first_molecule)
                {
                    optimal = *it;
                    optimal.pe = 10000000.0;
                }
                else if(it->pe < optimal.pe)
                {
                    //optimal = *it;
                }
                //myfile << "Molecule number "<< fi+1 << " PE is :  "<< it->pe << endl;

            }
            
            run();
             
            
            treeWriter.write(*optimal.structure, "OptimalTree.dnd");
            //myfile << "Optimal tree MP " << optimal.mp <<  "  ML " << optimal.ml<< endl;
            myfile  << "Optimal tree MP " << optimal.mp <<  "  ML " << optimal.ml<< endl;
            myfile << "For dataset " << st << endl;
            myfile << "Limit : " << limit << " Mole coll : " << mole_coll   << " Alpha : " <<  alpha <<  "Beta : " <<  beta << " KELossRate : " <<KELossRate;
            myfile << "Population size " << init_pop << endl;

            cout << "Optimal tree MP " << optimal.mp <<  "  ML " << optimal.ml<< endl;

            time_t tac = clock();

            myfile << "Elaspsed: " << (double)(tac - tic) / CLOCKS_PER_SEC << "Seconds" << endl;
            reform2();

            time_t toe = clock();
           	myfile << "Elaspsed with repair 2: " << (double)(toe - tac) / CLOCKS_PER_SEC << "Seconds" << endl;
             
        }

        

        CRO::~CRO()
        {
            //dtor
        }
        // calling this function from constructor is crushing
        
        

        float CRO::updatePE(float mp, float ml)
        {
            return ml;
        }

        
        float CRO::random1()
        {
            return static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        }

        int CRO::random_num(int start, int end)
        {
            int range = (end-start)+1;
            int random_int = start+(rand()%range);
            return random_int;
        }

        void CRO::run()
        {
            
            
            while(limit > cnt)
            {
                cnt += 1;
                printf("Iteration : %d\n", cnt);
                
                int ps = pop.size();
                
                
                if(random1() > mole_coll or ps < 2)
                {
                    int t3;
                    
                    t3 = rand()%ps;
                    printf("m1 : %d\n",t3);
                    //pop[t3].seeMolecule();
                    //decomposition(t3);
                    if(pop[t3].num_of_hits >= alpha )
                    {
                        
                        printf("___________Decomposition : %d\n",cnt);
                        nD++;
                        decomposition(t3);
                    }
                    else
                    {
                        
                        printf("_____________On-wall : %d\n", cnt);
                        nO++;
                        on_wall(t3);
                    }
                }else
                {
                    int t1, t2;
                
                    t1 = rand()%ps;
                    t2 = rand()%ps;
                    //pop[t1].seeMolecule();
                    //pop[t2].seeMolecule();
                    printf("____________m1 : %d, m2 : %d\n",t1, t2);
                    //synthesisNew(t1, t2);
                    if( pop[t1].num_of_hits > beta and pop[t2].num_of_hits > beta )
                    {    
                        printf("___________Synthesis : %d\n",cnt);
                        nS++;
                        //synthesis(t1, t2);
                        synthesisNew(t1, t2);
                    }
                    else
                    {
                           
                        printf("____________Intermolecular : %d\n", cnt);

                        nI++;
                        interaction(t1, t2);
                    }


                }
           
                
                
            }

            	printf("______________________________________\n");
                printf("___Minimum value %.2lf  ____ \n" , optimal.pe);
                printf("______________________________________\n"); 
                
            //printPopulation();
            cout << "Decomposition :" <<nD << ", On-wall :" << nO<< ", Synthesis :"<<nS<< "Inter-molecular :"<<nI<<endl;
            //branch_optimization();
            //branch_optimize_parameter();
            //timer = time(NULL);    printf("Final CRO Optimization Ends %s\n", ctime(&timer));
            //lastHope();
            //timer = time(NULL);    printf("Final BLO Optimization Ends %s\n", ctime(&timer));
        }

        void CRO::update(int m)
        {
            if(pop[m].pe < optimal.pe)
            {
                //printf("-*-*-*-{-}-}-{-->Optimal value updated from %.2lf to %.2lf\n",optimal.pe, pop[m].pe);
                optimal = pop[m];
                optimal.pe = updatePE(optimal.mp, optimal.ml);
            }

        }

        

        void CRO::on_wall(int m)
        {

            double peW1;
            
            
            Tree* w1 = pop[m].structure->clone();
            w1 = NNI(w1);

            peW1 = ML_value_from_tree(w1);

            pop[m].num_of_hits += 1;

            
            if(w1->getNumberOfLeaves() == nos )
            {
                
                if(peW1 < pop[m].pe)
                {
                    delete pop[m].structure;
                    pop[m].structure = w1;
                    pop[m].pe = peW1;
                    pop[m].ml = peW1;
                    //pop[m].mp = MP_value(pop[m]);
                    pop[m].mp = 0;

                }else
                {
                    delete w1;
                    pop[m].ke--;
                } 
                update(m);
            }
                                

        }

        void CRO::decomposition(int m)
        {
            
            Molecule new1 , new2;
            double peW1, peW2;

            // TreeTemplate<Node>* dt1 = TreeTemplateTools::getRandomTree(species, true);
            // dt1->setBranchLengths(0.05);
            // Tree* td1 = (Tree*) dt1;

            // TreeTemplate<Node>* dt2 = TreeTemplateTools::getRandomTree(species, true);
            // dt2->setBranchLengths(0.05);
            // Tree* td2 = (Tree*) dt2;

            // Tree * w1 = PDG(pop[m].structure,td1);
            // Tree * w2 = PDG(pop[m].structure,td2);

            Tree * w1 = decomposition_with_contribution(pop[m].structure);
            Tree * w2 = decomposition_with_contribution(pop[m].structure);


            peW1 = ML_value_from_tree(w1);
            peW2 = ML_value_from_tree(w2);

            if(w1->getNumberOfLeaves() == nos )
            {
                if(peW1 < pop[m].pe)
                {
                        new1.structure = w1;
                        new1.pe = peW1;
                        new1.ke = init_ke;
                        new1.ml = peW1;
                        //new1.mp = MP_value(new1);
                        new1.mp = 0;
                        new1.num_of_hits = 0;
                        pop.push_back(new1);
                        
                }else
                {
                        //delete w1;
                        pop[m].ke--;
                        new1.~Molecule();
                }
            }
            
            if(w2->getNumberOfLeaves() == nos )
            {


                if(peW2 < pop[m].pe)
                {
                        new2.structure = w2;
                        new2.pe = peW2;
                        new2.ke = init_ke;
                        new2.ml = peW2;
                        //new2.mp = MP_value(new2);
                        new2.mp = 0;
                        new2.num_of_hits = 0;
                        pop.push_back(new2);
                        
                }else
                {
                        //delete w2;
                        pop[m].ke--;
                        new2.~Molecule();
                }

            }    //delete td1,td2, dt1, dt2;
                  
                     
        }

                    
        void CRO::interaction(int m1, int m2)
        {
            
            double peW1, keW1, peW2, keW2;
            
            Tree *  w1 = SPR(pop[m1].structure);

            Tree * w2 = SPR(pop[m2].structure);
               

            peW1 = ML_value_from_tree(w1);
            peW2 = ML_value_from_tree(w2);
            

            pop[m1].num_of_hits += 1;
            pop[m2].num_of_hits += 1;

            
            if(w1->getNumberOfLeaves() == nos && w2->getNumberOfLeaves() == nos)
            {
                
                

                if(peW1 < pop[m1].pe)
                {
                    delete pop[m1].structure;
                    pop[m1].structure = w1;
                    pop[m1].pe = peW1;
                    pop[m1].ml = peW1;
                    //pop[m1].mp = MP_value(pop[m1]);
                    pop[m1].mp = 0;
                    
                }else
                {
                    delete w1;
                    pop[m1].ke--;

                }

                if(peW2 < pop[m2].pe)
                {
                    delete pop[m2].structure;
                    pop[m2].structure = w2;
                    pop[m2].pe = peW1;
                    pop[m2].ml = peW1;
                    //pop[m2].mp = MP_value(pop[m2]);
                    pop[m2].mp = 0;
                    
                }else
                {
                    delete w2;
                    pop[m2].ke--;
                }
                update(m1);
                update(m2);
            }

        }

        void CRO::synthesis(int m1, int m2)
        {
           if(pop.size() > init_pop)
           {
                if(pop[m1].pe < pop[m2].pe)
                {
                    erase_molecule(m2);
                }else
                {
                    erase_molecule(m1);
                }
           } 
            
            
        }

        void CRO::synthesisNew(int m1, int m2)
        {
            double peW1, keW1, peW2, keW2;
            Molecule n1;
            
            Tree *  new1 = synthesis_with_contribution(pop[m1].structure, pop[m2].structure);

            peW1 = ML_value_from_tree(new1);
            
            pop[m1].num_of_hits += 1;
            pop[m2].num_of_hits += 1;

            n1.structure = new1;
            n1.pe = peW1;
            n1.ke = init_ke;
            n1.ml = peW1;
            n1.mp = 0;
            n1.num_of_hits = 0;
            pop.push_back(n1);
                        

        }

        void CRO::erase_molecule( int m)
        {
            int count = 0;
            for(vector<Molecule>::iterator it = pop.begin(); it != pop.end(); it++)
            {
                if(count  == m)
                {
                    it->~Molecule();
                    pop.erase(it);
                    break;
                }
                count++;
            }
            
        }

        void CRO::branch_optimization()
        {
            NNIHomogeneousTreeLikelihood* bo = new NNIHomogeneousTreeLikelihood(*optimal.structure, *sites, model, rateDist); 
            bo->initialize();
            cout << "Log likelihood before: " << bo->getLogLikelihood() << endl;
            
            //bo = OptimizationTools::optimizeTreeNNI(bo, bo->getParameters(), true, 100, 0.000001, 300, 1, 0, 0, false, 3);
            bo = OptimizationTools::optimizeTreeNNI(bo, bo->getParameters(), true, 100, 0.0001, 50, 1, 0, 0, false, 0,"newton", 3);
            cout << "Value from optimal tree " << ML_value_from_tree(optimal.structure) << endl;
            cout << "Log likelihood after: " << bo->getLogLikelihood() << endl;
            myfile << "Log likelihood after: " << bo->getLogLikelihood() << endl;
            bo->getParameters().printParameters(cout);
            TreeTemplate<Node> mlTree(bo->getTree());
            treeWriter.write(mlTree, "branchOptimizedTree.dnd");
            myfile << "Final mp,ml value after branchlength optimization using our mp, ml value function" << endl;
            myfile << MP_value(optimal) << "is mp value and " << ML_value(optimal) << "is ml value" << endl;

            delete bo;
        }

        void CRO::branch_optimize_parameter()
        {
            cout << "Optimize Numerical Parameters" << endl;
            DRHomogeneousTreeLikelihood* likFunction =
            new DRHomogeneousTreeLikelihood(
                *optimal.structure, *sites, model, rateDist, true);
                likFunction->initialize();
                cout << "Before : " << likFunction->getLogLikelihood() << endl;
                // Save messages to separate files:
                ofstream* profiler = new ofstream("profile.txt", ios::out);
                ofstream* messenger = new ofstream("messages.txt", ios::out);
                OptimizationTools::optimizeNumericalParameters(likFunction,
            likFunction->getParameters(),
            0,
            5,
            0.000001, //Precision on log likelihood
            1000000,  //Max # of evaluations
            new StlOutputStreamWrapper(messenger),
            new StlOutputStreamWrapper(profiler),
            false,    //No reparametrization
            3);
            cout << "Done." << endl;
            cout <<"After : " << likFunction->getLogLikelihood() << endl;

            delete likFunction;
        }
        

        Molecule CRO::syn(Molecule & m1, Molecule & m2)
        {
                Molecule new1  = m1;

                Tree* t = PDG(m1.structure, m2.structure);
                new1.structure = t;
                float p1, p2, p3;
                p1 = m1.pe;
                p2 = m2.pe;
                p3 = MP_value(new1);

                if(max(p1, max(p2, p3)) == p1)
                {
                    
                    return m1;
                }
                else if(max(p1, max(p2, p3)) == p2)
                {
                    
                    return m2;
                }
                else{

                    return new1;
                }

        }

        pair<Molecule, Molecule> CRO::inter(Molecule & m1, Molecule & m2)
        {
                Molecule new1 = m1;
                Molecule new2 = m2;
                
                Tree* t1 = new1.structure;

                for(int i = 0; i < 5; i++)
                {
                    t1 = NNI(t1);
                    if(m1.ml < ML_value_from_tree(t1))
                        break;
                }
                new1.structure = t1;

                Tree* t2 = new2.structure;

                for(int i = 0; i < 5; i++)
                {
                    t2 = NNI(t2);
                    if(m2.ml < ML_value_from_tree(t2))
                        break;
                }
                new2.structure = t2;

            
                pair<Molecule, Molecule> pm;
                pm = make_pair(new1, new2);
                return pm;
        }

        Tree* CRO::wall( Molecule & m)
        {
            Tree* t1 = m.structure->clone();
            
            t1 = NNI(t1);
            
            return t1;
        }


        pair<Molecule,Molecule> CRO::dec(Molecule & m)
        {
                pair<Molecule, Molecule> pm;

                Tree* t1 = m.structure->clone();
                
                Tree* t2 = m.structure->clone();
                

                t1 = TBR(t1);
                t2 = SPR(t2);
                

                Molecule new1(t1);
                Molecule new2(t2);

                
                
                pm = make_pair(new1, new2);
                return pm;
        }

        float CRO::MP_value(Molecule & m){
            
            DRTreeParsimonyScore pars(*m.structure, *sites, true, true);
            

            return pars.getScore();
            

        }

        float CRO::ML_value(Molecule & m){
            
            double ml_score = 0;
            //NNIHomogeneousTreeLikelihood* tl = new NNIHomogeneousTreeLikelihood(*m.structure, *sites, model, rateDist); 
            RHomogeneousTreeLikelihood * tLK = new RHomogeneousTreeLikelihood(*m.structure, *sites, model, rateDist, true, false, true); 
            tLK->initialize();
            ml_score = tLK->getLogLikelihood();
            
            delete tLK;

            return -ml_score;
            

        }

        float CRO::ML_value_from_tree(Tree * t){
            
            double ml_score = 0;
            //NNIHomogeneousTreeLikelihood* tl = new NNIHomogeneousTreeLikelihood(*t, *sites, model, rateDist);
            RHomogeneousTreeLikelihood * tLK = new RHomogeneousTreeLikelihood(*t, *sites, model, rateDist, true, false, true); 
            tLK->initialize();
            ml_score = tLK->getLogLikelihood();
            delete tLK;

            return -ml_score;
            

        }

        float CRO::ofv(Tree* t){
            
            DRTreeParsimonyScore pars(*t, *sites, true, true);
            cout << "In operator function" <<pars.getScore();

            return pars.getScore();
            

        }

        

        bool CRO::NNIValidate(Node * Nodo){

        if (Nodo->getNumberOfSons()>1){
              if (!Nodo->getSon(0)->isLeaf() and !Nodo->getSon(1)->isLeaf()) return true;    
        }
        return false;

        }

        Tree* CRO::NNI(Tree * tr){

            
            tl = new NNIHomogeneousTreeLikelihood(*tr, *sites, model, rateDist); 
            tl->initialize();
            TreeTemplate<Node>  tre(tl->getTree());
            TreeTemplate<Node>* tree = new TreeTemplate<Node>(tre);
            

            Node * NodoSel;
            vector<Node *> nodes = tree->getNodes();
            
            do{
                  
                  NodoSel =  nodes[rand()%(nodes.size() - 1)];
                  
                  
            }while(!NNIValidate(NodoSel));
            
            Node * Nodo1;
            Node * Nodo2;
            int Pos1, Pos2;

            

            cout << NodoSel->getSon(0)->getNumberOfSons()-1 << endl;
            cout << NodoSel->getSon(1)->getNumberOfSons()-1 << endl;
            Pos1 = rand()%(NodoSel->getSon(0)->getNumberOfSons()-1);
            Pos2 = rand()%(NodoSel->getSon(1)->getNumberOfSons()-1);
            

            Nodo1=  NodoSel->getSon(0)->getSon(Pos1);
            Nodo2=  NodoSel->getSon(1)->getSon(Pos2);

            NodoSel->getSon(0)->setSon(Pos1,Nodo2);
            NodoSel->getSon(1)->setSon(Pos2,Nodo1);

            tr = (Tree*) tree;

            delete tl;

            return tr;

        }

        Tree* CRO::random_tree(Tree * tr){

            
            tl = new NNIHomogeneousTreeLikelihood(*tr, *sites, model, rateDist); 
            tl->initialize();
            TreeTemplate<Node>  tre(tl->getTree());
            TreeTemplate<Node>* tree = new TreeTemplate<Node>(tre);
            

            Node * NodoSel;

            
            vector<Node *> nodes = tree->getNodes();
            //experiment
            // vector<int> nodeIds = tree->getNodesId();
            // for(int i = 0; i < nodes.size(); i++)
            // {
            //         cout << "ID " << nodeIds[i] << "number of leafs" << TreeTemplateTools::getNumberOfLeaves(*nodes[i]) << endl;
            //         if(nodes[i]->hasName())
            //             cout << nodes[i]->getName() << endl ;
            //         for(int j = 0; j < nodes[i]->getNumberOfSons(); j++)
            //             cout <<  nodes[i]->getSon(j)->getId() << "  ";
            //         cout << endl;                
                
            // }
            
            do{
                  
                  NodoSel =  nodes[rand()%(nodes.size() - 1)];
                  
                  
            }while(!NNIValidate(NodoSel));
            
            Node * Nodo1;
            Node * Nodo2;
            int Pos1, Pos2;

            cout << NodoSel->getSon(0)->getNumberOfSons()-1 << endl;
            cout << NodoSel->getSon(1)->getNumberOfSons()-1 << endl;
            Pos1 = rand()%(NodoSel->getSon(0)->getNumberOfSons()-1);
            Pos2 = rand()%(NodoSel->getSon(1)->getNumberOfSons()-1);
            

            Nodo1=  NodoSel->getSon(0)->getSon(Pos1);
            Nodo2=  NodoSel->getSon(1)->getSon(Pos2);

            cout << "Node 1 : "<< TreeTemplateTools::getNumberOfLeaves(*Nodo1)<< " Node 2: "<< TreeTemplateTools::getNumberOfLeaves(*Nodo2);

            NodoSel->getSon(0)->setSon(Pos1,Nodo2);
            NodoSel->getSon(1)->setSon(Pos2,Nodo1);

            do{
                  
                  NodoSel =  nodes[rand()%(nodes.size() - 1)];
                  
                  
            }while(!NNIValidate(NodoSel));
           

            cout << NodoSel->getSon(0)->getNumberOfSons()-1 << endl;
            cout << NodoSel->getSon(1)->getNumberOfSons()-1 << endl;
            Pos1 = rand()%(NodoSel->getSon(0)->getNumberOfSons()-1);
            Pos2 = rand()%(NodoSel->getSon(1)->getNumberOfSons()-1);
            

            Nodo1=  NodoSel->getSon(0)->getSon(Pos1);
            Nodo2=  NodoSel->getSon(1)->getSon(Pos2);

            cout << "Node 1 : "<< TreeTemplateTools::getNumberOfLeaves(*Nodo1)<< " Node 2: "<< TreeTemplateTools::getNumberOfLeaves(*Nodo2);

            NodoSel->getSon(0)->setSon(Pos1,Nodo2);
            NodoSel->getSon(1)->setSon(Pos2,Nodo1);

            tr = (Tree*) tree;

            delete tl;

            return tr;

        }


        Tree* CRO::SPR(Tree* tr){

        
        tl = new NNIHomogeneousTreeLikelihood(*tr, *sites, model, rateDist); 
        tl->initialize();
        TreeTemplate<Node>  tre(tl->getTree());
        TreeTemplate<Node>* tree = new TreeTemplate<Node>(tre);

        
        int NextIDNode = tree->getNextId();
       
        bool b;
        Node* Nodo1;
        Node* Nodo2;
        
         vector<Node*> nodes = tree->getNodes();
         do{
            b=true;
            do {
                   Nodo1 =  nodes[rand()%nodes.size()];
                   if(Nodo1->hasFather()){
                       if(Nodo1->getFather()->hasFather()) b=false;
                    }
            }while(b);
         
            Nodo2 =  nodes[rand()%nodes.size()];
            
        }while(!SPRvalide (Nodo1,Nodo2));     
        
        
        int PosNodo;
        double distancetofather=0;
        Node * father;
        Node * father2;
        Node * GP;
        Node * brother;
      
        father=Nodo1->getFather();
        if(father->getNumberOfSons()==2){ 
           PosNodo= father->getSonPosition(Nodo1);
           brother = father->getSon(PosNodo==0?1:0);

           if (brother->hasDistanceToFather()) {
               distancetofather = brother->getDistanceToFather();
           }
           
          
           father->removeSon(brother);
           GP = father->getFather();
           GP->setSon(GP->getSonPosition(father),brother);
           
           if(father->hasDistanceToFather()) {
                   distancetofather+=father->getDistanceToFather();
           }
               
           brother->setDistanceToFather(distancetofather);

         }else{ 

           PosNodo= father->getSonPosition(Nodo1);
           brother = father->getSon(PosNodo==0?1:0);

           father->removeSon(Nodo1); 

           father = new Node(NextIDNode++);
           father->addSon(Nodo1);
           
         }

         distancetofather=0;
         
         father2 = Nodo2->getFather();
         
         if(Nodo2->hasDistanceToFather()) {
             distancetofather = Nodo2->getDistanceToFather();
         }

         father2->setSon(father2->getSonPosition(Nodo2),father);
         father->setDistanceToFather(distancetofather/2);
         
         
         father->addSon(Nodo2);
         
         Nodo2->setDistanceToFather(distancetofather/2);

         
        
        tr = tree;

        tr = (Tree*) tree;

        delete tl;

        return tr;
        
    }


    int  CRO::SPRvalide (Node* N1, Node* N2) {
        if (!N2->hasFather()) return 0;
        if (N1->getFather()==N2->getFather()) return 0;
        if (N1->getFather()==N2) return 0;
        if (N1 == N2) return 0;
        
        return 1;
    }

    Tree* CRO::TBR(Tree* tr)
    {

        vector<int> nodosIDs;
        tl = new NNIHomogeneousTreeLikelihood(*tr, *sites, model, rateDist);    
        
        tl->initialize();
        TreeTemplate<Node>  tre(tl->getTree());
        TreeTemplate<Node>* tree = new TreeTemplate<Node>(tre);
        
        Node * nodo;
        Node * nodoi;
        Node * nodoj;
        Node * parentNode;
        Node * nodoSubTree;

        nodosIDs = tree->getNodesId();
        
        do{
                  do{
                          nodo = tree->getNode(RandomTools::pickOne(nodosIDs,true));
                  }while(nodo->isLeaf() || nodo->getNumberOfSons() < 2 || !nodo->hasFather());

                  if (RandomTools::flipCoin()) {
                          nodoi= nodo->getSon(0); nodoj= nodo->getSon(1);
                  }else{
                          nodoi= nodo->getSon(1); nodoj= nodo->getSon(0);
                  }
                  
         }while(nodoi->isLeaf());

         TreeTemplate<Node> * subtree = new TreeTemplate<Node>(nodoi);
         subtree->resetNodesId();
         nodosIDs = subtree->getNodesId();
         do{
                nodoSubTree = subtree->getNode(RandomTools::pickOne(nodosIDs,true));
                
         }while(nodoSubTree->isLeaf());

         subtree->rootAt(nodoSubTree->getId());

         parentNode = nodo->getFather();
         parentNode->setSon(parentNode->getSonPosition(nodo),nodoj);
         tree->resetNodesId();
         nodosIDs = tree->getNodesId();
         do{
                  nodo = tree->getNode(RandomTools::pickOne(nodosIDs,true));
         }while(nodo->isLeaf());


         int posSon;
         if (RandomTools::flipCoin()) posSon=0; else posSon=1;

         Node * newNode = new Node();
         newNode->addSon(nodo->getSon(posSon));
         newNode->addSon(subtree->getRootNode());

         nodo->setSon(posSon,newNode);
         tree->resetNodesId();

         
         tr = (Tree*) tree;
         delete tl;
         
        return tr;

        }

        Tree* CRO::PDG(Tree * t1, Tree * t2 ) {

        tl = new NNIHomogeneousTreeLikelihood(*t1, *sites, model, rateDist); 
        tl->initialize();
        TreeTemplate<Node>  tre1(tl->getTree());
        TreeTemplate<Node>* tree1 = new TreeTemplate<Node>(tre1);
        //TreeTemplate<Node>* tree1 = TreeTemplateTools::parenthesisToTree(TreeTemplateTools::treeToParenthesis(t1));

        tl = new NNIHomogeneousTreeLikelihood(*t2, *sites, model, rateDist); 
        tl->initialize();
        TreeTemplate<Node>  tre2(tl->getTree());
        TreeTemplate<Node>* tree2 = new TreeTemplate<Node>(tre2);

        //TreeTemplate<Node>* tree2 = TreeTemplateTools::parenthesisToTree(TreeTemplateTools::treeToParenthesis(t2));


            
        Node *Nodo1;
        Node *Nodo2;
            Node *SubTree;
            vector<string> leaves;
            vector<int> nodosIDs;

        int numberOfLeavesT1 = tree1->getNumberOfLeaves();
            nodosIDs = tree1->getNodesId();

             bool b=true;
         do{
              Nodo1 = selectNodeToCross(tree1, nodosIDs );
              if(TreeTemplateTools::getNumberOfLeaves(*Nodo1) < numberOfLeavesT1-4){

                              SubTree=TreeTemplateTools::cloneSubtree<Node>(*Nodo1);
                  b=false;
              }
             }while(b);


           
        leaves = TreeTemplateTools::getLeavesNames(*SubTree);
            
    
            
        for (unsigned int i = 0; i < leaves.size(); i++)
                TreeTemplateTools::dropLeaf(*tree2,leaves[i]);


        Nodo2 = selectNodeToCross(tree2, tree2->getNodesId()); //Extrae PUNTERO refrencia a NODO, NO es copia


        Node * father = Nodo2->getFather();
        double distancetofather = Nodo2->getDistanceToFather();
        int PosNodo2= father->getSonPosition(Nodo2);


        Node * nodo = new Node();
        
        nodo->addSon(SubTree); nodo->addSon(Nodo2);
        Nodo2->setDistanceToFather(distancetofather/2);

        father->setSon(PosNodo2, nodo);
        nodo->setDistanceToFather(distancetofather/2);



        tree2->resetNodesId();
        delete tl, tree1;

        return (Tree*)tree2;
    }

    Node * CRO::selectNodeToCross(TreeTemplate<Node> * tree_, vector<int> nodosIDs ){
    Node * nodoSel;
        do{
            nodoSel =  tree_->getNode( nodosIDs[rand()% (nodosIDs.size() - 1)]);
        }while(!nodoSel->hasFather() || nodoSel->isLeaf() );

    return  nodoSel;
    }

    

    void CRO::init_pop_with_bionj()
        {
            cout << "Computing tree..." << endl;
            DistanceEstimation distanceMethod(model, rateDist, sites);
            DistanceMatrix* distances = distanceMethod.getMatrix();
 
            BioNJ bionj(*distances);
            Tree* tre = bionj.getTree();
            //Molecule m(tre);
            //Newick treeWriter;
            //treeWriter.write(*tre, "rbcl_BioNJ.dnd");
            //Tree* t;
            for(int i = 0; i < init_pop; i++)
            {
                
               Tree* t_new = tre->clone();
               t_new = SPR(t_new);
               Molecule mn(t_new);
               if(mn.structure->getNumberOfLeaves() == nos)
                {
                        pop.push_back(mn);
                    
                        update(i);    
                }
               

                
                //myfile << "Random tree : "<< i << " " << ofv(t) << endl;
                

            }

            for(vector<Molecule>::iterator it = pop.begin(); it != pop.end(); ++it)
            {
                it->mp = this->MP_value(*it);
                it->ml = this->ML_value(*it);

                it->pe = updatePE(it->mp,it->ml);
                it->ke = 100;
                it->update();
                if(!is_first_molecule)
                {
                    optimal = *it;
                }
                else if(it->pe < optimal.pe)
                {
                    optimal = *it;
                }

            }

        }

        void CRO::init_pop_with_rand_tree()
        {
            for(int i = 0; i < init_pop; i++)
            {
                
                TreeTemplate<Node>* tree = TreeTemplateTools::getRandomTree(species, true);
                tree->setBranchLengths(0.05);
                Tree* t = (Tree*) tree;
                //t = NNI(t);
                myfile << "Random tree : "<< i << " " << ofv(t) << endl;
                Molecule mn(t);
                pop.push_back(mn);

            }


            for(vector<Molecule>::iterator it = pop.begin(); it != pop.end(); ++it)
            {
                
                it->mp = this->MP_value(*it);
                it->ml = this->ML_value(*it);

                it->pe = it->mp;

                it->ke = 100;
                it->update();
                if(!is_first_molecule)
                {
                    optimal = *it;
                }
                else if(it->pe < optimal.pe)
                {
                    optimal = *it;
                }

            }
        }

        void CRO::readingTreesFromAFile()
	    {
	        //each tree is seperated with $ sign
            vector<string> treeString;
            bool isTaken[1000];
            int randIndex;
	        int tree_count = 0;
            printf("Enter tree file name: ");
            string s;
            cin >> s;
	        ifstream t(s);
	            stringstream buffer;
	            buffer << t.rdbuf();

	            string inTree = buffer.str();

	            stringstream ss( inTree);
	            vector<string> result1;
	            // TreeTemplate<Node>* tree19 = TreeTemplateTools::parenthesisToTree("((Zygnema_pe,((((Lycopodium,Anthoceros),(((((Larix,Picea_pun),Ginkgobil),Podocarpus),((Victoria,Austrobai),(((Saururus,((Avena,Carex),(Aloe,Iris))),(Tasmannia,Hedycarya)),(Platanus,(((Amaranthu,Heuchera),(Gunnera,((Nicotiana,Gentiana),(Dipsacus,Cornuscan)))),(Floerkea,(NothofBal,Celtis))))))),(Angiopter,(Osmunda_ci,(Marsilea,((Cheiropleu,Loxoma),(Monachosor,((Xiphopteri,Nephrolepi),Asplenium)))))))),((((Chara_conn,(Hookeria_2,((Hedwigia_7,Leucodon_5),(Fissidens_,Ptychomitr)))),Tetraphis_),Sphagnum_J),((Porella_4,Bazzania_J),(Metzgeria_,((Dumortiera,Marchantia),Conocephal))))),Sirogonium)),(Volvox_ro,(Chlor_ell,Pyrami_pa)),Anabaena_s);");
	            // //cout << TreeTemplateTools::treeToParenthesis(*tree9) << endl;
	            // Molecule fs(tree19);
	            //ss.good
	            while( ss.good() )
	            {
                    if(tree_count == 160)
                         break;
	                string substr;
	                getline( ss, substr,';');
	                substr = substr+";";
	                //cout << substr <<" " << substr.size() << endl;
	                //result1.push_back( substr);
	                 //cout << tree_count<< endl;
	                string::iterator end_pos = remove(substr.begin(), substr.end(), ' ');
	                substr.erase(end_pos, substr.end());

                    treeString.push_back(substr);
                    tree_count++;
	                //TreeTemplate<Node>* tree9 = TreeTemplateTools::parenthesisToTree(substr);
	                //cout << TreeTemplateTools::treeToParenthesis(*tree9) << endl;
	                //Tree * tree10 = (Tree*) tree9;

	                //Molecule mn(tree10);
	                //pop.push_back(mn);
	                //update(tree_count++);
	                

	            }
                cout << "Tree string " << treeString.size() << tree_count << endl;
                for(int i = 0; i < treeString.size(); i++)
                {
                    isTaken[i] = false;
                }
                tree_count = 0;
                while(1)
                {
                    randIndex = rand()%treeString.size() - 1;
                    if(!isTaken[randIndex])
                    {
                        //cout << treeString[randIndex] << endl;
                        TreeTemplate<Node>* tree9 = TreeTemplateTools::parenthesisToTree(treeString[randIndex]);
                        //cout << TreeTemplateTools::treeToParenthesis(*tree9) << endl;
                        Tree * tree10 = (Tree*) tree9;

                        Molecule mn(tree10);
                        pop.push_back(mn);
                        update(tree_count++);
                        isTaken[randIndex] = true;
                        tree_count++;
                    }

                    if(tree_count == 30)
                        break;
                    
                }
	    }

        void CRO::printPopulation()
        {
            double min = 1000000000.0;
            double mp, ml;


            myfile << "Molecules are : " << endl;
            for(vector<Molecule>::iterator it = pop.begin(); it != pop.end(); ++it)
            {
                if(it->pe < min)
                {
                    mp = it->mp;
                    ml = it->ml;
                    min = ml;
                }
               
                result.push_back(make_pair(it->mp, it->ml));
                numOfHitsCollection.push_back(it->num_of_hits);
                //myfile << "Molecule number is : " << it->ID << endl;
            }
            cout << "Optimal MP " << mp << " and Ml " << ml;
            
        }

    	  // DiscreteRatesAcrossSitesTreeLikelihood* tl,
		  // const ParameterList& parameters,
		  // OptimizationListener* listener,
		  // double tolerance,
		  // unsigned int tlEvalMax,
		  // OutputStream* messageHandler,
		  // OutputStream* profiler,
		  // unsigned int verbose,
		  // const std::string& optMethodDeriv
		  // NNIHomogeneousTreeLikelihood* tl,
		  // const ParameterList& parameters,
		  // bool optimizeNumFirst,
		  // double tolBefore,
		  // double tolDuring,
		  // unsigned int tlEvalMax,
		  // unsigned int numStep,
		  // OutputStream* messageHandler,
		  // OutputStream* profiler,
		  // bool reparametrization,
		  // unsigned int verbose,
		  // const std::string& optMethodDeriv,
		  // unsigned int nStep,
		  // const std::string& nniMethod

        void CRO::branch(Tree * t,int id)
        {
            NNIHomogeneousTreeLikelihood* bo = new NNIHomogeneousTreeLikelihood(*t, *sites, model, rateDist); 
            bo->initialize();
            myfile << "Tree : " << blo_tree << " Log likelihood before: " << bo->getLogLikelihood() << endl;
            
            //bo = OptimizationTools::optimizeTreeNNI(bo, bo->getParameters(), true, 100, 0.000001, 300, 1, 0, 0, false, 3);
            bo = OptimizationTools::optimizeTreeNNI(bo, bo->getParameters(), true, 100, 0.001, 100, 1, 0, 0, false, 0,"newton", 3);
            //bo = OptimizationTools::otimizeBranchlengthParameters(bo, bo->getParameters(), );
            //cout << "Value from optimal tree " << ML_value_from_tree(optimal.structure) << endl;
            myfile << "Log likelihood after: " << bo->getLogLikelihood() << endl;
            myfile << endl;
            //myfile << "Log likelihood after: " << bo->getLogLikelihood() << endl;
            bo->getParameters().printParameters(cout);
            TreeTemplate<Node> mlTree(bo->getTree());
            treeWriter.write(mlTree, "branchOptimizedTree"+to_string(id)+".dnd");
            //myfile << "Final mp,ml value after branchlength optimization using our mp, ml value function" << endl;
            //myfile << MP_value(optimal) << "is mp value and " << ML_value(optimal) << "is ml value" << endl;

            delete bo;
        }


        void CRO::reform2()
        {
            int toOpt = 0;
            sort(pop.begin(), pop.end());

            cout << "Enter number of trees to optimize" << endl;
            cin >> toOpt;
            //pop.resize(toOpt);
            bool track[pop.size()+1];
            memset(track, false , sizeof(track));
            
            for(int i = 0; i < toOpt; i++)
            {
            	branch(pop[i].structure, pop[i].ID);
            }

        }


        Tree *  CRO::synthesis_with_contribution(Tree * tr1, Tree * tr2)
        {  
            

            tl = new NNIHomogeneousTreeLikelihood(*tr1, *sites, model, rateDist); 
            tl->initialize();
            TreeTemplate<Node>  tre(tl->getTree());
            TreeTemplate<Node>* tree = new TreeTemplate<Node>(tre);
            
            delete tl;

            //printf("Hello Bio++!! How are you?\n");           
            //TreeTemplate<Node> * tree2 = TreeTemplateTools::parenthesisToTree("(((A, B), C), ((D, E), F));");
            //printf("Number of leafs %lu\n", tree->getNumberOfLeaves() );


            tl = new NNIHomogeneousTreeLikelihood(*tr2, *sites, model, rateDist); 
            tl->initialize();
            TreeTemplate<Node>  tre2(tl->getTree());
            TreeTemplate<Node>* tree2 = new TreeTemplate<Node>(tre);
            

            //vector<string> speciesForSecondTree = tree->getLeavesNames();
            //TreeTemplate<Node> * tree2 = TreeTemplateTools::getRandomTree(speciesForSecondTree,true);
            //tree2->setBranchLengths(0.05);
            Node * root = tree->getRootNode();
            cout << root->getId() << endl;
            
            vector<int> nodeIds = tree->getNodesId();
            Node * splitNode;
            bool found = false;
            double leafNumber;
            int firstSonLeaves, secondSonLeaves, length, firstId, secondId;
            length = nodeIds.size();
            leafNumber = tree->getNumberOfLeaves() * 0.25;
            Node * b1;
            Node * b2;
            for(int i = 0; i < length; i++)
            {
                firstSonLeaves = TreeTemplateTools::getNumberOfLeaves(*tree->getNode(nodeIds[i]));
                if(firstSonLeaves >= leafNumber)
                {
                    //b1 = tree->getNode(nodeIds[i]);
                    firstId = nodeIds[i];
                    break;
                }
                
            }
            nodeIds = tree2->getNodesId();
            leafNumber = tree2->getNumberOfLeaves() * 0.25;
            nodeIds = tree2->getNodesId();
            length = nodeIds.size();
            for(int i = 0; i <= length; i++)
            {
                secondSonLeaves = TreeTemplateTools::getNumberOfLeaves(*tree2->getNode(nodeIds[i]));
                if(secondSonLeaves >= leafNumber)
                {
                    //b2 = tree2->getNode(nodeIds[i]);
                    secondId = nodeIds[i+1];
                    break;
                }
                
            }

                cout << firstId << " " << secondId << endl;
            
                    b1 = TreeTemplateTools::cloneSubtree<Node>(*tree->getNode(firstId));
                    b2 = TreeTemplateTools::cloneSubtree<Node>(*tree2->getNode(secondId));
                    Node * rightTree = new Node();
                    rightTree->addSon(b1);
                    rightTree->addSon(b2);

                    TreeTemplate<Node> * contaminatedTree = new TreeTemplate<Node>(rightTree);
                    contaminatedTree->resetNodesId();
                    vector<string> contLeaves = contaminatedTree->getLeavesNames();

                    // for(int z = 0; z < contLeaves.size(); z++)
                    // {
                    // 	cout << contLeaves[z] << endl;
                    // }
                    // cout << "---> updated" << endl;

                    //before
                    //reform1
                    map<string,int> randomSpecies;
                    vector<string> leaves = tree->getLeavesNames();
                    vector<string> forRandom;
                    for(int j = 0; j < leaves.size(); j++)
                    {
                        randomSpecies[leaves[j]] = 0;
                    }
                    
                    for(int j = 0; j < contLeaves.size(); j++)
                    {
                        randomSpecies[contLeaves[j]]++;
                    }
                    vector<Node*> allNodes = contaminatedTree->getNodes();
                    
                    for (map<string,int>::iterator it=randomSpecies.begin(); it!=randomSpecies.end(); ++it)
                    {
                    	bool done = false;
                        if(it->second == 2)
                        {
                        	for (map<string,int>::iterator it1=randomSpecies.begin(); it1!=randomSpecies.end(); ++it1)
                        	{
                        		if(it1->second == 0)
                        		{
                        			for(int f = 0; f < allNodes.size(); f++)
                        			{
                        				if(allNodes[f]->hasName())
                        				{
                        					if(allNodes[f]->getName() == it->first)
	                        				{
	                        					allNodes[f]->setName(it1->first);
	                        					randomSpecies[it1->first] = 1;
	                        					done = true;
	                        					break;
	                        				}
                        				}
                        				
                        			}
                        		}
                        		if(done)
                        			break;

                        	}
                           
                        }
                    }
                    for (map<string,int>::iterator it=randomSpecies.begin(); it!=randomSpecies.end(); ++it)
                    {
                    	if(it->second == 0)
                    	{
                    		forRandom.push_back(it->first);
                    	}
                    }

                    TreeTemplate<Node>* rt = TreeTemplateTools::getRandomTree(forRandom, true);
                    rt->setBranchLengths(0.05);

                    Node * finalRoot = new Node();
                    finalRoot->addSon(rightTree);
                    finalRoot->addSon(rt->getRootNode());
                    rightTree->setDistanceToFather(0.05);
                    rt->getRootNode()->setDistanceToFather(0.05);

                    TreeTemplate<Node> * ft = new TreeTemplate<Node>(finalRoot);
                    ft->rootAt(finalRoot->getId());
                    ft->resetNodesId();
                    //cout << TreeTemplateTools::treeToParenthesis(*ft) << endl;
                    // leaves = ft->getLeavesNames();
                    // for(int z = 0; z < leaves.size(); z++)
                    // {
                    // 	cout << leaves[z] << endl;
                    // }
                    
                    //Delete duplicate species and add absent species
                    Tree * tr = (Tree*) ft;
                    delete tl;
         
                    return tr;
                      
         }
         Tree*  CRO::decomposition_with_contribution(Tree * tr1)
         {  
        

            tl = new NNIHomogeneousTreeLikelihood(*tr1, *sites, model, rateDist); 
            tl->initialize();
            TreeTemplate<Node>  tre(tl->getTree());
            TreeTemplate<Node>* tree = new TreeTemplate<Node>(tre);
        

            //printf("Hello Bio++!! How are you?I am decomposition\n");           
            //TreeTemplate<Node> * tree2 = TreeTemplateTools::parenthesisToTree("(((A, B), C), ((D, E), F));");
            //printf("Number of leafs %lu\n", tree->getNumberOfLeaves() );


            Node * root = tree->getRootNode();
            cout << root->getId() << endl;
            vector<Node*> allNodes = tree->getNodes();
            vector<int> nodeIds = tree->getNodesId();
            Node * splitNode;
            bool found = false;
            double leafNumber;
            int firstSonLeaves, secondSonLeaves, length, firstId, secondId;
            length = nodeIds.size();
            leafNumber = tree->getNumberOfLeaves() * 0.25;
            Node * b1;
            Node * b2;
            for(int i = 0; i < length; i++)
            {
                firstSonLeaves = TreeTemplateTools::getNumberOfLeaves(*tree->getNode(nodeIds[i]));
                if(firstSonLeaves >= leafNumber)
                {
                    //b1 = tree->getNode(nodeIds[i]);
                    firstId = nodeIds[i];
                    break;
                }
                
            }
            
                cout << firstId  << endl;
            
                    b1 = TreeTemplateTools::cloneSubtree<Node>(*tree->getNode(firstId));
                    
                    
                    map<string,bool> randomSpecies;
                    vector<string> leaves = tree->getLeavesNames();
                    vector<string> forRandom;
                    for(int j = 0; j < leaves.size(); j++)
                    {
                        randomSpecies[leaves[j]] = false;
                    }
                    vector<string> firstLeaves = TreeTemplateTools::getLeavesNames(*b1);
                    
                    for(int j = 0; j < firstLeaves.size(); j++)
                    {
                        randomSpecies[firstLeaves[j]] =  true;
                    }

                    for (map<string,bool>::iterator it=randomSpecies.begin(); it!=randomSpecies.end(); ++it)
                    {
                        if(it->second == false)
                            forRandom.push_back(it->first);
                    }

                    TreeTemplate<Node>* rt = TreeTemplateTools::getRandomTree(forRandom, true);
                    rt->setBranchLengths(0.05);

                    Node * finalRoot = new Node();
                    finalRoot->addSon(b1);
                    finalRoot->addSon(rt->getRootNode());
                    b1->setDistanceToFather(0.05);
                    rt->getRootNode()->setDistanceToFather(0.05);

                    TreeTemplate<Node> * ft = new TreeTemplate<Node>(finalRoot);
                    ft->rootAt(finalRoot->getId());
                    //cout << TreeTemplateTools::treeToParenthesis(*ft) << endl;
                    ft->resetNodesId();
                    Tree * tr = (Tree*) ft;
                    delete tl;
         
                    return tr;

        }




        


