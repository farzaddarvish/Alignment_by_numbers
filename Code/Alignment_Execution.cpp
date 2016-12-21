
//
//  Alignment_Execution.cpp
//  Al_by_num_v2_array
//
//  Created by Avraam Tapinos on 20/02/2016.
//  Copyright (c) 2016 Avraam_Tapinos. All rights reserved.
//

#include "Alignment_Execution.h"


//////////////////////////////////////////////////////////////////////////////////////////////
//                                         Private:                                         //
//////////////////////////////////////////////////////////////////////////////////////////////
void Alignment_Execution::executeFastaParser(std::vector<std::string> fafilepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> * RefReads, int kmer){
    FastaParser2 fa(fafilepath, repmeth, tranmeth, tranlev, RefReads, kmer);
    reference = fa.returnRead();
};

//Call FastqParser to extract information from fats a file
void Alignment_Execution::executeFastqParser(std::string fqfilepath, std::string repmeth, std::string tranmeth, int tranlev, std::vector<NGS> *NuSequences, int kmer){
    FastqParser fq(fqfilepath, repmeth, tranmeth, tranlev, NuSequences, kmer);
};


//////////////////////////////////////////////////////////////////////////////////////////////
//                                         Public:                                          //
//////////////////////////////////////////////////////////////////////////////////////////////
//Destractor
Alignment_Execution::~Alignment_Execution(){};


//Constractor
//Initate alignemnt process
//Alignment_Execution::Alignment_Execution(std::string fafilepath, std::string fqfilepath, std::string outdir, std::string repmeth ,std::string tranmeth, int tranlev, int kmer, int knn, bool sens, std::string repdir, bool sr){
Alignment_Execution::Alignment_Execution(std::vector<std::string> fafilepath, std::string fqfilepath, std::string outdir, std::string repmeth ,std::string tranmeth, int tranlev, int kmer, int knn, bool sens, std::string repdir, bool sr){
    //Variables
    //Vector containing subsequence NGS reads from reference
    std::vector<NGS> RefReads;
    //Vector containing  short NGS reads from fastq file
    std::vector<NGS> NuSequences;
    
    
    std::size_t found =fqfilepath.find_last_of(std::string("/"));
    std::size_t found1 =fqfilepath.find_last_of(std::string("."));
    found++;
    std::string sns;
    
    std::size_t foundsl = fqfilepath.find_last_of(std::string("/"));
    std::string rep_dir = repdir;
    
    
    struct stat sb;
    srand (time(nullptr));
    
    if (!(stat(rep_dir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)))
    {
        mkdir(rep_dir.c_str(), 0777);
    }
    
    //sns variable to use to the alignment file name (Tag if an exhaustive or not search performaance is used)
    if(sens == true){
        sns = "exhaustive_settings";
    }
    else{
        sns = "non_exhaustive_settings";
        std::cout << "> Exhaustive comparison: OFF\n";
    }
    
    if(sr == true){
        sns+="_shuffle_reads";
    }
    else{
        sns+="_non_shuffle_reads";
    }
    
    
    std::string rkmer, rtrans;
    
    if(kmer < 100){
        rkmer = "0"+ std::to_string(kmer);
    }
    else{
        rkmer = std::to_string(kmer);
    }
    
    if(tranlev < 100){
        if(tranlev < 10){
            rtrans = "00" + std::to_string(tranlev);
        }
        else{
            rtrans = "0" + std::to_string(tranlev);
        }
    }
    else{
        rtrans = std::to_string(tranlev);
    }
    
    //Sam file directory
    
    rep_dir+=fqfilepath.substr(foundsl+1, fqfilepath.size()-foundsl+1)+"_"+repmeth+"_"+rkmer+ "kmer_"+tranmeth+"_"+ rtrans + "complvl_" + sns + "_KNN_"+std::to_string(knn)+ ".txt";
    
    std::ofstream o(rep_dir);
    
    
    //Write info to output file
    o << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    o << "%%%%%%%%%%                                  Alignment by numbers                                  %%%%%%%%%%\n";
    o << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    
    
    o << "> Alignment process\n";// << fqfilepath.substr(found, (fqfilepath.size()-found)-(fqfilepath.size()-found1) )<< " reads dataset.\n";
    o << "> Parameters:\n";
    o << "> Representation method: " << repmeth << " || Kmer size: " << kmer << "\n";
    o << "> Transformation method: " << tranmeth << " || Transformation level: " << tranlev << "\n";
    o << "> Number of K-NN neighbours considered: " << knn << "\n";
    
    
    //Print alignment info
    std::cout << "> Alignment process\n";// << fqfilepath.substr(found, (fqfilepath.size()-found)-(fqfilepath.size()-found1) )<< " reads dataset.\n";
    std::cout << "> Parameters:\n";
    std::cout << "> Representation method: " << repmeth << " || Kmer size: " << kmer << "\n";
    std::cout << "> Transformation method: " << tranmeth << " || Transformation level: " << tranlev << "\n";
    std::cout << "> Number of K-NN neighbours considered: " << knn << "\n";
    if(sens == true){
        
        o << "> Exhaustive comparison: ON\n";
        std::cout << "> Exhaustive comparison: ON\n";
    }
    else{
        
        o << "> Exhaustive comparison: OFF\n";
        std::cout << "> Exhaustive comparison: OFF\n";
    }
    
    
    if(sr == true){
        o << "> Randomly shuffle reads: ON\n";
        std::cout << "> Randomly shuffle reads: ON\n";
    }
    else{
        o << "> Non randomly shuffle reads: ON\n";
        std::cout << "> Non randomly shuffle reads: ON\n";
    }
    
    
    
    o << "\n> Process:\n";
    std::cout << "\n> Process:\n";
    //Staring time for alignment process
    clock_t Alstart = clock();
    
    //Staring time for fastaparser process
    clock_t fastart = clock();
    executeFastaParser(fafilepath, repmeth, tranmeth, tranlev, &RefReads, kmer);
    //Ending time for fastaparser process
    clock_t faend = clock();
    int readssize = (int)RefReads.size();
    
    o << "> Importing reference genome and creating transformations for "<< readssize << " kmers elapse time: " << (faend - fastart)/CLOCKS_PER_SEC << " seconds\n";
    std::cout << "> Importing reference genome and creating transformations for "<< readssize << " kmers elapse time: " << (faend - fastart)/CLOCKS_PER_SEC << " seconds\n";
    
    //Build the vector of vectors base don the blocksize variable
    if(sr == true){
        //Randomly shuffled the RefRead vector
        std::random_shuffle (RefReads.begin(), RefReads.end());
    }
    
    //Staring time for building the VP-Forest
    clock_t trstart = clock();
    
    
    std::cout << "> Start creating VPtree for " << RefReads.size() <<" reads.\n";

    //MVpTree2 <NGS,  Fwdist> VPTree;
    //VPTree.create(RefReads, 0);
    //VPTree.create(RefReads, 20);
    
    
    
    //MVpTree <NGS,  Fwdist> VPTree;
    //VpTree <NGS,  Fwdist> VPTree;
    //vptree <NGS,  Fwdist> VPTree;
    //int num_reads = static_cast<int>(RefReads.size());
    ngstree VP;
    VP.addMany(RefReads.begin(), RefReads.end());
    //VP.update();

    
    

    //VpTree1 <NGS,  Fwdist1> VPTree;
    //VpTree2 <NGS,  Fwdist1> VPTree;
    //VPTree.create(RefReads);
    
    
    //Ending time for building the VP-Forest
    clock_t trend = clock();
    
    o << "> Creating tree for  " << readssize <<" reads elapse time: " << (trend - trstart)/CLOCKS_PER_SEC << " seconds.\n";
    std::cout << "> Creating tree for  " << readssize <<" reads elapse time: " << (trend - trstart)/CLOCKS_PER_SEC << " seconds.\n";
    
    o << "> Start importing reads from the fastq file.\n";
    std::cout << "> Start importing reads the from fastq file.\n";
    clock_t fqstart = clock();
    executeFastqParser(fqfilepath, repmeth, tranmeth, tranlev, &NuSequences, kmer);
    clock_t fqend = clock();
    o << "> Importing "<< NuSequences.size() <<" reads from the fastq file and creating transformations elapse time: " << (fqend - fqstart)/CLOCKS_PER_SEC << " seconds.\n";
    std::cout << "> Importing "<< NuSequences.size() <<" reads from the fastq file and creating transformations elapse time: " << (fqend - fqstart)/CLOCKS_PER_SEC << " seconds.\n";
    
    
    Process2dArray<double> CA;
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    Representation refrep, ngsrep;

    
    

    clock_t start = clock();
    int readcount= 0, countper, percentageval= 0;
    int maxsearch = 0, minsearch = std::numeric_limits<int>::max();
    double  searchcountv = 0.0 ;
    countper = (int)NuSequences.size()/10;
    std::cout << "> " << percentageval << "%\n";
    double abansearch = 0.0, avsearch = 0.0;
    //int vln = 0, minvln = std::numeric_limits<int>::max(), maxvln = 0 , minvlnptimes = 0, maxvlnptimes = 0;
    //int totalsearch = 0, mintotalsearch = std::numeric_limits<int>::max(), maxtotalsearch = 0, mintotalsearchptimes = 0, maxtotalsearchptimes = 0;
    //int leavecheck = 0, minleavecheck = std::numeric_limits<int>::max(), maxleavecheck = 0, minleavecheckptimes = 0, maxleavecheckptimes = 0;
    
    double At = 0.0;
    std::vector<std::string> minsearchreads, maxsearchreads;
    for (int i = 0; i < NuSequences.size();i++){
        
        readcount++;
        if(readcount == countper){
            percentageval += 10;
            std::cout << "> " << percentageval << "%\n";
            
            readcount = 0;
        }
        
        //std::vector <NGS> neighbours;
        //std::vector <double> distances;
        
        //clock_t t = clock();
        //VPTree.search(NuSequences.at(i), knn, &neighbours, &distances, tau);
        const NGS &query = NuSequences.at(i);
        int sc = 0, absc = 0;
        std::vector<const NGS *> neighbours = VP.nearestNeighbors(query, knn, &sc, &absc);
        //std::vector<const NGS *> neighbours = VP.approxNearestNeighbors(query, knn);
        
        //t = clock() - t;
        abansearch += (((double)absc)/((double)NuSequences.size()));
        //std::cout << i <<") KNN search visited " << sc << " vantage points took " << ((float)t)/CLOCKS_PER_SEC << " seconds\n";
       // At +=(((double)t)/CLOCKS_PER_SEC);
        
        
        
        if(maxsearch < sc){
            maxsearch = sc;
        }
        if(minsearch > sc){
            minsearch = sc;
        }
        
        searchcountv += (((double)sc)/((double)NuSequences.size()));
        
        
        
        
        
        
        if(sens == false){
            if(neighbours.size() > 1){
                
                std::string srseq;
                srseq = NuSequences.at(i).returnRead();
                double **ngsrepvec;
                ngsrep.createRepresentation(repmeth, srseq, ngsrepvec);
                
                
                double bestdist = std::numeric_limits<double>::infinity();
                int ydim = (int)NuSequences.at(i).returnY();
                std::vector<const NGS *> newneighbours;
                
                for (int j = 0; j < neighbours.size(); j++){
                    if(neighbours.at(j) == NULL){
                        break;
                    }
                    int xdim;
                    std::string refseq;
                    refseq = neighbours.at(j)->returnRead();
                    double **refrepvec;
                    refrep.createRepresentation(repmeth, refseq, refrepvec);
                    
                    srseq = NuSequences.at(i).returnRead();
                    
                    if(refseq.size() <= srseq.size()){
                        xdim = (int)refseq.size();
                    }
                    else{
                        xdim = (int)srseq.size();
                    }
                    double tempdist = Nudist(refrepvec, ngsrepvec, ydim, xdim);
                    
                    CA.DeleteArray(refrepvec);
                    if(tempdist <= bestdist){
                        if(tempdist < bestdist){
                            bestdist = tempdist;
                            std::vector<const NGS *>().swap(newneighbours);
                        }
                        newneighbours . push_back(neighbours.at(j));
                    }
                }
                
                
                CA.DeleteArray(ngsrepvec);
                
                std::vector<const NGS *>().swap(neighbours);
                neighbours = newneighbours;
                std::vector<const NGS *>().swap(newneighbours);
                
            }
        }
        
        
        
        std::string cigarstring;
        uint16_t bestscore = 0;
        int bestalpos = 0;
        bool bestdir = false;
        NuSequences.at(i).clearAlVec();
        
        for (int j = 0; j < neighbours.size(); j++){
            std::string refcre;
            int startbinsize = 10;
            int pos;
            if(neighbours.at(j)->returnDir() == false){
                pos = neighbours.at(j)->returnIdx();
            }
            else{
                pos = neighbours.at(j)->returnIdx()-(NuSequences.at(i).returnFullX()-kmer);
            }
            int start = pos-startbinsize;
            int endpos;
            int endbinsize = 10;
            if(start < 0){
                start = 0;
                startbinsize = abs(pos-0);
            }
            if(pos + (int) NuSequences.at(i).returnRead().size() + endbinsize >= (int) reference.size())
            {
                endbinsize = ((int)reference.size() - endpos);
            }
            
            int binsize = startbinsize + endbinsize + (int) (NuSequences.at(i).returnRead().size());
            refcre = reference[neighbours.at(j)->returnID()].substr(start, binsize);
            if(neighbours.at(j)->returnDir()==false){
                aligner.Align(NuSequences.at(i).returnRead().c_str(), refcre.c_str(), (int)refcre.size(), filter, &alignment);
            }
            else{
                aligner.Align(NuSequences.at(i).returnRvRead().c_str(), refcre.c_str(), (int)refcre.size(), filter, &alignment);
            }
            if(bestscore <= alignment.sw_score){
                if(bestscore < alignment.sw_score){
                    bestscore = alignment.sw_score;
                    NuSequences.at(i).clearAlVec();
                }
                NuSequences.at(i).setAlDirection(neighbours.at(j)->returnDir());
                NuSequences.at(i).setCigar(alignment.cigar_string);
                NuSequences.at(i).setAliPos(start + alignment.ref_begin);
                NuSequences.at(i).setAlRef(neighbours.at(j)->returnID());
            }
        }
        NuSequences.at(i).UniAl();
        
    }
    
    clock_t end = clock();
    
    o << "> Performing " << knn << "-NN neighbours search for " << NuSequences.size() <<" reads elapse time: " << ((double)(end - start))/CLOCKS_PER_SEC << " seconds.\n";
    std::cout << "> Performing " << knn << "-NN neighbours search for " << NuSequences.size() <<" reads elapse time: " << (end - start)/CLOCKS_PER_SEC << " seconds.\n";
    
    o << "> Minimum query search comparisons performed on the VP-tree: " <<  minsearch << "\n";
    std::cout << "> Minimum query search comparisons performed on the VP-tree: " <<  minsearch <<"\n";
    
    o << "> Maximum query search comparisons performed on the VP-tree: " << maxsearch << "\n";
    std::cout << "> Maximum query search comparisons performed on the VP-tree: " <<  maxsearch <<"\n";
    
    o << "> Average query search comparisons performed on the VP-tree: " << round(searchcountv) << ".\n";
    std::cout << "> Average query search comparisons performed on the VP-tree: " << round(searchcountv) << ".\n";
    
    o<< "> Average abandon search: " << round(abansearch) << "\n";
    std::cout << "> Average abandon search: " << round(abansearch) << "\n";
    
    o << "> Average search time: " << At/((double)NuSequences.size()) << " seconds\n";
    std::cout<< "> Average search time: " << At/((double)NuSequences.size()) << " seconds\n";
    std::cout<< "\n";
    
    
    std::string knnstr = std::to_string(knn);
    std::string kmerlendir, transintstr;
    if(kmer < 100){
        kmerlendir = "0"+std::to_string(kmer);
    }
    else{
        kmerlendir = std::to_string(kmer);
    }
    
    if(tranlev < 100){
        if(tranlev < 10){
            transintstr = "00" +std::to_string(tranlev);
        }
        else{
            transintstr = "0"+std::to_string(tranlev);
        }
    }
    else{
        transintstr = std::to_string(tranlev);
    }
    
    start = clock();
    //Initiate SamCreator class to generate a sam file for the alignent
    SamCreator samgenerator (fqfilepath,repmeth, tranmeth, kmerlendir, transintstr, &NuSequences, &reference, outdir, sns, knnstr);
    end =clock();
    o << "> Creating SAM file elapse time " << (end-start)/CLOCKS_PER_SEC << ".\n";
    std::cout << "> Creating SAM file elapse time " << ((double)(end-start))/CLOCKS_PER_SEC << ".\n";
    clock_t Alend = clock();
    
    o << "> Overall alignment time " <<  (Alend - Alstart)/CLOCKS_PER_SEC << " in seconds.\n";
    o << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
    o.close();
    std::cout << "> Overall alignment time " <<  (Alend - Alstart)/CLOCKS_PER_SEC << " in seconds.\n";
};
