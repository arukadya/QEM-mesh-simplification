#ifndef UnionFind_hpp
#define UnionFind_hpp

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

struct UnionFind{
    std::vector<int> parents;
    UnionFind(unsigned int size)
    {
        parents.resize(size);
        for(auto &x:parents)x = -1;
    }
    unsigned int find(unsigned int id,int &rank)
    {
        while(parents[id] >= 0)
        {
            id = parents[id];
        }
        rank = -parents[id];
        return id;
    }
    void union_ab(unsigned int id_a,unsigned int id_b)
    {
        int rank_a,rank_b;
        unsigned int root_a,root_b;
        root_a = find(id_a,rank_a);
        root_b = find(id_b,rank_b);
        if(root_a != root_b)
        {
            if(rank_a > rank_b)parents[id_b] = root_a;
            else if(rank_a < rank_b)parents[id_a] = root_b;
            else
            {
                parents[id_b] = root_a;
                --parents[root_a];
            }
        }
    }
    void print_parents()
    {
        for(int i=0;i<parents.size();++i)std::cout << "i,parents = " << i << "," << parents[i] <<std::endl;
    }
    void print_ranks()
    {
        for(int i=0;i<parents.size();++i)
        {
            int rank = 0;
            find(i, rank);
            std::cout << "i,ranks = " << i << "," << rank <<std::endl;
        }
    }
};

#endif /* UnionFind_hpp */
