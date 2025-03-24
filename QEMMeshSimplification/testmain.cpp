/*N 頂点 0 辺の無向グラフに Q 個のクエリが飛んできます。処理してください。

0 u v: 辺(u,v)を追加する。
1 u v: u,v が連結ならば 1、そうでないなら 0 を出力する。

 input
 4 7
 1 0 1
 0 0 1
 0 2 3
 1 0 1
 1 1 2
 0 0 2
 1 1 3
 
 output
 0
 1
 0
 1
*/

#include "UnionFind.hpp"
int main(int argc, const char * argv[])
{
    
    const char* InputFileName = "data.txt";
    std::ifstream Inputfile(InputFileName);
    if (!Inputfile.is_open()) {
        std::cerr << "Could not open the file - '"
             << InputFileName << "'" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::string line;
    std::vector<std::string> words(3);
    std::getline(Inputfile,line);
    std::stringstream ss_nums{line};
    std::getline(ss_nums, words[0], ' ');
    std::getline(ss_nums, words[1], ' ');
    int N,Q;
    N = stoi(words[0]);
    Q = stoi(words[1]);
    UnionFind uf(N);
//    std::cout << "input" <<std::endl;
    for(int i=0;i<Q;++i)
    {
        std::getline(Inputfile,line);
        std::stringstream ss_nums{line};
        for(int i=0;i<3;++i)
        {
            std::getline(ss_nums, words[i], ' ');
//            std::cout << words[i] << " ";
        }
//        std::cout << std::endl;
        unsigned int u = stoi(words[1]);
        unsigned int v = stoi(words[2]);
        if(words[0] == "0")
        {
//            std::cout<< "0command" << std::endl;;
            uf.union_ab(u,v);
        }
        if(words[0] == "1")
        {
            int dummy = 0;
            unsigned int u_p = uf.find(u, dummy);
            unsigned int v_p = uf.find(v, dummy);
            bool flg = (uf.find(u, dummy) == uf.find(v, dummy));
            uf.print_parents();
            uf.print_ranks();
            std::cout << "answer ";
            if(flg)std::cout << 1 << std::endl;
            else std::cout << 0 << std::endl;
        }
    }
    return 0;
}
