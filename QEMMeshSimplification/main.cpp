//レポートA 2024/11/16 須之内俊樹 スノウチトシキ 23N8100018B
#include "Mesh.hpp"
#include <algorithm>
#include <tuple>
#include <queue>
#include <set>
#include "UnionFind.hpp"
//時間計測用
struct TimeDisplayer{
    std::chrono::system_clock::time_point startTime;
    std::chrono::system_clock::time_point endTime;
    const char* str;
    void startTimer(const char* s);
    double endTimer();
};
void TimeDisplayer::startTimer(const char* s){
    startTime = std::chrono::system_clock::now();
    str = s;
}
double TimeDisplayer::endTimer(){
    endTime = std::chrono::system_clock::now();
    double time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count());
    std::cout << str << ": " << time << "ms" << std::endl;
    return time;
}
//3次元配列の添字(x,y,z)を1次元配列の添字に変換する．
unsigned int volume_to_array_id(unsigned int x,unsigned int y,unsigned int z,
                                 unsigned int ny,unsigned int nz)
{
    return nz*ny*x + nz*y + z;
}
//頂点座標から1次元配列の添字に変換する．
unsigned int position_to_array_id(Eigen::Vector3d &position,Eigen::Vector3d &origin,double lx,double ly, double lz)
{
    //structured gridの空間解像度はnx*nx*nx
    unsigned int nx = 100;
    //structured gridの一辺の長さ．AABBの最長辺の1/nx
    double dx = std::max(std::max(lx,ly),lz) / nx;
    Eigen::Vector3d delta = {1e-6,1e-6,1e-6};
    Eigen::Vector3d L = position - origin - delta;
    unsigned int id_x = (unsigned int)(L.x() / dx);
    unsigned int id_y = (unsigned int)(L.y() / dx);
    unsigned int id_z = (unsigned int)(L.z() / dx);
    return volume_to_array_id(id_x, id_y, id_z, nx, nx);
}
//QEMの構造体
struct QEM{
    //QEMの係数
    Eigen::Matrix3d A;
    Eigen::Vector3d b;
    double c;
    //QEMを最小にする位置
    Eigen::Vector3d minimizer;
    //コンストラクタ
    QEM()
    {
        A = Eigen::Matrix3d::Zero(3, 3);
        b = Eigen::Vector3d::Zero();
        c = 0.0;
    }
    //デバッグ用標準出力
    void print()
    {
        std::cout << "A=" << std::endl;
        std::cout << A << std::endl;
        std::cout << "b=";
        std::cout << b.transpose() << std::endl;
        std::cout << "c=";
        std::cout << c << std::endl;
    }
    //QEMを最小化するminimizerの位置を計算し，その位置でのQEMの値を返す．
    //引数にminimizerを近づけたい位置を渡す．
    double min_value(Eigen::Vector3d &edge_pos)
    {
        //minimizer
        Eigen::Vector3d x;
        //SVDのソルバー
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,Eigen::ComputeFullU | Eigen::ComputeFullV);
        //特異値
        Eigen::VectorXd singular_values = svd.singularValues();
        //0.0へ置き換えた後の特異値
        Eigen::VectorXd fixed_singular_values = svd.singularValues();
        //特異値の逆数
        Eigen::VectorXd fixed_singular_values_inverse = svd.singularValues();
        double max = singular_values(0);
        //特異性の計算と，0.0への置き換え
        for(int i=0;i<3;++i)
        {
            //特異値が非常に小さいとき
            if(singular_values(i) < 1e-8)
            {
                fixed_singular_values(i) = 0.0;
                fixed_singular_values_inverse(i) = 0.0;
                
            }
            //条件数による特異性の判定
            else if(max / singular_values(i) > 1e03)
            {
                fixed_singular_values(i) = 0.0;
                fixed_singular_values_inverse(i) = 0.0;
            }
            else fixed_singular_values_inverse(i) = 1 / fixed_singular_values_inverse(i);
        }
        //特異値行列
        Eigen::MatrixXd fixed_singular_value_matrix = fixed_singular_values.asDiagonal();
        //特異値の逆数の行列
        Eigen::MatrixXd fixed_singular_value_matrix_inverce = fixed_singular_values_inverse.asDiagonal();
        //特異値を置き換えた後の行列A
        Eigen::MatrixXd edited_A = svd.matrixU() * fixed_singular_value_matrix * svd.matrixV().transpose();
        //Moore_Penrose形一般逆行列
        Eigen::MatrixXd Moore_Penrose = svd.matrixV() * fixed_singular_value_matrix_inverce * svd.matrixU().transpose();
        //minimizer
        x = edge_pos + Moore_Penrose * (b - edited_A * edge_pos);
        minimizer = x;
        //QEMの最小値を計算して返す
        double s = x.transpose()*A*x;
        double t = 2*b.transpose()*x;
        double cost = s - t + c;
        return cost;
    }
    //QEMの係数を計算する．
    void calQEMCoefficient(Eigen::Vector3d &f0_pos,Eigen::Vector3d &f1_pos,Eigen::Vector3d &f2_pos)
    {
        Eigen::Vector3d face_normal = (f1_pos - f0_pos).cross(f2_pos - f0_pos);
        Eigen::Vector3d onface_point = f0_pos;
        double area = face_normal.norm();
        face_normal.normalize();
        A = area*face_normal*(face_normal.transpose());
        b = A * onface_point;
        c = b.transpose()*onface_point;
    }
    QEM operator +(const QEM qem)
    {
        QEM tmp;
        tmp.A = A + qem.A;
        tmp.b = b + qem.b;
        tmp.c = c + qem.c;
        return tmp;
    }
};

//QEMメッシュ簡略化
void QEMmesh_simplificatiton(Mesh &input_mesh,Mesh &output_mesh,unsigned int terget_num_vertex)
{
    //output_mesh = input_mesh;
    //簡略化中の頂点数
    int n = input_mesh.vertices.size();
    int dummy = 0;
    using edge = std::pair<int,int>;
    using CostEdgeTime = std::tuple<double,edge,int>;
    //簡略化後の頂点リスト
    std::vector<Eigen::Vector3d> coarse_v_list;
    //簡略化メッシュの面リスト
    std::vector<std::vector<unsigned int>> coarse_f_list;
    //簡略化前後の頂点番号のmap
    std::vector<int> coarse_v_idx(n,-1);
    //頂点のQEM
    std::vector<QEM> Vertices_QEM(input_mesh.vertices.size(),QEM());
    //UnionFindDS
    UnionFind uf(input_mesh.vertices.size());
    //timestamp
    int t = 0;
    //heap
    std::priority_queue<CostEdgeTime> pq;
    //辺の重複確認用のset．keyは辺の数が頂点数の2乗以下とする．
    std::set<int>st;
    //全頂点のQEMを計算する
    for(int i=0;i<input_mesh.faces.size();++i)
    {
        //面を作る頂点の，入力メッシュでの頂点番号
        unsigned int p0_id = input_mesh.faces[i][0];
        unsigned int p1_id = input_mesh.faces[i][1];
        unsigned int p2_id = input_mesh.faces[i][2];
        //面を作る頂点の，入力メッシュでの頂点座標
        Eigen::Vector3d p0_pos = input_mesh.vertices[p0_id];
        Eigen::Vector3d p1_pos = input_mesh.vertices[p1_id];
        Eigen::Vector3d p2_pos = input_mesh.vertices[p2_id];
        Eigen::Vector3d face_normal = input_mesh.face_normals[i];
        QEM qem;
        //面のQEMの係数を計算
        qem.calQEMCoefficient(p0_pos, p1_pos,p2_pos);
        Vertices_QEM[p0_id] = Vertices_QEM[p0_id] + qem;
        Vertices_QEM[p1_id] = Vertices_QEM[p1_id] + qem;
        Vertices_QEM[p2_id] = Vertices_QEM[p2_id] + qem;
    }
//    for(int i=0;i<Vertices_QEM.size();++i)Vertices_QEM[i].print();
    //ヒープを構成する．コストは辺の両端のQEMの和の最小値
    for(int i=0;i<input_mesh.faces.size();++i)
    {
        std::vector<edge> f_edges;
        for(int j=0;j<3;++j)
        {
            edge e = std::make_pair(input_mesh.faces[i][j], input_mesh.faces[i][(j+1)%3]);
            if(input_mesh.faces[i][j] > input_mesh.faces[i][(j+1)%3])std::swap(e.first,e.second);
            int key = e.first*input_mesh.vertices.size()*input_mesh.vertices.size() + e.second;
            if(st.find(key) == st.end())
            {
                st.insert(key);
                Eigen::Vector3d onface_pos = input_mesh.vertices[e.first];
                //std::priority_queueはデフォルトだと降順なので，QEMの昇順にするためにコストはQEMを-1倍したものにする．
                double cost = -(Vertices_QEM[e.first] + Vertices_QEM[e.second]).min_value(onface_pos);
                pq.push({cost,e,t});
            }
        }
    }
    //簡略化ループ
    //目標の頂点数以上の間
    while(n > terget_num_vertex)
    {
        auto[cost,e_min,timestamp] = pq.top();
        pq.pop();
        int r0 = uf.find(e_min.first, dummy);
        int r1 = uf.find(e_min.second, dummy);
        Eigen::Vector3d onface_vertex = input_mesh.vertices[r0];
        //辺は縮約ずみ
        if(r0 == r1)continue;
        //コストが古いので更新する
        double c = -(Vertices_QEM[r0] + Vertices_QEM[r1]).min_value(onface_vertex);
        if(timestamp != t)pq.push({c,e_min,t});
        else//辺を縮約する
        {
            uf.union_ab(r0, r1);
            Vertices_QEM[uf.find(r0,dummy)] = Vertices_QEM[r0] + Vertices_QEM[r1];
            ++t;
            --n;
        }
    }
    //簡略化メッシュの抽出
    //頂点の登録
    int new_vertex_idx = 0;
    for(int i=0;i<input_mesh.vertices.size();++i)
    {
        if(uf.find(i, dummy) == i)
        {
            coarse_v_idx[i] = new_vertex_idx;
            ++new_vertex_idx;
            Vertices_QEM[i].min_value(input_mesh.vertices[i]);
            coarse_v_list.push_back(Vertices_QEM[i].minimizer);
        }
    }
    //面リストの作成
    for(auto &f:input_mesh.faces)
    {
        unsigned int v0 = coarse_v_idx[uf.find(f[0], dummy)];
        unsigned int v1 = coarse_v_idx[uf.find(f[1], dummy)];
        unsigned int v2 = coarse_v_idx[uf.find(f[2], dummy)];
        if(v0 != v1 && v1 != v2 && v2 != v0)coarse_f_list.push_back({v0,v1,v2});
    }
//    出力メッシュの頂点リストと面リストを更新
    output_mesh.set_vertices(coarse_v_list);
    output_mesh.set_faces(coarse_f_list);
}

int main(int argc, const char * argv[]) {
//    std::string inputFileName = "bun_zipper.obj";
//    int terget_v_num = 24333;
    std::string inputFileName = "lucy.obj";
    //頂点クラスタリングによる簡略化メッシュの頂点数まで削減する
    int terget_v_num = 13824;
    std::string QEM_outputFileName = "QEM_simplificate_" + inputFileName;
    Mesh in_mesh,centroid_out_mesh,QEM_out_mesh;
    TimeDisplayer td;
    double time = 0;
    //読み込み
    td.startTimer("input");
    in_mesh.input(inputFileName);
    in_mesh.cal_face_normals(true);
    time += td.endTimer();
    //メッシュ簡略化
    td.startTimer("simplification");
    QEMmesh_simplificatiton(in_mesh, QEM_out_mesh,terget_v_num);
    time += td.endTimer();
    //簡略化前後の頂点数，面数出力
    std::cout << "before: vertices.size(), faces.size() = "
    << in_mesh.vertices.size() << "," << in_mesh.faces.size() << std::endl;
    std::cout << "centroid after: vertices.size(), faces.size() = "
    << centroid_out_mesh.vertices.size() << "," << centroid_out_mesh.faces.size() << std::endl;
    std::cout << "QEM after: vertices.size(), faces.size() = "
    << QEM_out_mesh.vertices.size() << "," << QEM_out_mesh.faces.size() << std::endl;
    //簡略化メッシュの出力
    td.startTimer("output");
    QEM_out_mesh.outputOBJ(QEM_outputFileName);
    time += td.endTimer();
    std::cout << "total: " << time << "ms"<< std::endl;
    return 0;
}
