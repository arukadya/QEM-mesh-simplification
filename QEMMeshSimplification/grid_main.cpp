//演習04A 2024/10/21 須之内俊樹 スノウチトシキ 23N8100018B
#include "Mesh.hpp"
#include <algorithm>
struct QEM{
    Eigen::Matrix3d A;
    Eigen::Vector3d b;
    Eigen::Vector3d minimizer;
    double c;
    QEM()
    {
        A = Eigen::Matrix3d::Zero(3, 3);
        b = Eigen::Vector3d::Zero();
        c = 0.0;
    }
    void print()
    {
        std::cout << "A=" << std::endl;
        std::cout << A << std::endl;
        std::cout << "b=";
        std::cout << b.transpose() << std::endl;
        std::cout << "c=";
        std::cout << c << std::endl;
    }
    double min_value(Eigen::Vector3d &edge_pos)
    {
        //minimizer
        Eigen::Vector3d x;
        //SVDを利用
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A,Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::VectorXd singular_values = svd.singularValues();
        Eigen::VectorXd e_singular_values = svd.singularValues();
        Eigen::VectorXd e_singular_values_inverse = svd.singularValues();
        double max = singular_values(0);
        for(int i=0;i<3;++i)
        {
            if(singular_values(i) < 1e-8)
            {
                e_singular_values(i) = 0.0;
                e_singular_values_inverse(i) = 0.0;
                
            }
            else if(max / singular_values(i) > 1e03)
            {
                e_singular_values(i) = 0.0;
                e_singular_values_inverse(i) = 0.0;
//                std::cout << "rankochi" << std::endl;
            }
            else e_singular_values_inverse(i) = 1 / e_singular_values_inverse(i);
        }
        Eigen::MatrixXd e_singular_value_matrix = e_singular_values.asDiagonal();
        Eigen::MatrixXd e_singular_value_matrix_inverce = e_singular_values_inverse.asDiagonal();
        Eigen::MatrixXd edited_A = svd.matrixU() * e_singular_value_matrix * svd.matrixV().transpose();
        Eigen::MatrixXd Moore_Penrose = svd.matrixV() * e_singular_value_matrix_inverce * svd.matrixU().transpose();
        x = edge_pos + Moore_Penrose * (b - edited_A * edge_pos);
        minimizer = x;
        double s = x.transpose()*A*x;
        double t = 2*b.transpose()*x;
        double cost = s - t + c;
        return cost;
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
//grid内の頂点の重心を計算する．
bool get_grid_centroid(std::vector<Eigen::Vector3d>&grid_vertices,Eigen::Vector3d &centroid)
{
    if(!grid_vertices.empty())
    {
        centroid = Eigen::Vector3d::Zero();
        for(auto &v:grid_vertices)
        {
            centroid += v;
        }
        centroid /= grid_vertices.size();
        return true;
    }
    else return false;
}
//重複面の削除
std::vector<std::vector<unsigned int>> remove_deplicate_face
 (std::vector<std::vector<unsigned int>> &faces,
  std::vector<Eigen::Vector3d> &vertex_pos,
  std::vector<Eigen::Vector3d> &vertex_normals
  )
{
    std::vector<std::vector<unsigned int>> deplicated_faces;
    std::map<unsigned int, std::vector<unsigned int>>map;
    int size = (int)faces.size();
    for(int i=0;i<size;++i)
    {
        std::vector<unsigned int>face_index;
        for(int j=0;j<3;++j)face_index.push_back(faces[i][j]);
        std::sort(face_index.begin(),face_index.end());
        unsigned int key = volume_to_array_id(face_index[0], face_index[1], face_index[2], size, size);
        if(map.count(key) == 0)map.insert(std::make_pair(key, faces[i]));
    }
    for(auto &x:map)
    {
        std::vector<unsigned int>virtual_face = x.second;
        deplicated_faces.push_back(virtual_face);
    }
    return deplicated_faces;
}

//重心座標によるメッシュ簡略化
void mesh_simplificatiton(Mesh &input_mesh,Mesh &output_mesh)
{
    //入力メッシュのAABBの，x,y,z座標最小の点
    Eigen::Vector3d origin;
    //入力メッシュのAABBの，x,y,z軸方向の長さ
    double lx,ly,lz;
    input_mesh.getAABB(origin, lx, ly, lz);
    //structured gridの空間解像度はnx*nx*nx
    unsigned int nx = 100;
    //簡略化前後の頂点idのmap
    std::vector<int> simplificate_vertices_map(input_mesh.vertices.size(),-1);
    //structure_gridに属する頂点のid集合
    std::vector<std::vector<unsigned int>> classtaring_vertices_id(nx*nx*nx);
    //structure_gridに対するQEMminimizerの行列A
    std::vector<Eigen::Matrix3d> grid_sum_QEM_A(nx*nx*nx);
    //structure_gridに対するQEMminimizerのベクトルb
    std::vector<Eigen::Vector3d> grid_sum_QEM_b(nx*nx*nx);
    //簡略化後の頂点リスト
    std::vector<Eigen::Vector3d> simplificate_vertices;
    //簡略化メッシュの面リスト
    std::vector<std::vector<unsigned int>>simplificate_faces;
    //面リストで使う頂点のみマークする
    for(auto &f:input_mesh.faces)
    {
        for(int i=0;i<3;++i)simplificate_vertices_map[f[i]] = 0;
    }
    //入力メッシュの頂点idをstructure_gridに振り分ける
    for(int i=0; i<input_mesh.vertices.size(); ++i)
    {
        if(simplificate_vertices_map[i] == -1)continue;
        classtaring_vertices_id[position_to_array_id(input_mesh.vertices[i],origin,lx,ly,lz)].push_back(i);
    }
    //structure_gridごとに重心を計算し，頂点ごとに簡略化後の頂点番号を登録
    int simplificate_vertex_id = 0;
    for(auto &id_contenna : classtaring_vertices_id)
    {
        Eigen::Vector3d centroid;
        //各strcture_gridに属する頂点座標リスト
        std::vector<Eigen::Vector3d>pos_contenna;
        for(auto &id : id_contenna)pos_contenna.push_back(input_mesh.vertices[id]);
        //重心を計算
        if(get_grid_centroid(pos_contenna, centroid))
        {
            simplificate_vertices.push_back(centroid);
            //簡略化前後の頂点idをmapに登録
            for(auto &id : id_contenna)simplificate_vertices_map[id] = simplificate_vertex_id;
            ++simplificate_vertex_id;
        }
    }
    //面リストの作成
    for(auto &f:input_mesh.faces)
    {
        //面を構成する頂点の簡略化後のid
        unsigned int simplificate_f0 = simplificate_vertices_map[f[0]];
        unsigned int simplificate_f1 = simplificate_vertices_map[f[1]];
        unsigned int simplificate_f2 = simplificate_vertices_map[f[2]];
        //簡略化後のidが全て異なるとき，簡略化後の面リストに面を追加
        if(simplificate_f0 != simplificate_f1
           && simplificate_f1 != simplificate_f2
           && simplificate_f2 != simplificate_f0)
        {
            std::vector<unsigned int>simplificate_face = 
            {
                simplificate_f0,
                simplificate_f1,
                simplificate_f2
            };
            simplificate_faces.push_back(simplificate_face);
        }
    }
    simplificate_faces = remove_deplicate_face(simplificate_faces, input_mesh.vertices, input_mesh.vertex_normals);
    //出力メッシュの頂点リストと面リストを更新
    output_mesh.set_vertices(simplificate_vertices);
    output_mesh.set_faces(simplificate_faces);
}

//QEMの係数を計算する．
void calQEMCoefficient(Eigen::Vector3d &f0_pos,Eigen::Vector3d &f1_pos,Eigen::Vector3d &f2_pos,QEM &qem)
{
    Eigen::Vector3d face_normal = (f1_pos - f0_pos).cross(f2_pos - f0_pos);
    Eigen::Vector3d onface_point = f0_pos;
    double area = face_normal.norm();
    face_normal.normalize();
    qem.A = area*face_normal*(face_normal.transpose());
    qem.b = qem.A * onface_point;
    qem.c = qem.b.transpose()*onface_point;
//    qem.print();
}
//QEMメッシュ簡略化
void QEMmesh_simplificatiton(Mesh &input_mesh,Mesh &output_mesh)
{
    //入力メッシュのAABBの，x,y,z座標最小の点
    Eigen::Vector3d origin;
    //入力メッシュのAABBの，x,y,z軸方向の長さ
    double lx,ly,lz;
    input_mesh.getAABB(origin, lx, ly, lz);
    //structured gridの空間解像度はnx*nx*nx
    unsigned int nx = 100;
    //頂点のQEM
    std::vector<QEM> grid_QEM(nx*nx*nx,QEM());
    //簡略化前後の頂点idのmap
    std::vector<unsigned int> simplificate_vertices_map(input_mesh.vertices.size());
    //structure_gridに属する頂点のid集合
    std::vector<std::vector<unsigned int>> classtaring_vertices_id(nx*nx*nx);
    //簡略化後の頂点リスト
    std::vector<Eigen::Vector3d> simplificate_vertices;
    //簡略化メッシュの面リスト
    std::vector<std::vector<unsigned int>>simplificate_faces;
    //面リストで使う頂点のみマークする
    
    input_mesh.faces = remove_deplicate_face(
    input_mesh.faces, input_mesh.vertices, input_mesh.vertex_normals);
    //structured gridごとのQEMMinimizerの線型方程式の係数を計算する
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
        //面を作る頂点の，grid番号
        unsigned int p0_grid_id = position_to_array_id(p0_pos, origin, lx, ly, lz);
        unsigned int p1_grid_id = position_to_array_id(p1_pos, origin, lx, ly, lz);
        unsigned int p2_grid_id = position_to_array_id(p2_pos, origin, lx, ly, lz);
        
        Eigen::Vector3d face_normal = input_mesh.face_normals[i];
        QEM qem;
        //面のQEMの係数を計算
        calQEMCoefficient(p0_pos, p1_pos,p2_pos, qem);
        //頂点が属するgridにQEMの係数を加算
        grid_QEM[p0_grid_id] = grid_QEM[p0_grid_id] + qem;
        grid_QEM[p1_grid_id] = grid_QEM[p1_grid_id] + qem;
        grid_QEM[p1_grid_id] = grid_QEM[p1_grid_id] + qem;
        //gridに所属する頂点番号を記録
        classtaring_vertices_id[p0_grid_id].push_back(p0_id);
        classtaring_vertices_id[p1_grid_id].push_back(p1_id);
        classtaring_vertices_id[p2_grid_id].push_back(p2_id);
    }
    //structure_gridごとにQEMMinimizerを計算し，頂点ごとに簡略化後の頂点番号を登録
    int simplificate_vertex_id = 0;
    for(int i=0;i<nx*nx*nx;++i)
    {
        std::vector<unsigned int>id_contenna = classtaring_vertices_id[i];
        if(id_contenna.empty())continue;
        Eigen::Vector3d grid_center = input_mesh.vertices[id_contenna[0]];
        grid_QEM[i].min_value(grid_center);
        Eigen::Vector3d simplificate_vertex = grid_QEM[i].minimizer;
        simplificate_vertices.push_back(simplificate_vertex);
        for(auto &id : id_contenna)simplificate_vertices_map[id] = simplificate_vertex_id;
        ++simplificate_vertex_id;
    }
    //面リストの作成
    for(auto &f:input_mesh.faces)
    {
        //面を構成する頂点の簡略化後のid
        unsigned int simplificate_f0 = simplificate_vertices_map[f[0]];
        unsigned int simplificate_f1 = simplificate_vertices_map[f[1]];
        unsigned int simplificate_f2 = simplificate_vertices_map[f[2]];
//        簡略化後のidが全て異なるとき，簡略化後の面リストに面を追加
        if(simplificate_f0 != simplificate_f1
           && simplificate_f1 != simplificate_f2
           && simplificate_f2 != simplificate_f0)
        {
            std::vector<unsigned int>simplificate_face =
            {
                simplificate_f0,
                simplificate_f1,
                simplificate_f2
            };
            simplificate_faces.push_back(simplificate_face);
        }
    }
    simplificate_faces = remove_deplicate_face(simplificate_faces, input_mesh.vertices, input_mesh.vertex_normals);
//    出力メッシュの頂点リストと面リストを更新
    output_mesh.set_vertices(simplificate_vertices);
    output_mesh.set_faces(simplificate_faces);
}
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

int main(int argc, const char * argv[]) {
    std::string source_dir = "resources/";
    std::string output_dir = "output/";
//    std::string inputFileName = "bun_zipper.obj";
    std::string inputFileName = "bun_zipper.obj";
//    int terget_v_num = 23333;
    std::string centroid_outputFileName = output_dir + "grid_centroid_simplificate_" + inputFileName;
    std::string QEM_outputFileName = output_dir + "grid_QEM_simplificate_" + inputFileName;
    Mesh in_mesh,centroid_out_mesh,QEM_out_mesh;
    TimeDisplayer td;
    double time = 0;
    //読み込み
    td.startTimer("input");
    in_mesh.input(source_dir + inputFileName);
    in_mesh.cal_face_normals(true);
//    in_mesh.cal_face_normals(false);
    time += td.endTimer();
    //メッシュ簡略化
    td.startTimer("simplification");
    QEMmesh_simplificatiton(in_mesh, QEM_out_mesh);
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
    centroid_out_mesh.outputOBJ(centroid_outputFileName);
    QEM_out_mesh.outputOBJ(QEM_outputFileName);
    time += td.endTimer();
    std::cout << "total: " << time << "ms"<< std::endl;
}

