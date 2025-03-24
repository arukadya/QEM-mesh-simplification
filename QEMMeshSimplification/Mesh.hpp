#ifndef Mesh_hpp
#define Mesh_hpp

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include "Eigen/Core"
#include "Eigen/Dense"
struct Mesh
{
    //頂点リスト
    std::vector<Eigen::Vector3d>vertices;
    //頂点法線リスト
    std::vector<Eigen::Vector3d>vertex_normals;
    //面リスト
    std::vector<std::vector<unsigned int>>faces;
    //面法線リスト
    std::vector<Eigen::Vector3d>face_normals;
    //ハーフエッジデータ構造
    
    //入出力用関数
    void input(std::string InputFlieName);
    int inputOBJ(std::string InputFlieName);
    void inputOFF(const char* InputFileName);
    void outputOBJ(std::string OutputFileName);
    //AABBのx,y,z座標最小の点と，x,y,z軸方向の長さを計算
    void getAABB(Eigen::Vector3d &origin, double &lx, double &ly, double &lz);
    Mesh(){}
    //setter
    void set_vertices(std::vector<Eigen::Vector3d>&vertex_list);
    void set_faces(std::vector<std::vector<unsigned int>>&face_list);
    void cal_face_normals(bool is_normalize);
};
#endif /* Mesh_hpp */
