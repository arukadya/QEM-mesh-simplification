#include "Mesh.hpp"
void Mesh::cal_face_normals(bool is_normalize)
{
    for(int i=0;i<faces.size();++i)
    {
        Eigen::Vector3d vn0 = vertex_normals[faces[i][0]];
        Eigen::Vector3d vn1 = vertex_normals[faces[i][1]];
        Eigen::Vector3d vn2 = vertex_normals[faces[i][2]];
        Eigen::Vector3d face_normal = (vn0 + vn1 + vn2)/3;
        if(is_normalize)face_normal.normalize();
//        std::cout << face_normal.transpose() << std::endl;
        face_normals.push_back(face_normal);
    }
}
void Mesh::input(std::string InputFileName){
    std::cout << "InputFileName:" << InputFileName << std::endl;
    if(std::filesystem::path(InputFileName).extension() == ".obj")inputOBJ(InputFileName);
    if(std::filesystem::path(InputFileName).extension() == ".off")inputOFF(InputFileName.c_str());
}
int Mesh::inputOBJ(std::string InputFileName){
    std::ifstream Inputfile(InputFileName);
    if (!Inputfile.is_open()) {
        std::cerr << "Could not open the file - '"
             << InputFileName << "'" << std::endl;
        return EXIT_FAILURE;
    }
    std::string line;
    std::string word;
    
    while(std::getline(Inputfile,line)){
        std::stringstream ss_nums{line};
        getline(ss_nums,word,' ');
        if(word[0] == '#' || word[0] == '\n')continue;
        if(word[0] == 'v' && word.size() == 1){
            //std::cout << word << ",";
            Eigen::Vector3d v0;
            do{
                getline(ss_nums,word,' ');
            }while(word.length() <= 0);
            //std::cout << word << ",";
            v0.x() = std::stof(word);
            
            do{
                getline(ss_nums,word,' ');
            }while(word.length() <= 0);
            v0.y() = std::stof(word);
            //std::cout << word << ",";
            do{
                getline(ss_nums,word,' ');
            }while(word.length() <= 0);
            v0.z() = std::stof(word);
            //std::cout << word << ",";
            vertices.push_back(v0);
            //std::cout << std::endl;
        }
        if(word[0] == 'v' && word[1] == 'n' && word.size() == 2){
            //std::cout << word << ",";
            Eigen::Vector3d n;
            do{
                getline(ss_nums,word,' ');
            }while(word.length() <= 0);
            //std::cout << word << ",";
            n.x() = std::stof(word);
            
            do{
                getline(ss_nums,word,' ');
            }while(word.length() <= 0);
            n.y() = std::stof(word);
            //std::cout << word << ",";
            do{
                getline(ss_nums,word,' ');
            }while(word.length() <= 0);
            n.z() = std::stof(word);
            //std::cout << word << ",";
            vertex_normals.push_back(n);
            //std::cout << std::endl;
        }
        if(word[0] == 'f'){
            std::vector<unsigned int> f;
            while(getline(ss_nums,word,' ')){
                f.push_back(std::stoi(word)-1);
            };
            faces.push_back(f);
        }
    };
    Inputfile.close();
    return EXIT_SUCCESS;
}
void Mesh::inputOFF(const char* InputFileName){
    FILE *ifp = fopen(InputFileName,"r");
    int num_vertices, num_faces,dummy;
    fscanf(ifp, "OFF %d %d %d", &num_vertices, &num_faces, &dummy);
    for(int i=0;i<num_vertices;i++){
        double x,y,z;//点の入力
        fscanf(ifp, "%lf %lf %lf", &x, &y, &z);
        Eigen::Vector3d v(x,y,z);
        vertices.push_back(v);
    }
    for(int i=0;i<num_faces;i++){//面の入力
        unsigned int num_size, v0, v1, v2;
        fscanf(ifp, "%d %d %d %d", &num_size, &v0, &v1, &v2);
        faces.push_back({v0,v1,v2});
    }
    fclose(ifp);
}
void Mesh::outputOBJ(std::string OutputFileName)
{
    std::ofstream outputfile(OutputFileName);
    for(int i=0;i<vertices.size();++i)
    {
        if(!vertex_normals.empty())
        {
            Eigen::Vector3d vn = vertex_normals[i];
            outputfile << "vn " << vn.x() << " " << vn.y() << " " << vn.z() << std::endl;
        }
        Eigen::Vector3d v = vertices[i];
        outputfile << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
    }
    for(int i=0;i<faces.size();++i)
    {
        if(!vertex_normals.empty())
        {
            std::vector<unsigned int>f = faces[i];
            outputfile << "f "
            << f[0] + 1 << "//" << f[0] + 1 << " "
            << f[1] + 1 << "//" << f[1] + 1 << " "
            << f[2] + 1 << "//" << f[2] + 1 << " "
            << std::endl;
        }
        else{
            std::vector<unsigned int>f = faces[i];
            outputfile << "f " << f[0] + 1 << " " << f[1] + 1 << " " << f[2] + 1 << " " << std::endl;
        }
    }
    outputfile.close();
}
void Mesh::getAABB(Eigen::Vector3d &origin, double &lx, double &ly, double &lz)
{
    double min_x = vertices[0].x(),min_y = vertices[0].y(),min_z = vertices[0].z();
    double max_x = vertices[0].x(),max_y = vertices[0].y(),max_z = vertices[0].z();
    for(auto &v : vertices)
    {
        if(v.x() < min_x)min_x = v.x();
        if(v.x() > max_x)max_x = v.x();
        if(v.y() < min_y)min_y = v.y();
        if(v.y() > max_y)max_y = v.y();
        if(v.z() < min_z)min_z = v.z();
        if(v.z() > max_z)max_z = v.z();
    }
    origin.x() = min_x;origin.y() = min_y;origin.z() = min_z;
    lx = max_x - min_x;ly = max_y - min_y;lz = max_z-min_z;
}
void Mesh::set_vertices(std::vector<Eigen::Vector3d>&vertex_list)
{
    vertices = vertex_list;
}
void Mesh::set_faces(std::vector<std::vector<unsigned int>>&face_list)
{
    faces = face_list;
}

