#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/function_output_iterator.hpp>
#include <fstream>
#include <vector>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>

#include <igl/readOFF.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;
struct halfedge2edge
{
  halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
    : m_mesh(m), m_edges(edges)
  {}
  void operator()(const halfedge_descriptor& h) const
  {
    m_edges.push_back(edge(h, m_mesh));
  }
  const Mesh& m_mesh;
  std::vector<edge_descriptor>& m_edges;
};



template <typename P>
bool read_off(const Eigen::MatrixXd& V,
              const Eigen::MatrixXi& F,
              CGAL::Surface_mesh<P>& sm)
{
  using namespace CGAL;
  typedef Surface_mesh<P> Mesh;
  typedef typename Kernel_traits<P>::Kernel K;
  typedef typename K::Vector_3 Vector_3;
  typedef typename Mesh::Face_index Face_index;
  typedef typename Mesh::Vertex_index Vertex_index;
  typedef typename Mesh::size_type size_type;
  int n, f, e = 0;
  std::string off;

  n = V.rows();
  f = F.rows();

  sm.reserve(V.rows(), F.rows()*2, F.rows());
  std::vector<Vertex_index> vertexmap(n);
  P p;
  Vector_3 v;
  typename Mesh::template Property_map<Vertex_index,CGAL::Color> vcolor;
  typename Mesh::template Property_map<Vertex_index,Vector_3> vnormal;
  bool vcolored = false, v_has_normals = false;

  char ci;

  for(int i=0; i < n; i++){
    Vertex_index vi = sm.add_vertex(P(V(i,0), V(i,1), V(i,2)));
    vertexmap[i] = vi;
  }
  std::vector<Vertex_index> vr;
  size_type d, vi;
  bool fcolored = false;
  typename Mesh::template Property_map<Face_index,CGAL::Color> fcolor;

  for(int i=0; i < f; i++){
    d = 3;
    vr.resize(3);
    for(std::size_t j=0; j<d; j++){
      vi = F(i,j);
      vr[j] = vertexmap[vi];
    }
    Face_index fi = sm.add_face(vr);
    if(fi == sm.null_face())
    {
      sm.clear();
      return false;
    }
   
  }
  return true;
}

  template <typename P>
  bool write_off(const CGAL::Surface_mesh<P>& sm, Eigen::MatrixXd& V,
              Eigen::MatrixXi& F) {
    using namespace CGAL;
    typedef Surface_mesh<P> Mesh;
    typedef typename Mesh::Vertex_index Vertex_index;
    typedef typename Mesh::Face_index Face_index;

    V.resize(sm.number_of_vertices(), 3);
    F.resize(sm.number_of_faces(), 3);
    std::vector<int> reindex;
    reindex.resize(sm.num_vertices());
    int n = 0;
    BOOST_FOREACH(Vertex_index v, sm.vertices()){
      auto p = sm.point(v);
      V.row(n) << p.x(), p.y(), p.z();
      reindex[v]=n++;
    }

    int n_f = 0;
    BOOST_FOREACH(Face_index f, sm.faces()){
      // os << sm.degree(f);
      int fj = 0;
      BOOST_FOREACH(Vertex_index v, CGAL::vertices_around_face(sm.halfedge(f),sm)){
        F(n_f, fj) = reindex[v];
        fj++;
      }
      n_f ++;
    }
    return true;
  }


int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "build/bumpy.off";
  std::ifstream input(filename);

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::readOFF(filename,V,F);
  Mesh mesh;
  if (!input || !::read_off(V,F, mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }
  // for(int i=0; i<V.rows(); i++)
  double target_edge_length = 0.4;
  unsigned int nb_iter = 3;
  std::cout << "#Split border...";
    std::vector<edge_descriptor> border;
    PMP::border_halfedges(faces(mesh),
      mesh,
      boost::make_function_output_iterator(halfedge2edge(mesh, border)));
    PMP::split_long_edges(border, target_edge_length, mesh);
  std::cout << "#done." << std::endl;
  std::cout << "#Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;
  PMP::isotropic_remeshing(
      faces(mesh),
      target_edge_length,
      mesh,
      PMP::parameters::number_of_iterations(nb_iter)
      .protect_constraints(true)//i.e. protect border, here
      );
  V.setZero();
  F.setZero();
  ::write_off(mesh, V, F);
  igl::writeOBJ("bumpy.obj",V,F);
  std::cout << "#Remeshing done." << std::endl;
  return 0;
}
