#pragma once

#include <string>
#include <array>
#include <memory>
#include <vector>
#include <sstream>
#include <algorithm>
#include <limits>
#include <unordered_map>

#include "Mollifier.hpp"

class FemProjection;
class Hex27Projection;
class Tet15Projection;
class Wedge21Projection;
class Quad9Projection;
class Tri7Projection;
class Line3Projection;
class Mesh;
class Field;

struct InitialMeshData {
    std::vector<unsigned> elLevel;
    std::vector<unsigned> elType;
    std::vector<std::vector<unsigned>> elTplgy;
    std::vector<std::vector<double>> X;
};

// struct CellMarkersData {
//   std::vector<std::vector<double>> markers;  // globale: [dim][marker]
//   std::unordered_map<unsigned, std::vector<unsigned>> own;
//   std::unordered_map<unsigned, std::vector<unsigned>> combined;
// };

struct CellMarkersData;

namespace {
  constexpr unsigned INVALID_ID = std::numeric_limits<unsigned>::max();

  inline long long qcoord(const double x, const double eps = 1.0e-12) {
    return static_cast<long long>(std::llround(x / eps));
  }

  inline std::string makeNodeKey(const std::vector<std::vector<double>>& X,
                                 const unsigned dim,
                                 const unsigned inode,
                                 const double eps = 1.0e-12) {
    std::ostringstream os;
    for (unsigned d = 0; d < dim; ++d) {
      os << qcoord(X[d][inode], eps) << '#';
    }
    return os.str();
  }

  inline std::string makeElemKey(std::vector<unsigned> conn) {
    std::sort(conn.begin(), conn.end());
    std::ostringstream os;
    for (const auto v : conn) {
      os << v << '#';
    }
    return os.str();
  }
}

class Ex35 {
private:
    int _levelS;
    int _levelN;
    int _dim;
    unsigned _neighMode;
    std::string _filename;
    std::string _filenamefixed;
    std::vector<std::vector<double>> _markers;
    std::vector<std::vector<double>> _nodeCoords;
    std::vector<std::vector<double>> _markerVelocity;
    std::vector<std::vector<double>> _nodeVelocity;
    std::unordered_map<unsigned, std::vector<unsigned>> _elemMarkers;
    Mollifier _mollifier;
    unsigned _matchedLevel = 0;

    std::vector<unsigned> _origToThisNode;
    std::vector<unsigned> _thisToOrigNode;

    std::vector<unsigned> _origToThisElem;
    std::vector<unsigned> _thisToOrigElem;

    std::array<std::unique_ptr<FemProjection>, 6> _elProj;

    std::vector<std::vector<std::vector<double>>> _Xvel;
    std::vector<std::vector<unsigned>> _elTypevel;
    std::vector<std::vector<std::vector<unsigned>>> _elTplgyvel;
    std::vector<std::vector<unsigned>> _elLevelvel;
    std::vector<std::vector<int>> _AMRvel;
    std::vector<std::vector<unsigned>> _elFathervel;
    std::vector<std::vector<std::vector<unsigned>>> _elChildernvel;
    std::vector<std::vector<std::vector<unsigned>>> _nodeElementsvel;
    std::vector<std::vector<std::vector<unsigned>>> _neighborsvel;

    std::vector<std::vector<std::vector<double>>> _X0;
    std::vector<std::vector<unsigned>> _elType0;
    std::vector<std::vector<std::vector<unsigned>>> _elTplgy0;
    std::vector<std::vector<unsigned>> _elLevel0;
    std::vector<std::vector<int>> _AMR0;
    std::vector<std::vector<unsigned>> _elFather0;
    std::vector<std::vector<std::vector<unsigned>>> _elChildern0;
    std::vector<std::vector<std::vector<unsigned>>> _nodeElements0;
    std::vector<std::vector<std::vector<unsigned>>> _neighbors0;

    std::vector<std::vector<std::vector<double>>> _X1;
    std::vector<std::vector<unsigned>> _elType1;
    std::vector<std::vector<std::vector<unsigned>>> _elTplgy1;
    std::vector<std::vector<unsigned>> _elLevel1;
    std::vector<std::vector<int>> _AMR1;
    std::vector<std::vector<unsigned>> _elFather1;
    std::vector<std::vector<std::vector<unsigned>>> _elChildern1;
    std::vector<std::vector<std::vector<unsigned>>> _nodeElements1;
    std::vector<std::vector<std::vector<unsigned>>> _neighbors1;

    std::vector<std::vector<std::vector<double>>> _X2;
    std::vector<std::vector<unsigned>> _elType2;
    std::vector<std::vector<std::vector<unsigned>>> _elTplgy2;
    std::vector<std::vector<unsigned>> _elLevel2;
    std::vector<std::vector<int>> _AMR2;
    std::vector<std::vector<unsigned>> _elFather2;
    std::vector<std::vector<std::vector<unsigned>>> _elChildern2;
    std::vector<std::vector<std::vector<unsigned>>> _nodeElements2;
    std::vector<std::vector<std::vector<unsigned>>> _neighbors2;

    std::vector<Mesh> _meshvel;
    std::vector<Field> _fieldvel;
    std::vector<Mesh> _mesh0;
    std::vector<Field> _field0;
    std::vector<Mesh> _mesh1;
    std::vector<Field> _field1;
    std::vector<Mesh> _mesh2;
    std::vector<Field> _field2;

public:
    explicit Ex35(int levelS, int levelN);
    void initializeVELMesh(const InitialMeshData& mesh_data);
    void initializeLSMesh();
    InitialMeshData buildSeed1x2() const;
    InitialMeshData buildSeed5x10() const;
    void buildMapsToOriginalMesh(const InitialMeshData& mesh_data, const unsigned level);
    void initializeField_Ball(std::vector<double> xc, double r);
    void initMarkers();
    void advectField(const std::vector<std::vector<double>>& vel, const int it, const double dt);
    CellMarkersData getCellMarkers();
    std::vector<double> interpolatePsiOnOriginalNodes() const;
    void interpolatePsiOnNodes(const std::vector<std::vector<double>>& node_coords, std::vector<double>& psi) const;
    void reinit();
    ~Ex35();
};
