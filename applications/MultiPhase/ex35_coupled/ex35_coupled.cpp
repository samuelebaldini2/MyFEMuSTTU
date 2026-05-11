#include "ex35_coupled.hpp"

#include <vector>
#include <cmath>
#include <memory>
#include <numeric>
#include <ctime>
#include <tuple>
#include <array>
#include <stdexcept>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <sstream>

#include <gperftools/profiler.h>

#include "Mesh.hpp"
#include "Field.hpp"
#include "FemProjection.hpp"   // FemProjection, Quad9Projection, Tri7Projection
#include "RefineMesh.hpp"
#include "VtkOutput.hpp"
#include "VtuOutput.hpp"
#include "Mollifier.hpp"
#include "Psi.hpp"
#include "RungeKutta.hpp"
#include "PointLocator.hpp"
#include "Errors.hpp"
#include "MeshSeed.hpp"
#include "AdvectMarkers.hpp"

#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"

#include "SimConfig.hpp"
#include "Reinit.hpp"

struct CellMarkersData {
  std::vector<std::vector<double>> markers;  // globale: [dim][marker]
  std::unordered_map<unsigned, std::vector<unsigned>> own;
  std::unordered_map<unsigned, std::vector<unsigned>> combined;
};

Ex35::Ex35(int levelS, int levelN)
  : _levelS(levelS),
    _levelN(levelS + levelN),
    _Xvel(levelS + levelN + 1),
    _elTypevel(levelS + levelN + 1),
    _elTplgyvel(levelS + levelN + 1),
    _elLevelvel(levelS + levelN + 1),
    _AMRvel(levelS + levelN + 1),
    _elFathervel(levelS + levelN + 1),
    _elChildernvel(levelS + levelN + 1),
    _nodeElementsvel(levelS + levelN + 1),
    _neighborsvel(levelS + levelN + 1),
    _X0(levelS + levelN + 1),
    _elType0(levelS + levelN + 1),
    _elTplgy0(levelS + levelN + 1),
    _elLevel0(levelS + levelN + 1),
    _AMR0(levelS + levelN + 1),
    _elFather0(levelS + levelN + 1),
    _elChildern0(levelS + levelN + 1),
    _nodeElements0(levelS + levelN + 1),
    _neighbors0(levelS + levelN + 1),
    _X1(levelS + levelN + 1),
    _elType1(levelS + levelN + 1),
    _elTplgy1(levelS + levelN + 1),
    _elLevel1(levelS + levelN + 1),
    _AMR1(levelS + levelN + 1),
    _elFather1(levelS + levelN + 1),
    _elChildern1(levelS + levelN + 1),
    _nodeElements1(levelS + levelN + 1),
    _neighbors1(levelS + levelN + 1),
    _X2(levelS + levelN + 1),
    _elType2(levelS + levelN + 1),
    _elTplgy2(levelS + levelN + 1),
    _elLevel2(levelS + levelN + 1),
    _AMR2(levelS + levelN + 1),
    _elFather2(levelS + levelN + 1),
    _elChildern2(levelS + levelN + 1),
    _nodeElements2(levelS + levelN + 1),
    _neighbors2(levelS + levelN + 1)
{
  _elProj[0] = std::make_unique<Hex27Projection>();
  _elProj[1] = std::make_unique<Tet15Projection>();
  _elProj[2] = std::make_unique<Wedge21Projection>();
  _elProj[3] = std::make_unique<Quad9Projection>();
  _elProj[4] = std::make_unique<Tri7Projection>();
  _elProj[5] = std::make_unique<Line3Projection>();

  _meshvel.reserve(levelS );
  _fieldvel.reserve(levelS );
  _mesh0.reserve(levelS + levelN + 1);
  _field0.reserve(levelS + levelN + 1);
  _mesh1.reserve(levelS + levelN + 1);
  _field1.reserve(levelS + levelN + 1);
  _mesh2.reserve(levelS + levelN + 1);
  _field2.reserve(levelS + levelN + 1);

  for (unsigned l = 0; l <= levelS + levelN; l++) {
    _meshvel.emplace_back(_elTplgyvel[l], _elTypevel[l], _elLevelvel[l], _Xvel[l], _AMRvel[l], _elFathervel[l], _elChildernvel[l], _nodeElementsvel[l], _neighborsvel[l]);
    _fieldvel.emplace_back(_meshvel[l]);

    _mesh0.emplace_back(_elTplgy0[l], _elType0[l], _elLevel0[l], _X0[l], _AMR0[l], _elFather0[l], _elChildern0[l], _nodeElements0[l], _neighbors0[l]);
    _field0.emplace_back(_mesh0[l]);

    _mesh1.emplace_back(_elTplgy1[l], _elType1[l], _elLevel1[l], _X1[l], _AMR1[l], _elFather1[l], _elChildern1[l], _nodeElements1[l], _neighbors1[l]);
    _field1.emplace_back(_mesh1[l]);

    _mesh2.emplace_back(_elTplgy2[l], _elType2[l], _elLevel2[l], _X2[l], _AMR2[l], _elFather2[l], _elChildern2[l], _nodeElements2[l], _neighbors2[l]);
    _field2.emplace_back(_mesh2[l]);
  }

}

void Ex35::printHeader() {
  const char* C_RST = "\033[0m";
  const char* C_BLD = "\033[1m";
  const char* C_CYN = "\033[36m";
  const char* C_YEL = "\033[33m";
  const char* C_GRN = "\033[32m";

  std::cout << "\n"
              << C_BLD << C_CYN
              << " --------------------------------------------------- \n"
              << "         AMR Level-set\n"
              << " --------------------------------------------------- "
              << C_RST << std::endl;

  return;
}

void Ex35::initializeVELMesh(const InitialMeshData& mesh_data) {

  auto t_start = std::chrono::high_resolution_clock::now();

  InitialMeshData seed;
  if (mesh_data.X.size() == 2 && mesh_data.elType[0] == 3)
    seed = buildSeed1x2();
  else if (mesh_data.X.size() == 3 && mesh_data.elType[0] == 0)
    seed = buildSeed1x1x2();

  _elLevelvel[0] = seed.elLevel;
  _elTypevel[0]  = seed.elType;
  _elTplgyvel[0] = seed.elTplgy;
  _Xvel[0]       = seed.X;
  _dim         = _Xvel[0].size();

  std::vector<unsigned> candidateIndices(_meshvel[0].numNodes());
  std::iota(candidateIndices.begin(), candidateIndices.end(), 0u);
  dedupNodesHash(_meshvel[0], candidateIndices);

  _meshvel[0].resetAllFathersToNoFather();
  _meshvel[0].buildNodeToElementAdjacency();
  _meshvel[0].buildFaceNeighborsFromNodeToElement();

  const unsigned targetNodes = mesh_data.X[0].size();

  std::cout << "target nodes = " << targetNodes << std::endl;

  bool found = false;

  if (_meshvel[0].numNodes() == targetNodes) {
    _matchedLevel = 0;
    found = true;
  }

  for (unsigned l = 1; l <= _levelS && !found; ++l) {
    _meshvel[l - 1].setUniformRefinement();
    refineAndProjectMesh(_elProj, _meshvel[l - 1], _meshvel[l]);
    _meshvel[l].buildNodeToElementAdjacency();
    _meshvel[l].buildFaceNeighborsFromNodeToElement();

    const unsigned nNodes = _meshvel[l].numNodes();

    if (nNodes == targetNodes) {
      std::cout << "level " << l << " nodes = " << nNodes << std::endl;
      _matchedLevel = l;
      found = true;
      break;
    }

    if (nNodes > targetNodes) {
      throw std::runtime_error("Ex35::initializeMesh: target mesh not compatible with uniform refinement from 2x1 seed");
    }
  }

  if (!found) {
    throw std::runtime_error("Ex35::initializeMesh: could not match target mesh with available refinement levels");
  }

  std::cout << "matched level = " << _matchedLevel << std::endl;

  buildMapsToOriginalMesh(mesh_data, _matchedLevel);

  _fieldvel[_levelS].addField("U", Field::Location::Nodal, 0.);
  _fieldvel[_levelS].addField("V", Field::Location::Nodal, 0.);
  _fieldvel[_levelS].addField("W", Field::Location::Nodal, 0.);

  auto t_end = std::chrono::high_resolution_clock::now();

  const double elapsed =
    std::chrono::duration<double>(t_end - t_start).count();

  std::cout << "Ex35::initializeVELMesh time = "
            << elapsed << " s" << std::endl;

}

void Ex35::initializeLSMesh() {

  auto t_start = std::chrono::high_resolution_clock::now();

  InitialMeshData seed;
  if (_dim == 2 && _elTypevel[0][0] == 3)
    seed = buildSeed1x2();
  else if (_dim == 3 && _elTypevel[0][0] == 0)
    seed = buildSeed1x1x2();

  _elLevel0[0] = seed.elLevel;
  _elType0[0]  = seed.elType;
  _elTplgy0[0] = seed.elTplgy;
  _X0[0]       = seed.X;
  _dim         = _X0[0].size();

  std::vector<unsigned> candidateIndices(_mesh0[0].numNodes());
  std::iota(candidateIndices.begin(), candidateIndices.end(), 0u);
  dedupNodesHash(_mesh0[0], candidateIndices);

  _mesh0[0].resetAllFathersToNoFather();
  _mesh0[0].buildNodeToElementAdjacency();
  _mesh0[0].buildFaceNeighborsFromNodeToElement();

  _mesh1[0] = _mesh0[0];

  auto t_end = std::chrono::high_resolution_clock::now();

  const double elapsed =
    std::chrono::duration<double>(t_end - t_start).count();

  std::cout << "Ex35::initializeLSMesh time = "
            << elapsed << " s" << std::endl;
}

void Ex35::initializeField_Ball(std::vector<double> xc, double r) {
  if (xc.size() != _dim)
      throw std::runtime_error("WRONG DIMENSION");

  auto t_start = std::chrono::high_resolution_clock::now();

  _field0[0].clear();

  _filename = "../output/refined_mesh" + std::to_string(_dim) + "D_level" + std::to_string(_levelN) + ".";
  _filenamefixed = "../output/fixed_mesh" + std::to_string(_dim) + "D_level" + std::to_string(_levelN) + ".";

  double eps = (_dim == 2)
    ? 1. / pow(2, std::max(_levelN  - 5u, 1u))
    : 1. / pow(2, std::max(_levelN  - 2u, 1u));

  _neighMode = 3; // 0=no-ring, 1=vertices, 2=faces, 3=hybrid
  for (unsigned l = 1; l <= _levelN; l++) {
    _mesh0[l - 1].setRefinementFromBallLevelSetCrossing_OneRing(xc, r, _neighMode);
    _mesh0[l - 1].adjustAMRForOneLevelDiscontinuity();
    refineAndProjectMesh(_elProj, _mesh0[l - 1], _mesh0[l]);
    _field0[l].clear();
  }

  const unsigned topLevel = _levelN;

  PsiBall psi(xc, r, eps);
  _mollifier = psi._m;

  _field0[topLevel].addField("Psi", Field::Location::Nodal, 1.);
  const unsigned psiId0_init = _field0[topLevel].id("Psi");
  auto& Psi0 = _field0[topLevel].getById(psiId0_init);
  for (std::size_t k = 0; k < Psi0.size(); ++k) {
    const std::vector<double> x = _field0[topLevel].dofCoordById(psiId0_init, k);
    Psi0[k] = psi(x);
  }

  std::cout << "Iteration = " << 0 << " Number of Points = " << _mesh0[topLevel].X()[0].size() << std::endl;
  writeMeshFieldVTU(_filename + std::to_string(0) + ".vtu", _field0[topLevel]);

  for (unsigned l = 0; l <= _levelN; l++) _field2[l] = _field0[l];

  auto t_end = std::chrono::high_resolution_clock::now();

  const double elapsed =
    std::chrono::duration<double>(t_end - t_start).count();

  std::cout << "Ex35::initializeField time = "
            << elapsed << " s" << std::endl;
}

void Ex35::initMarkers() {
  auto t_start = std::chrono::high_resolution_clock::now();

  Reinit reinit(_elProj, _field0, _field0[_levelN].id("Psi"), _mollifier);
  reinit.computeMarkersAdvection(_markers);

  auto t_end = std::chrono::high_resolution_clock::now();

  const double elapsed =
    std::chrono::duration<double>(t_end - t_start).count();

  std::cout << "Ex35::initMarkers time = "
            << elapsed << " s" << std::endl;
}

void Ex35::advectField(const std::vector<std::vector<double>>& vel, const int it, const double dt) {
  const char* C_RST = "\033[0m";
  const char* C_BLD = "\033[1m";
  const char* C_CYN = "\033[36m";
  const char* C_YEL = "\033[33m";
  const char* C_GRN = "\033[32m";

  auto print_reinit_header = [&]() {
    std::cout << "\n"
              << C_BLD << C_CYN
              << " --------------------------------------------------- \n"
              << "         AMR Advection\n"
              << " --------------------------------------------------- "
              << C_RST << std::endl;
  };
  print_reinit_header();
  if (vel.size() < _dim) {
    throw std::runtime_error("Ex35::setVelocity: wrong number of velocity components");
  }

  // -----------------------------
  // 1) set coarse velocity on level 0
  // -----------------------------
  const unsigned uId0 = _fieldvel[_levelS].id("U");
  const unsigned vId0 = _fieldvel[_levelS].id("V");

  auto& U0 = _fieldvel[_levelS].getById(uId0);
  auto& V0 = _fieldvel[_levelS].getById(vId0);

  if (vel[0].size() != U0.size() || vel[1].size() != V0.size()) {
    throw std::runtime_error("Ex35::setVelocity: wrong coarse velocity size");
  }

  for (std::size_t k = 0; k < U0.size(); ++k) U0[_origToThisNode[k]] = vel[0][k];
  for (std::size_t k = 0; k < V0.size(); ++k) V0[_origToThisNode[k]] = vel[1][k];

  unsigned wId0 = 0;
  if (_dim == 3) {
    wId0 = _fieldvel[_levelS].id("W");   
    auto& W0 = _fieldvel[_levelS].getById(wId0);

    if (vel.size() < 3 || vel[2].size() != W0.size()) {
      throw std::runtime_error("Ex35::setVelocity: wrong coarse W size");
    }

    for (std::size_t k = 0; k < W0.size(); ++k) W0[_origToThisNode[k]] = vel[2][k];
  }

  if (it % 10 == 0)
  writeMeshFieldVTU(_filenamefixed + std::to_string(it) + ".vtu", _fieldvel[_levelS]);

  // -----------------------------
  // 2) interpolate velocity on markers
  // -----------------------------
  if (_markers.empty()) {
    return;
  }

  if (_markers.size() != _dim) {
    throw std::runtime_error("Ex35::setVelocity: _markers has wrong dimension");
  }

  const std::size_t nMarkers = _markers[0].size();
  for (unsigned d = 1; d < _dim; ++d) {
    if (_markers[d].size() != nMarkers) {
      throw std::runtime_error("Ex35::setVelocity: inconsistent _markers sizes");
    }
  }

  std::vector<PointLocatorResult> in, out;

  PointLocator pl1(_meshvel[0], 0.1);
  out.clear();
  pl1.locateAll(out, _markers);

  for (unsigned l = 1; l <= _levelS; ++l) {
    std::swap(in, out);
    _meshvel[l - 1].projectPointLocatorResultsToNextLevel(_meshvel[l], in, out);
  }

  if (out.size() != nMarkers) {
    throw std::runtime_error("Ex35::setVelocity: wrong number of located marker points after projection");
  }

  _markerVelocity.clear();
  _markerVelocity.resize(_dim);
  _markerVelocity[0].assign(nMarkers, 0.0);
  _markerVelocity[1].assign(nMarkers, 0.0);
  if (_dim == 3) {
    _markerVelocity[2].assign(nMarkers, 0.0);
  }

  _fieldvel[_levelS].evalNodalAtLocatedPointsById(uId0, out, _elProj, _markerVelocity[0], 0.0);
  _fieldvel[_levelS].evalNodalAtLocatedPointsById(vId0, out, _elProj, _markerVelocity[1], 0.0);

  if (_dim == 3) {
    _fieldvel[_levelS].evalNodalAtLocatedPointsById(wId0, out, _elProj, _markerVelocity[2], 0.0);
  }

  // -----------------------------
  // 3) advect forward and update grid
  // -----------------------------

  AdvectMarkers::forward(_markers, _markerVelocity, dt);
  PointLocator pl2(_mesh0[0], 0.1);
  pl2.locateAll(out, _markers);

  _mesh1[0].setRefinementFromLocatedPoints_OneRing(out, _neighMode);
  _mesh1[0].adjustAMRForOneLevelDiscontinuity();
  _field1[0].clear();

  for (unsigned l = 1; l <= _levelN; l++) {
    _mesh1[l].clearAllData();
    refineAndProjectMesh(_elProj, _mesh1[l - 1], _mesh1[l]);

    std::swap(in, out);
    _mesh1[l - 1].projectPointLocatorResultsToNextLevel(_mesh1[l], in, out);
    _mesh1[l].setRefinementFromLocatedPoints_OneRing(out, _neighMode);
    _mesh1[l].adjustAMRForOneLevelDiscontinuity();

    _field1[l].clear();
  }

  // -----------------------------
  // 4) advect backward and update field
  // -----------------------------

  _nodeCoords = _mesh1[_levelN].X();

  out.clear();
  pl1.locateAll(out, _nodeCoords);

  for (unsigned l = 1; l <= _levelS; ++l) {
    std::swap(in, out);
    _meshvel[l - 1].projectPointLocatorResultsToNextLevel(_meshvel[l], in, out);
  }

  if (out.size() != _mesh1[_levelN].X()[0].size()) {
    throw std::runtime_error("Ex35::setVelocity: wrong number of located marker points after projection");
  }

  _nodeVelocity.clear();
  _nodeVelocity.resize(_dim);
  _nodeVelocity[0].assign(nMarkers, 0.0);
  _nodeVelocity[1].assign(nMarkers, 0.0);
  if (_dim == 3) {
    _nodeVelocity[2].assign(nMarkers, 0.0);
  }

  _fieldvel[_levelS].evalNodalAtLocatedPointsById(uId0, out, _elProj, _nodeVelocity[0], 0.0);
  _fieldvel[_levelS].evalNodalAtLocatedPointsById(vId0, out, _elProj, _nodeVelocity[1], 0.0);

  if (_dim == 3) {
    _fieldvel[_levelS].evalNodalAtLocatedPointsById(wId0, out, _elProj, _nodeVelocity[2], 0.0);
  }

  AdvectMarkers::backward(_nodeCoords, _nodeVelocity, dt);
  out.clear();
  pl2.locateAll(out, _nodeCoords);
  for (unsigned l = 1; l <= _levelN; l++) {
    std::swap(in, out);
    _mesh0[l - 1].projectPointLocatorResultsToNextLevel(_mesh0[l], in, out);
  }

  _field1[_levelN].addField("Psi", Field::Location::Nodal, 1.);
  const unsigned psiId1 = _field1[_levelN].id("Psi");
  auto& Psi1 = _field1[_levelN].getById(psiId1);
  const unsigned psiId0 = _field0[_levelN].id("Psi");
  _field0[_levelN].evalNodalAtLocatedPointsById(psiId0, out, _elProj, Psi1, 1.0);

  for (unsigned l = 0; l <= _levelN; l++) swap(_field0[l], _field1[l]);

  if (it % 100 == 0)
  writeMeshFieldVTU(_filename + std::to_string(it) + ".vtu", _field0[_levelN]);

}

void Ex35::interpolatePsiOnNodes(const std::vector<std::vector<double>>& node_coords, std::vector<double>& psi) const {

  auto t_start = std::chrono::high_resolution_clock::now();

  std::vector<PointLocatorResult> in, out;
  PointLocator pl(_mesh0[0], 0.1);
  out.clear();

  pl.locateAll(out, node_coords);
  for (unsigned l = 1; l <= _levelN; l++) {
    std::swap(in, out);
    _mesh0[l - 1].projectPointLocatorResultsToNextLevel(_mesh0[l], in, out);
  }

  psi.clear();
  psi.assign(node_coords[0].size(), 0.0);
  const unsigned psiId0 = _field0[_levelN].id("Psi");
  _field0[_levelN].evalNodalAtLocatedPointsById(psiId0, out, _elProj, psi, 1.0);

  auto t_end = std::chrono::high_resolution_clock::now();

  const double elapsed =
    std::chrono::duration<double>(t_end - t_start).count();

  std::cout << "Ex35::interpolatePsiOnNodes time = "
            << elapsed << " s" << std::endl;

}

CellMarkersData Ex35::getCellMarkers(const InitialMeshData & mesh_data) {

  auto t_start = std::chrono::high_resolution_clock::now();

  const std::size_t nMarkers = _markers[0].size();
  for (unsigned d = 1; d < _dim; ++d) {
    if (_markers[d].size() != nMarkers) {
      throw std::runtime_error("Ex35::getCellMarkers: inconsistent _markers sizes");
    }
  }

  std::vector<PointLocatorResult> in, out;

  PointLocator pl(_meshvel[0], 0.1);
  out.clear();
  pl.locateAll(out, _markers);

  for (unsigned l = 1; l <= _matchedLevel; ++l) {
    std::swap(in, out);
    _meshvel[l - 1].projectPointLocatorResultsToNextLevel(_meshvel[l], in, out);
  }

  if (out.size() != nMarkers) {
    throw std::runtime_error("Ex35::getCellMarkers: wrong number of located marker points after projection");
  }

  CellMarkersData data;
  data.markers = _markers;
  data.own.reserve(nMarkers);

  // --------------------------------------------------
  // 1) own markers only
  // --------------------------------------------------
  for (std::size_t im = 0; im < nMarkers; ++im) {
    if (!out[im].ok) continue;

    const unsigned thisElem = static_cast<unsigned>(out[im].elem);

    if (thisElem >= _thisToOrigElem.size()) {
      throw std::runtime_error("Ex35::getCellMarkers: thisElem out of _thisToOrigElem range");
    }

    const unsigned origElem = _thisToOrigElem[thisElem];
    if (origElem == INVALID_ID) {
      throw std::runtime_error("Ex35::getCellMarkers: invalid origElem in _thisToOrigElem");
    }

    data.own[mesh_data.elID[origElem]].push_back(static_cast<unsigned>(im));
  }

  // --------------------------------------------------
  // 2) combined markers = own + cut neighbours
  // --------------------------------------------------
  // data.combined = data.own;

  // const auto& neighAll = _meshvel[_matchedLevel].neighbors();

  // if (neighAll.size() != _meshvel[_matchedLevel].numElements()) {
  //   throw std::runtime_error("Ex35::getCellMarkers: neighbors not built on matched level");
  // }

  // for (const auto& kv : data.own) {
  //   const unsigned origElem = kv.first;

  //   if (origElem >= _origToThisElem.size() || _origToThisElem[origElem] == INVALID_ID) {
  //     throw std::runtime_error("Ex35::getCellMarkers: invalid _origToThisElem map");
  //   }

  //   const unsigned thisElem = _origToThisElem[origElem];

  //   if (thisElem >= neighAll.size()) {
  //     throw std::runtime_error("Ex35::getCellMarkers: thisElem out of neighbours range");
  //   }

  //   auto& markersComb = data.combined[origElem];

  //   for (unsigned ineigh = 0; ineigh < neighAll[thisElem].size(); ++ineigh) {
  //     const unsigned thisNeigh = neighAll[thisElem][ineigh];

  //     if (thisNeigh == INVALID_ID) continue;
  //     if (thisNeigh >= _thisToOrigElem.size()) continue;
  //     if (_thisToOrigElem[thisNeigh] == INVALID_ID) continue;

  //     const unsigned origNeigh = _thisToOrigElem[thisNeigh];

  //     // add only neighbours that are cut
  //     auto itCutNeigh = data.own.find(origNeigh);
  //     if (itCutNeigh == data.own.end()) continue;

  //     const auto& neighMarkers = itCutNeigh->second;
  //     markersComb.insert(markersComb.end(), neighMarkers.begin(), neighMarkers.end());
  //   }

  //   std::sort(markersComb.begin(), markersComb.end());
  //   markersComb.erase(std::unique(markersComb.begin(), markersComb.end()), markersComb.end());
  // }

  auto t_end = std::chrono::high_resolution_clock::now();

  const double elapsed =
    std::chrono::duration<double>(t_end - t_start).count();

  std::cout << "Ex35::getCellMarkers time = "
            << elapsed << " s" << std::endl;

  return data;
}

void Ex35::reinit() {
  const char* C_RST = "\033[0m";
  const char* C_BLD = "\033[1m";
  const char* C_CYN = "\033[36m";
  const char* C_YEL = "\033[33m";
  const char* C_GRN = "\033[32m";

  auto print_reinit_header = [&]() {
    std::cout << "\n"
              << C_BLD << C_CYN
              << " --------------------------------------------------- \n"
              << "         AMR Reinit\n"
              << " --------------------------------------------------- "
              << C_RST << std::endl;
      
  };
  print_reinit_header();
  const unsigned psiId0 = _field0[_levelN].id("Psi");
  Reinit reinit(_elProj, _field0, psiId0, _mollifier);
  reinit.reinitializeSignedDistance();
}

Ex35::~Ex35() = default;


InitialMeshData Ex35::buildSeed1x2() const {
  InitialMeshData seed;

  const double xmin = 0;
  const double xmax = 1;
  const double ymin = 0;
  const double ymax = 2;

  const double xmid = 0.5 * (xmin + xmax);
  const double ymid = 0.5 * (ymin + ymax);

  const double yq1  = 0.5 * (ymin + ymid);
  const double yq3  = 0.5 * (ymid + ymax);

  // -----------------------------
  // nodi globali della mesh 1x2 Quad9
  //
  // griglia 3 x 5:
  //
  // 12 -- 13 -- 14
  //  9 -- 10 -- 11
  //  6 --  7 --  8
  //  3 --  4 --  5
  //  0 --  1 --  2
  // -----------------------------
  seed.X.resize(2);

  seed.X[0] = {
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    xmin, xmid, xmax
  };

  seed.X[1] = {
    ymin, ymin, ymin,
    yq1,  yq1,  yq1,
    ymid, ymid, ymid,
    yq3,  yq3,  yq3,
    ymax, ymax, ymax
  };

  seed.elLevel = {0, 0};
  seed.elType  = {3, 3}; // Quad9

  // convenzione:
  // [v0, v1, v2, v3, e01, e12, e23, e30, center]

  seed.elTplgy.resize(2);

  // elemento basso
  seed.elTplgy[0] = {
    0,  // bottom-left
    2,  // bottom-right
    8,  // top-right
    6,  // top-left
    1,  // mid bottom
    5,  // mid right
    7,  // mid top
    3,  // mid left
    4   // center
  };

  // elemento alto
  seed.elTplgy[1] = {
    6,   // bottom-left
    8,   // bottom-right
    14,  // top-right
    12,  // top-left
    7,   // mid bottom
    11,  // mid right
    13,  // mid top
    9,   // mid left
    10   // center
  };

  return seed;
}

InitialMeshData Ex35::buildSeed1x1x2() const {
  InitialMeshData seed;

  const double xmin = 0;
  const double xmax = 1;
  const double ymin = 0;
  const double ymax = 1;
  const double zmin = 0;
  const double zmax = 2;

  const double xmid = 0.5 * (xmin + xmax);
  const double ymid = 0.5 * (ymin + ymax);
  const double zmid = 0.5 * (zmin + zmax);

  const double zq1  = 0.5 * (zmin + zmid);
  const double zq3  = 0.5 * (zmid + zmax);

  // -----------------------------
  // nodi globali della mesh 1x1x2 Hex27
  //
  // griglia 9 x 5:
  //

  seed.X.resize(3);

  seed.X[0] = {
    // z =0
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    // z = 0.5
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    // z = 1
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    // z = 1.5
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    // z = 2
    xmin, xmid, xmax,
    xmin, xmid, xmax,
    xmin, xmid, xmax
  };

  seed.X[1] = {
    // z = 0
    ymin, ymin, ymin,
    ymid, ymid, ymid,
    ymax, ymax, ymax,

    // z = 0.5
    ymin, ymin, ymin,
    ymid, ymid, ymid,
    ymax, ymax, ymax,

    // z = 1
    ymin, ymin, ymin,
    ymid, ymid, ymid,
    ymax, ymax, ymax,

    // z = 1.5
    ymin, ymin, ymin,
    ymid, ymid, ymid,
    ymax, ymax, ymax,

    // z = 2
    ymin, ymin, ymin,
    ymid, ymid, ymid,
    ymax, ymax, ymax
  };

  seed.X[2] = {
    zmin, zmin, zmin,
    zmin, zmin, zmin,
    zmin, zmin, zmin,

    zq1,  zq1,  zq1,
    zq1,  zq1,  zq1,
    zq1,  zq1,  zq1,

    zmid, zmid, zmid,
    zmid, zmid, zmid,
    zmid, zmid, zmid,

    zq3,  zq3,  zq3,
    zq3,  zq3,  zq3,
    zq3,  zq3,  zq3,

    zmax, zmax, zmax,
    zmax, zmax, zmax,
    zmax, zmax, zmax
  };

  seed.elLevel = {0, 0};
  seed.elType  = {0, 0}; // Quad9

  // convenzione:
  // [v0, v1, v2, v3, e01, e12, e23, e30, center]

  seed.elTplgy.resize(2);

  // elemento basso
  seed.elTplgy[0] = {
    0, 2, 8, 6, 18, 20, 26, 24,
    1, 5, 7, 3,
    19, 23, 25, 21,
    9, 11, 17, 15,
    10, 14, 16, 12, 4, 22,
    13
  };

  // elemento alto
  seed.elTplgy[1] = {
    18, 20, 26, 24, 36, 38, 44, 42,
    19, 23, 25, 21,
    37, 41, 43, 39,
    27, 29, 35, 33,
    28, 32, 34, 30, 22, 40,
    31
  };

  return seed;
}


void Ex35::buildMapsToOriginalMesh(const InitialMeshData& mesh_data, const unsigned level) {
  const auto& Xnew      = _Xvel[level];
  const auto& elTplgyNew = _elTplgyvel[level];
  const auto& elTypeNew  = _elTypevel[level];

  const unsigned nOrigNodes = mesh_data.X[0].size();
  const unsigned nNewNodes  = Xnew[0].size();

  _origToThisNode.assign(nOrigNodes, INVALID_ID);
  _thisToOrigNode.assign(nNewNodes, INVALID_ID);

  // ---------------------------------
  // 1) node map by coordinates
  // ---------------------------------
  std::unordered_map<std::string, unsigned> newNodeFromKey;
  newNodeFromKey.reserve(nNewNodes);

  for (unsigned j = 0; j < nNewNodes; ++j) {
    const std::string key = makeNodeKey(Xnew, _dim, j);
    const auto ok = newNodeFromKey.emplace(key, j).second;
    if (!ok) {
      throw std::runtime_error("buildMapsToOriginalMesh: duplicated node key in reconstructed mesh");
    }
  }

  for (unsigned i = 0; i < nOrigNodes; ++i) {
    const std::string key = makeNodeKey(mesh_data.X, _dim, i);
    const auto it = newNodeFromKey.find(key);
    if (it == newNodeFromKey.end()) {
      throw std::runtime_error("buildMapsToOriginalMesh: original node not found in reconstructed mesh");
    }

    const unsigned j = it->second;
    _origToThisNode[i] = j;
    _thisToOrigNode[j] = i;
  }

  // ---------------------------------
  // 2) element map by connectivity
  // ---------------------------------
  const unsigned nOrigElem = mesh_data.elTplgy.size();
  const unsigned nNewElem  = elTplgyNew.size();

  _origToThisElem.assign(nOrigElem, INVALID_ID);
  _thisToOrigElem.assign(nNewElem, INVALID_ID);

  std::unordered_map<std::string, unsigned> newElemFromKey;
  newElemFromKey.reserve(nNewElem);

  for (unsigned e = 0; e < nNewElem; ++e) {
    const std::string key = makeElemKey(elTplgyNew[e]);
    const auto ok = newElemFromKey.emplace(key, e).second;
    if (!ok) {
      throw std::runtime_error("buildMapsToOriginalMesh: duplicated element key in reconstructed mesh");
    }
  }

  for (unsigned e0 = 0; e0 < nOrigElem; ++e0) {
    std::vector<unsigned> mappedConn;
    mappedConn.reserve(mesh_data.elTplgy[e0].size());

    for (const auto inodeOrig : mesh_data.elTplgy[e0]) {
      if (inodeOrig >= _origToThisNode.size() || _origToThisNode[inodeOrig] == INVALID_ID) {
        throw std::runtime_error("buildMapsToOriginalMesh: invalid node map while building element map");
      }
      mappedConn.push_back(_origToThisNode[inodeOrig]);
    }

    const std::string key = makeElemKey(mappedConn);
    const auto it = newElemFromKey.find(key);
    if (it == newElemFromKey.end()) {
      throw std::runtime_error("buildMapsToOriginalMesh: original element not found in reconstructed mesh");
    }

    const unsigned e1 = it->second;

    if (mesh_data.elType[e0] != elTypeNew[e1]) {
      throw std::runtime_error("buildMapsToOriginalMesh: matched element has wrong type");
    }

    _origToThisElem[e0] = e1;
    _thisToOrigElem[e1] = e0;
  }

  std::cout << "node map built: " << _origToThisNode.size() << " nodes\n";
  std::cout << "elem map built: " << _origToThisElem.size() << " elements\n";
}

std::vector<double> Ex35::interpolatePsiOnOriginalNodes() const {

  if (_origToThisNode.empty()) {
    throw std::runtime_error("Ex35::interpolatePsiOnOriginalNodes: node map not built");
  }

  auto t_start = std::chrono::high_resolution_clock::now();

  const unsigned nOrigNodes = _origToThisNode.size();

  const auto& Xvel = _meshvel[_matchedLevel].X();

  if (Xvel.size() != _dim) {
    throw std::runtime_error("Ex35::interpolatePsiOnOriginalNodes: wrong Xvel dimension");
  }

  std::vector<std::vector<double>> Xquery(_dim, std::vector<double>(nOrigNodes, 0.0));

  for (unsigned inodeOrig = 0; inodeOrig < nOrigNodes; ++inodeOrig) {
    const unsigned inodeThis = _origToThisNode[inodeOrig];

    if (inodeThis == INVALID_ID) {
      throw std::runtime_error("Ex35::interpolatePsiOnOriginalNodes: invalid mapped node");
    }

    for (unsigned d = 0; d < _dim; ++d) {
      if (inodeThis >= Xvel[d].size()) {
        throw std::runtime_error("Ex35::interpolatePsiOnOriginalNodes: mapped node out of range");
      }
      Xquery[d][inodeOrig] = Xvel[d][inodeThis];
    }
  }

  std::vector<PointLocatorResult> in, out;
  PointLocator pl(_mesh0[0], 0.1);
  pl.locateAll(out, Xquery);

  for (unsigned l = 1; l <= _levelN; ++l) {
    std::swap(in, out);
    _mesh0[l - 1].projectPointLocatorResultsToNextLevel(_mesh0[l], in, out);
  }

  if (out.size() != nOrigNodes) {
    throw std::runtime_error("Ex35::interpolatePsiOnOriginalNodes: wrong number of located points");
  }

  const unsigned psiId = _field0[_levelN].id("Psi");

  std::vector<double> psiOrig;
  _field0[_levelN].evalNodalAtLocatedPointsById(psiId, out, _elProj, psiOrig,
                                                std::numeric_limits<double>::quiet_NaN());

  if (psiOrig.size() != nOrigNodes) {
    throw std::runtime_error("Ex35::interpolatePsiOnOriginalNodes: wrong psi vector size");
  }

  for (unsigned i = 0; i < nOrigNodes; ++i) {
    if (!std::isfinite(psiOrig[i])) {
      throw std::runtime_error("Ex35::interpolatePsiOnOriginalNodes: non-finite Psi on queried node");
    }
  }

  auto t_end = std::chrono::high_resolution_clock::now();

  const double elapsed =
    std::chrono::duration<double>(t_end - t_start).count();

  std::cout << "Ex35::interpolatePsiOnOriginalNodes time = "
            << elapsed << " s" << std::endl;

  return psiOrig;
}




























