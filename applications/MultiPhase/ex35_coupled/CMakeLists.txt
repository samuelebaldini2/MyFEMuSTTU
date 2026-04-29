cmake_minimum_required(VERSION 3.20)
project(ex35core LANGUAGES C CXX)

set(TTU_ROOT "/home/samuele/femus_ttu/MyFEMuSTTU")

if (NOT TARGET Eigen3::Eigen)
  find_package(Eigen3 REQUIRED NO_MODULE)
endif()

# =========================================================
# jsoncpp locale
# =========================================================
file(GLOB JSONCPP_SRC CONFIGURE_DEPENDS
  "${TTU_ROOT}/external/jsoncpp/jsoncpp-src-0.5.0/src/lib_json/*.cpp"
  "${TTU_ROOT}/external/jsoncpp/jsoncpp-src-0.5.0/src/*.cpp"
)

add_library(jsoncpp_ttu STATIC ${JSONCPP_SRC})

set_target_properties(jsoncpp_ttu PROPERTIES
  POSITION_INDEPENDENT_CODE ON
)

target_include_directories(jsoncpp_ttu PUBLIC
  "${TTU_ROOT}/external/jsoncpp/jsoncpp-src-0.5.0/include"
)

target_compile_features(jsoncpp_ttu PUBLIC cxx_std_17)

# =========================================================
# adept locale
# =========================================================
file(GLOB ADEPT_SRC CONFIGURE_DEPENDS
  "${TTU_ROOT}/external/adept/adept-1.1/adept/*.cpp"
)

add_library(adept_ttu STATIC ${ADEPT_SRC})

set_target_properties(adept_ttu PROPERTIES
  POSITION_INDEPENDENT_CODE ON
)

target_include_directories(adept_ttu PUBLIC
  "${TTU_ROOT}/external/adept/adept-1.1/include"
)

target_compile_features(adept_ttu PUBLIC cxx_std_17)

# =========================================================
# b64 locale
# =========================================================
file(GLOB B64_SRC CONFIGURE_DEPENDS
  "${TTU_ROOT}/external/b64/b64-1.4.2/src/*.c"
  "${TTU_ROOT}/external/b64/b64-1.4.2/src/*.cpp"
)

add_library(b64_ttu STATIC ${B64_SRC})

set_target_properties(b64_ttu PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE C
)

target_include_directories(b64_ttu PUBLIC
  "${TTU_ROOT}/external/b64/b64-1.4.2/include"
)

# =========================================================
# backend TTU locale
# =========================================================
file(GLOB_RECURSE FEMUS_TTU_LOCAL_SRC CONFIGURE_DEPENDS
  "${TTU_ROOT}/src/*.cpp"
)

# Escludi file sperimentali/non presenti nel CMake originale di src
list(REMOVE_ITEM FEMUS_TTU_LOCAL_SRC
  "${TTU_ROOT}/src/mesh/refinement/MeshRefinementNew.cpp"
)

add_library(femus_ttu_local SHARED ${FEMUS_TTU_LOCAL_SRC})

target_compile_features(femus_ttu_local PUBLIC cxx_std_17)

target_include_directories(femus_ttu_local
  PUBLIC
    ${TTU_ROOT}/include
    ${TTU_ROOT}/src
    ${TTU_ROOT}/src/algebra
    ${TTU_ROOT}/src/enums
    ${TTU_ROOT}/src/equations
    ${TTU_ROOT}/src/equations/assemble
    ${TTU_ROOT}/src/fe
    ${TTU_ROOT}/src/ism
    ${TTU_ROOT}/src/mesh
    ${TTU_ROOT}/src/mesh/gencase
    ${TTU_ROOT}/src/mesh/geom_elements
    ${TTU_ROOT}/src/mesh/input
    ${TTU_ROOT}/src/mesh/refinement
    ${TTU_ROOT}/src/mesh/partitioning
    ${TTU_ROOT}/src/parallel
    ${TTU_ROOT}/src/physics
    ${TTU_ROOT}/src/quadrature
    ${TTU_ROOT}/src/quadrature/cutfem
    ${TTU_ROOT}/src/solution
    ${TTU_ROOT}/src/solution/output
    ${TTU_ROOT}/src/utils
    ${TTU_ROOT}/src/utils/input_parser
    ${TTU_ROOT}/src/uq
    ${TTU_ROOT}/external/adept/adept-1.1/include
    ${TTU_ROOT}/external/b64/b64-1.4.2/include
    ${TTU_ROOT}/external/jsoncpp/jsoncpp-src-0.5.0/include
)

target_link_libraries(femus_ttu_local
  PUBLIC
    jsoncpp_ttu
    adept_ttu
    b64_ttu
    femus
    Eigen3::Eigen
)

# =========================================================
# ex35
# =========================================================
add_library(ex35core STATIC
  ex35_coupled.cpp
)

target_compile_features(ex35core PUBLIC cxx_std_17)

target_include_directories(ex35core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(ex35core PUBLIC
  femus_ttu_local
)

# =========================================================
# gperftools (optional)
# =========================================================
find_library(GPERFTOOLS_LIB profiler)

if (GPERFTOOLS_LIB)
  target_link_libraries(ex35core PUBLIC ${GPERFTOOLS_LIB})
endif()
