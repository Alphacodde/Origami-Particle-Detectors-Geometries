cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(OrigamiDetector)

# --- GEANT4 ---
# ui_all / vis_all pull in Qt/OpenGL visualization drivers if present on your build.
# If your GEANT4 was built without Qt, drop ui_all/vis_all and this still builds fine
# (you just lose the interactive viewer; macros/batch mode still work).
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# --- CADMesh ---
# CADMesh is a header-mostly third-party library for importing STL/PLY/OBJ meshes
# into GEANT4 as G4VSolid. Clone it alongside this project:
#   git clone https://github.com/christopherpoole/CADMesh.git
# and point CADMESH_DIR at it, e.g.:
#   cmake .. -DCADMESH_DIR=/path/to/CADMesh
if(NOT DEFINED CADMESH_DIR)
  set(CADMESH_DIR "${CMAKE_SOURCE_DIR}/../CADMesh")
endif()
# CADMesh's own CMakeLists.txt builds a full shared library and requires
# tetgen/assimp, which we don't want or have on Windows. Instead of
# add_subdirectory(), we replicate just the one step we actually need: the
# configure_file() call that turns Configuration.hh.in into a real
# Configuration.hh (this mirrors CADMesh's own CMakeLists.txt line for line;
# DEFAULT_READER defaults to "BuiltIn", which is what we want since we don't
# need the assimp-based OBJ/PLY readers, only the built-in STL reader).
set(DEFAULT_READER "BuiltIn" CACHE STRING "The CAD file reader to use by default.")
configure_file(${CADMESH_DIR}/include/Configuration.hh.in "${CMAKE_BINARY_DIR}/include/Configuration.hh" @ONLY)
include_directories(${CADMESH_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/include)  # generated Configuration.hh lands here
# CADMesh vendors its own tiny STL/OBJ parsers; no extra linking needed for the
# STL path (the OBJ/PLY path in some CADMesh versions wants assimp - we don't need it).

include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# CADMesh is NOT purely header-only: CADMeshTemplate.hh, TessellatedMesh.hh etc.
# only declare things like FromSTL() - the implementations live in CADMesh's
# own src/*.cc, which normally get compiled into a separate "cadmesh" shared
# library by CADMesh's own CMakeLists.txt. We deliberately did NOT
# add_subdirectory() that (it hard-requires tetgen + assimp, which we don't
# have on Windows), so instead we compile CADMesh's needed .cc files directly
# into our own executable. We EXCLUDE:
#   - ASSIMPReader.cc: the only one of the reader .cc files that actually
#     needs the assimp library, which we don't have installed. OBJReader.cc
#     and PLYReader.cc are CADMesh's own lightweight parsers (unrelated to
#     ASSIMPReader.cc despite similar names) and DON'T need assimp - they
#     must be included because BuiltInReader.cc references them directly as
#     part of its reader dispatch, even though we only invoke the STL path
#     at runtime (link-time resolution needs all symbols present regardless).
#   - TetrahedralMesh.cc: CADMesh's tetgen-backed volume mesh path, unused by
#     this project (we only use CADMeshTemplate<TessellatedMesh>::FromSTL()).
# If a future LNK2019/LNK2001 surfaces for a symbol not covered here, that
# means another CADMesh .cc needs to be added back in - check which reader/
# component it belongs to before blanket-including everything (that's how we
# ended up needing to exclude ASSIMPReader.cc specifically, to avoid an
# assimp dependency we don't have).
set(CADMESH_EXCLUDE_SOURCES
  "${CADMESH_DIR}/src/ASSIMPReader.cc"
  "${CADMESH_DIR}/src/TetrahedralMesh.cc"
)
file(GLOB cadmesh_all_sources ${CADMESH_DIR}/src/*.cc)
list(REMOVE_ITEM cadmesh_all_sources ${CADMESH_EXCLUDE_SOURCES})
list(APPEND sources ${cadmesh_all_sources})

add_executable(origamiDet origamiDet.cc ${sources} ${headers})
target_link_libraries(origamiDet ${Geant4_LIBRARIES})

# Copy macros/geometry into the build dir so you can run from wherever you build
file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/geometry DESTINATION ${PROJECT_BINARY_DIR})

# Also copy them directly to the configuration folder (Release/Debug) containing the executable
add_custom_command(TARGET origamiDet POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${PROJECT_SOURCE_DIR}/macros
    $<TARGET_FILE_DIR:origamiDet>/macros
)
add_custom_command(TARGET origamiDet POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${PROJECT_SOURCE_DIR}/geometry
    $<TARGET_FILE_DIR:origamiDet>/geometry
)