cmake_minimum_required(VERSION 3.10)
project(pycvs)

set(PYBIND11_FINDPYTHON ON)
set(Python3_ROOT_DIR "C:/Anaconda3" CACHE PATH "Python 3 installation root")
set(Python3_EXECUTABLE "C:/Anaconda3/python.exe" CACHE FILEPATH "Path to Python 3 interpreter")

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)


pybind11_add_module(pycvs extractCurves.cpp PixelGraph.cpp UnionFind.cpp)

target_include_directories(pycvs PRIVATE
    ${PYTHON_INCLUDE_DIR}
)

target_link_libraries(pycvs PRIVATE pybind11::module)
