cmake_minimum_required(VERSION 3.31)
project(sbh)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_DEFAULT_PROJECT_FILE src/main.cpp)

set(SOURCE_FILES src/main.cpp
		src/config.h
		src/instance_runner.h
		src/instance_runner.cpp
		src/sbh.h
		src/sbh.cpp
		src/levenshtein.h
		src/levenshtein.cpp
		src/ant.cpp
		src/ant.h
)

add_executable(sbh ${SOURCE_FILES})
