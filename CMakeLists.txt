cmake_minimum_required(VERSION 3.31)
project(sbh)

set(CMAKE_CXX_STANDARD 23)

set(CMAKE_DEFAULT_PROJECT_FILE src/main.cpp)

set(SOURCE_FILES src/main.cpp
		src/config.h
        src/instance_parser.h
        src/instance_parser.cpp
		src/sbh.h
		src/sbh.cpp
		src/levenshtein.h
		src/levenshtein.cpp
        src/aco.cpp
        src/aco.h
		src/settings.h
)

add_executable(sbh ${SOURCE_FILES})
