md build
cd build
cmake .. -G "Unix Makefiles" -D use_OpenMP=OFF -D CMAKE_BUILD_TYPE=Debug
make
cd ..