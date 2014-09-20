debug:
	g++ -std=c++11 -fno-strict-aliasing -DONLINE_JUDGE -D_LINUX -DDEBUG -lm -g -O0 -Wall -Wno-unknown-pragmas -o MyStrategy `./file-list.sh`

release:
	g++ -std=c++11 -fno-optimize-sibling-calls -fno-strict-aliasing -fno-omit-frame-pointer -DONLINE_JUDGE -D_LINUX -g -O2 -Wall -Wno-unknown-pragmas -o MyStrategy `./file-list.sh`

map: utils/map.cpp
	g++ -std=c++11 -march=native -g -Ofast utils/map.cpp -lpnglite -lz -o map

sweep: utils/sweep.cpp
	g++ -std=c++11 -march=native -g -Ofast utils/sweep.cpp -o sweep

opt: opt/MyStrategy opt/opt

opt/MyStrategy: MyStrategy.cpp
	g++ -std=c++11 -DONLINE_JUDGE -D_LINUX -DOPTIMIZING -march=native -Ofast -o opt/MyStrategy `./file-list.sh`

opt/opt: opt/opt.cpp
	g++ opt/opt.cpp -g -O0 -Wall -o opt/opt
