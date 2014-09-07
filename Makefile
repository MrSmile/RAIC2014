debug:
	g++ -std=c++11 -fno-optimize-sibling-calls -fno-strict-aliasing -DONLINE_JUDGE -D_LINUX -DDEBUG -lm -g -O0 -Wall -Wno-unknown-pragmas -o MyStrategy `./file-list.sh`

release:
	g++ -std=c++11 -fno-optimize-sibling-calls -fno-strict-aliasing -DONLINE_JUDGE -D_LINUX -lm -s -x c++ -O2 -Wall -Wno-unknown-pragmas -o MyStrategy `./file-list.sh`

opt: opt/MyStrategy opt/opt

opt/MyStrategy: MyStrategy.cpp
	g++ -std=c++11 -DONLINE_JUDGE -D_LINUX -DOPTIMIZING -march=native -Ofast -o opt/MyStrategy `./file-list.sh`

opt/opt: opt/opt.cpp
	g++ opt/opt.cpp -g -O0 -Wall -o opt/opt
