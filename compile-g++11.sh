name="MyStrategy"

if [ ! -f $name.cpp ]
then
    echo Unable to find $name.cpp > compilation.log
    exit 1
fi

rm -f $name

files=""

for i in *.cpp
do
    files="$files $i"
done

for i in model/*.cpp
do
    files="$files $i"
done

for i in csimplesocket/*.cpp
do
    files="$files $i"
done

g++ -std=c++11 -static -fno-optimize-sibling-calls -fno-strict-aliasing -DONLINE_JUDGE -D_LINUX -lm -s -x c++ -O2 -Wall -Wno-unknown-pragmas -o $name $files 2>compilation.log
