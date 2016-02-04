ROOT=`git rev-parse --show-toplevel`
cd $ROOT
FILE=$ROOT/toulbar2/src/ToulbarVersion.hpp
VERSION=`git describe --abbrev=0 --tags`
BRANCH=`git rev-parse --abbrev-ref HEAD`
HASH=`git log -1 --format=%h`

git diff --quiet HEAD > /dev/null 2>&1
if [ $? -ne 0 ]; then
  TAINTED="-tainted"
fi

echo "// Cmake generated version" > $FILE
echo "#define Toulbar_VERSION \"$VERSION-$BRANCH-$HASH$TAINTED\"" >> $FILE

