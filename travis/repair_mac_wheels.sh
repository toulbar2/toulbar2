#!\bin/bash

cd $1
for wheel in *.whl
do
  mkdir $wheel.tmp
  cd $wheel.tmp
  unzip ../$wheel
  cd pytoulbar2/.dylibs
  libicuuc=`ls libicuuc*.dylib`
  libicui18n=`ls libicui18n*.dylib`
  libicudata=`ls libicudata*.dylib`
  slibicuuc=${libicuuc%.*}
  slibicuuc=${slibicuuc%.*}.dylib
  slibicui18n=${libicui18n%.*}
  slibicui18n=${slibicui18n%.*}.dylib
  slibicudata=${libicudata%.*}
  slibicudata=${slibicudata%.*}.dylib
  
  chmod u+w $libicuuc
  chmod u+w $libicui18n

  install_name_tool -change "@loader_path/$slibicudata" "@loader_path/$libicudata" $libicuuc
  install_name_tool -change "@loader_path/$slibicudata" "@loader_path/$libicudata" $libicui18n
  install_name_tool -change "@loader_path/$slibicuuc" "@loader_path/$libicuuc" $libicui18n

  cd ../..
  rm ../$wheel
  zip -r ../$wheel *
  cd ..
  rm -rf $wheel.tmp

done

