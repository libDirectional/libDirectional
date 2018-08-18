if ispc
    mex tessellate_S3.cpp hypersphere.cpp tetramesh.cpp octetramesh.cpp util.cpp
elseif ismac
    %todo test
    mex CFLAGS='-I/usr/include/malloc' tessellate_S3.cpp hypersphere.cpp tetramesh.cpp octetramesh.cpp util.cpp
else
    %linux
    mex CFLAGS='$CFLAGS' tessellate_S3.cpp hypersphere.cpp tetramesh.cpp octetramesh.cpp util.cpp
end