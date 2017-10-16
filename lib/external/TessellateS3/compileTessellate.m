if ispc
    mex tessellate_S3.cpp hypersphere.cpp tetramesh.cpp octetramesh.cpp util.cpp
elseif ismac
    %todo test
    mex tessellate_S3.cpp hypersphere.cpp tetramesh.cpp octetramesh.cpp util.cpp
else
    %linux
    mex CFLAGS='$CFLAGS -std=c99' tessellate_S3.cpp hypersphere.cpp tetramesh.cpp octetramesh.cpp util.cpp
end