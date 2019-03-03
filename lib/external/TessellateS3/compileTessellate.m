if ispc
    mex tessellate_S3.cpp hypersphere.cpp tetramesh.cpp octetramesh.cpp util.cpp
elseif ismac
    %todo test
    %on macOS Mojave, calling "xcode-select --install" is required prior to
    %this working.
    mex CFLAGS='$CFLAGS -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/malloc/ -I/usr/include/malloc' tessellate_S3.cpp hypersphere.cpp tetramesh.cpp octetramesh.cpp util.cpp
else
    %linux
    mex CFLAGS='$CFLAGS' tessellate_S3.cpp hypersphere.cpp tetramesh.cpp octetramesh.cpp util.cpp
end