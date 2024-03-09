@ECHO OFF
SET ZLIB_VERSION=1.3

IF NOT EXIST build_ci\libs (
  MKDIR build_ci\libs
)
CD build_ci\libs
IF NOT EXIST zlib-%ZLIB_VERSION%.zip (
  ECHO Downloading https://github.com/libarchive/zlib/archive/v%ZLIB_VERSION%.zip
  curl -L -o zlib-%ZLIB_VERSION%.zip https://github.com/libarchive/zlib/archive/v%ZLIB_VERSION%.zip || EXIT /b 1
)
IF NOT EXIST zlib-%ZLIB_VERSION% (
  ECHO Unpacking zlib-%ZLIB_VERSION%.zip
  C:\windows\system32\tar.exe -x -f zlib-%ZLIB_VERSION%.zip || EXIT /b 1
)
CD zlib-%ZLIB_VERSION%
cmake -G "Visual Studio 17 2022" . || EXIT /b 1
cmake --build . --target ALL_BUILD --config Release || EXIT /b 1
cmake --build . --target RUN_TESTS --config Release || EXIT /b 1
cmake --build . --target INSTALL --config Release || EXIT /b 1

@EXIT /b 0
