@ECHO OFF
SET ZLIB_VERSION=1.3
SET PCRE2_VERSION=10.45

IF NOT EXIST build_ci\libs (
  MKDIR build_ci\libs
)
CD build_ci\libs

@REM Install ZLib
IF NOT EXIST zlib-%ZLIB_VERSION%.zip (
  ECHO Downloading https://github.com/libarchive/zlib/archive/v%ZLIB_VERSION%.zip
  curl -L -o zlib-%ZLIB_VERSION%.zip https://github.com/libarchive/zlib/archive/v%ZLIB_VERSION%.zip || EXIT /b 1
)
IF NOT EXIST zlib-%ZLIB_VERSION% (
  ECHO Unpacking zlib-%ZLIB_VERSION%.zip
  C:\windows\system32\tar.exe -x -f zlib-%ZLIB_VERSION%.zip || EXIT /b 1
)
CD zlib-%ZLIB_VERSION%
cmake -B build || EXIT /b 1
cmake --build build --target ALL_BUILD --config Release || EXIT /b 1
cmake --build build --target RUN_TESTS --config Release || EXIT /b 1
cmake --build build --target INSTALL --config Release || EXIT /b 1

@REM Install PCRE2
IF NOT EXIST pcre2-%PCRE2_VERSION%.zip (
  ECHO Downloading https://github.com/PCRE2Project/pcre2/releases/download/pcre2-%PCRE2_VERSION%/pcre2-%PCRE2_VERSION%.zip
  curl -L -o pcre2-%PCRE2_VERSION%.zip https://github.com/PCRE2Project/pcre2/releases/download/pcre2-%PCRE2_VERSION%/pcre2-%PCRE2_VERSION%.zip || EXIT /b 1
)
IF NOT EXIST pcre2-%PCRE2_VERSION% (
  ECHO Unpacking pcre2-%PCRE2_VERSION%.zip
  C:\windows\system32\tar.exe -x -f pcre2-%PCRE2_VERSION%.zip || EXIT /b 1
)
CD pcre2-%PCRE2_VERSION%
cmake -B build -DPCRE2_BUILD_PCRE2GREP=OFF || EXIT /b 1
cmake --build build --target ALL_BUILD --config Release || EXIT /b 1
@REM cmake --build build --target RUN_TESTS --config Release || EXIT /b 1
cmake --build build --target INSTALL --config Release || EXIT /b 1

@EXIT /b 0
