// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
// $Authors: OpenMS contributors $

#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <string>

int main(int argc, char** argv)
{
  if (argc < 2) { return 2; }
  // A reused PID can inherit a zero-byte file from a crashed old libomp process.
  // libomp maps that existing file without checking its size, before main().
  // Remove only this process's old, empty, regular file owned by this user.
  const auto launched = std::time(nullptr);
  const std::string path = "/tmp/__KMP_REGISTERED_LIB_" + std::to_string(getpid());
  struct stat state {};
  if (lstat(path.c_str(), &state) == 0 && S_ISREG(state.st_mode) &&
      state.st_uid == getuid() && state.st_size == 0 && state.st_mtime < launched)
  {
    if (unlink(path.c_str()) == 0)
    {
      std::fprintf(stderr, "XML_BENCHMARK_REMOVED_STALE_KMP\n");
    }
  }
  // exec preserves the PID and avoids inheriting a Python parent's peak RSS.
  execv(argv[1], argv + 1);
  std::fprintf(stderr, "Cannot execute benchmark: %s\n", std::strerror(errno));
  return 127;
}
