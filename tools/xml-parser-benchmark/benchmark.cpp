// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Timo Sachsenberg $
// $Authors: OpenMS contributors $

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

#include <barrier>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#if defined(__APPLE__)
#include <mach-o/dyld.h>
#elif defined(__linux__)
#include <link.h>
#endif

namespace
{
  using Clock = std::chrono::steady_clock;

  std::string loadedLibrary()
  {
    std::string result;
#if defined(__APPLE__)
    for (std::uint32_t i = 0; i < _dyld_image_count(); ++i)
    {
      const std::string path(_dyld_get_image_name(i));
      if (std::filesystem::path(path).filename() == "libOpenMS.dylib") { result = path; }
    }
#elif defined(__linux__)
    dl_iterate_phdr([](dl_phdr_info* info, std::size_t, void* data)
    {
      const std::string path(info->dlpi_name);
      if (std::filesystem::path(path).filename().string().starts_with("libOpenMS.so"))
      {
        *static_cast<std::string*>(data) = path;
      }
      return 0;
    }, &result);
#endif
    return result;
  }

  // Hash every decoded coordinate/intensity, independently of XML text formatting.
  struct Digest final : OpenMS::Interfaces::IMSDataConsumer
  {
    std::uint64_t hash = 14695981039346656037ULL;
    std::uint64_t spectra = 0, chromatograms = 0, points = 0;

    void word(std::uint64_t value)
    {
      for (unsigned i = 0; i < 8; ++i)
      {
        hash = (hash ^ ((value >> (8 * i)) & 255)) * 1099511628211ULL;
      }
    }

    void text(const std::string& value)
    {
      word(value.size());
      for (const unsigned char character : value) { hash = (hash ^ character) * 1099511628211ULL; }
    }

    void consumeSpectrum(OpenMS::MSSpectrum& spectrum) override
    {
      ++spectra;
      points += spectrum.size();
      word(1);
      word(spectrum.size());
      text(spectrum.getNativeID());
      for (const auto& peak : spectrum)
      {
        word(std::bit_cast<std::uint64_t>(peak.getMZ()));
        word(std::bit_cast<std::uint64_t>(static_cast<double>(peak.getIntensity())));
      }
    }

    void consumeChromatogram(OpenMS::MSChromatogram& chromatogram) override
    {
      ++chromatograms;
      points += chromatogram.size();
      word(2);
      word(chromatogram.size());
      text(chromatogram.getNativeID());
      for (const auto& peak : chromatogram)
      {
        word(std::bit_cast<std::uint64_t>(peak.getRT()));
        word(std::bit_cast<std::uint64_t>(static_cast<double>(peak.getIntensity())));
      }
    }

    void setExpectedSize(std::size_t, std::size_t) override {}
    void setExperimentalSettings(const OpenMS::ExperimentalSettings&) override {}
    bool operator==(const Digest& other) const
    {
      return hash == other.hash && spectra == other.spectra &&
        chromatograms == other.chromatograms && points == other.points;
    }
  };

  struct Result
  {
    Digest digest;
    double open_seconds = 0;
    double read_seconds = 0;
    std::size_t index_spectra = 0, index_chromatograms = 0;
  };

  Result run(const std::string& mode, const std::string& path, std::size_t queries)
  {
    Result result;
    const auto begin = Clock::now();
    if (mode == "stream")
    {
      OpenMS::MzMLFile().transform(path, &result.digest, true, true);
    }
    else if (mode == "indexed" || mode == "index-open")
    {
      OpenMS::OnDiscMSExperiment experiment;
      if (!experiment.openFile(path, true))
      {
        throw std::runtime_error("Input has no usable mzML index");
      }
      const auto opened = Clock::now();
      result.open_seconds = std::chrono::duration<double>(opened - begin).count();
      result.index_spectra = experiment.getNrSpectra();
      result.index_chromatograms = experiment.getNrChromatograms();
      if (mode == "indexed")
      {
        if (experiment.getNrSpectra() == 0)
        {
          throw std::runtime_error("Indexed input has no spectra");
        }
        std::mt19937_64 random(20261019);
        for (std::size_t i = 0; i < queries; ++i)
        {
          auto spectrum = experiment.getSpectrum(random() % experiment.getNrSpectra());
          result.digest.consumeSpectrum(spectrum);
        }
      }
      result.read_seconds = std::chrono::duration<double>(Clock::now() - opened).count();
    }
    else
    {
      throw std::runtime_error("Expected stream, indexed, or index-open");
    }
    if (mode == "stream")
    {
      result.read_seconds = std::chrono::duration<double>(Clock::now() - begin).count();
    }
    return result;
  }

  void generate(const std::string& path, std::size_t count, std::size_t points, bool compressed)
  {
    OpenMS::PlainMSDataWritingConsumer writer(path);
    OpenMS::PeakFileOptions options;
    options.setCompression(compressed);
    options.setWriteIndex(true);
    writer.setOptions(options);
    writer.setExpectedSize(count, 1);
    writer.setExperimentalSettings(OpenMS::ExperimentalSettings());
    std::mt19937 random(20261019);
    for (std::size_t i = 0; i < count; ++i)
    {
      OpenMS::MSSpectrum spectrum;
      spectrum.setNativeID("scan=" + std::to_string(i + 1));
      spectrum.setRT(static_cast<double>(i) * 0.5);
      spectrum.setMSLevel(1);
      spectrum.setType(OpenMS::SpectrumSettings::SpectrumType::CENTROID);
      spectrum.setMetaValue("xml_poc_note", std::string("Unicode \xCE\xB2 & <escaped> \"value\""));
      spectrum.reserve(points);
      for (std::size_t p = 0; p < points; ++p)
      {
        OpenMS::Peak1D peak;
        peak.setMZ(100.0 + static_cast<double>(p) * 0.01 + static_cast<double>(i % 17) * 0.00001);
        peak.setIntensity(static_cast<float>(random() % 1000000) / 8.0F);
        spectrum.push_back(peak);
      }
      writer.consumeSpectrum(spectrum);
    }
    OpenMS::MSChromatogram chromatogram;
    chromatogram.setNativeID("TIC");
    for (std::size_t i = 0; i < count; ++i)
    {
      OpenMS::ChromatogramPeak peak;
      peak.setRT(static_cast<double>(i) * 0.5);
      peak.setIntensity(static_cast<float>(i % 101));
      chromatogram.push_back(peak);
    }
    writer.consumeChromatogram(chromatogram);
  }

  void verify(const std::string& path)
  {
    OpenMS::MSExperiment loaded;
    OpenMS::MzMLFile().load(path, loaded);
    Digest expected;
    for (auto& spectrum : loaded) { expected.consumeSpectrum(spectrum); }
    for (auto& chromatogram : loaded.getChromatograms()) { expected.consumeChromatogram(chromatogram); }
    const auto streamed = run("stream", path, 0).digest;
    OpenMS::OnDiscMSExperiment indexed;
    if (!indexed.openFile(path)) { throw std::runtime_error("Index verification failed"); }
    Digest actual;
    for (std::size_t i = 0; i < indexed.getNrSpectra(); ++i)
    {
      auto spectrum = indexed.getSpectrum(i);
      const auto& original = loaded.getSpectra().at(i);
      if (spectrum.getRT() != original.getRT() || spectrum.getMSLevel() != original.getMSLevel() ||
          spectrum.getFloatDataArrays() != original.getFloatDataArrays() ||
          spectrum.getIntegerDataArrays() != original.getIntegerDataArrays() ||
          spectrum.getStringDataArrays() != original.getStringDataArrays() ||
          spectrum.getMetaValue("xml_poc_note") != original.getMetaValue("xml_poc_note"))
      {
        throw std::runtime_error("Indexed spectrum metadata differ from full load");
      }
      actual.consumeSpectrum(spectrum);
    }
    for (std::size_t i = 0; i < indexed.getNrChromatograms(); ++i)
    {
      auto chromatogram = indexed.getChromatogram(i);
      actual.consumeChromatogram(chromatogram);
    }
    if (!(expected == streamed) || !(expected == actual))
    {
      throw std::runtime_error("Decoded data differ between load, streaming and indexed access");
    }
    std::cout << "{\"verified\":true,\"spectra\":" << expected.spectra
              << ",\"points\":" << expected.points << ",\"hash\":\"" << expected.hash << "\"}\n";
  }
}

int main(int argc, char** argv)
{
  std::cerr << "XML_BENCHMARK_MAIN_STARTED\n";
  try
  {
    if (argc < 3) { throw std::runtime_error("Usage: xml_parser_benchmark MODE FILE [WORKERS|SPECTRA] [QUERIES|POINTS] [ZLIB]"); }
    const std::string mode(argv[1]), path(argv[2]);
    const std::size_t workers = argc > 3 ? std::stoull(argv[3]) : 1;
    const std::size_t queries = argc > 4 ? std::stoull(argv[4]) : 1000;
    if (workers == 0) { throw std::runtime_error("Worker/spectrum count must be positive"); }
    if (mode == "generate")
    {
      generate(path, workers, queries, argc > 5 && std::string(argv[5]) == "zlib");
      return 0;
    }
    if (mode == "verify") { verify(path); return 0; }

    // Warm shared CV/parser initialization before independent worker jobs.
    OpenMS::MzMLFile warmup;
    OpenMS::MSExperiment empty;
    std::string warmup_xml;
    warmup.storeBuffer(warmup_xml, empty);
    warmup.loadBuffer(warmup_xml, empty);
    std::vector<Result> results(workers);
    std::vector<std::exception_ptr> errors(workers);
    std::barrier start(static_cast<std::ptrdiff_t>(workers + 1));
    std::vector<std::jthread> threads;
    for (std::size_t i = 0; i < workers; ++i)
    {
      threads.emplace_back([&, i]
      {
        start.arrive_and_wait();
        try { results[i] = run(mode, path, queries); }
        catch (...) { errors[i] = std::current_exception(); }
      });
    }
    const auto begin = Clock::now();
    start.arrive_and_wait();
    threads.clear();
    const double seconds = std::chrono::duration<double>(Clock::now() - begin).count();
    for (const auto& error : errors) { if (error) { std::rethrow_exception(error); } }
    for (const auto& result : results)
    {
      if (!(result.digest == results.front().digest)) { throw std::runtime_error("Worker digest mismatch"); }
      if (result.index_spectra != results.front().index_spectra ||
          result.index_chromatograms != results.front().index_chromatograms)
      {
        throw std::runtime_error("Worker index count mismatch");
      }
    }
    std::cout << std::setprecision(12) << "{\"mode\":\"" << mode << "\",\"workers\":" << workers
              << ",\"openms_library\":" << std::quoted(loadedLibrary())
              << ",\"seconds\":" << seconds << ",\"file_bytes\":" << std::filesystem::file_size(path)
              << ",\"spectra_per_worker\":" << results.front().digest.spectra
              << ",\"chromatograms_per_worker\":" << results.front().digest.chromatograms
              << ",\"points_per_worker\":" << results.front().digest.points
              << ",\"index_spectra\":" << results.front().index_spectra
              << ",\"index_chromatograms\":" << results.front().index_chromatograms
              << ",\"hash_per_worker\":\"" << results.front().digest.hash << "\",\"jobs\":[";
    for (std::size_t i = 0; i < workers; ++i)
    {
      if (i) { std::cout << ','; }
      std::cout << "{\"open_seconds\":" << results[i].open_seconds
                << ",\"read_seconds\":" << results[i].read_seconds << '}';
    }
    std::cout << "]}\n";
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
