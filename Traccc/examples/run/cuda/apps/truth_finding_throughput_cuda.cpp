// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/// @file
///
/// Throughput benchmark of the CUDA track finding on events with perfect
/// hit-to-track knowledge.
///
/// Every input event provides measurements, one spacepoint per measurement
/// (the usual @c -hits.csv / @c -measurements.csv /
/// @c -measurement-simhit-map.csv triplet), and a @c -seeds.csv file listing
/// the bottom/middle/top spacepoint indices of one truth seed per particle.
/// All events are loaded into device memory up front. The timed loop then
/// runs only the GPU seed parameter estimation and the combinatorial Kalman
/// filter, on several host threads with one CUDA stream each, and reports the
/// number of tracks found per second.

// Project include(s).
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/cuda/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/cuda/utils/make_magnetic_field.hpp"
#include "traccc/cuda/utils/stream_wrapper.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/examples/make_magnetic_field.hpp"
#include "traccc/examples/utils/printable.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_spacepoints.hpp"
#include "traccc/io/utils.hpp"
#include "traccc/options/detector.hpp"
#include "traccc/options/details/interface.hpp"
#include "traccc/options/input_data.hpp"
#include "traccc/options/magnetic_field.hpp"
#include "traccc/options/program_options.hpp"
#include "traccc/options/track_finding.hpp"
#include "traccc/options/track_propagation.hpp"
#include "traccc/seeding/detail/track_params_estimation_config.hpp"
#include "traccc/utils/logging.hpp"

// VecMem include(s).
#include <vecmem/memory/binary_page_memory_resource.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/copy.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>
#include <vecmem/utils/cuda/stream_wrapper.hpp>

// System include(s).
#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

/// Options specific to this benchmark
class throughput_options : public traccc::opts::interface {
 public:
  /// Number of host threads (each with its own CUDA stream)
  std::size_t threads = 1;
  /// Number of untimed warm-up event iterations
  std::size_t cold_run_events = 10;
  /// Number of timed event iterations
  std::size_t processed_events = 100;
  /// Minimum acceptable found-tracks / seeds ratio
  double min_track_fraction = 0.95;

  throughput_options() : interface("Truth Finding Throughput Options") {
    m_desc.add_options()(
        "cpu-threads",
        boost::program_options::value(&threads)->default_value(threads),
        "Number of host threads, each driving its own CUDA stream")(
        "cold-run-events",
        boost::program_options::value(&cold_run_events)
            ->default_value(cold_run_events),
        "Number of untimed warm-up event iterations")(
        "processed-events",
        boost::program_options::value(&processed_events)
            ->default_value(processed_events),
        "Number of timed event iterations (events are reused round-robin)")(
        "min-track-fraction",
        boost::program_options::value(&min_track_fraction)
            ->default_value(min_track_fraction),
        "Warn if found tracks / seeds falls below this value");
  }

  void read(const boost::program_options::variables_map&) override {
    if (threads == 0) {
      throw std::invalid_argument{"Must use cpu-threads > 0"};
    }
  }

  std::unique_ptr<traccc::configuration_printable> as_printable()
      const override {
    auto cat = std::make_unique<traccc::configuration_category>(m_description);
    cat->add_child(std::make_unique<traccc::configuration_kv_pair>(
        "CPU threads / CUDA streams", std::to_string(threads)));
    cat->add_child(std::make_unique<traccc::configuration_kv_pair>(
        "Cold run events", std::to_string(cold_run_events)));
    cat->add_child(std::make_unique<traccc::configuration_kv_pair>(
        "Processed events", std::to_string(processed_events)));
    cat->add_child(std::make_unique<traccc::configuration_kv_pair>(
        "Minimum track fraction", std::to_string(min_track_fraction)));
    return cat;
  }
};

/// Read the truth seed triplets of one event
void read_seeds(traccc::edm::seed_collection::host& seeds,
                const std::string& filename, std::size_t n_spacepoints) {
  std::ifstream in(filename);
  if (!in) {
    throw std::runtime_error("Could not open " + filename);
  }
  std::string line;
  std::getline(in, line);  // header
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    std::istringstream ss(line);
    unsigned int b = 0, m = 0, t = 0;
    char c1 = 0, c2 = 0;
    if (!(ss >> b >> c1 >> m >> c2 >> t) || c1 != ',' || c2 != ',') {
      throw std::runtime_error("Malformed seed line in " + filename + ": " +
                               line);
    }
    if (b >= n_spacepoints || m >= n_spacepoints || t >= n_spacepoints) {
      throw std::runtime_error("Seed spacepoint index out of range in " +
                               filename);
    }
    seeds.push_back({b, m, t, 0.f});
  }
}

/// One event, resident in device memory
struct device_event {
  traccc::edm::measurement_collection::buffer measurements;
  traccc::edm::spacepoint_collection::buffer spacepoints;
  traccc::edm::seed_collection::buffer seeds;
  unsigned int n_seeds = 0;
};

/// Per-thread state: stream, memory caches and algorithm instances
struct worker {
  worker(vecmem::memory_resource& device_mr,
         vecmem::memory_resource& pinned_host_mr,
         const traccc::track_params_estimation_config& est_cfg,
         const traccc::finding_config& find_cfg, const traccc::Logger& logger)
      : cached_device_mr(device_mr),
        cached_pinned_host_mr(pinned_host_mr),
        stream(vecmem_stream.stream()),
        copy(stream.cudaStream()),
        mr{cached_device_mr, &cached_pinned_host_mr},
        estimation(est_cfg, mr, copy, stream,
                   logger.clone("CudaSeedParamEstimation")),
        finding(find_cfg, mr, copy, stream, logger.clone("CudaFindingAlg")) {}

  vecmem::binary_page_memory_resource cached_device_mr;
  vecmem::binary_page_memory_resource cached_pinned_host_mr;
  vecmem::cuda::stream_wrapper vecmem_stream;
  traccc::cuda::stream_wrapper stream;
  vecmem::cuda::async_copy copy;
  traccc::memory_resource mr;
  traccc::cuda::seed_parameter_estimation_algorithm estimation;
  traccc::cuda::combinatorial_kalman_filter_algorithm finding;
};

}  // namespace

int main(int argc, char* argv[]) {
  std::unique_ptr<const traccc::Logger> ilogger = traccc::getDefaultLogger(
      "TracccTruthFindingThroughputCuda", traccc::Logging::Level::INFO);
  TRACCC_LOCAL_LOGGER(ilogger->clone());

  // Program options.
  traccc::opts::detector detector_opts;
  traccc::opts::magnetic_field bfield_opts;
  traccc::opts::input_data input_opts;
  traccc::opts::track_finding finding_opts;
  traccc::opts::track_propagation propagation_opts;
  throughput_options throughput_opts;
  traccc::opts::program_options program_opts{
      "Truth-Seeded Track Finding Throughput Using CUDA",
      {detector_opts, bfield_opts, input_opts, finding_opts, propagation_opts,
       throughput_opts},
      argc,
      argv,
      ilogger->cloneWithSuffix("Options")};

  // Memory resources.
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::host_memory_resource pinned_host_mr;
  vecmem::cuda::device_memory_resource device_mr;
  vecmem::cuda::copy sync_copy;

  // Detector and magnetic field.
  traccc::host_detector detector;
  traccc::io::read_detector(detector, host_mr, detector_opts.detector_file,
                            detector_opts.material_file,
                            detector_opts.grid_file);
  const traccc::detector_buffer detector_buffer =
      traccc::buffer_from_host_detector(detector, device_mr, sync_copy);
  const auto host_field = traccc::details::make_magnetic_field(bfield_opts);
  const auto device_field = traccc::cuda::make_magnetic_field(host_field);

  // Read all events and copy them to the device.
  std::vector<device_event> events;
  events.reserve(input_opts.events);
  std::uint64_t total_seeds = 0;
  std::uint64_t total_measurements = 0;
  for (std::size_t event = input_opts.skip;
       event < input_opts.skip + input_opts.events; ++event) {
    traccc::edm::spacepoint_collection::host spacepoints{host_mr};
    traccc::edm::measurement_collection::host measurements{host_mr};
    traccc::io::read_spacepoints(
        spacepoints, measurements, event, input_opts.directory,
        (input_opts.use_acts_geom_source ? &detector : nullptr), nullptr,
        nullptr, input_opts.format);

    traccc::edm::seed_collection::host seeds{host_mr};
    read_seeds(seeds,
               traccc::io::get_absolute_path(
                   (std::filesystem::path(input_opts.directory) /
                    traccc::io::get_event_filename(event, "-seeds.csv"))
                       .native()),
               spacepoints.size());

    device_event& dev = events.emplace_back();
    dev.n_seeds = static_cast<unsigned int>(seeds.size());
    dev.measurements = traccc::edm::measurement_collection::buffer(
        static_cast<unsigned int>(measurements.size()), device_mr);
    sync_copy.setup(dev.measurements)->wait();
    sync_copy(vecmem::get_data(measurements), dev.measurements)->wait();
    dev.spacepoints = traccc::edm::spacepoint_collection::buffer(
        static_cast<unsigned int>(spacepoints.size()), device_mr);
    sync_copy.setup(dev.spacepoints)->wait();
    sync_copy(vecmem::get_data(spacepoints), dev.spacepoints)->wait();
    dev.seeds = traccc::edm::seed_collection::buffer(dev.n_seeds, device_mr);
    sync_copy.setup(dev.seeds)->wait();
    sync_copy(vecmem::get_data(seeds), dev.seeds)->wait();

    total_seeds += dev.n_seeds;
    total_measurements += measurements.size();
  }
  if (events.empty()) {
    TRACCC_ERROR("No input events");
    return 1;
  }
  TRACCC_INFO("Loaded " << events.size() << " events, "
                        << total_measurements << " measurements, "
                        << total_seeds << " seeds ("
                        << static_cast<double>(total_seeds) /
                               static_cast<double>(events.size())
                        << " per event)");

  // Algorithm configuration.
  const traccc::track_params_estimation_config est_cfg;
  traccc::finding_config find_cfg(finding_opts);
  find_cfg.propagation = detray::propagation::config(propagation_opts);

  std::vector<std::unique_ptr<worker>> workers;
  for (std::size_t i = 0; i < throughput_opts.threads; ++i) {
    workers.push_back(std::make_unique<worker>(device_mr, pinned_host_mr,
                                               est_cfg, find_cfg, logger()));
  }

  // Run @c n_iterations event iterations spread over all workers. Iteration
  // @c i processes event @c i % n_events.
  auto run = [&](std::size_t n_iterations, std::uint64_t& n_tracks,
                 std::uint64_t& n_seeds) {
    std::atomic<std::size_t> next{0};
    std::atomic<std::uint64_t> tracks{0}, seeds{0};
    std::vector<std::thread> threads;
    std::vector<std::exception_ptr> errors(workers.size());
    for (std::size_t w = 0; w < workers.size(); ++w) {
      threads.emplace_back([&, w]() {
        try {
          worker& wk = *workers[w];
          for (std::size_t i = next++; i < n_iterations; i = next++) {
            const device_event& ev = events[i % events.size()];
            const auto params = wk.estimation(device_field, ev.measurements,
                                              ev.spacepoints, ev.seeds);
            const auto found = wk.finding(detector_buffer, device_field,
                                          ev.measurements, params);
            // Synchronises the stream, so that the output can be released.
            tracks += wk.copy.get_size(found.tracks);
            seeds += ev.n_seeds;
          }
        } catch (...) {
          errors[w] = std::current_exception();
        }
      });
    }
    for (auto& t : threads) {
      t.join();
    }
    for (const auto& e : errors) {
      if (e) {
        std::rethrow_exception(e);
      }
    }
    n_tracks = tracks;
    n_seeds = seeds;
  };

  std::uint64_t warm_tracks = 0, warm_seeds = 0;
  run(throughput_opts.cold_run_events, warm_tracks, warm_seeds);

  std::uint64_t n_tracks = 0, n_seeds = 0;
  const auto start = std::chrono::steady_clock::now();
  run(throughput_opts.processed_events, n_tracks, n_seeds);
  const auto stop = std::chrono::steady_clock::now();
  const double seconds = std::chrono::duration<double>(stop - start).count();

  const double fraction =
      n_seeds ? static_cast<double>(n_tracks) / static_cast<double>(n_seeds)
              : 0.;
  std::ostringstream report;
  report << std::fixed << std::setprecision(4);
  report << "\n  Events processed (timed):  " << throughput_opts.processed_events
         << "\n  Seeds (expected tracks):   " << n_seeds
         << "\n  Tracks found:              " << n_tracks
         << "\n  Found / expected:          " << fraction
         << "\n  Wall time:                 " << seconds << " s"
         << std::setprecision(1)
         << "\n  Throughput:                "
         << static_cast<double>(n_tracks) / seconds << " found tracks/s, "
         << static_cast<double>(n_seeds) / seconds << " seeds/s, "
         << static_cast<double>(throughput_opts.processed_events) / seconds
         << " events/s";
  TRACCC_INFO("Results:" << report.str());
  if (fraction < throughput_opts.min_track_fraction) {
    TRACCC_WARNING("Found / expected track ratio "
                   << fraction << " is below " <<
                   throughput_opts.min_track_fraction);
  }

  return 0;
}
