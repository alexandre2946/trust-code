// unique_token_scratch.cpp
#include <Kokkos_Core.hpp>
#include <iostream>

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  {
    using ExecSpace = Kokkos::DefaultExecutionSpace;
    using Token     = Kokkos::Experimental::UniqueToken<ExecSpace>;

    const int N = 1'000'000;

    // UniqueToken can provide either unique-per-thread or unique-per-team IDs
    // depending on how you use it. Here we use it "per thread".
    Token token(ExecSpace{});

    // How many possible concurrent token IDs?
    const int S = token.size();

    // Per-token scratch (e.g., each token has its own accumulator)
    Kokkos::View<double*, ExecSpace> scratch("scratch", S);
    Kokkos::deep_copy(scratch, 0.0);

    Kokkos::parallel_for(
      "use_unique_token_scratch",
      Kokkos::RangePolicy<ExecSpace>(0, N),
      KOKKOS_LAMBDA(const int i) {
        const int tid = token.acquire();      // get my unique slot
        scratch(tid) += 1.0;                  // no atomics needed
        token.release(tid);                   // give it back
      });

    // Reduce scratch on host
    auto scratch_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), scratch);
    double sum = 0.0;
    for (int t = 0; t < S; ++t) sum += scratch_h(t);

    std::cout << "N=" << N << " token.size()=" << S << " sum=" << sum << "\n";
  }
  Kokkos::finalize();
  return 0;
}
