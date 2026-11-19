// M0: proves the test framework is wired into CTest.
#include "doctest/doctest.h"

TEST_CASE("smoke: doctest runs") {
    CHECK(1 + 1 == 2);
}
